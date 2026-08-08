"""
    FlatOperator

A [`CompiledOperator`](@ref) with its nesting removed: the same data, in a handful of
rectangular arrays.

`CompiledOperator` is a vector of terms, each a vector of factors, each holding three more
vectors. That is the right shape for a CPU kernel and an impossible one for an accelerator —
`isbits` is false, so it cannot be uploaded, and chasing three levels of pointers is exactly
what a GPU is worst at. Flattening replaces the nesting with offsets:

- terms `1:n_diagonal` are diagonal, the rest off-diagonal, so the kernel needs no predicate;
- term `t` owns factors `term_start[t]:term_start[t+1]-1`;
- factor `f` owns column pointers `factor_col_start[f]:factor_col_start[f+1]-1`, and those
  pointers index `outs` and `vals` **absolutely**, so no second offset is needed inside the
  innermost loop.

Nothing is duplicated and nothing is padded; this is a change of layout, not of content. Every
field is a plain vector of `Int32` or of the matrix-element type, which makes the whole thing
`Adapt`-able to a device in one step.

Build one with [`flatten`](@ref). [`connected_padded!`](@ref) accepts it wherever it accepts a
`CompiledOperator`, and returns the same answer.
"""
struct FlatOperator{T,VI<:AbstractVector{Int32},VT<:AbstractVector{T}}
    n_diagonal::Int
    n_terms::Int
    term_start::VI          # n_terms + 1
    factor_position::VI     # n_factors
    factor_col_start::VI    # n_factors + 1
    colptr::VI              # absolute 1-based indices into outs/vals
    outs::VI                # zero-based output digits
    vals::VT
    max_conn::Int
    max_branch::Int
end

Base.eltype(::FlatOperator{T}) where {T} = T
max_conn_size(op::FlatOperator) = op.max_conn

function Base.show(io::IO, op::FlatOperator{T}) where {T}
    print(io, "FlatOperator{", T, "}(", op.n_terms, " terms, ", op.n_diagonal, " diagonal, ",
          "max_conn=", op.max_conn, ")")
end
Base.show(io::IO, ::MIME"text/plain", op::FlatOperator) = show(io, op)

"""
    flatten(operator) -> FlatOperator

Rewrite a compiled operator into flat arrays. Compiles first if it has not been compiled.

The diagonal terms are placed first so that a kernel can walk `1:n_diagonal` and then
`n_diagonal+1:n_terms` without asking which kind each term is.
"""
flatten(operator) = flatten(compile(operator))

function flatten(op::CompiledOperator{T}) where {T}
    terms = vcat(op.diagonal, op.offdiagonal)

    term_start = Int32[1]
    factor_position = Int32[]
    factor_col_start = Int32[1]
    colptr = Int32[]
    outs = Int32[]
    vals = T[]

    for term in terms
        for f in term.factors
            push!(factor_position, Int32(f.position))
            # `f.colptr` indexes this factor's own `outs`/`vals`; shifting by what is already
            # written makes it index the shared arrays instead, which removes an addition from
            # the innermost loop.
            offset = length(vals)
            for c in f.colptr
                push!(colptr, Int32(c + offset))
            end
            append!(outs, Int32.(f.outs))
            append!(vals, f.vals)
            push!(factor_col_start, Int32(length(colptr) + 1))
        end
        push!(term_start, Int32(length(factor_position) + 1))
    end

    return FlatOperator{T,Vector{Int32},Vector{T}}(
        length(op.diagonal), length(terms), term_start, factor_position, factor_col_start,
        colptr, outs, vals, op.max_conn, op.max_branch
    )
end

"""
Value of diagonal term `t` on `state`, or zero if any factor's column is empty on the diagonal.
"""
@inline function _flat_diagonal_value(op::FlatOperator{T}, t::Integer, state) where {T}
    v = one(T)
    @inbounds for fi in op.term_start[t]:(op.term_start[t+1]-1)
        d = read_digit(state, op.factor_position[fi])
        cs = op.factor_col_start[fi]
        lo = op.colptr[cs+d]
        lo < op.colptr[cs+d+1] || return zero(T)
        v *= op.vals[lo]
    end
    return v
end

"""
Expand off-diagonal term `t` on `state` into `dst_*`, returning how many entries it produced.

The two scratch pairs alternate as the term's factors are applied, exactly as the nested kernel
does. They are passed in rather than allocated so that a caller in a loop — or a device thread
with its own slice — owns the memory.
"""
function _flat_run_term!(
    op::FlatOperator{T}, t::Integer, state::S,
    a_states::AbstractVector{S}, a_vals::AbstractVector{T},
    b_states::AbstractVector{S}, b_vals::AbstractVector{T},
) where {S,T}
    cur_s, cur_v, alt_s, alt_v = a_states, a_vals, b_states, b_vals
    n = 1
    @inbounds cur_s[1] = state
    @inbounds cur_v[1] = one(T)

    @inbounds for fi in op.term_start[t]:(op.term_start[t+1]-1)
        pos = op.factor_position[fi]
        cs = op.factor_col_start[fi]
        m = 0
        for k in 1:n
            s = cur_s[k]
            amp = cur_v[k]
            d = read_digit(s, pos)
            for p in op.colptr[cs+d]:(op.colptr[cs+d+1]-1)
                m += 1
                alt_s[m] = write_digit(s, pos, op.outs[p])
                alt_v[m] = amp * op.vals[p]
            end
        end
        # An empty column annihilates the configuration, and with it the whole term.
        m == 0 && return (0, cur_s, cur_v)
        n = m
        cur_s, alt_s = alt_s, cur_s
        cur_v, alt_v = alt_v, cur_v
    end
    return (n, cur_s, cur_v)
end

"""
Fill column `b` of the output from `state`, returning how many entries it holds.

Identical in behaviour to the nested kernel's `_fill_column!`, including reserving slot one for
the diagonal and closing the gap when the diagonal turns out to be zero — every retained row
costs one full evaluation of the wavefunction downstream.
"""
function _flat_fill_column!(
    configs::AbstractMatrix{S}, mels::AbstractMatrix{T}, b::Integer,
    op::FlatOperator{T}, state::S,
    a_states, a_vals, b_states, b_vals,
) where {S,T}
    diagonal = zero(T)
    for t in 1:op.n_diagonal
        diagonal += _flat_diagonal_value(op, t, state)
    end

    k = 1
    @inbounds for t in (op.n_diagonal+1):op.n_terms
        n, res_s, res_v = _flat_run_term!(op, t, state, a_states, a_vals, b_states, b_vals)
        for j in 1:n
            v = res_v[j]
            iszero(v) && continue
            s′ = res_s[j]
            # A non-diagonal matrix can still map a digit to itself for a particular input.
            if s′ == state
                diagonal += v
            else
                k += 1
                configs[k, b] = s′
                mels[k, b] = v
            end
        end
    end

    if iszero(diagonal)
        @inbounds for j in 2:k
            configs[j-1, b] = configs[j, b]
            mels[j-1, b] = mels[j, b]
        end
        return k - 1
    end

    @inbounds configs[1, b] = state
    @inbounds mels[1, b] = diagonal
    return k
end

function connected_padded!(
    configs::AbstractArray{S}, mels::AbstractArray{T}, counts::AbstractArray{Int},
    op::FlatOperator{T}, states::AbstractArray{S}
) where {S,T}
    height = op.max_conn
    expected = (height, size(states)...)
    size(configs) == expected || throw(DimensionMismatch(
        "configs is $(size(configs)); for these states and operator it must be $expected"
    ))
    size(mels) == expected || throw(DimensionMismatch(
        "mels is $(size(mels)); for these states and operator it must be $expected"
    ))
    size(counts) == size(states) || throw(DimensionMismatch(
        "counts is $(size(counts)); it must have the shape of states, $(size(states))"
    ))

    flat = vec(states)
    n = length(flat)
    flat_configs = reshape(configs, height, n)
    flat_mels = reshape(mels, height, n)
    flat_counts = vec(counts)

    width = max(op.max_branch, 1)
    a_states = Vector{S}(undef, width)
    a_vals = Vector{T}(undef, width)
    b_states = Vector{S}(undef, width)
    b_vals = Vector{T}(undef, width)

    for b in 1:n
        s = @inbounds flat[b]
        k = _flat_fill_column!(
            flat_configs, flat_mels, b, op, s, a_states, a_vals, b_states, b_vals
        )
        @inbounds flat_counts[b] = k
        _pad_column!(flat_configs, flat_mels, b, k + 1, height, s)
    end

    return (; configs=configs, mels=mels, counts=counts)
end

function connected_padded(op::FlatOperator{T}, states::AbstractArray{S}) where {S,T}
    height = op.max_conn
    configs = Array{S}(undef, height, size(states)...)
    mels = Array{T}(undef, height, size(states)...)
    counts = Array{Int}(undef, size(states))
    connected_padded!(configs, mels, counts, op, states)

    flat_counts = vec(counts)
    max_conn = isempty(flat_counts) ? 0 : maximum(flat_counts)
    return _reshape_result(
        reshape(configs, height, length(flat_counts)),
        reshape(mels, height, length(flat_counts)),
        flat_counts, max_conn, size(states)
    )
end
