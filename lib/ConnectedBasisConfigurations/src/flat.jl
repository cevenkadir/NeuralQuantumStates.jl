"""
    FlatOperator

A [`CompiledOperator`](@ref) with its nesting removed: the same data, in a handful of
rectangular arrays.

A `CompiledOperator` is a vector of terms of factors of vectors, which `isbits` cannot describe
and so cannot be uploaded to a device. Flattening replaces the nesting with offsets:

- terms `1:n_diagonal` are diagonal, the rest off-diagonal, so the kernel needs no predicate;
- term `t` owns factors `term_start[t]:term_start[t+1]-1`;
- factor `f` owns column pointers `factor_col_start[f]:factor_col_start[f+1]-1`, and those
  pointers index `outs` and `vals` **absolutely**, saving an addition in the innermost loop.

Nothing is duplicated or padded; this is a change of layout, not of content. Build one with
[`flatten`](@ref); [`connected_padded!`](@ref) accepts it wherever it accepts a
`CompiledOperator` and returns the same answer.
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

Diagonal terms are placed first so a kernel can walk `1:n_diagonal` and then
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
            # written makes it index the shared arrays instead.
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

    # Two buffers to expand one term into, alternating as its factors are applied, reused
    # across the batch. From here down this is the same code in the same order as
    # `connected_kernel!` in the KernelAbstractions extension, which gives each thread its own
    # slice instead; the two are meant to be read side by side.
    width = max(op.max_branch, 1)
    cur_s, alt_s = Vector{S}(undef, width), Vector{S}(undef, width)
    cur_v, alt_v = Vector{T}(undef, width), Vector{T}(undef, width)

    @inbounds for b in 1:n
        state = flat[b]

        # Diagonal terms. A factor whose column is empty on the diagonal kills its whole term.
        diagonal = zero(T)
        for t in 1:op.n_diagonal
            v = one(T)
            alive = true
            for fi in op.term_start[t]:(op.term_start[t+1]-1)
                d = read_digit(state, op.factor_position[fi])
                cs = op.factor_col_start[fi]
                lo = op.colptr[cs+d]
                if lo >= op.colptr[cs+d+1]
                    alive = false
                    break
                end
                v *= op.vals[lo]
            end
            alive && (diagonal += v)
        end

        # Off-diagonal terms. Slot one is held for the diagonal and filled last.
        k = 1
        for t in (op.n_diagonal+1):op.n_terms
            src_s, src_v, dst_s, dst_v = cur_s, cur_v, alt_s, alt_v
            m = 1
            src_s[1] = state
            src_v[1] = one(T)

            alive = true
            for fi in op.term_start[t]:(op.term_start[t+1]-1)
                pos = op.factor_position[fi]
                cs = op.factor_col_start[fi]
                w = 0
                for j in 1:m
                    s = src_s[j]
                    amp = src_v[j]
                    d = read_digit(s, pos)
                    for q in op.colptr[cs+d]:(op.colptr[cs+d+1]-1)
                        w += 1
                        dst_s[w] = write_digit(s, pos, op.outs[q])
                        dst_v[w] = amp * op.vals[q]
                    end
                end
                if w == 0
                    alive = false
                    break
                end
                m = w
                src_s, dst_s = dst_s, src_s
                src_v, dst_v = dst_v, src_v
            end

            if alive
                for j in 1:m
                    v = src_v[j]
                    iszero(v) && continue
                    s′ = src_s[j]
                    if s′ == state
                        diagonal += v
                    else
                        k += 1
                        flat_configs[k, b] = s′
                        flat_mels[k, b] = v
                    end
                end
            end
        end

        # Close the gap: every retained row costs one wavefunction evaluation downstream.
        if iszero(diagonal)
            for j in 2:k
                flat_configs[j-1, b] = flat_configs[j, b]
                flat_mels[j-1, b] = flat_mels[j, b]
            end
            k -= 1
        else
            flat_configs[1, b] = state
            flat_mels[1, b] = diagonal
        end
        flat_counts[b] = k

        # Inert padding: the sample itself with a zero matrix element.
        for j in (k+1):height
            flat_configs[j, b] = state
            flat_mels[j, b] = zero(T)
        end
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

"""
    to_backend(operator, backend) -> FlatOperator

Move a flat operator onto a KernelAbstractions backend.

Defined by the extension that KernelAbstractions activates; without it there is no backend to
move to, and this says so rather than failing later and less clearly.
"""
function to_backend end
