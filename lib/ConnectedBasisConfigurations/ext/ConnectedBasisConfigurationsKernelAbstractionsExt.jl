"""
The connected-configuration kernel, written once and run on any backend KernelAbstractions
supports.

It is the same algorithm as the CPU kernel in `flat.jl`, with the sample index coming from
`@index` instead of a loop — deliberately so, since the two are checked against each other slot
for slot. Two things had to exist first: [`FlatOperator`](@ref), which removed the nesting that
made a compiled operator impossible to upload, and the unchecked accessors below.

The work is parallel over samples — one thread owns one column of the output — but each thread
needs somewhere to expand a term into, so scratch is allocated per sample and passed in.

A second, smaller kernel ([`ConnectedBasisConfigurations.configurations!`](@ref)) unpacks the
packed states into the numeric array a network consumes. The two belong together: computing the
connections here and unpacking them on the host would put the larger array back on the wire.
"""
module ConnectedBasisConfigurationsKernelAbstractionsExt

using ConnectedBasisConfigurations
using ConnectedBasisConfigurations: FlatOperator
using KernelAbstractions
using SymBasis.DigitBase: BaseInt

# ------------------------------------------------------------------ digit access, unchecked

"""
    unchecked_read(state, position) -> Int
    unchecked_read(value, position, Val(base)) -> Int

The digit at `position`, without the bounds checks of the ordinary accessor.

Those checks `throw` with an interpolated message, which cannot be compiled into a GPU kernel —
it would allocate a string. They are also unnecessary: positions come from a compiled operator's
factors and digits from the columns of a `d × d` local matrix, so both are in range by
construction.

The arithmetic is written on the raw integer and `BaseInt` forwards to it, so the kernel can run
over either. That matters for XLA, whose tensors carry primitive element types only: a
`BaseInt` array is not one, and `Reactant.to_rarray` hands it back unconverted rather than
refusing it. The base travels as a `Val` because a raw integer does not carry it the way
`BaseInt{V,Ti,B}` does, and it has to be known at compile time for the `ispow2` branch to fold.
"""
@inline function unchecked_read(value::V, position::Integer, ::Val{B}) where {V<:Integer,B}
    if ispow2(B)
        bits = trailing_zeros(B)
        return Int((value >> ((position - 1) * bits)) & V(B - 1))
    else
        return Int(rem(div(value, V(B)^(position - 1)), V(B)))
    end
end

@inline unchecked_read(state::BaseInt{V,Ti,B}, position::Integer) where {V,Ti,B} =
    unchecked_read(state.value, position, Val(B))
@inline unchecked_read(state::BaseInt, position::Integer, ::Val) =
    unchecked_read(state, position)

"""
    unchecked_write(state, position, digit) -> BaseInt
    unchecked_write(value, position, digit, Val(base)) -> Integer

`state` with `position` set to `digit`. See [`unchecked_read`](@ref) for why it is unchecked, and
for why the raw-integer form is the one carrying the arithmetic.
"""
@inline function unchecked_write(
    value::V, position::Integer, digit::Integer, ::Val{B}
) where {V<:Integer,B}
    if ispow2(B)
        bits = trailing_zeros(B)
        shift = (position - 1) * bits
        mask = V(B - 1) << shift
        return (value & ~mask) | (V(digit) << shift)
    else
        p = V(B)^(position - 1)
        old = rem(div(value, p), V(B))
        # Modular arithmetic: the difference may wrap for an unsigned type, and adding a
        # wrapped difference still lands on the right value.
        return value + (V(digit) - old) * p
    end
end

@inline unchecked_write(
    state::BaseInt{V,Ti,B}, position::Integer, digit::Integer
) where {V,Ti,B} = BaseInt{V,Ti,B}(unchecked_write(state.value, position, digit, Val(B)))
@inline unchecked_write(state::BaseInt, position::Integer, digit::Integer, ::Val) =
    unchecked_write(state, position, digit)

"""
    digit_base(S) -> Val

The base an element type carries, for the entry points that can work it out from the states they
were handed. A raw integer carries no base, so a caller passing one has to say.
"""
digit_base(::Type{BaseInt{V,Ti,B}}) where {V,Ti,B} = Val(B)
digit_base(::Type{S}) where {S} = throw(ArgumentError(
    "states of type $S carry no digit base; pass `base=Val(b)` explicitly"
))

# ------------------------------------------------------------------------------- the kernel

@kernel function connected_kernel!(
    configs, mels, counts,
    @Const(term_start), @Const(factor_position), @Const(factor_col_start),
    @Const(colptr), @Const(outs), @Const(vals), @Const(states),
    scratch_states, scratch_vals,
    n_diagonal::Int, n_terms::Int, height::Int, base::Val,
)
    b = @index(Global)
    T = eltype(mels)

    @inbounds begin
        state = states[b]

        # Diagonal terms. A factor whose column is empty on the diagonal kills its whole term.
        diagonal = zero(T)
        for t in 1:n_diagonal
            v = one(T)
            alive = true
            for fi in term_start[t]:(term_start[t+1]-1)
                d = unchecked_read(state, factor_position[fi], base)
                cs = factor_col_start[fi]
                lo = colptr[cs+d]
                if lo >= colptr[cs+d+1]
                    alive = false
                    break
                end
                v *= vals[lo]
            end
            alive && (diagonal += v)
        end

        # Off-diagonal terms. Slot one is held for the diagonal and filled last.
        k = 1
        for t in (n_diagonal+1):n_terms
            n = 1
            src, dst = 1, 2
            scratch_states[1, src, b] = state
            scratch_vals[1, src, b] = one(T)

            alive = true
            for fi in term_start[t]:(term_start[t+1]-1)
                pos = factor_position[fi]
                cs = factor_col_start[fi]
                m = 0
                for j in 1:n
                    s = scratch_states[j, src, b]
                    amp = scratch_vals[j, src, b]
                    d = unchecked_read(s, pos, base)
                    for p in colptr[cs+d]:(colptr[cs+d+1]-1)
                        m += 1
                        scratch_states[m, dst, b] = unchecked_write(s, pos, outs[p], base)
                        scratch_vals[m, dst, b] = amp * vals[p]
                    end
                end
                if m == 0
                    alive = false
                    break
                end
                n = m
                src, dst = dst, src
            end

            if alive
                for j in 1:n
                    v = scratch_vals[j, src, b]
                    iszero(v) && continue
                    s′ = scratch_states[j, src, b]
                    if s′ == state
                        diagonal += v
                    else
                        k += 1
                        configs[k, b] = s′
                        mels[k, b] = v
                    end
                end
            end
        end

        # Close the gap: every retained row costs one wavefunction evaluation downstream.
        if iszero(diagonal)
            for j in 2:k
                configs[j-1, b] = configs[j, b]
                mels[j-1, b] = mels[j, b]
            end
            k -= 1
        else
            configs[1, b] = state
            mels[1, b] = diagonal
        end
        counts[b] = k

        # Inert padding: the sample itself with a zero matrix element.
        for j in (k+1):height
            configs[j, b] = state
            mels[j, b] = zero(T)
        end
    end
end

# ------------------------------------------------------------------------- unpacking states

"""One thread per entry: `out[i, b]` is the local value of digit `i` of state `b`."""
@kernel function configurations_kernel!(out, @Const(values), @Const(states), base::Val)
    i, b = @index(Global, NTuple)
    @inbounds out[i, b] = values[unchecked_read(states[b], i, base)+1]
end

function ConnectedBasisConfigurations.configurations!(
    out::AbstractMatrix, values::AbstractVector, states::AbstractArray{S},
    nsites::Integer, backend; base::Val=digit_base(S),
) where {S}
    n = length(states)
    size(out) == (nsites, n) || throw(DimensionMismatch(
        "out is $(size(out)); for $n states of $nsites sites it must be $((Int(nsites), n))"
    ))
    n == 0 && return out

    kernel = configurations_kernel!(backend)
    kernel(out, values, reshape(states, n), base; ndrange=(Int(nsites), n))
    KernelAbstractions.synchronize(backend)
    return out
end

# ---------------------------------------------------------------------------- entry points

"""
`op` may be a [`FlatOperator`](@ref) or a `NamedTuple` carrying the same fields, and the element
types are left open, because a compiled region can hold neither the struct nor the eltypes.

`FlatOperator{T,VI<:AbstractVector{Int32},VT<:AbstractVector{T}}` rejects a traced array outright:
`TracedRArray{Int32,1}` has element type `TracedRNumber{Int32}`, so it is not an
`AbstractVector{Int32}`. `Int` and `T` on `counts` and `mels` fail the same way. None of it is
load-bearing here — this launches a kernel, and the kernel reads fields and writes slots — so the
signature asks for the shape it checks and nothing more. The operator's own type keeps its
guarantees for everyone who builds and inspects one.
"""
function ConnectedBasisConfigurations.connected_padded!(
    configs::AbstractArray{S}, mels::AbstractArray, counts::AbstractArray,
    op::Union{FlatOperator,NamedTuple}, states::AbstractArray{S}, backend;
    base::Val=digit_base(S), workgroupsize::Integer=64,
) where {S}
    T = eltype(op.vals)
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

    n = length(states)
    n == 0 && return (; configs=configs, mels=mels, counts=counts)

    width = max(op.max_branch, 1)
    # `similar(states, ...)` and not `KernelAbstractions.allocate(backend, ...)`. They agree on
    # every backend KernelAbstractions knows about, and this one also works where the arrays are
    # traced rather than real — inside a Reactant compiled region, `states` is a tensor being
    # built and `similar` is how the documented examples ask for another one.
    scratch_states = similar(states, S, width, 2, n)
    scratch_vals = similar(mels, T, width, 2, n)

    kernel = connected_kernel!(backend, workgroupsize)
    kernel(
        reshape(configs, height, n), reshape(mels, height, n), reshape(counts, n),
        op.term_start, op.factor_position, op.factor_col_start,
        op.colptr, op.outs, op.vals, reshape(states, n),
        scratch_states, scratch_vals,
        op.n_diagonal, op.n_terms, height, base;
        ndrange=n,
    )
    KernelAbstractions.synchronize(backend)

    return (; configs=configs, mels=mels, counts=counts)
end

"""
    ConnectedBasisConfigurations.to_backend(op, backend) -> FlatOperator

Move a flat operator's arrays onto `backend`, so a kernel there can read them.

This is what [`FlatOperator`](@ref) exists to make possible: every field is a plain vector of
numbers, so moving the operator is moving seven arrays and nothing else.
"""
function ConnectedBasisConfigurations.to_backend(op::FlatOperator{T}, backend) where {T}
    move(v) = copyto!(KernelAbstractions.allocate(backend, eltype(v), length(v)), v)
    return FlatOperator(
        op.n_diagonal, op.n_terms,
        move(op.term_start), move(op.factor_position), move(op.factor_col_start),
        move(op.colptr), move(op.outs), move(op.vals),
        op.max_conn, op.max_branch,
    )
end

end # module ConnectedBasisConfigurationsKernelAbstractionsExt
