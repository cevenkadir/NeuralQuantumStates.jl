"""
The connected-configuration kernel, written once and run on any backend KernelAbstractions
supports — a CPU, CUDA, ROCm, Metal or oneAPI.

Two things had to change before the kernel could exist at all, and both live elsewhere:
[`FlatOperator`](@ref) removed the nesting that made a compiled operator impossible to upload,
and the accessors below removed the bounds checks that make a digit read impossible to compile.
What is left here is the same algorithm the CPU kernel runs, with the sample index coming from
`@index` instead of a loop.

The work is embarrassingly parallel over samples — one thread owns one column of the output and
touches nothing else — but each thread needs somewhere to expand a term into. That scratch is
allocated once, per sample rather than per thread, and passed in.

A second, much smaller kernel unpacks the packed states into the numeric array a network
consumes ([`ConnectedBasisConfigurations.configurations!`](@ref)). The two belong together
because they are the pair that keeps a batch on the device from end to end: computing the
connections there and then unpacking them on the host would put the larger of the two arrays
back on the wire.

# Measured

A transverse-field Ising chain of twelve sites, 4096 configurations, on a Quadro GV100: the host
kernel takes 1.217 ms and moving its results to the device a further 2.716 ms, against 48.9 µs
here. That is 24.9× against the host kernel alone and **80.4× against the host kernel plus the
transfer it replaces**, with output identical to the host kernel's in every slot.

Unpacking costs a further 21.4 µs on the device against 1.331 ms on the host, so the pair is
54× the host path it replaces. Wired into `NQSCore.local_energy`, that took a device `expect`
on the same model from 5.897 ms to 1.753 ms.
"""
module ConnectedBasisConfigurationsKernelAbstractionsExt

using ConnectedBasisConfigurations
using ConnectedBasisConfigurations: FlatOperator
using KernelAbstractions
using SymBasis.DigitBase: BaseInt

# ------------------------------------------------------------------ digit access, unchecked

"""
    unchecked_read(state, position) -> Int

The digit at `position`, without the bounds checks of the ordinary accessor.

Those checks `throw`, and a `throw` carrying an interpolated message cannot be compiled into a
GPU kernel — it needs to allocate a string, which is exactly what device code may not do. They
are also unnecessary here: every position comes from a compiled operator's factor, and every
digit written comes from a column of a `d × d` local matrix, so both are in range by
construction rather than by inspection.

Kept private to this extension. The checked accessors remain what everything else uses.
"""
@inline function unchecked_read(state::BaseInt{V,Ti,B}, position::Integer) where {V,Ti,B}
    if ispow2(B)
        bits = trailing_zeros(B)
        return Int((state.value >> ((position - 1) * bits)) & V(B - 1))
    else
        return Int(rem(div(state.value, V(B)^(position - 1)), V(B)))
    end
end

"""
    unchecked_write(state, position, digit) -> BaseInt

`state` with `position` set to `digit`. See [`unchecked_read`](@ref) for why it is unchecked.
"""
@inline function unchecked_write(
    state::BaseInt{V,Ti,B}, position::Integer, digit::Integer
) where {V,Ti,B}
    if ispow2(B)
        bits = trailing_zeros(B)
        shift = (position - 1) * bits
        mask = V(B - 1) << shift
        return BaseInt{V,Ti,B}((state.value & ~mask) | (V(digit) << shift))
    else
        p = V(B)^(position - 1)
        old = rem(div(state.value, p), V(B))
        # Modular arithmetic: the difference may wrap for an unsigned type, and adding a
        # wrapped difference still lands on the right value.
        return BaseInt{V,Ti,B}(state.value + (V(digit) - old) * p)
    end
end

# ------------------------------------------------------------------------------- the kernel

@kernel function connected_kernel!(
    configs, mels, counts,
    @Const(term_start), @Const(factor_position), @Const(factor_col_start),
    @Const(colptr), @Const(outs), @Const(vals), @Const(states),
    scratch_states, scratch_vals,
    n_diagonal::Int, n_terms::Int, height::Int,
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
                d = unchecked_read(state, factor_position[fi])
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
                    d = unchecked_read(s, pos)
                    for p in colptr[cs+d]:(colptr[cs+d+1]-1)
                        m += 1
                        scratch_states[m, dst, b] = unchecked_write(s, pos, outs[p])
                        scratch_vals[m, dst, b] = amp * vals[p]
                    end
                end
                if m == 0
                    alive = false
                    break
                end
                n = m
                # Ping-pong by index rather than by branching on a flag.
                src, dst = dst, src
            end

            if alive
                for j in 1:n
                    v = scratch_vals[j, src, b]
                    iszero(v) && continue
                    s′ = scratch_states[j, src, b]
                    # A non-diagonal matrix can still map a digit to itself for a given input.
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

        # Close the gap rather than leave an inert row: every retained row costs one full
        # evaluation of the wavefunction downstream.
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

        # Pad inert: the sample itself with a zero matrix element, so a consumer that reduces
        # over the whole column gets exactly zero from these and never `0 * Inf`.
        for j in (k+1):height
            configs[j, b] = state
            mels[j, b] = zero(T)
        end
    end
end

# ------------------------------------------------------------------------- unpacking states

"""
Write the physical local value of every digit of every state into a `(nsites, batch)` array.

One thread owns one entry, which is the whole of the parallelism here: no thread reads what
another writes, and the digit extraction is a shift and a mask. `values` is the `d`-element
table of local values, resident on the same backend, so the kernel indexes it rather than
knowing anything about degrees of freedom.
"""
@kernel function configurations_kernel!(out, @Const(values), @Const(states))
    i, b = @index(Global, NTuple)
    @inbounds out[i, b] = values[unchecked_read(states[b], i)+1]
end

function ConnectedBasisConfigurations.configurations!(
    out::AbstractMatrix{T}, values::AbstractVector{T}, states::AbstractArray,
    nsites::Integer, backend,
) where {T}
    n = length(states)
    size(out) == (nsites, n) || throw(DimensionMismatch(
        "out is $(size(out)); for $n states of $nsites sites it must be $((Int(nsites), n))"
    ))
    n == 0 && return out

    kernel = configurations_kernel!(backend)
    kernel(out, values, reshape(states, n); ndrange=(Int(nsites), n))
    KernelAbstractions.synchronize(backend)
    return out
end

# ---------------------------------------------------------------------------- entry points

function ConnectedBasisConfigurations.connected_padded!(
    configs::AbstractArray{S}, mels::AbstractArray{T}, counts::AbstractArray{Int},
    op::FlatOperator{T}, states::AbstractArray{S}, backend;
    workgroupsize::Integer=64,
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

    n = length(states)
    n == 0 && return (; configs=configs, mels=mels, counts=counts)

    width = max(op.max_branch, 1)
    scratch_states = KernelAbstractions.allocate(backend, S, width, 2, n)
    scratch_vals = KernelAbstractions.allocate(backend, T, width, 2, n)

    kernel = connected_kernel!(backend, workgroupsize)
    kernel(
        reshape(configs, height, n), reshape(mels, height, n), reshape(counts, n),
        op.term_start, op.factor_position, op.factor_col_start,
        op.colptr, op.outs, op.vals, reshape(states, n),
        scratch_states, scratch_vals,
        op.n_diagonal, op.n_terms, height;
        ndrange=n,
    )
    KernelAbstractions.synchronize(backend)

    return (; configs=configs, mels=mels, counts=counts)
end

"""
    ConnectedBasisConfigurations.to_backend(op, backend) -> FlatOperator

Move a flat operator's arrays onto `backend`, so a kernel there can read them.

This is the step [`FlatOperator`](@ref) exists to make possible: every field is a plain vector
of numbers, so moving the operator is moving seven arrays and nothing else.
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
