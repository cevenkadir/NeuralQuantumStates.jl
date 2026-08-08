"""
    NQSCore.LocalEnergyKernel(source, compiled, states)

A compiled operator together with the buffers its connected-configuration kernel writes into.

Two costs disappear when this is kept across calls rather than rebuilt. The operator is compiled
**once** instead of on every `connected_padded` call, and the `(max_conn, batch)` configuration
and matrix-element arrays are written in place instead of freshly allocated — which is what
`ConnectedBasisConfigurations.connected_padded!` exists for and what its docstring means by "the
form to use inside an optimization loop".

`source` is the operator as the caller handed it over, kept only so that a later call can ask
"same operator?" by identity. An optimization loop passes the same `OpSum` object every step, so
that test succeeds and nothing is rebuilt; a caller that builds a fresh operator each time gets
a fresh kernel and is no worse off than before.

Held by the variational states, and rebuilt whenever the operator or the batch size changes. A
`CompiledSector` has no in-place kernel in ConnectedBasisConfigurations, so a sector operator
falls back to the allocating path instead of caching anything.
"""
struct LocalEnergyKernel{R,O,S,T}
    source::R
    operator::O
    configs::Matrix{S}
    mels::Matrix{T}
    counts::Vector{Int}
end

function LocalEnergyKernel(
    source, operator::ConnectedBasisConfigurations.CompiledOperator,
    states::AbstractVector{S}
) where {S}
    height = ConnectedBasisConfigurations.max_conn_size(operator)
    n = length(states)
    T = eltype(operator)
    return LocalEnergyKernel(
        source, operator,
        Matrix{S}(undef, height, n), Matrix{T}(undef, height, n), Vector{Int}(undef, n)
    )
end

"""
The compiled form to build a kernel around, or `nothing` when there is no in-place kernel for it.
"""
_kernel_operator(op::ConnectedBasisConfigurations.CompiledOperator) = op
_kernel_operator(::ConnectedBasisConfigurations.CompiledSector) = nothing
_kernel_operator(op) = ConnectedBasisConfigurations.compile(op)

"""Whether `k` was built for this operator and this batch size."""
_fits(k::LocalEnergyKernel, source, n::Integer) = k.source === source && length(k.counts) == n

"""
Fill the kernel's buffers for `states` and return the filled region.

`connected_padded!` always pads to the operator's static `max_conn`, while the allocating form
trims to the batch's actual maximum. The views returned here restore that trimming, because
every retained row costs one full evaluation of the wavefunction downstream — the rows above the
batch maximum are inert, but paying for them would hand back the allocation saving as ansatz
time.
"""
function _fill!(k::LocalEnergyKernel, states::AbstractVector)
    ConnectedBasisConfigurations.connected_padded!(
        k.configs, k.mels, k.counts, k.operator, states
    )
    m = isempty(k.counts) ? 0 : maximum(k.counts)
    return (;
        configs=view(k.configs, 1:m, :),
        mels=view(k.mels, 1:m, :),
        counts=k.counts,
    )
end

"""
    energy_kernel(state) / energy_kernel!(state, kernel)

Where a variational state keeps its [`LocalEnergyKernel`](@ref NQSCore.LocalEnergyKernel).

The generic methods decline to cache, so a state type defined elsewhere keeps working and simply
allocates per call. The two shipped here override them.
"""
energy_kernel(::AbstractVariationalState) = nothing
energy_kernel!(::AbstractVariationalState, kernel) = kernel

"""
    connections(state, operator, states) -> (; configs, mels, counts)

The connected configurations of every sample, reusing the state's buffers where possible.
"""
function connections(vs::AbstractVariationalState, operator, states::AbstractVector)
    cached = energy_kernel(vs)
    cached !== nothing && _fits(cached, operator, length(states)) && return _fill!(cached, states)

    compiled = _kernel_operator(operator)
    # A sector operator has no in-place kernel; nothing to cache, so do not try.
    compiled === nothing &&
        return ConnectedBasisConfigurations.connected_padded(operator, states)

    kernel = LocalEnergyKernel(operator, compiled, states)
    energy_kernel!(vs, kernel)
    return _fill!(kernel, states)
end

"""
Put `mels` on the same device as `like`, leaving it alone when both are already in host memory.

On the CPU nothing needs to happen: broadcasting a real matrix-element array against complex
log-amplitudes promotes elementwise for free, and materializing a converted copy would be a
full `(max_conn, batch)` allocation bought for nothing. When the log-amplitudes live on a device
— because the ansatz put them there — the matrix elements have to follow, and `copyto!` into a
`similar` of the log-amplitudes does that using nothing but Base, which is why this file needs
no GPU dependency to be GPU-ready.
"""
_colocate(::Array, mels::AbstractArray) = mels
_colocate(like::AbstractArray, mels::AbstractArray) =
    copyto!(similar(like, eltype(mels), size(mels)), mels)
