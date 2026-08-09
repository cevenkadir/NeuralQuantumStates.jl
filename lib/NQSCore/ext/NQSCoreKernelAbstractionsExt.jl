"""
The device half of [`NQSCore.connections`](@ref) and [`NQSCore.configurations_of`](@ref):
connected configurations computed where the ansatz already is, rather than on the host and then
shipped. Samples go up as packed integers, eight bytes each, and the `nsites × max_conn` floats
they expand into never cross the bus.

KernelAbstractions is a weak dependency, so a user who never touches an accelerator loads
nothing extra and gets the host path, which is the faster one for a CPU. Loading it alone does
not divert anything either; see `NQSCore.device_backend` below.

The reduction is not here. It lives in `NQSCore` and already runs unchanged on a device array,
so all this extension supplies is the pair of arrays that feed it.
"""
module NQSCoreKernelAbstractionsExt

using ConnectedBasisConfigurations: ConnectedBasisConfigurations, FlatOperator
using KernelAbstractions
using NQSCore

"""
Whether `x` names a device, and which one.

An `Array` is answered `nothing` on purpose. `get_backend` would call it `CPU()`, and acting on
that would let an unrelated `using` swap every host run onto a launch-per-batch kernel — which
exists to be portable, not to beat a serial loop over a few thousand samples.
"""
NQSCore.device_backend(::Array) = nothing
NQSCore.device_backend(x::AbstractArray) = KernelAbstractions.get_backend(x)

"""
Move `x` onto `backend`, or leave it if it is already somewhere other than host memory. The same
`Array`/`AbstractArray` split as above: "not a plain `Array`" is as close as this gets to
"already on the device" without naming a GPU package.
"""
_resident(backend, x::Array) =
    copyto!(KernelAbstractions.allocate(backend, eltype(x), size(x)...), x)
_resident(::Any, x::AbstractArray) = x

"""
A flat operator resident on `backend`.

One already on the device is handed back untouched, so a loop can upload once with
`to_backend(flatten(H), backend)`; anything else is flattened and uploaded here, at seven
transfers per call.
"""
_device_operator(op::FlatOperator{<:Any,<:Array}, backend) =
    ConnectedBasisConfigurations.to_backend(op, backend)
_device_operator(op::FlatOperator, ::Any) = op
_device_operator(op, backend) = ConnectedBasisConfigurations.to_backend(
    ConnectedBasisConfigurations.flatten(op), backend
)

function NQSCore.device_configurations(
    vs::NQSCore.AbstractVariationalState, states::AbstractVector, T,
    reference::AbstractArray, backend,
)
    a = NQSCore.ansatz(vs)
    nsites = NQSCore.n_sites(a)
    # The kernel writes the requested type straight out, so an ansatz that wants a complex batch
    # never pays a second pass to widen a real one.
    S = T === nothing ? real(eltype(reference)) : T

    values = _resident(backend, collect(S, ConnectedBasisConfigurations.local_values(
        NQSCore.dof(a)
    )))
    x = KernelAbstractions.allocate(backend, S, nsites, length(states))
    ConnectedBasisConfigurations.configurations!(
        x, values, _resident(backend, states), nsites, backend
    )
    return x
end

function NQSCore.device_connections(
    vs::NQSCore.AbstractVariationalState, operator, states::AbstractVector,
    like::AbstractArray, backend,
)
    a = NQSCore.ansatz(vs)
    nsites = NQSCore.n_sites(a)

    op = _device_operator(operator, backend)
    n = length(states)
    height = ConnectedBasisConfigurations.max_conn_size(op)

    # Uploaded per call rather than cached on the state: packed samples are eight bytes each,
    # which is small beside the kernel that consumes them.
    dstates = _resident(backend, states)
    configs = KernelAbstractions.allocate(backend, eltype(states), height, n)
    mels = KernelAbstractions.allocate(backend, eltype(op), height, n)
    counts = KernelAbstractions.allocate(backend, Int, n)
    ConnectedBasisConfigurations.connected_padded!(
        configs, mels, counts, op, dstates, backend
    )

    # Full `max_conn` height, where the host path trims to the largest count it saw. Trimming
    # here would need a device-side `maximum` — a launch and a synchronising read — and would
    # then hand the network a non-contiguous view. The extra rows are inert padding.
    T = real(eltype(like))
    values = _resident(backend, collect(T, ConnectedBasisConfigurations.local_values(
        NQSCore.dof(a)
    )))
    x = KernelAbstractions.allocate(backend, T, nsites, height * n)
    ConnectedBasisConfigurations.configurations!(x, values, configs, nsites, backend)

    return x, mels
end

end # module NQSCoreKernelAbstractionsExt
