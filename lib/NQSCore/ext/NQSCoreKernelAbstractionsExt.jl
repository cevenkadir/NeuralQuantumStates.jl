"""
The device half of [`NQSCore.connections`](@ref): connected configurations computed where the
ansatz already is, instead of on the host and then shipped.

# Why this is worth an extension

Measured on a Quadro GV100, twelve sites and 4096 configurations, an `RBM(12, 4)`. Before this
path a device `expect` cost 5.897 ms and was made of the host connected-configuration kernel
(1.217 ms), unpacking its output and moving it to the device (2.716 ms), and the network and
reduction that were already running there. Two thirds of a "GPU" `expect` was host work and a
transfer.

With the connections computed on the device from packed states that are 8 bytes each, and
unpacked into the network's input array without ever being materialised as `Rational`s or
crossing the bus, the same `expect` costs **1.753 ms** — 132× the CPU's 231.6 ms, where it was
39×. `expect_and_grad` went from 26× to 44×. The transfer this removed, 2.099 ms per step, is
larger than the whole `expect` that is left. Agreement with the host is unchanged at a relative
6.6e-16.

The two kernels together take 70.6 µs of that 1.753 ms — 49.2 µs for the connections and
21.4 µs to unpack them — against 1.240 ms and a further 1.331 ms on the host. What dominates
now is the network, which is where the time ought to be.

# What stays out of it

`NQSCore` keeps its five hard dependencies. KernelAbstractions is a weak dependency, so a user
who never touches an accelerator loads nothing extra and the host path — which is the tested,
faster one for a CPU — is what runs. Loading KernelAbstractions alone does not divert anything
either; see [`NQSCore._device_backend`](@ref).

The reduction itself is not here. It lives in `NQSCore` and is already written in terms that
run unchanged on a device array, so the only thing this extension supplies is the pair of arrays
that feed it.
"""
module NQSCoreKernelAbstractionsExt

using ConnectedBasisConfigurations: ConnectedBasisConfigurations, FlatOperator
using KernelAbstractions
using NQSCore

"""
Whether `x` names a device, and which one.

An `Array` is answered `nothing` on purpose: `KernelAbstractions.get_backend` would call it
`CPU()`, and taking that as an instruction would mean that merely loading KernelAbstractions
swapped every host run onto a launch-per-batch kernel. That kernel exists to be portable, not
to beat a serial loop over a few thousand samples, and the choice belongs to the caller who put
their parameters somewhere rather than to an unrelated `using`.
"""
NQSCore._device_backend(::Array) = nothing
NQSCore._device_backend(x::AbstractArray) = KernelAbstractions.get_backend(x)

"""
Move `x` onto `backend`, or leave it if it is already somewhere other than host memory.

The same `Array`/`AbstractArray` split as above, and for the same reason: "not a plain `Array`"
is as close as this can get to "already on the device" without naming a GPU package.
"""
_resident(backend, x::Array) =
    copyto!(KernelAbstractions.allocate(backend, eltype(x), size(x)...), x)
_resident(::Any, x::AbstractArray) = x

"""
A flat operator resident on `backend`.

Three cases, and the middle one is the point: an operator that is *already* on the device is
handed back untouched, so a loop can upload once with `to_backend(flatten(H), backend)` and
stop paying for it. Anything else is flattened and uploaded here, which is correct but is seven
transfers per call.
"""
_device_operator(op::FlatOperator{<:Any,<:Array}, backend) =
    ConnectedBasisConfigurations.to_backend(op, backend)
_device_operator(op::FlatOperator, ::Any) = op
_device_operator(op, backend) = ConnectedBasisConfigurations.to_backend(
    ConnectedBasisConfigurations.flatten(op), backend
)

"""The table of local values, in the element type the network works in, resident on `backend`."""
_local_value_table(::Type{T}, spec, backend) where {T} =
    _resident(backend, collect(T, ConnectedBasisConfigurations.local_values(spec)))

function NQSCore._configurations(
    vs::NQSCore.AbstractVariationalState, states::AbstractVector, T,
    reference::AbstractArray, backend,
)
    a = NQSCore.ansatz(vs)
    nsites = NQSCore.n_sites(a)
    # The ansatz's answer when it has one, and the narrowest faithful type otherwise. The
    # kernel writes it straight out, so an ansatz that wants a complex batch never pays for a
    # second pass to widen a real one.
    S = T === nothing ? real(eltype(reference)) : T

    values = _local_value_table(S, NQSCore.dof(a), backend)
    x = KernelAbstractions.allocate(backend, S, nsites, length(states))
    ConnectedBasisConfigurations.configurations!(
        x, values, _resident(backend, states), nsites, backend
    )
    return x
end

function NQSCore._connections(
    vs::NQSCore.AbstractVariationalState, operator, states::AbstractVector,
    like::AbstractArray, backend,
)
    a = NQSCore.ansatz(vs)
    nsites = NQSCore.n_sites(a)

    op = _device_operator(operator, backend)
    n = length(states)
    height = ConnectedBasisConfigurations.max_conn_size(op)

    # The samples go up as packed integers — eight bytes each, against the `nsites × max_conn`
    # floats they expand into. Uploading them per call rather than caching them on the state is
    # deliberate: 4096 of them is 32 kB against a 48.9 µs kernel, and an earlier cache added on
    # the same reasoning measured as a loss and was reverted. Add one when a benchmark asks.
    dstates = _resident(backend, states)
    configs = KernelAbstractions.allocate(backend, eltype(states), height, n)
    mels = KernelAbstractions.allocate(backend, eltype(op), height, n)
    counts = KernelAbstractions.allocate(backend, Int, n)
    ConnectedBasisConfigurations.connected_padded!(
        configs, mels, counts, op, dstates, backend
    )

    # Full `max_conn` height, where the host path trims to the largest count it actually saw.
    # Trimming here would need a device-side `maximum` — a launch and a synchronising read —
    # and would then hand the network a non-contiguous view. The padded rows cost a wasted
    # column of the forward pass each; they contribute exactly zero, because a padded slot is
    # the sample itself with a zero matrix element.
    T = real(eltype(like))
    values = _local_value_table(T, NQSCore.dof(a), backend)
    x = KernelAbstractions.allocate(backend, T, nsites, height * n)
    ConnectedBasisConfigurations.configurations!(x, values, configs, nsites, backend)

    return x, mels
end

end # module NQSCoreKernelAbstractionsExt
