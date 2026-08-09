```@meta
CurrentModule = ConnectedBasisConfigurations
```

# Running the kernel on a device

A [`CompiledOperator`](@ref) is a vector of terms holding vectors of factors holding three more
vectors. That is the right shape for a CPU kernel and an impossible one for an accelerator:
`isbits` is false, so it cannot be uploaded at all, and chasing three levels of pointers is what
a GPU is worst at.

[`flatten`](@ref) rewrites it as a [`FlatOperator`](@ref) — the same content, in seven plain
vectors of `Int32` and matrix elements, with the nesting replaced by offsets. [`to_backend`](@ref)
then moves those seven arrays onto a KernelAbstractions backend, and the kernel there reads them
exactly as the host kernel reads the nested form.

Two methods appear once `KernelAbstractions` is loaded:

- [`connected_padded!`](@ref) with a trailing `backend` argument, which computes the connections
  on that backend;
- [`configurations!`](@ref), which unpacks the packed results into the numeric array a network
  consumes — writing the float type directly, since the host [`configurations`](@ref) returns
  `Rational` for a spin and no accelerator can hold one.

The pair is what keeps a batch on the device from end to end. Computing the connections there
and unpacking them on the host would put the larger of the two arrays back on the wire, which is
the transfer the whole exercise exists to remove.

```julia
using ConnectedBasisConfigurations, KernelAbstractions, CUDA

op = to_backend(flatten(H), CUDABackend())     # once, outside the loop
states = CuArray(basis(dof_object(spec), nsites).states)

h, n = max_conn_size(op), length(states)
configs = similar(states, h, n)
mels = CuArray{Float64}(undef, h, n)
counts = CuArray{Int}(undef, n)
connected_padded!(configs, mels, counts, op, states, CUDABackend())

values = CuArray(collect(Float64, local_values(spec)))
x = CuArray{Float64}(undef, nsites, h * n)
configurations!(x, values, configs, nsites, CUDABackend())
```

Uploading the operator is seven transfers, so a variational loop should do it once and reuse the
result rather than flattening on every step — the device counterpart of hoisting
[`compile`](@ref) out of the loop.

!!! note "The device path keeps the full static bound"
    [`connected_padded`](@ref) on the host trims its output to the largest connection count it
    actually saw. The device kernel fills `max_conn` rows and leaves it at that: trimming would
    need a device-side `maximum` — a launch and a synchronising read — and would then hand the
    consumer a non-contiguous view. The extra rows are inert, being the sample itself with a
    zero matrix element, so a consumer that reduces over the whole column is unaffected.

## Index

```@index
Pages = ["device.md"]
```

```@autodocs
Modules = [ConnectedBasisConfigurations]
Pages = ["flat.jl"]
```
