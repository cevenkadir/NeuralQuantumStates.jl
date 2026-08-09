```@meta
CurrentModule = NQSCore
```

# Variational states

A variational state is an [`AbstractAnsatz`](@ref), a set of parameters, and a way of averaging.
Two ship here, and they differ only in the last of those.

| Type | Averages by | Error bar |
|---|---|---|
| [`FullSumState`](@ref) | exact summation over the basis | none — the answer is exact |
| [`MCState`](@ref) | Monte Carlo sampling | standard error, autocorrelation-corrected |

Both satisfy the same [`expect`](@ref) interface, which is what makes the development loop work:
build a model against exact summation, where any disagreement is a bug rather than a
fluctuation, then scale up by swapping the state type and changing nothing else.

## What a state owes you

```julia
expect(state, operator)            # -> Stats
expect_and_grad(state, operator)   # -> (Stats, gradient shaped like the parameters)
local_energy(state, operator)      # -> the per-configuration energies behind that mean
local_estimators(state, operator)  # -> (; E, O, weights) for a preconditioner
samples(state)                     # -> the configurations currently being averaged over
```

`operator` is anything `ConnectedBasisConfigurations.connected_padded` accepts, so an `OpSum`
works directly. Passing a `compile`d operator avoids recompiling it on each call; inside an
optimization loop that is free to do and worth doing, though for a single large batch the
connected-configuration kernel dominates and the difference is not measurable.

## Which basis is summed over

A `FullSumState` needs to know what "the whole basis" means, and that belongs to the state
rather than to the ansatz — an ansatz is a functional form and has no opinion about which
configurations exist. Pass one explicitly to sum over a symmetry sector:

```julia
sector = basis(dof_object(spec), nsites, sym(TotalMagnetization(0 // 1, nsites), dof_object(spec)))
vs = FullSumState(a, θ; backend=AutoForwardDiff(), basis=sector)
```

Otherwise [`default_basis`](@ref) is consulted. There is deliberately no generic fallback: a
package that guessed here would be wrong for every ansatz defined on a sector, and would have to
depend on a basis library to do the guessing. `LogStateVector` knows its own basis, `NQSAnsatze`
defines the default for a `LuxAnsatz`, and anything else says so explicitly.

## Sampling, and why the chain stays warm

An `MCState` draws its samples when it is constructed and reuses them until something
invalidates them. Reuse is not an optimization — it is what makes [`expect`](@ref) and
[`expect_and_grad`](@ref) consistent with one another, since a gradient computed from different
samples than the energy does not correspond to it.

Changing the parameters marks the samples stale, and the next call to [`samples`](@ref) redraws.
The *sampler's* state survives that, which is the point of the contract in [`sample`](@ref):

```julia
drawn, state = sample(sampler, ansatz, θ, rng, previous_state)
```

A Metropolis sampler returns the configuration each chain finished on, and handing it back
resumes those chains instead of restarting them from scratch and re-paying the burn-in. A
parameter update moves the distribution only slightly, so the previous chain is already very
nearly equilibrated for the new parameters. Over a run this turns equilibration from a
per-iteration cost into a one-off one. A sampler that draws independently — [`ExactSampler`](@ref)
is the one here — has nothing to carry and returns `nothing`.

## Chain layout survives

Samples may arrive as a `(steps, chains)` matrix, and [`local_energy`](@ref) keeps that shape.
[`statistics`](@ref) needs it: with more than one chain the error of the mean comes from the
spread *between* chain means, which assumes only that the chains are independent, and split-R̂
compares them. Flattening happens inside the local-energy kernel rather than in the sampler,
leaving the chain layout the sampler's business.

## Running on a device

Put the parameters on a device and everything downstream follows them. `NQSCore` never names a
GPU package to do it: [`log_amplitude`](@ref) is the ansatz's business, and the local-energy
reduction is written so that it runs unchanged on whatever array type comes back.

```julia
using CUDA, cuDNN, Functors, KernelAbstractions

θ = fmap(CuArray, θ)                    # not gpu_device(), which demotes ComplexF64
vs = FullSumState(a, θ; backend=AutoZygote(), basis=b)
expect(vs, H)
```

With `KernelAbstractions` loaded, the connected configurations are computed on the device too,
by the kernel in `ConnectedBasisConfigurations` — the samples go up as packed integers, eight
bytes each, and the array `max_conn` times larger never crosses the bus at all. On a Quadro
GV100 with twelve sites and 4096 configurations, that host work and its transfer were two thirds
of a device `expect`.

Without `KernelAbstractions` everything still works: the connections are computed on the host
and moved, which is what the extra `_colocate` step in the kernel is for. Loading
`KernelAbstractions` on a machine with no accelerator changes nothing either — an `Array` is
deliberately not treated as a backend, because a portable kernel is not the way to beat a serial
loop over a few thousand samples.

The one thing worth doing by hand is the operator. The device path uploads it on every call
unless it is already resident, so a loop should hoist that exactly as it hoists `compile`:

```julia
H_dev = to_backend(flatten(H), CUDABackend())    # once
expect(vs, H_dev)                                # per step
```
