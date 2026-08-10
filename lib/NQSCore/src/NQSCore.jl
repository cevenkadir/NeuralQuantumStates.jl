"""
    NQSCore

The interface backbone of the neural-quantum-states stack, plus the two variational state types
and the log-derivative machinery they share.

This corresponds to NetKet's `netket/vqs` together with `netket/stats`. It exists so that
`NQSAnsatze`, `NQSSamplers`, and `NQSOptimisers` can be written against one agreed set of
abstract types and function names without ever depending on each other — the role `SciMLBase`
plays for SciML, or `ChainRulesCore` for the autodiff ecosystem.

# What lives here

- **Interfaces**: [`AbstractAnsatz`](@ref), [`AbstractVariationalState`](@ref),
  [`AbstractSampler`](@ref), [`AbstractPreconditioner`](@ref), and the verbs
  [`log_amplitude`](@ref), [`local_energy`](@ref), [`log_derivatives`](@ref),
  [`parameters`](@ref), [`setparameters!`](@ref), [`expect`](@ref),
  [`expect_and_grad`](@ref), [`sample`](@ref).
- **Two state types**, mirroring NetKet's `vqs/mc/` and `vqs/full_summ/`:
  [`MCState`](@ref) and [`FullSumState`](@ref). Both satisfy the same `expect` interface, so a
  model can be validated by exact summation before any Markov chain is involved.
- **[`Stats`](@ref)**: the single return type of every `expect` — mean, error of the mean,
  variance, integrated autocorrelation time, and split-R̂.
- **[`log_derivatives`](@ref)**: the `O_k = ∂ log ψ(s) / ∂θ_k` matrix, including the complex
  arithmetic that variational Monte Carlo needs and that is easy to get subtly wrong.
- **Reference implementations** with no approximations in them — [`LogStateVector`](@ref) and
  [`ExactSampler`](@ref) — so that the machinery can be tested against exact answers rather than
  against itself.

# Backends

Automatic differentiation is reached through DifferentiationInterface.jl, so no concrete backend
— ForwardDiff, Zygote, Enzyme, Reactant — is a dependency here: you load the one you want and
pass it as `backend`. `using NQSCore` pulls in no autodiff, no neural-network library and no GPU
code, and its test suite asserts as much.

Accelerators work the same way. Put the parameters on a device and the local energy follows them
there, because the reduction is written in terms that run unchanged on any array type. Loading
KernelAbstractions — a *weak* dependency — additionally moves the connected-configuration kernel
onto the device, so the batch never crosses the bus. Neither is required, and neither is named
in `[deps]`.

# Example

```julia
using NQSCore, ConnectedBasisConfigurations, SymBasis, DifferentiationInterface, ForwardDiff

spec, nsites = Spin(1 // 2), 4
b = basis(dof_object(spec), nsites)

a = LogStateVector(spec, nsites, b)
vs = FullSumState(a, init_parameters(a); backend=AutoForwardDiff())

expect(vs, H)                      # exact, zero error bar
E, ∇ = expect_and_grad(vs, H)      # ...and its gradient
```
"""
module NQSCore

using Random
using Random: AbstractRNG
using Statistics: mean, var

using ComponentArrays: ComponentArray, getaxes, getdata
using DifferentiationInterface

using ConnectedBasisConfigurations

include("interface.jl")
include("stats.jl")
include("log_derivatives.jl")
include("ansatz.jl")
include("sampler.jl")
include("states.jl")

# interfaces
export AbstractAnsatz, AbstractVariationalState, AbstractSampler, AbstractPreconditioner
export log_amplitude, local_energy, parameters, setparameters!, ansatz, samples
export expect, expect_and_grad, sample, precondition, Compiled

# statistics
export Stats, statistics, weighted_statistics, exact_stats
export integrated_autocorrelation, split_rhat

# log-derivatives
export log_derivatives, centered, flatten_parameters, match_parameter_shape

# reference implementations
export LogStateVector, init_parameters, n_parameters
export ExactSampler

# variational states
export FullSumState, MCState, probabilities, resample!, default_basis
export local_estimators, sample_weights, sampler_state

end # module NQSCore
