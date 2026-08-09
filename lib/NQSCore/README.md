<div align="center">

# NQSCore.jl

*The interface backbone of the neural-quantum-states stack for Julia*

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://cevenkadir.github.io/NeuralQuantumStates.jl/NQSCore/stable/) [![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://cevenkadir.github.io/NeuralQuantumStates.jl/NQSCore/dev/) [![Build Status](https://github.com/cevenkadir/NeuralQuantumStates.jl/actions/workflows/CI-NQSCore.yml/badge.svg?branch=main)](https://github.com/cevenkadir/NeuralQuantumStates.jl/actions/workflows/CI-NQSCore.yml?query=branch%3Amain) [![Coverage](https://codecov.io/gh/cevenkadir/NeuralQuantumStates.jl/branch/main/graph/badge.svg?flag=NQSCore)](https://codecov.io/gh/cevenkadir/NeuralQuantumStates.jl/flags?flag=NQSCore) [![Downloads](https://img.shields.io/badge/dynamic/json?url=http%3A%2F%2Fjuliapkgstats.com%2Fapi%2Fv1%2Fmonthly_downloads%2FNQSCore&query=total_requests&suffix=%2Fmonth&label=Downloads)](https://juliapkgstats.com/pkg/NQSCore) [![Total Downloads](https://img.shields.io/badge/dynamic/json?url=http%3A%2F%2Fjuliapkgstats.com%2Fapi%2Fv1%2Ftotal_downloads%2FNQSCore&query=total_requests&label=Total%20Downloads)](https://juliapkgstats.com/pkg/NQSCore)
</div>

A variational Monte Carlo run is assembled from four independent pieces — an ansatz, a sampler, an estimator, a preconditioner — and the usual way of writing one couples all four together, so that swapping the network means rewriting the optimizer. NQSCore.jl is the agreement that stops that happening: a handful of abstract types, the verbs that act on them, and the two variational state types everything else is written against. An ansatz is anything that answers `log_amplitude`; a sampler is anything that answers `sample`. Nothing here knows what a neural network is. This is the role `SciMLBase` plays for SciML or `ChainRulesCore` for the autodiff ecosystem, and it is why `NQSAnsatze`, `NQSSamplers` and `NQSOptimisers` can each be written and tested without depending on one another. Its dependencies are [ConnectedBasisConfigurations.jl](https://github.com/cevenkadir/NeuralQuantumStates.jl/tree/main/lib/ConnectedBasisConfigurations), [DifferentiationInterface.jl](https://github.com/JuliaDiff/DifferentiationInterface.jl) and [ComponentArrays.jl](https://github.com/jonniedie/ComponentArrays.jl) — no autodiff backend, no neural-network library, and no GPU code. A run on an accelerator needs none of that either: parameters on a device carry the local energy there by themselves, and loading [KernelAbstractions.jl](https://github.com/JuliaGPU/KernelAbstractions.jl) — a weak dependency, so it costs nothing to anyone who does not — moves the connected-configuration kernel there too, which is where two thirds of a device `expect` used to go.

## Key features
- **Two states, one contract**: [`FullSumState`](https://cevenkadir.github.io/NeuralQuantumStates.jl/NQSCore/stable/) sums exactly over the basis and [`MCState`](https://cevenkadir.github.io/NeuralQuantumStates.jl/NQSCore/stable/) samples, but both answer `expect` and `expect_and_grad` identically. A model is developed against exact summation, where any disagreement is a bug rather than a fluctuation, and then scaled up by changing one type.
- **Error bars you can believe**: every `expect` returns a `Stats` carrying the autocorrelation-corrected error of the mean, the variance, the integrated autocorrelation time and split-R̂. Reporting `std/sqrt(n)` on correlated samples is the standard way to converge confidently to a wrong answer, so the correction is applied here rather than left to the caller.
- **Complex arithmetic done once**: `log_derivatives` handles a complex log-amplitude with real parameters, and complex parameters in both the holomorphic and the general convention, in a single differentiation pass. Getting this subtly wrong is the classic source of gradients that almost work.
- **Gradients without the Jacobian**: the energy gradient is one contraction of the log-derivative matrix, so it is computed as the gradient of a scalar rather than by building the matrix and contracting it. On a reverse-mode backend that is one pass instead of one per sample — 17× faster for a small RBM under Zygote, and the gap widens with the sample count.
- **Warm Markov chains**: `sample` returns the sampler's own state alongside the draw, and an `MCState` keeps it across a parameter update. Burn-in is paid once per run instead of once per optimization step.
- **Any autodiff backend, none of them a dependency**: differentiation goes through DifferentiationInterface.jl, so `backend` is `AutoForwardDiff()`, `AutoZygote()`, `AutoEnzyme()` or whatever else you load yourself.

## The interface in one table
| You provide | By defining | And you get |
| ----------- | ----------- | ----------- |
| An ansatz | `log_amplitude(ansatz, θ, x)` | `expect`, `expect_and_grad`, `log_derivatives`, sampling |
| A sampler | `sample(sampler, ansatz, θ, rng, state)` | `MCState` and everything built on it |
| A preconditioner | `precondition(p, state, operator)` | A drop-in replacement for plain descent |

`ConnectedBasisConfigurations.connected_padded` supplies the connected configurations, so any operator it accepts works here without further arrangement.

## Installation
**Requirements**: Julia 1.11 or later.

To install the latest stable version, use the Julia package manager. Either use the Julia REPL package mode (by pressing `]`):
```julia
pkg> add NQSCore
```
or open the Julia REPL and run the following command:
```julia
julia> import Pkg; Pkg.add("NQSCore")
```

## Documentation
For detailed information on using this package, check out the [stable documentation](https://cevenkadir.github.io/NeuralQuantumStates.jl/NQSCore/stable/).

## Manual
- [Variational states](https://cevenkadir.github.io/NeuralQuantumStates.jl/NQSCore/stable/manual/states/)
- [Statistics](https://cevenkadir.github.io/NeuralQuantumStates.jl/NQSCore/stable/manual/statistics/)
- [Log-derivatives and gradients](https://cevenkadir.github.io/NeuralQuantumStates.jl/NQSCore/stable/manual/log_derivatives/)
- [Extending the interfaces](https://cevenkadir.github.io/NeuralQuantumStates.jl/NQSCore/stable/manual/interfaces/)

## Quick example
Build the exact ansatz on a four-site transverse-field Ising chain and ask for its energy. `LogStateVector` carries one parameter per basis state, so it can represent any state exactly — which is what makes it the reference the machinery is tested against:
```julia
julia> using NQSCore, ConnectedBasisConfigurations, SymBasis, OperatorAlgebra

julia> using DifferentiationInterface, ForwardDiff, Random

julia> spec, nsites = Spin(1 // 2), 4;

julia> b = basis(dof_object(spec), nsites);

julia> ops = local_operators(spec);

julia> H = OpSum(vcat([Op(2ops.sz, i) * Op(2ops.sz, mod1(i + 1, nsites)) for i in 1:nsites],
                      [2.0 * Op(2ops.sx, i) for i in 1:nsites]));

julia> a = LogStateVector(spec, nsites, b);

julia> vs = FullSumState(a, init_parameters(a, Xoshiro(0); scale=0.3); backend=AutoForwardDiff());

julia> expect(vs, H)
7.7427896 + -4.68375e-17im (variance=7.857, τ=1.0, R̂=1.0)
```
The error bar is exactly zero because nothing was sampled. The gradient comes back in the same shape as the parameters, whatever that shape is, so a descent step is a subtraction:
```julia
julia> for _ in 1:2000
           _, ∇ = expect_and_grad(vs, H)
           setparameters!(vs, parameters(vs) .- 0.05 .* ∇)
       end

julia> expect(vs, H)
-8.5431168 + -1.38051e-30im (variance=3.691e-28, τ=1.0, R̂=1.0)
```
Exact diagonalization gives `-8.543116820279426`, and the vanishing variance is the sharper statement: the energy variance is zero exactly on an eigenstate. Swapping in Monte Carlo changes the state type and nothing else, and now the error bar is real:
```julia
julia> mc = MCState(a, init_parameters(a, Xoshiro(0); scale=0.3),
                    ExactSampler(b, 100_000); backend=AutoForwardDiff(), rng=Xoshiro(1));

julia> expect(mc, H)
7.7258267 + 0.0141259im ± 0.00889 (variance=7.848, τ=1.01, R̂=0.99999)
```
`τ ≈ 1` says the samples are independent, which they are — `ExactSampler` enumerates the basis and draws from it directly, so it isolates Monte Carlo noise from Markov chain pathology. `NQSSamplers` provides the Metropolis and autoregressive samplers meant for real use.

## Scope
Ansatze, samplers and preconditioners are deliberately not here: this package defines what they must satisfy and provides one reference implementation of each — [`LogStateVector`](https://cevenkadir.github.io/NeuralQuantumStates.jl/NQSCore/stable/) and [`ExactSampler`](https://cevenkadir.github.io/NeuralQuantumStates.jl/NQSCore/stable/) — chosen because they are exact rather than because they are useful. Their job is to give the machinery something to be tested against that is not itself. Neural networks live in `NQSAnsatze`, Markov chains in `NQSSamplers`, stochastic reconfiguration in `NQSOptimisers`.

## Part of a larger ecosystem
NQSCore.jl is developed in the [NeuralQuantumStates.jl](https://github.com/cevenkadir/NeuralQuantumStates.jl) monorepo. It is the layer the rest of the stack agrees on, and it is usable on its own — it has no idea any particular neural network library exists, and its test suite asserts as much.

## Important notice
This project is still under active development. While it includes an extensive test suite and is developed with high scientific rigor, you should always benchmark your own code. Please report any issues you encounter via the [GitHub issue tracker](https://github.com/cevenkadir/NeuralQuantumStates.jl/issues/new).

## Citation
If you use this package in your work, we would appreciate the following reference as in [CITATION.bib](https://github.com/cevenkadir/NeuralQuantumStates.jl/blob/main/CITATION.bib).
