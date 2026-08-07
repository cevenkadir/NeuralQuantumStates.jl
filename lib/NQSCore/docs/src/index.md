```@meta
CurrentModule = NQSCore
```

# NQSCore.jl

*The interface backbone of the neural-quantum-states stack.*

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://cevenkadir.github.io/NeuralQuantumStates.jl/NQSCore/stable/) [![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://cevenkadir.github.io/NeuralQuantumStates.jl/NQSCore/dev/) [![Build Status](https://github.com/cevenkadir/NeuralQuantumStates.jl/actions/workflows/CI-NQSCore.yml/badge.svg?branch=main)](https://github.com/cevenkadir/NeuralQuantumStates.jl/actions/workflows/CI-NQSCore.yml?query=branch%3Amain) [![Coverage](https://codecov.io/gh/cevenkadir/NeuralQuantumStates.jl/branch/main/graph/badge.svg?flag=NQSCore)](https://codecov.io/gh/cevenkadir/NeuralQuantumStates.jl/flags?flag=NQSCore)

A variational Monte Carlo run is assembled from four independent pieces — an ansatz, a sampler,
an estimator, a preconditioner — and the usual way of writing one couples all four together, so
that swapping the network means rewriting the optimizer. NQSCore is the agreement that stops
that happening: a handful of abstract types, the verbs that act on them, and the two variational
state types everything else is written against.

An ansatz is anything that answers [`log_amplitude`](@ref); a sampler is anything that answers
[`sample`](@ref). Nothing here knows what a neural network is. This is the role `SciMLBase`
plays for SciML, or `ChainRulesCore` for the autodiff ecosystem, and it is why `NQSAnsatze`,
`NQSSamplers` and `NQSOptimisers` can each be written and tested without depending on one
another.

## Key features

- **Two states, one contract** — [`FullSumState`](@ref) sums exactly over the basis and
  [`MCState`](@ref) samples, but both answer [`expect`](@ref) and [`expect_and_grad`](@ref)
  identically. A model is developed against exact summation, where any disagreement is a bug
  rather than a fluctuation, and then scaled up by changing one type.
- **Error bars you can believe** — every `expect` returns a [`Stats`](@ref) carrying the
  autocorrelation-corrected error of the mean, the variance, the integrated autocorrelation time
  and split-R̂.
- **Complex arithmetic done once** — [`log_derivatives`](@ref) handles a complex log-amplitude
  with real parameters, and complex parameters in both conventions, in a single differentiation
  pass.
- **Gradients without the Jacobian** — the energy gradient is one contraction of the
  log-derivative matrix, so it is computed as the gradient of a scalar rather than by building
  the matrix and contracting it. On a reverse-mode backend that is one pass instead of one per
  sample.
- **Warm Markov chains** — [`sample`](@ref) returns the sampler's own state alongside the draw,
  and an `MCState` keeps it across a parameter update, so burn-in is paid once per run rather
  than once per optimization step.
- **Any autodiff backend, none of them a dependency** — differentiation goes through
  [DifferentiationInterface.jl](https://github.com/JuliaDiff/DifferentiationInterface.jl).

## Installation

**Requirements**: Julia 1.11 or later.

```julia
julia> import Pkg; Pkg.add("NQSCore")
```

## Quick example

```@example nqscore
using NQSCore, ConnectedBasisConfigurations, SymBasis, OperatorAlgebra
using DifferentiationInterface, ForwardDiff, Random

spec, nsites = Spin(1 // 2), 4
b = basis(dof_object(spec), nsites)

ops = local_operators(spec)
H = OpSum(vcat(
    [Op(2 .* ops.sz, i) * Op(2 .* ops.sz, mod1(i + 1, nsites)) for i in 1:nsites],
    [2.0 * Op(2 .* ops.sx, i) for i in 1:nsites],
))

a = LogStateVector(spec, nsites, b)
vs = FullSumState(a, init_parameters(a, Xoshiro(0); scale=0.3); backend=AutoForwardDiff())

expect(vs, H)
```

The error bar is exactly zero because nothing was sampled. Descending on the gradient reaches
the exact ground state, and the vanishing variance is the sharper statement — the energy
variance is zero precisely on an eigenstate:

```@example nqscore
for _ in 1:2000
    _, ∇ = expect_and_grad(vs, H)
    setparameters!(vs, parameters(vs) .- 0.05 .* ∇)
end

expect(vs, H)
```

Swapping in Monte Carlo changes the state type and nothing else:

```@example nqscore
mc = MCState(a, init_parameters(a, Xoshiro(0); scale=0.3), ExactSampler(b, 100_000);
             backend=AutoForwardDiff(), rng=Xoshiro(1))
expect(mc, H)
```

## Scope

Ansatze, samplers and preconditioners are deliberately not here: this package defines what they
must satisfy and provides one reference implementation of each — [`LogStateVector`](@ref) and
[`ExactSampler`](@ref) — chosen because they are exact rather than because they are useful.
Neural networks live in `NQSAnsatze`, Markov chains in `NQSSamplers`, stochastic reconfiguration
in `NQSOptimisers`.

## Module

What `?NQSCore` shows at the REPL:

```@docs
NQSCore
```
