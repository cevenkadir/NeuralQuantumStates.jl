```@meta
CurrentModule = NQSCore
```

# NQSCore.jl

*The interface backbone of the neural-quantum-states stack.*

NQSCore defines the abstract types and function names that `NQSAnsatze`, `NQSSamplers`, and
`NQSOptimisers` are all written against — the role `SciMLBase` plays for SciML, or
`ChainRulesCore` for the autodiff ecosystem. Because they share this vocabulary, none of them
has to depend on any other.

It also provides the two variational state types, the statistics, and the log-derivative
machinery they share.

## Two states, one interface

| Type | Averages by | Error bar |
|---|---|---|
| [`FullSumState`](@ref) | exact summation over the basis | none — the answer is exact |
| [`MCState`](@ref) | Monte Carlo sampling | standard error, autocorrelation-corrected |

Both satisfy the same [`expect`](@ref) interface, which is what makes the development loop
work: build a model against exact summation, where any disagreement is a bug rather than a
fluctuation, then scale up by swapping the state type and changing nothing else.

## Reference implementations

[`LogStateVector`](@ref) and [`ExactSampler`](@ref) contain no approximation at all. The first
carries one parameter per basis state, so it can represent *any* state exactly; the second
enumerates the basis and draws genuinely independent samples. Both are useless for anything
large and indispensable for testing — they let the machinery be checked against exact answers
rather than against itself.

`LogStateVector` also has a closed-form Jacobian (`log ψ(s) = θ_s`, so `O` is a one-hot
indicator matrix), which is how [`log_derivatives`](@ref) itself is verified.

## Statistics

Every [`expect`](@ref) returns a [`Stats`](@ref): mean, error of the mean, variance, integrated
autocorrelation time, and split-R̂. Downstream code never branches on where a number came from.

The error bar is **corrected for autocorrelation**. Markov chain samples are correlated, so the
naive `std/sqrt(n)` understates it by roughly `sqrt(tau_corr)` — reporting that unadjusted is
the most common way to converge confidently to a wrong answer. The energy variance is worth
watching in its own right: it vanishes at an exact eigenstate, which makes it the sharpest
convergence diagnostic available.

## Log-derivatives

[`log_derivatives`](@ref) computes `O[s, k] = ∂ log ψ(x_s) / ∂θ_k`, the object behind both the
energy gradient (a covariance between `O` and the local energies) and the quantum geometric
tensor that stochastic reconfiguration inverts.

Two independent things can be complex, and conflating them is the usual source of wrong
gradients:

- **A complex log-amplitude with real parameters** — the common case for a network emitting a
  modulus and a phase. `O` is complex even though `θ` is real; both parts come from one
  autodiff pass over `[real(log ψ); imag(log ψ)]`.
- **Complex parameters**, which are ambiguous until a convention is fixed. `holomorphic=true`
  assumes `ψ` is holomorphic in `θ`; the default treats real and imaginary parts as `2n`
  independent real parameters, which is always valid. The holomorphic form is a shortcut that
  is silently wrong when the ansatz does not satisfy it, so it must be asked for explicitly.

Whatever the parameterization, the gradient handed back always matches the *structure of the
parameters*, so an optimizer never has to ask how the ansatz was parameterized.

## Backends

Differentiation goes through
[DifferentiationInterface.jl](https://github.com/JuliaDiff/DifferentiationInterface.jl), so any
of its backends works — `AutoForwardDiff()`, `AutoZygote()`, `AutoEnzyme()`. None of those
packages is a dependency here; load the one you want and pass it.

## Quick example

```@example nqscore
using NQSCore, ConnectedBasisConfigurations, SymBasis, OperatorAlgebra
using DifferentiationInterface, ForwardDiff, Random

spec, nsites = Spin(1 // 2), 4
b = basis(dof_object(spec), nsites)

ops = local_operators(spec)
H = OpSum([Op(2 .* ops.sz, i) * Op(2 .* ops.sz, mod1(i + 1, nsites)) for i in 1:nsites])

a = LogStateVector(spec, nsites, b)
vs = FullSumState(a, init_parameters(a, Xoshiro(0); scale=0.3), AutoForwardDiff())

expect(vs, H)
```

```@example nqscore
stats, gradient = expect_and_grad(vs, H)
stats
```
