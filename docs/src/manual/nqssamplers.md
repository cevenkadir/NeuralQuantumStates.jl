```@meta
CurrentModule = NQSSamplers
```

# NQSSamplers.jl

*Monte Carlo sampling over discrete quantum configurations.*

[`MetropolisSampler`](@ref) owns the Markov chain machinery; an [`AbstractRule`](@ref) owns the
physics of which moves to propose. That separation follows NetKet, and it matters because the
rule is where all the problem knowledge lives.

## Choosing a rule

This is the decision that determines whether the chain works at all.

| Situation | Rule |
|---|---|
| No conserved quantity | [`LocalRule`](@ref) |
| Fixed magnetization or particle number | [`ExchangeRule`](@ref) |
| Sparse Hamiltonian, hard-to-satisfy constraints | [`HamiltonianRule`](@ref) |

!!! warning "A local rule cannot move inside a conserved sector"
    Every single-site change alters the total, so every proposal leaves the sector and is
    rejected. The chain then returns one configuration forever — while reporting perfectly
    reasonable-looking error bars, because a constant sequence has no variance. Use
    [`ExchangeRule`](@ref), which swaps two sites and therefore conserves any quantity that is
    a sum over sites.

[`HamiltonianRule`](@ref) proposes only configurations the Hamiltonian actually connects to, so
it respects whatever the Hamiltonian conserves without being told what that is. It is
**asymmetric** — a configuration and its image generally have different numbers of connections
— so it carries a proposal correction. Dropping that correction biases the sampled distribution
while leaving everything else looking fine.

## Multiple chains

Not for speed. Several chains started from dispersed configurations give a convergence
diagnostic no single chain can: split-R̂ compares them, and a chain stuck in one region is
invisible from the inside. The error bar is also estimated from the spread *between* chain
means, which assumes only that the chains are independent rather than relying on an
autocorrelation model.

All chains advance together, so each step evaluates the ansatz once on a batch rather than once
per chain — the ansatz being the expensive part of sampling.

## Diagnostics

Check two numbers before believing a result:

- `r_hat` from the returned `Stats`. Anything much above `1.01` means the chains disagree and
  the error bar cannot be trusted.
- [`ACCEPTANCE`](@ref) after a run. Near zero means a barely-moving chain; near one usually
  means proposals too timid to explore. Both produce confident-looking error bars from samples
  that carry no information.

## What is deliberately elsewhere

Exact summation is a variational *state* (`NQSCore.FullSumState`), not a sampler, because it
involves no sampling. `NQSCore.ExactSampler` likewise stays there as the reference
implementation these samplers are validated against — a Metropolis chain must reproduce the
exact `|ψ|²`, and the test suite checks the total-variation distance between them.

## Quick example

```@example nqssamplers
using NQSSamplers, NQSCore, ConnectedConfigs, SymBasis, OperatorAlgebra
using DifferentiationInterface, ForwardDiff, Random

spec, nsites = Spin(1 // 2), 6
b = basis(dof_object(spec), nsites)
a = LogStateVector(spec, nsites, b)
θ = init_parameters(a, Xoshiro(0); scale=0.4)

ops = local_operators(spec)
H = OpSum(vcat(
    [Op(2 .* ops.sz, i) * Op(2 .* ops.sz, mod1(i + 1, nsites)) for i in 1:nsites],
    [Op(2 .* ops.sx, i) for i in 1:nsites],
))

starts = random_configurations(spec, nsites, 8, Xoshiro(1))
sampler = MetropolisSampler(LocalRule(), starts;
    n_chains=8, n_samples=5_000, burn_in=1_000)

vs = MCState(a, θ, sampler; backend=AutoForwardDiff(), rng=Xoshiro(2))
expect(vs, H)
```

Compare against the exact answer for the same parameters:

```@example nqssamplers
expect(FullSumState(a, θ, AutoForwardDiff()), H)
```
