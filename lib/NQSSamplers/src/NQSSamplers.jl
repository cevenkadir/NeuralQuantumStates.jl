"""
    NQSSamplers

Monte Carlo samplers over discrete quantum configurations.

This corresponds to NetKet's `netket/sampler`, including its `rules/` submodule, and follows the
same separation: the [`MetropolisSampler`](@ref) owns the Markov chain machinery, while an
[`AbstractRule`](@ref) owns the physics of which moves to propose.

# What is here
- [`MetropolisSampler`](@ref) — multi-chain Metropolis–Hastings with burn-in and thinning.
- Rules: [`LocalRule`](@ref) for unconstrained systems, [`ExchangeRule`](@ref) for conserved
  sectors, and [`HamiltonianRule`](@ref) for proposals drawn from the Hamiltonian's own
  connectivity.

Exact summation lives elsewhere on purpose: following NetKet, `NQSCore.FullSumState` is a
variational *state* rather than a sampler, because it involves no sampling at all.
`NQSCore.ExactSampler` likewise stays there as the reference implementation that these samplers
are validated against.

# Choosing a rule

The rule is the part that has to match the problem:

| Situation | Rule |
|---|---|
| No conserved quantity | [`LocalRule`](@ref) |
| Fixed magnetization or particle number | [`ExchangeRule`](@ref) |
| Sparse Hamiltonian, hard-to-satisfy constraints | [`HamiltonianRule`](@ref) |

Using `LocalRule` inside a conserved sector is the classic mistake: every proposal leaves the
sector, everything is rejected, and the chain returns the same configuration forever with
error bars that look perfectly reasonable.

# Example
```julia
using NQSCore, NQSSamplers

starts = random_configurations(dof, nsites, 8, rng)
sampler = MetropolisSampler(LocalRule(), starts; n_chains=8, n_samples=2000, burn_in=500)
vs = MCState(ansatz, θ, sampler; backend=AutoForwardDiff(), rng=rng)
expect(vs, H)
```
"""
module NQSSamplers

using Random: AbstractRNG, rand

using NQSCore
using NQSCore: AbstractAnsatz, AbstractSampler, log_amplitude

using ConnectedConfigs
using ConnectedConfigs: connected, local_dimension, local_values
using SymBasis.DigitBase: read, write

include("rules.jl")
include("metropolis.jl")

export AbstractRule, LocalRule, ExchangeRule, HamiltonianRule, propose
export MetropolisSampler, ACCEPTANCE
export random_configuration, random_configurations

end # module NQSSamplers
