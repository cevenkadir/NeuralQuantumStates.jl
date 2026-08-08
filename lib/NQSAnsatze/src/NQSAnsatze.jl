"""
    NQSAnsatze

Neural-network ansätze for quantum many-body wavefunctions, built on Lux.jl.

This corresponds to NetKet's `netket/nn` and `netket/models`: [`LuxAnsatz`](@ref) is the bridge
that makes any Lux model usable as a wavefunction, and the layers here are the wavefunction
architectures themselves.

This is the one package in the stack allowed to depend on Lux. Everything else reaches neural
networks through the `NQSCore` interface, which is why a sampler or an optimizer never has to
know that a network is involved at all.

# Ansätze
- [`RBM`](@ref) — the canonical neural quantum state, with its hidden units summed analytically.
- [`SymmetricRBM`](@ref) — weights shared across a symmetry group, so `|ψ|` is invariant by
  construction rather than by training.
- [`Jastrow`](@ref) — a pair-correlation factor.

# Complex wavefunctions

Parameters default to `ComplexF64`. A wavefunction has a phase, and a real-parameter network can
only represent a positive one — which silently fails for any model with a sign structure, such
as a frustrated magnet. A real network can still be used by having it return two rows, read as
the real and imaginary parts of `log ψ`; see [`LuxAnsatz`](@ref).

# Example
```julia
using NQSAnsatze, NQSCore

model = RBM(nsites, 2)
a = LuxAnsatz(model, dof, nsites; rng)
vs = FullSumState(a, init_parameters(a, rng); backend=AutoZygote())
```

Loading `LatticeSpaceGroups` alongside activates an extension that builds symmetric ansätze
straight from a lattice, deriving the permutations from its geometry.
"""
module NQSAnsatze

using LinearAlgebra
using Random: AbstractRNG, default_rng, randn

using ChainRulesCore
using Lux
using SymBasis

using NQSCore
using NQSCore: AbstractAnsatz

include("layers.jl")
include("lux_ansatz.jl")

export LuxAnsatz
export RBM, Jastrow, SymmetricRBM, logtwocosh

end # module NQSAnsatze
