```@meta
CurrentModule = NQSAnsatze
```

# NQSAnsatze.jl

*Neural-network ansätze for quantum many-body wavefunctions.*

[`LuxAnsatz`](@ref) makes any [Lux.jl](https://github.com/LuxDL/Lux.jl) model usable as a
variational wavefunction, and the layers here are wavefunction architectures in their own right.

This is the one package in the stack that depends on Lux. Everything else reaches neural
networks through the `NQSCore` interface, which is why a sampler or an optimizer never has to
know a network is involved.

## The ansätze

| Ansatz | Captures |
|---|---|
| [`RBM`](@ref) | The canonical neural quantum state; hidden units summed analytically |
| [`SymmetricRBM`](@ref) | An RBM whose weights are shared across a symmetry group |
| [`Jastrow`](@ref) | Pair correlations exactly, and nothing beyond them |

## Complex parameters by default

Parameters default to `ComplexF64`. A wavefunction has a phase, and a real-parameter network can
only represent a positive one — which fails silently for any model with a sign structure, such
as a frustrated magnet. A real network can still be used by having it emit two rows, read as the
real and imaginary parts of `log ψ`; see [`LuxAnsatz`](@ref).

## Symmetry by construction

[`SymmetricRBM`](@ref) applies each hidden filter to every permuted copy of the input, so `|ψ|`
is invariant under the group *by construction* rather than by training. The parameter count
drops by roughly the order of the group, and — more importantly — the ansatz cannot spend
capacity representing states the ground state is known not to occupy.

Loading `LatticeSpaceGroups` activates an extension that derives the permutations from lattice
geometry, so a translation-invariant ansatz on a kagome torus is no harder to write than one on
a chain:

```julia
using NQSAnsatze, LatticeSpaceGroups

lat = build(Kagome([3, 3], 1.0; periodic=true))
model = SymmetricRBM(lat, 2)      # invariant under every lattice translation
```

## Numerical care

`log(2 cosh z)` is evaluated as [`logtwocosh`](@ref), which factors out the dominant exponential
rather than computing the expression literally. A literal `log(2cosh(z))` overflows to `Inf`
once `|Re z|` passes about 710 — which an RBM reaches routinely as its weights grow — and one
`Inf` poisons the whole batch rather than a single sample.

## Quick example

```@example nqsansatze
using NQSAnsatze, NQSCore, ConnectedBasisConfigurations, SymBasis, OperatorAlgebra
using DifferentiationInterface, ForwardDiff, Random

spec, nsites = Spin(1 // 2), 10
b = basis(dof_object(spec), nsites)

a = LuxAnsatz(RBM(nsites, 2), spec, nsites; rng=Xoshiro(0))
NQSCore.n_parameters(a), length(b.states)
```

An RBM's parameter count grows *quadratically* in the system size while the Hilbert space grows
exponentially, so the two cross over quickly — at six sites an `alpha = 2` RBM actually has more
parameters than there are basis states, and only past roughly eight does it start compressing.
That crossover is why an RBM is worth nothing on a toy system and everything on a real one.

```@example nqsansatze
ops = local_operators(spec)
H = OpSum(vcat(
    [Op(2 .* ops.sz, i) * Op(2 .* ops.sz, mod1(i + 1, nsites)) for i in 1:nsites],
    [Op(2 .* ops.sx, i) for i in 1:nsites],
))

vs = FullSumState(a, init_parameters(a, Xoshiro(0)), AutoForwardDiff())
expect(vs, H)
```
