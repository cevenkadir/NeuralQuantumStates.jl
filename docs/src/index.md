```@meta
CurrentModule = NeuralQuantumStates
```

# NeuralQuantumStates.jl

*Neural quantum states in Julia.*

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://cevenkadir.github.io/NeuralQuantumStates.jl/dev/) [![Build Status](https://github.com/cevenkadir/NeuralQuantumStates.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/cevenkadir/NeuralQuantumStates.jl/actions/workflows/CI.yml?query=branch%3Amain) [![Coverage](https://codecov.io/gh/cevenkadir/NeuralQuantumStates.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/cevenkadir/NeuralQuantumStates.jl)

NeuralQuantumStates.jl facilitates the training of neural quantum states (NQS) by variational
Monte Carlo (VMC). It aims to provide an efficient and extensible environment for simulating
closed many-body quantum systems, taking inspiration from [NetKet](https://github.com/netket/netket)
and [jVMC](https://github.com/markusschmitt/vmc_jax).

!!! warning "Work in progress"
    This package is under active development and is being reorganized (see below). Most
    functionality is still being implemented, performance is not yet optimized for CPU or GPU,
    and the API is unstable.

## The package ecosystem

NeuralQuantumStates.jl is being split from a single monolithic package into a set of focused
packages, so that each piece can be used on its own without paying for the whole stack. In
Julia this matters more than it does in Python: depending on a package means paying its compile
latency, so a user who only wants lattice geometry should not have to load an autodiff engine
and a GPU backend.

The packages are organized in two tiers.

**Tier 1 — no machine-learning dependencies.** Usable on their own for plain exact
diagonalization.

| Package | Role |
|---|---|
| [SymBasis.jl](https://github.com/cevenkadir/SymBasis.jl) | States, degrees of freedom, symmetry groups, symmetry-reduced bases |
| [OperatorAlgebra.jl](https://github.com/h-mnzlr/OperatorAlgebra.jl) | Operator algebra: `Op`, `OpChain`, `OpSum`, sparse/dense conversion, fermionic sites |
| `LatticeSpaceGroups` | Lattice geometry, bonds, and the site permutations symmetry groups need |
| `ConnectedConfigs` | The batched local-energy kernel: connected configurations and their matrix elements |

**Tier 2 — the neural-network stack.**

| Package | Role |
|---|---|
| `NQSCore` | Interfaces, the `MCState` and `FullSumState` variational states, log-derivatives, statistics |
| `NQSAnsatze` | Ansätze built on Lux: RBM, symmetric RBM, Jastrow |
| `NQSSamplers` | Metropolis sampling with local, exchange and Hamiltonian transition rules |
| `NQSOptimisers` | Stochastic reconfiguration and its kernel-trick (MinSR) form |

`NeuralQuantumStates.jl` itself is the umbrella: it re-exports all of the above and adds the
predefined models and the [`VMC`](@ref) driver with callbacks and logging. Installing it gives
you the whole stack.

## Installation

```julia
import Pkg; Pkg.add(url="https://github.com/cevenkadir/NeuralQuantumStates.jl")
```

## Getting started

See [Basics](@ref) for a worked example that builds a lattice, a Hilbert space, and a
Hamiltonian, and then inspects its connected basis configurations.

`LatticeSpaceGroups` and `ConnectedConfigs` are registered as packages in their own right and
carry their own documentation sites. The Manual has a chapter for each —
[Lattices and symmetries](manual/latticespacegroups.md) and
[Connected configurations](manual/connectedconfigs.md) — covering what you need to use them
from here, and linking on to the full sites for the rest.

## Bug reports and feature requests

Please open an [issue](https://github.com/cevenkadir/NeuralQuantumStates.jl/issues/new).

## Citation

If you use this package in your work, we would appreciate the reference in
[CITATION.bib](https://github.com/cevenkadir/NeuralQuantumStates.jl/blob/main/CITATION.bib).
