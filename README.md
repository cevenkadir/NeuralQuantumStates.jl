<p align="center">
    <img width="200px" src="docs/src/assets/logo.svg#gh-light-mode-only"/>
    <img width="200px" src="docs/src/assets/logo-dark.svg#gh-dark-mode-only"/>
</p>
<div align="center">

# NeuralQuantumStates.jl

*Neural quantum states in Julia*

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://cevenkadir.github.io/NeuralQuantumStates.jl/dev/) [![Build Status](https://github.com/cevenkadir/NeuralQuantumStates.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/cevenkadir/NeuralQuantumStates.jl/actions/workflows/CI.yml?query=branch%3Amain) [![Coverage](https://codecov.io/gh/cevenkadir/NeuralQuantumStates.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/cevenkadir/NeuralQuantumStates.jl) [![Code Style: Blue](https://img.shields.io/badge/code%20style-blue-4495d1.svg)](https://github.com/invenia/BlueStyle)
</div>

NeuralQuantumStates.jl is a Julia package under development to facilitate the training of neural quantum states (NQS) by variational Monte Carlo (VMC).

The package aims to provide an efficient and extensible environment for the simulation of closed many-body quantum systems by exploiting the power of neural networks and modern computational resources. Inspired by established Python libraries such as [NetKet](https://github.com/netket/netket) and [jVMC](https://github.com/markusschmitt/vmc_jax), NeuralQuantumStates.jl focuses on providing a machine learning toolbox for quantum many-body systems in a Julia-based environment.

## Project status
**This package is a <ins>work in progress</ins>.** The full variational Monte Carlo pipeline is
in place — lattices and space groups, symmetry-reduced bases, the local-energy kernel, ansätze,
samplers, and stochastic reconfiguration — and is verified end to end against exact
diagonalization. What is not yet done: performance work for CPU and GPU, distributed execution,
and time evolution (TDVP). The API is still unstable.

## Installation
If you still want to try it out, you can install it from the Julia REPL by entering:
```julia
import Pkg; Pkg.add(url="https://github.com/cevenkadir/NeuralQuantumStates.jl")
```

## Documentation
For information on using this package, check out the [in-development documentation](https://cevenkadir.github.io/NeuralQuantumStates.jl/dev/).

## The package ecosystem

This project is being reorganized from a single monolithic package into a set of focused
packages developed together in this repository under [`lib/`](lib/), following the layout
[Lux.jl](https://github.com/LuxDL/Lux.jl) uses. In Julia, depending on a package means paying
its compile latency, so someone who wants only lattice geometry should not have to load an
autodiff engine and a GPU backend.

**Tier 1 — no machine-learning dependencies.** Usable on their own for exact diagonalization.

| Package | Role | Status |
|---|---|---|
| [SymBasis.jl](https://github.com/cevenkadir/SymBasis.jl) | States, degrees of freedom, symmetry groups, symmetry-reduced bases | released |
| [OperatorAlgebra.jl](https://github.com/h-mnzlr/OperatorAlgebra.jl) | Operator algebra: `Op`/`OpChain`/`OpSum`, sparse and dense conversion, fermionic sites | released |
| `LatticeSpaceGroups` | Lattice geometry, neighbour graphs, and the site permutations symmetry groups need | working |
| `ConnectedConfigs` | Batched local-energy kernel: connected configurations and their matrix elements | working |

**Tier 2 — the neural-network stack.**

| Package | Role | Status |
|---|---|---|
| `NQSCore` | Interfaces, `MCState` and `FullSumState`, log-derivatives, Monte Carlo statistics | working |
| `NQSAnsatze` | Ansätze on [Lux.jl](https://github.com/LuxDL/Lux.jl): RBM, symmetric RBM, Jastrow | working |
| `NQSSamplers` | Metropolis sampler with local, exchange and Hamiltonian transition rules | working |
| `NQSOptimisers` | Stochastic reconfiguration, MinSR, natural gradient | working |

`NeuralQuantumStates.jl` itself is the umbrella: it re-exports the above and provides the
predefined models and the `VMC` driver with callbacks and logging. Installing it gives you the
whole stack.

## Quick example

```julia
using NeuralQuantumStates, DifferentiationInterface, ForwardDiff, Random

lat   = build(Hypercube([6], 1.0; periodic=[true]))
model = build(TransverseFieldIsing(lat; J=1.0, h_x=1.0))

ansatz = LuxAnsatz(RBM(6, 2), model.dof, 6; rng=Xoshiro(0))
state  = FullSumState(ansatz, init_parameters(ansatz, Xoshiro(0)), AutoForwardDiff())

log = run!(VMC(state, model.hamiltonian;
        preconditioner=StochasticReconfiguration(; diag_shift=1e-2),
        optimizer=Descent(0.05));
    iterations=1000, callbacks=(InvalidLossStopping(),))

final_energy(log)
```

The pre-split `Lattices`, `Hilberts`, and `Operators` modules have been retired. `Lattices`
became `LatticeSpaceGroups`, `Hilberts` is replaced by SymBasis.jl, and `Operators` by
OperatorAlgebra.jl together with `ConnectedConfigs` — which was proved to reproduce the old
matrix elements exactly against the reference data in [`test/golden/`](test/golden/) before the
old code was removed. The retired sources are kept under
[`archive/pre-split/`](archive/pre-split/) as the provenance of that data.

### Migrating

| Was | Now |
|---|---|
| `Lattices.build(:Hypercube, [8], 1.0; periodic=[true])` | `build(Hypercube([8], 1.0; periodic=[true]))` |
| `Hilberts.build(:Spin, 1//2, N)` | `basis(dof_object(Spin(1//2)), N)` |
| `Operators.build(:TransverseFieldIsing, h, l; ...)` | `build(TransverseFieldIsing(l; ...))` |
| `Operators.connected_basis_configs(H, samples)` | `connected_padded(H, states)` |

Two behavioural changes are deliberate: connected configurations are padded with a **zero**
matrix element rather than `missing`, and the degree-of-freedom axis is always **first** (the
old code put it last for Ising and first for Bose-Hubbard).

### Further goals
- [ ] Support for distributed and parallel computing via [MPI.jl](https://github.com/JuliaParallel/MPI.jl/tree/master).
- [ ] GPU support via [CUDA.jl](https://github.com/JuliaGPU/CUDA.jl), [AMDGPU.jl](https://github.com/JuliaGPU/AMDGPU.jl), and [Metal.jl](https://github.com/JuliaGPU/Metal.jl).

## Bugs report and feature requests
If you think you have found a bug or have a feature request, you can open an [issue](https://github.com/cevenkadir/NeuralQuantumStates.jl/issues/new).

## Citation
If you use this package in your work, 
we would appreciate the following reference as in [CITATION.bib](https://github.com/cevenkadir/NeuralQuantumStates.jl/blob/main/CITATION.bib).
