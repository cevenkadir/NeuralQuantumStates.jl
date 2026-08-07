<p align="center">
    <img width="200px" src="docs/src/assets/logo.svg#gh-light-mode-only"/>
    <img width="200px" src="docs/src/assets/logo-dark.svg#gh-dark-mode-only"/>
</p>
<div align="center">

# NeuralQuantumStates.jl

*Neural quantum states in Julia*

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://cevenkadir.github.io/NeuralQuantumStates.jl/dev/) [![Build Status](https://github.com/cevenkadir/NeuralQuantumStates.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/cevenkadir/NeuralQuantumStates.jl/actions/workflows/CI.yml?query=branch%3Amain) [![Coverage](https://codecov.io/gh/cevenkadir/NeuralQuantumStates.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/cevenkadir/NeuralQuantumStates.jl)
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

Two external packages carry the foundations and are not developed here.
[SymBasis.jl](https://github.com/cevenkadir/SymBasis.jl) provides states, degrees of freedom,
symmetry groups and symmetry-reduced bases;
[OperatorAlgebra.jl](https://github.com/h-mnzlr/OperatorAlgebra.jl) provides `Op`/`OpChain`/`OpSum`,
sparse and dense conversion, and fermionic sites. Everything below lives in this repository.

**Tier 1 — no machine-learning dependencies.** Usable on their own for exact diagonalization.

| Package | Role | Stable Version | Monthly Downloads | Total Downloads | Build Status |
| :------ | :--- | :------------- | :---------------- | :-------------- | :----------- |
| [LatticeSpaceGroups.jl](lib/LatticeSpaceGroups) | Lattice geometry, bonds, and the site permutations symmetry groups need | [![][lsg-version]][lsg-juliahub] | [![][lsg-monthly]][lsg-downloads] | [![][lsg-total]][lsg-downloads] | [![][lsg-ci]][lsg-ci-url] |
| [ConnectedConfigs.jl](lib/ConnectedConfigs) | Batched local-energy kernel: connected configurations and their matrix elements | [![][cc-version]][cc-juliahub] | [![][cc-monthly]][cc-downloads] | [![][cc-total]][cc-downloads] | [![][cc-ci]][cc-ci-url] |

**Tier 2 — the neural-network stack.**

| Package | Role | Stable Version | Monthly Downloads | Total Downloads | Build Status |
| :------ | :--- | :------------- | :---------------- | :-------------- | :----------- |
| [NQSCore.jl](lib/NQSCore) | Interfaces, `MCState` and `FullSumState`, log-derivatives, statistics | — | [![][core-monthly]][core-downloads] | [![][core-total]][core-downloads] | [![][core-ci]][core-ci-url] |
| [NQSAnsatze.jl](lib/NQSAnsatze) | Ansätze on [Lux.jl](https://github.com/LuxDL/Lux.jl): RBM, symmetric RBM, Jastrow | — | [![][ansatze-monthly]][ansatze-downloads] | [![][ansatze-total]][ansatze-downloads] | [![][ansatze-ci]][ansatze-ci-url] |
| [NQSSamplers.jl](lib/NQSSamplers) | Metropolis sampling with local, exchange and Hamiltonian rules | — | [![][samplers-monthly]][samplers-downloads] | [![][samplers-total]][samplers-downloads] | [![][samplers-ci]][samplers-ci-url] |
| [NQSOptimisers.jl](lib/NQSOptimisers) | Stochastic reconfiguration, MinSR, natural gradient | — | [![][optimisers-monthly]][optimisers-downloads] | [![][optimisers-total]][optimisers-downloads] | [![][optimisers-ci]][optimisers-ci-url] |

An em dash in **Stable Version** means the package has not been registered in General yet, so
it has no released version to point at.

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
old code was removed. The retired sources remain in git history — recover them with
`git show 0774c25:archive/pre-split/` — as does the pre-split scratch work under
`archive/scratch/`.

### Migrating

| Was                                                     | Now                                         |
| ------------------------------------------------------- | ------------------------------------------- |
| `Lattices.build(:Hypercube, [8], 1.0; periodic=[true])` | `build(Hypercube([8], 1.0; periodic=true))` |
| `nv(lattice)`                                           | `n_sites(lattice)`                          |
| `Hilberts.build(:Spin, 1//2, N)`                        | `basis(dof_object(Spin(1//2)), N)`          |
| `Operators.build(:TransverseFieldIsing, h, l; ...)`     | `build(TransverseFieldIsing(l; ...))`       |
| `Operators.connected_basis_configs(H, samples)`         | `connected_padded(H, states)`               |

Three behavioural changes are deliberate: connected configurations are padded with a **zero**
matrix element rather than `missing`; the degree-of-freedom axis is always **first** (the old
code put it last for Ising and first for Bose-Hubbard); and two terms reaching the same
configuration now produce two rows rather than one summed row, which every consumer sums over
anyway.

`connected_padded` also accepts a batch of any shape, and `compile(H)` flattens a Hamiltonian
once so that the per-step cost stops including it — worth about two orders of magnitude on the
kernel. See the [ConnectedConfigs README](lib/ConnectedConfigs/README.md).

### Further goals
- [ ] Support for distributed and parallel computing via [MPI.jl](https://github.com/JuliaParallel/MPI.jl/tree/master).
- [ ] GPU support via [CUDA.jl](https://github.com/JuliaGPU/CUDA.jl), [AMDGPU.jl](https://github.com/JuliaGPU/AMDGPU.jl), and [Metal.jl](https://github.com/JuliaGPU/Metal.jl).

## License

MIT — see [LICENSE](LICENSE). Each package under [`lib/`](lib/) carries the same licence, so a
subdirectory can be registered on its own.

## Bugs report and feature requests
If you think you have found a bug or have a feature request, you can open an [issue](https://github.com/cevenkadir/NeuralQuantumStates.jl/issues/new).

## Citation
If you use this package in your work, 
we would appreciate the following reference as in [CITATION.bib](https://github.com/cevenkadir/NeuralQuantumStates.jl/blob/main/CITATION.bib).

<!-- Badge definitions for the package ecosystem tables. Kept out of the tables themselves:
     inlined, a single row would run past a thousand characters. -->

[lsg-version]: https://juliahub.com/docs/General/LatticeSpaceGroups/stable/version.svg
[lsg-juliahub]: https://juliahub.com/ui/Packages/General/LatticeSpaceGroups
[lsg-monthly]: https://img.shields.io/badge/dynamic/json?url=http%3A%2F%2Fjuliapkgstats.com%2Fapi%2Fv1%2Fmonthly_downloads%2FLatticeSpaceGroups&query=total_requests&suffix=%2Fmonth&label=Downloads
[lsg-total]: https://img.shields.io/badge/dynamic/json?url=http%3A%2F%2Fjuliapkgstats.com%2Fapi%2Fv1%2Ftotal_downloads%2FLatticeSpaceGroups&query=total_requests&label=Total%20Downloads
[lsg-downloads]: https://juliapkgstats.com/pkg/LatticeSpaceGroups
[lsg-ci]: https://github.com/cevenkadir/NeuralQuantumStates.jl/actions/workflows/CI-LatticeSpaceGroups.yml/badge.svg?branch=main
[lsg-ci-url]: https://github.com/cevenkadir/NeuralQuantumStates.jl/actions/workflows/CI-LatticeSpaceGroups.yml?query=branch%3Amain

[cc-version]: https://juliahub.com/docs/General/ConnectedConfigs/stable/version.svg
[cc-juliahub]: https://juliahub.com/ui/Packages/General/ConnectedConfigs
[cc-monthly]: https://img.shields.io/badge/dynamic/json?url=http%3A%2F%2Fjuliapkgstats.com%2Fapi%2Fv1%2Fmonthly_downloads%2FConnectedConfigs&query=total_requests&suffix=%2Fmonth&label=Downloads
[cc-total]: https://img.shields.io/badge/dynamic/json?url=http%3A%2F%2Fjuliapkgstats.com%2Fapi%2Fv1%2Ftotal_downloads%2FConnectedConfigs&query=total_requests&label=Total%20Downloads
[cc-downloads]: https://juliapkgstats.com/pkg/ConnectedConfigs
[cc-ci]: https://github.com/cevenkadir/NeuralQuantumStates.jl/actions/workflows/CI-ConnectedConfigs.yml/badge.svg?branch=main
[cc-ci-url]: https://github.com/cevenkadir/NeuralQuantumStates.jl/actions/workflows/CI-ConnectedConfigs.yml?query=branch%3Amain

[core-monthly]: https://img.shields.io/badge/dynamic/json?url=http%3A%2F%2Fjuliapkgstats.com%2Fapi%2Fv1%2Fmonthly_downloads%2FNQSCore&query=total_requests&suffix=%2Fmonth&label=Downloads
[core-total]: https://img.shields.io/badge/dynamic/json?url=http%3A%2F%2Fjuliapkgstats.com%2Fapi%2Fv1%2Ftotal_downloads%2FNQSCore&query=total_requests&label=Total%20Downloads
[core-downloads]: https://juliapkgstats.com/pkg/NQSCore
[core-ci]: https://github.com/cevenkadir/NeuralQuantumStates.jl/actions/workflows/CI-NQSCore.yml/badge.svg?branch=main
[core-ci-url]: https://github.com/cevenkadir/NeuralQuantumStates.jl/actions/workflows/CI-NQSCore.yml?query=branch%3Amain

[ansatze-monthly]: https://img.shields.io/badge/dynamic/json?url=http%3A%2F%2Fjuliapkgstats.com%2Fapi%2Fv1%2Fmonthly_downloads%2FNQSAnsatze&query=total_requests&suffix=%2Fmonth&label=Downloads
[ansatze-total]: https://img.shields.io/badge/dynamic/json?url=http%3A%2F%2Fjuliapkgstats.com%2Fapi%2Fv1%2Ftotal_downloads%2FNQSAnsatze&query=total_requests&label=Total%20Downloads
[ansatze-downloads]: https://juliapkgstats.com/pkg/NQSAnsatze
[ansatze-ci]: https://github.com/cevenkadir/NeuralQuantumStates.jl/actions/workflows/CI-NQSAnsatze.yml/badge.svg?branch=main
[ansatze-ci-url]: https://github.com/cevenkadir/NeuralQuantumStates.jl/actions/workflows/CI-NQSAnsatze.yml?query=branch%3Amain

[samplers-monthly]: https://img.shields.io/badge/dynamic/json?url=http%3A%2F%2Fjuliapkgstats.com%2Fapi%2Fv1%2Fmonthly_downloads%2FNQSSamplers&query=total_requests&suffix=%2Fmonth&label=Downloads
[samplers-total]: https://img.shields.io/badge/dynamic/json?url=http%3A%2F%2Fjuliapkgstats.com%2Fapi%2Fv1%2Ftotal_downloads%2FNQSSamplers&query=total_requests&label=Total%20Downloads
[samplers-downloads]: https://juliapkgstats.com/pkg/NQSSamplers
[samplers-ci]: https://github.com/cevenkadir/NeuralQuantumStates.jl/actions/workflows/CI-NQSSamplers.yml/badge.svg?branch=main
[samplers-ci-url]: https://github.com/cevenkadir/NeuralQuantumStates.jl/actions/workflows/CI-NQSSamplers.yml?query=branch%3Amain

[optimisers-monthly]: https://img.shields.io/badge/dynamic/json?url=http%3A%2F%2Fjuliapkgstats.com%2Fapi%2Fv1%2Fmonthly_downloads%2FNQSOptimisers&query=total_requests&suffix=%2Fmonth&label=Downloads
[optimisers-total]: https://img.shields.io/badge/dynamic/json?url=http%3A%2F%2Fjuliapkgstats.com%2Fapi%2Fv1%2Ftotal_downloads%2FNQSOptimisers&query=total_requests&label=Total%20Downloads
[optimisers-downloads]: https://juliapkgstats.com/pkg/NQSOptimisers
[optimisers-ci]: https://github.com/cevenkadir/NeuralQuantumStates.jl/actions/workflows/CI-NQSOptimisers.yml/badge.svg?branch=main
[optimisers-ci-url]: https://github.com/cevenkadir/NeuralQuantumStates.jl/actions/workflows/CI-NQSOptimisers.yml?query=branch%3Amain
