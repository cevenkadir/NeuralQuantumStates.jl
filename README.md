# NeuralQuantumStates.jl

*Neural quantum states in Julia*

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://cevenkadir.github.io/NeuralQuantumStates.jl/dev/)
[![Build Status](https://github.com/cevenkadir/NeuralQuantumStates.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/cevenkadir/NeuralQuantumStates.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/cevenkadir/NeuralQuantumStates.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/cevenkadir/NeuralQuantumStates.jl)
[![Code Style: Blue](https://img.shields.io/badge/code%20style-blue-4495d1.svg)](https://github.com/invenia/BlueStyle)

NeuralQuantumStates.jl aims to provide a Julia package for training neural quantum states (NQS) using the variational Monte Carlo (VMC) method.

## Project status
**This package is a <ins>work in progress</ins>.** Most of the functionality still needs to be implemented. The performance still needs to be optimized for both CPU and GPU. Its API is not yet stable.

## Installation
If you still want to try it out, you can install it from the Julia REPL by pressing `]` to enter `pkg` mode and then entering:
```julia
pkg> add https://github.com/cevenkadir/NeuralQuantumStates.jl
```

## Documentation
For information on using this package, check out the [stable documentation](https://kadirceven.com/NeuralQuantumStates.jl/stable/). Since this package still lacks many features, we strongly recommend you also check out the [in-development documentation](https://kadirceven.com/NeuralQuantumStates.jl/dev/) for recent updates in the package.

## Development goals
- [ ] `Lattices` submodule to generate any Bravais lattice.
- [ ] `Networks` submodule to generate canonical artificial neural networks (ANN) via [Flux.jl](https://github.com/FluxML/Flux.jl).
- [ ] `VarStates` submodule to define variational quantum states.
- [ ] `Hilberts` submodule to define Hilbert spaces.
- [ ] `Operators` submodule to define arbitrary quantum operators on a computational basis.
- [ ] `Samplers` submodule to sample variational quantum states with Markov chain Monte-Carlo (MCMC) methods.
- [ ] `OptiHandlers` submodule to optimize variational quantum states with gradient-based methods.
- [ ] Support for distributed and parallel computing via [MPI.jl](https://github.com/JuliaParallel/MPI.jl/tree/master).
- [ ] GPU support via [CUDA.jl](https://github.com/JuliaGPU/CUDA.jl), [AMDGPU.jl](https://github.com/JuliaGPU/AMDGPU.jl), and [Metal.jl](https://github.com/JuliaGPU/Metal.jl).

## Bugs report and feature requests
If you think you have found a bug or have a feature request, you can open an [issue](https://github.com/cevenkadir/NeuralQuantumStates.jl/issues/new).

## Citation
If you use this package in your work, 
we would appreciate the following reference as in [CITATION.bib](https://github.com/cevenkadir/NeuralQuantumStates.jl/main/CITATION.bib).

## Acknowledgements
This package is mainly inspired by the Python libraries of [NetKet](https://github.com/netket/netket) and [jVMC](https://github.com/markusschmitt/vmc_jax).


