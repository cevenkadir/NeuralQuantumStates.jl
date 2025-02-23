```@raw html
---
# https://vitepress.dev/reference/default-theme-home-page
layout: home

hero:
  name: "NeuralQuantumStates.jl"
  tagline: Neural quantum states in Julia
  actions:
    - theme: brand
      text: Getting Started
      link: /basics
    - theme: alt
      text: Public API
      link: /lib/public
    - theme: alt
      text: View on GitHub
      link: https://github.com/cevenkadir/NeuralQuantumStates.jl
  image:
    src: /logo.svg
    alt: NeuralQuantumStates.jl

features:
  - icon: 🤔
    title: What is NeuralQuantumStates.jl?
    details: NeuralQuantumStates.jl is a Julia package under development to facilitate the training of neural quantum states (NQS) by variational Monte Carlo (VMC). The package aims to provide an efficient and extensible environment for the simulation of closed many-body quantum systems by exploiting the power of neural networks and modern computational resources.
  
  - icon: ⚙️
    title: Is the package ready to use?
    details: This package is a work in progress. Most of the functionality still needs to be implemented. The performance still needs to be optimized for both CPU and GPU. The API for this package might still be unstable. However, you are still welcomed to try it. <em><strong>Click for the details!</em></strong>
    link: /basics
---

## Development goals
 - <input type="checkbox" disabled checked> `Lattices` module to generate any Bravais lattice.
 - <input type="checkbox" disabled> `Networks` module to generate canonical artificial neural networks (ANN) via [Flux.jl](https://github.com/FluxML/Flux.jl). (*work in progress*)
 - <input type="checkbox" disabled> `VarStates` module to define variational quantum states. (*work in progress*)
 - <input type="checkbox" disabled checked> `Hilberts` module to define Hilbert spaces. 
 - <input type="checkbox" disabled checked> `Operators` module to define arbitrary quantum operators on a computational basis.
 - <input type="checkbox" disabled> `Samplers` module to sample variational quantum states with Markov chain Monte-Carlo (MCMC) methods.
 - <input type="checkbox" disabled> `Handlers` module to optimize variational quantum states with gradient-based methods.
 - <input type="checkbox" disabled> Support for distributed and parallel computing via [MPI.jl](https://github.com/JuliaParallel/MPI.jl/tree/master).
 - <input type="checkbox" disabled> GPU support via [CUDA.jl](https://github.com/JuliaGPU/CUDA.jl), [AMDGPU.jl](https://github.com/JuliaGPU/AMDGPU.jl), and [Metal.jl](https://github.com/JuliaGPU/Metal.jl).
```
## Bugs report and feature requests
If you think you have found a bug or have a feature request, you can open an [issue](https://github.com/cevenkadir/NeuralQuantumStates.jl/issues/new).

## Citation
If you use this package in your work, 
we would appreciate the following reference as in [CITATION.bib](https://github.com/cevenkadir/NeuralQuantumStates.jl/blob/main/CITATION.bib).