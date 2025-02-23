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
    details: This package is a work in progress. Most of the functionality still needs to be implemented. The performance still needs to be optimized for both CPU and GPU. The API for this package might still be unstable. However, you are still welcomed to try it.
---
```