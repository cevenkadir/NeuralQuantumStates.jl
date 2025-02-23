## Installation
If you still want to try it out, you can install it from the Julia REPL by entering:
```julia
julia> import Pkg; Pkg.add("https://github.com/cevenkadir/NeuralQuantumStates.jl")
```

## Basics
You can start to use the package by entering:
```@example 1
using NeuralQuantumStates: Lattices, Hilberts, Operators
```
As an example, let's construct the Hamiltonian of a transverse-field Ising chain. For that, first initialize a chain lattice of length 8:
```@example 1
lat = Lattices.build(:Hypercube, [8], 1.0; periodic=[true])
```
Then, choose a spin-1/2 as the basis for the Hilbert space:
```@example 1
hil = Hilberts.build(:Spin, 1 // 2, Lattices.nv(lat))
```
With these two, construct the Hamiltonian operator of the model:
```@example 1
ham = Operators.build(:TransverseFieldIsing, hil, lat; J=1.0, h_x=1.0, h_z=1.0)
```
To test it, let's find out the connected basis configurations, `s_prime`, to the Hamiltonian for the following 2 basis configurations with their matrix elements, `mels`:
```@example 1
s_prime, mels = Operators.connected_basis_configs(
    ham,
    [
        1//2 1//2 -1//2 1//2 -1//2 -1//2 1//2 -1//2;
        -1//2 1//2 1//2 1//2 1//2 -1//2 -1//2 -1//2;
    ]
);
nothing # hide
```
```@example 1
s_prime
```
```@example 1
mels
```
