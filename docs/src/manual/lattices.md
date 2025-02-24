```@meta
CurrentModule = NeuralQuantumStates
```

# Lattices

🚧 **Under construction** 🚧

## Functions

### `Lattices.build`
[`Lattices.build`](@ref) is used to build a `Lattices.Lattice` instance for common lattice structures, such as a hypercube.

::::tabs

== Hypercube

For instance, to build a simple chain lattice of length ``7`` with the edge length ``2.0`` and periodic boundary conditions, run
```@example
using NeuralQuantumStates: Lattices # hide
Lattices.build(
    :Hypercube, [7,], 2.0;
    periodic=[true,]
)
```

For a cuboid lattice of shape ``[3, 4, 5]`` with the edge length ``2.0`` and periodic boundary conditions in the second dimension, run
```@example
using NeuralQuantumStates: Lattices # hide
Lattices.build(
    :Hypercube, [3, 4, 5], 2.0;
    periodic=[false, true, false]
)
```

== Triclinic

For instance, to build a triclinic lattice of shape ``[3, 2, 1]`` with the edge lengths ``[1.0, 1.5, 2.0]`` and the angles ``[40.0, 65.0, 90.0]``, run
```@example
using NeuralQuantumStates: Lattices # hide
Lattices.build(
    :Triclinic, [3, 2, 1], [1.0, 1.5, 2.0], [40.0, 65.0, 90.0];
    periodic=false
)
```

== Triangular

For instance, to build a 2D triangular lattice of shape ``[2, 4]`` with the edge length ``1.0`` and periodic boundary conditions, run
```@example
using NeuralQuantumStates: Lattices # hide
Lattices.build(
    :Triangular, [2, 4], 1.0;
    periodic=[true, true]
)
```

== Honeycomb

For instance, to build a 2D honeycomb lattice of shape ``[3, 6]`` with the edge length ``2.5`` and open boundary conditions, run
```@example
using NeuralQuantumStates: Lattices # hide
Lattices.build(
    :Honeycomb, [3, 6], 2.5;
    periodic=false,
)
```

== Kagome

For instance, to build a 2D kagome lattice of shape ``[4, 2]`` with the edge length ``3.0`` and periodic boundary conditions, run
```@example
using NeuralQuantumStates: Lattices # hide
Lattices.build(
    :Kagome, [4, 2], 3.0;
    periodic=[true, true]
)
```
::::

### `Lattices.vertices`
[`Lattices.vertices`](@ref) is used to get the labels and the label data from of all the vertices of a given a `Lattices.Lattice` instance:
```@example 1
using NeuralQuantumStates: Lattices # hide
lat = Lattices.build(
    :Hypercube, [3, 4, 5], 2.0;
    periodic=[false, true, false]
);
labels, label_data = Lattices.vertices(lat.shape, lat.basis);
nothing # hide
```
```@example 1
labels
```
```@example 1
label_data
```