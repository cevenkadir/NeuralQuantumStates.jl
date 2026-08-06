```@meta
CurrentModule = LatticeSpaceGroups
```

# LatticeSpaceGroups.jl

*Lattice geometry, and the space groups it induces.*

LatticeSpaceGroups.jl builds Bravais lattices and derives their symmetry operations as **site
permutations** — the form that symmetry-reduced basis construction actually consumes.

It has no machine-learning dependencies. It is part of the
[NeuralQuantumStates.jl](https://github.com/cevenkadir/NeuralQuantumStates.jl) ecosystem, but it
is equally usable on its own for exact diagonalization.

## Why not just a graph?

Graphs.jl gives you vertices and edges, and this package is built on it. But a quantum lattice
needs more than its connectivity:

- **Positions and a Bravais basis**, so that "next-nearest neighbour" means something.
- **A periodic-boundary-aware metric**, so neighbour orders are well defined on a torus rather
  than being cut by the boundary.
- **Space groups**, which are the whole point: `SymBasis.Translational`,
  `SymBasis.SpatialReflection`, and `SymBasis.Rotational` each need a site permutation, and
  writing one by hand only ever works for a one-dimensional chain.

## Quick example

```@example index
using LatticeSpaceGroups

lat = build(Hypercube([4], 1.0; periodic=[true]))
translation_permutation(lat, 1)
```

The point group is found from lattice geometry, including for lattices whose symmetry centre is
not the coordinate origin:

```@example index
length(point_group(build(Kagome([3, 3], 1.0; periodic=true))))    # D6
```

## Predefined lattices

`build` constructs hypercubic, triclinic, triangular, honeycomb, and kagome lattices; see
[Lattices](@ref). Anything else can be assembled from a [`LatticeBasis`](@ref) and a
[`Lattice`](@ref) directly, including with custom edges.

## Manual

- [Lattices](@ref) — building them
- [Symmetries](@ref) — translations, point groups, space groups, and the SymBasis bridge
