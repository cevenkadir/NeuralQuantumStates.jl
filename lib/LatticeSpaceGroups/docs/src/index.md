```@meta
CurrentModule = LatticeSpaceGroups
```

# LatticeSpaceGroups.jl

*Lattice geometry, and the space groups it induces.*

LatticeSpaceGroups.jl builds Bravais lattices and derives their symmetry operations as **site
permutations** — the form that symmetry-reduced basis construction actually consumes.

It depends on StaticArrays and nothing else. It is part of the
[NeuralQuantumStates.jl](https://github.com/cevenkadir/NeuralQuantumStates.jl) ecosystem, but it
is equally usable on its own for exact diagonalization.

## Why not just a graph?

A graph gives you vertices and edges. A quantum lattice needs more:

- **Positions and a Bravais basis**, so that "next-nearest neighbour" means something.
- **A periodic-boundary-aware metric**, so neighbour orders are well defined on a torus rather
  than being cut by the boundary.
- **Space groups**, which are the whole point: `SymBasis.Translational`,
  `SymBasis.SpatialReflection`, and `SymBasis.Rotational` each need a site permutation, and
  writing one by hand only ever works for a one-dimensional chain.

A [`Lattice`](@ref) is a plain struct holding its sites and bonds, so none of this costs a graph
dependency. If you *want* a graph, load Graphs.jl or MetaGraphsNext.jl and convert — see
[Optional integrations](@ref) below.

## Quick example

```@example index
using LatticeSpaceGroups

lat = build(Hypercube([4]; periodic=true))
translation_permutation(lat, 1)
```

The point group is found from lattice geometry, including for lattices whose symmetry centre is
not the coordinate origin:

```@example index
length(point_group(build(Kagome([3, 3], 1.0; periodic=true))))    # D6
```

Bonds come back as site-index pairs in the numbering [`site_positions`](@ref) uses, so they drop
straight into a Hamiltonian:

```@example index
bonds(build(Honeycomb([2, 2], 1.0; periodic=true)))
```

## Predefined lattices

[`build`](@ref) turns a spec into a [`Lattice`](@ref):

| Spec | Dim | Sites/cell | Coordination |
|---|---|---|---|
| [`Hypercube`](@ref), and [`Square`](@ref) / [`Cube`](@ref) | any | 1 | `2D` |
| [`Triclinic`](@ref) | 3 | 1 | 6 |
| [`Triangular`](@ref) | 2 | 1 | 6 |
| [`Honeycomb`](@ref) | 2 | 2 | 3 |
| [`Kagome`](@ref) | 2 | 3 | 4 |
| [`BCC`](@ref) | 3 | 1 | 8 |
| [`FCC`](@ref) | 3 | 1 | 12 |
| [`Diamond`](@ref) | 3 | 2 | 4 |
| [`Pyrochlore`](@ref) | 3 | 4 | 6 |

Anything else can be assembled from a [`LatticeBasis`](@ref) and a [`Lattice`](@ref) directly,
including with explicitly given bonds. See [Lattices](@ref).

## Optional integrations

None of these is a dependency; loading one activates an extension.

- **SymBasis.jl** — the symmetry constructors gain methods taking a lattice directly, so a
  momentum sector is one line. See [Symmetries](@ref).
- **Graphs.jl** — `SimpleGraph(lattice)`, plus `Graphs.nv` and `Graphs.ne` methods. (This
  package spells the site count [`n_sites`](@ref) precisely so it does not collide with
  `Graphs.nv` when both are loaded.)
- **MetaGraphsNext.jl** — `MetaGraph(lattice)`, carrying site labels as vertex labels,
  Cartesian positions as vertex data, and neighbour-shell orders as edge data.

## Manual

- [Lattices](@ref) — building them
- [Symmetries](@ref) — translations, point groups, space groups, and the SymBasis bridge
