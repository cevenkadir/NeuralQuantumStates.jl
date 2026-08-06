```@meta
CurrentModule = LatticeSpaceGroups
```

# Lattices

```@example lattices
using LatticeSpaceGroups
```

## Why not just a graph?

A graph gives you vertices and edges. A quantum lattice needs more:

- **Positions and a Bravais basis**, so that "next-nearest neighbour" means something.
- **A periodic-boundary-aware metric**, so neighbour orders are well defined on a torus rather
  than being cut by the boundary.
- **Space groups**, which are the whole point: `SymBasis.Translational`,
  `SymBasis.SpatialReflection`, and `SymBasis.Rotational` each need a site permutation, and
  writing one by hand only ever works for a one-dimensional chain. See [Symmetries](@ref).

A [`Lattice`](@ref) is a plain struct holding its sites and bonds, so none of this costs a graph
dependency. If you *want* a graph, see [Converting to a graph](@ref) below.

## Predefined lattices

A lattice is described by a **spec**, and [`build`](@ref) turns a spec into a
[`Lattice`](@ref). Specs validate on construction, so a mistake is caught where you wrote it
rather than deep inside the build.

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

### Hypercubic, and its named cases

[`Square`](@ref) and [`Cube`](@ref) are shorthands for the two commonest
[`Hypercube`](@ref) shapes. A one-dimensional chain is `Hypercube([n])` — there is
deliberately no exported `Chain`, because that name collides with `Lux.Chain` and
`Flux.Chain`, and this package is meant to be loaded alongside them.

```@example lattices
build(Hypercube([7], 2.0; periodic=true))
```

`periodic` may be one flag per dimension, which is how you get a cylinder or a slab:

```@example lattices
build(Hypercube([3, 4, 5], 2.0; periodic=[false, true, false]))
```

### Triclinic

The least symmetric Bravais lattice: three independent edge lengths at three independent
angles, given in degrees.

```@example lattices
build(Triclinic([3, 2, 1], [1.0, 1.5, 2.0], [40.0, 65.0, 90.0]; periodic=false))
```

### Two-dimensional lattices

```@example lattices
build(Triangular([2, 4], 1.0; periodic=true))
```

```@example lattices
build(Honeycomb([3, 6], 2.5; periodic=false))
```

```@example lattices
build(Kagome([4, 2], 3.0; periodic=true))
```

For [`Honeycomb`](@ref) and [`Kagome`](@ref), `edge_length` is the length of the *primitive
vectors*, not the bond length — the sublattice offsets place the bonds at `edge_length / √3`
and `edge_length / 2` respectively.

### Three-dimensional lattices

[`BCC`](@ref) and [`FCC`](@ref) are built in the **primitive** basis, so `shape` counts
primitive cells: an `FCC([3, 3, 3])` has 27 sites, not the 108 of 27 conventional cubic cells.
Their `edge_length` is the conventional cubic cell parameter `a`.

```@example lattices
build(BCC([3, 3, 3], 1.0; periodic=true))
```

```@example lattices
build(FCC([3, 3, 3], 1.0; periodic=true))
```

[`Diamond`](@ref) and [`Pyrochlore`](@ref) put a two- and four-site basis on the FCC lattice:

```@example lattices
build(Diamond([2, 2, 2], 1.0; periodic=true))
```

```@example lattices
build(Pyrochlore([2, 2, 2], 1.0; periodic=true))
```

## Sites and bonds

Sites are numbered with the sublattice index varying fastest, then the first cell coordinate,
and so on. That numbering is the lattice's interface to everything else — it is what
[`site_permutation`](@ref) permutes and what a Hilbert space assigns its degrees of freedom to.

```@example lattices
lat = build(Honeycomb([2, 2], 1.0; periodic=true))
site_labels(lat)
```

```@example lattices
site_positions(lat)
```

[`bonds`](@ref) reports pairs in that same numbering, each once with `i < j`, in lexicographic
order — so they can be used directly as the site identifiers of an operator term.

```@example lattices
bonds(lat)
```

## Further neighbour shells

`max_order` includes more than the nearest neighbours. [`bonds`](@ref) can then select a shell.

```@example lattices
square = build(Square(4, 1.0; periodic=true); max_order=2)
length(bonds(square; order=1)), length(bonds(square; order=2))
```

The second shell of a square lattice is its cell diagonals — which is why `max_order` is not
simply "more grid steps".

## Building a lattice by hand

For a lattice this package does not predefine, give a [`LatticeBasis`](@ref) its primitive
vectors (one per column) and its site offsets:

```@example lattices
basis_vectors = [
    1.0 0.25
    0.1 0.4
]
site_offsets = [
    0.0 0.0
    0.15 0.20
]
lat_basis = LatticeBasis(basis_vectors, site_offsets)
```

The primitive vectors must be linearly independent; a singular set describes no
two-dimensional lattice and is rejected.

### Bonds from neighbour shells

```@example lattices
Lattice([4, 3], lat_basis, [true, true]; max_order=1)
```

Boundary conditions are per dimension, so a cylinder is `[true, false]`:

```@example lattices
Lattice([4, 3], lat_basis, [true, false]; max_order=1)
```

### Bonds given explicitly

When the connectivity is not distance-derived, pass the bonds as **site-index pairs**. The
`orders` keyword records which shell each belongs to, if that distinction matters to you.

```@example lattices
Lattice([4, 3], lat_basis, [(1, 8), (2, 15)], [true, true])
```

## Converting to a graph

Loading Graphs.jl or MetaGraphsNext.jl activates an extension that converts a lattice.
A `MetaGraph` carries the geometry a `SimpleGraph` has nowhere to put:

```julia
using LatticeSpaceGroups, MetaGraphsNext

mg = MetaGraph(build(Honeycomb([3, 3], 1.0; periodic=true)))
mg[(1, 1, 1)]                # position of sublattice 1 in cell (1, 1)
mg[(1, 1, 1), (2, 1, 1)]     # 1, a nearest-neighbour bond
```
