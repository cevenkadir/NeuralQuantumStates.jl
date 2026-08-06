```@meta
CurrentModule = LatticeSpaceGroups
```

# Lattices

```@example lattices
using LatticeSpaceGroups
```

## How to build a lattice?

### via `build`
The most straightforward way to build a common lattice structure, such as a hypercube, is to utilize [`build`](@ref LatticeSpaceGroups.build). This function depending on arguments and keywords creates an instance of [`Lattice`](@ref LatticeSpaceGroups.Lattice) of the given lattice.


#### Hypercube

For instance, to build a simple chain lattice of length ``7`` with the edge length ``2.0`` and periodic boundary conditions, run
```@example lattices
build(Hypercube([7,], 2.0;
    periodic=[true,]
))
```

For a cuboid lattice of shape ``[3, 4, 5]`` with the edge length ``2.0`` and periodic boundary conditions in the second dimension, run
```@example lattices
build(Hypercube([3, 4, 5], 2.0;
    periodic=[false, true, false]
))
```

#### Triclinic

For instance, to build a triclinic lattice of shape ``[3, 2, 1]`` with the edge lengths ``[1.0, 1.5, 2.0]`` and the angles ``[40.0, 65.0, 90.0]``, run
```@example lattices
build(Triclinic([3, 2, 1], [1.0, 1.5, 2.0], [40.0, 65.0, 90.0];
    periodic=false
))
```

#### Triangular

For instance, to build a 2D triangular lattice of shape ``[2, 4]`` with the edge length ``1.0`` and periodic boundary conditions, run
```@example lattices
build(Triangular([2, 4], 1.0;
    periodic=[true, true]
))
```

#### Honeycomb

For instance, to build a 2D honeycomb lattice of shape ``[3, 6]`` with the edge length ``2.5`` and open boundary conditions, run
```@example lattices
build(Honeycomb([3, 6], 2.5;
    periodic=false,
))
```

#### Kagome

For instance, to build a 2D kagome lattice of shape ``[4, 2]`` with the edge length ``3.0`` and periodic boundary conditions, run
```@example lattices
build(Kagome([4, 2], 3.0;
    periodic=[true, true]
))
```

### via `LatticeBasis` and `Lattice`
If your lattice has a special structure or not already defined in this package, first, define your [`LatticeBasis`](@ref LatticeSpaceGroups.LatticeBasis) instance with the chosen basis vectors and lattice site offsets:
```@example lattices
basis_vectors = [
    1.0 0.25;
    0.1 0.4
];
site_offsets = [
    0.0 0.0;
    0.15 0.20
];
lat_basis = LatticeBasis(basis_vectors, site_offsets)
```
#### for edges from the ``k``-th nearest neighbors
If you want to consider a lattice of shape ``[4, 3]`` for the nearest neighbors, run


##### open boundary conditions

```@example lattices
Lattice(
    [4, 3], lat_basis;
    max_order=1
)
```

##### periodic boundary conditions

```@example lattices
Lattice(
    [4, 3], lat_basis, [true, true];
    max_order=1
)
```

##### custom boundary conditions

If the lattice has periodic boundary conditions only in the first dimension,
```@example lattices
Lattice(
    [4, 3], lat_basis, [true, false];
    max_order=1
)
```

#### for custom edges
If you want to define the custom edges for the lattice, first, define the custom edges

```@example lattices
custom_edges = [((:A, 1, 1), (:B, 3, 3)),], [1,];
nothing # hide
```


##### open boundary conditions

```@example lattices
Lattice(
    [4, 3], lat_basis, custom_edges
)
```

##### periodic boundary conditions

```@example lattices
Lattice(
    [4, 3], lat_basis, custom_edges, [true, true]
)
```

##### custom boundary conditions

If the lattice has periodic boundary conditions only in the first dimension,
```@example lattices
Lattice(
    [4, 3], lat_basis, custom_edges, [true, false]
)
```
