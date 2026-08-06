```@meta
CurrentModule = LatticeSpaceGroups
```

# Symmetries

A lattice is more than a graph: it has geometry, and that geometry has symmetries. This page
covers how to get at them, and — the part that matters most in practice — how to turn them into
the **site permutations** that symmetry-reduced basis construction needs.

```@example sym
using LatticeSpaceGroups
```

## Site indexing

Everything here is expressed in terms of site numbers. Sites are numbered with the sublattice
index varying fastest, then the first cell coordinate, and so on; [`site_positions`](@ref) and
[`site_labels`](@ref) give the Cartesian positions and the lattice labels in that order, and it
is the same order a Hilbert space built on the lattice will use for its degrees of freedom.

```@example sym
lat = build(Hypercube([4]; periodic=true))
site_positions(lat)
```

## Translations

[`translation_permutation`](@ref) gives the permutation induced by shifting the lattice by whole
primitive cells. On a periodic chain this is just the cyclic shift:

```@example sym
translation_permutation(lat, 1)
```

[`translation_generators`](@ref) returns one generator per periodic direction — these are what
you hand to a symmetry group, which derives the cyclic group from the single generator. An open
direction admits no translation symmetry and is skipped:

```@example sym
open_lat = build(Hypercube([4]; periodic=false))
translation_generators(open_lat)
```

## Point and space groups

[`point_group`](@ref) finds the orthogonal operations that are symmetries of the lattice. Rather
than guessing rotation angles, candidates are found from the condition that a point operation
maps a Bravais lattice onto itself exactly when it is an integer matrix in the basis of
primitive vectors that preserves the metric tensor.

```@example sym
square = build(Square(4; periodic=true))
length(point_group(square))    # D4
```

```@example sym
triangular = build(Triangular([3, 3], 1.0; periodic=true))
length(point_group(triangular))    # D6
```

!!! note "The symmetry centre is found for you"
    A lattice with more than one site per unit cell generally has its symmetry centre somewhere
    other than the coordinate origin — a honeycomb's six-fold axis runs through a hexagon
    centre, not through the origin — so its rotations are symmetries only when paired with a
    compensating fractional translation. [`compensating_translation`](@ref) finds that
    translation, which is why a honeycomb reports the full ``D_6`` rather than the smaller
    group you would get by insisting on operations about the origin. The same mechanism puts an
    open chain's mirror at its midpoint.

```@example sym
honeycomb = build(Honeycomb([3, 3], 1.0; periodic=true))
length(point_group(honeycomb))    # D6, about the hexagon centre
```

The same machinery works in three dimensions, where the multi-site bases are what make the
compensating translation indispensable:

```@example sym
length(point_group(build(Pyrochlore([2, 2, 2], 1.0; periodic=true))))    # O_h
```

[`space_group`](@ref) combines the point group with the translation group, returning each
distinct site permutation once.

## Checking an operation

[`is_symmetry`](@ref) tests an operation without throwing, and by default requires that the
lattice's **bond set** be preserved as well as its site set. That distinction is not pedantic: a
lattice built with explicit bonds, or one with mixed boundary conditions, can admit an operation
that permutes sites correctly while mapping a bond onto a non-bond, and such an operation is not
a symmetry of any Hamiltonian defined on those bonds.

```@example sym
mixed = build(Hypercube([4, 4], 1.0; periodic=[true, false]))
is_symmetry(mixed, Translation([1, 0])), is_symmetry(mixed, Translation([0, 1]))
```

## Using them with SymBasis.jl

Loading [SymBasis.jl](https://github.com/cevenkadir/SymBasis.jl) activates an extension that
lets its symmetry specifications be built straight from a lattice, instead of from a
hand-written permutation:

```julia
using LatticeSpaceGroups, SymBasis

lat  = build(Hypercube([8]; periodic=true))
dofo = dof_object(Spin(1 // 2))

# Momentum-zero sector of an 8-site spin-1/2 chain.
b = basis(dofo, 8, sym(Translational(0, lat), dofo))

# Combined with magnetization conservation.
sz = sym(TotalMagnetization(0 // 1, 8), dofo)
b = basis(dofo, 8, sz ∘ sym(Translational(0, lat), dofo))
```

`Translational`, `SpatialReflection`, and `Rotational` all gain a method taking a lattice.
Summing the dimensions of all momentum sectors recovers the full space, which is the sharpest
check that the generated permutations are correct.
