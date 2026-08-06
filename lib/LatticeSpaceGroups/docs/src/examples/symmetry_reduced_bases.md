```@meta
CurrentModule = LatticeSpaceGroups
```

# Symmetry-reduced bases

This is what the package is for. [SymBasis.jl](https://github.com/cevenkadir/SymBasis.jl) builds
a basis that conserves a symmetry, and it needs that symmetry as a **site permutation**. Loading
the two packages together activates an extension so the symmetry constructors take a lattice
directly.

```@example sectors
using LatticeSpaceGroups, SymBasis
```

## Momentum sectors of a chain

Take a periodic chain of eight spin-1/2 sites.

```@example sectors
nsites = 8
lat = build(Hypercube([nsites]; periodic=true))
dofo = dof_object(Spin(1 // 2))

full = length(basis(dofo, nsites).states)
```

Translational symmetry has one sector per momentum number `k`. `Translational(k, lat)` derives
its generator from the lattice rather than making you write it:

```@example sectors
dims = [length(basis(dofo, nsites, sym(Translational(k, lat), dofo)).states) for k in 0:(nsites-1)]
```

The decisive check is that these **partition** the full space — every state belongs to exactly
one momentum sector, so the dimensions must sum to `2^8`. If the permutation were wrong, this
sum would not come out:

```@example sectors
sum(dims), full
```

## Composing symmetries

Symmetries compose with `∘`, which is how you reach a small sector of a large space. Fixing the
total magnetization to zero first:

```@example sectors
sz = sym(TotalMagnetization(0 // 1, nsites), dofo)
length(basis(dofo, nsites, sz).states), binomial(nsites, nsites ÷ 2)
```

and then splitting *that* by momentum:

```@example sectors
sector_dims = [
    length(basis(dofo, nsites, sz ∘ sym(Translational(k, lat), dofo)).states)
    for k in 0:(nsites-1)
]
sum(sector_dims), binomial(nsites, nsites ÷ 2)
```

The largest sector is a fraction of the full ``2^8``-dimensional space, which is the entire
point of the exercise:

```@example sectors
maximum(sector_dims), full
```

## The same code in two dimensions

Nothing above was specific to a chain. On a 2 × 3 torus, `axis` selects which direction to
translate along — and this is the case where writing the permutation by hand stops being
realistic:

```@example sectors
square = build(Hypercube([2, 3]; periodic=true))
n = n_sites(square)

translation_permutation(square, 2)
```

```@example sectors
dims2 = [
    length(basis(dofo, n, sym(Translational(k, square; axis=2), dofo)).states)
    for k in 0:2
]
sum(dims2), 2^n
```

## Reflections and rotations

[`SpatialReflection`](https://cevenkadir.github.io/SymBasis.jl/stable/) takes a parity, and the
mirror plane is placed wherever the lattice actually admits one — not assumed to be at the
origin:

```@example sectors
open_chain = build(Hypercube([6]; periodic=false))
reflection_permutation(open_chain, 1)
```

That is the reversal `6, 5, 4, 3, 2, 1`, the mirror sitting at the chain's midpoint rather than
at site 1. See [Boundary conditions and symmetry](@ref) for why that matters.

For rotations, inspect what a lattice offers before asking for one — a chain has no non-trivial
proper rotation about a point, and `Rotational` throws rather than silently giving you the
identity:

```@example sectors
length(rotation_permutations(build(Square(3; periodic=true))))
```
