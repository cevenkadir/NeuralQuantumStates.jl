```@meta
CurrentModule = LatticeSpaceGroups
```

# Three-dimensional lattices

The three-dimensional specs are where the machinery earns its keep: BCC and FCC have primitive
cells that are not cubes, and diamond and pyrochlore add multi-site bases whose symmetry centres
are not at the coordinate origin. Their point groups come out right anyway.

```@example threed
using LatticeSpaceGroups
```

## Coordination numbers

The number of nearest neighbours per site is the quickest check that a lattice is what you
think it is. Build each on a periodic supercell and count:

```@example threed
coordination(lat) = 2 * length(bonds(lat)) / n_sites(lat)

specs = [
    "cube"       => Cube(3; periodic=true),
    "BCC"        => BCC([3, 3, 3], 1.0; periodic=true),
    "FCC"        => FCC([3, 3, 3], 1.0; periodic=true),
    "diamond"    => Diamond([3, 3, 3], 1.0; periodic=true),
    "pyrochlore" => Pyrochlore([3, 3, 3], 1.0; periodic=true),
]

for (name, spec) in specs
    lat = build(spec)
    println(rpad(name, 12), " sites: ", lpad(n_sites(lat), 4),
            "   coordination: ", Int(coordination(lat)))
end
```

Those are the textbook values — 6 for simple cubic, 8 for BCC, 12 for FCC, 4 for diamond and 6
for pyrochlore. Diamond is the interesting one: it is an FCC Bravais lattice with a two-site
basis, and its low coordination is what makes it a diamond rather than an FCC.

## Point groups

Every one of these lattices has the full cubic point group ``O_h``, of order 48:

```@example threed
for (name, spec) in specs
    println(rpad(name, 12), " |point group| = ", length(point_group(build(spec))))
end
```

That is a stronger statement than it looks for diamond and pyrochlore. Their symmetry axes do
not pass through the coordinate origin, so the bare orthogonal operation `R` maps sites to
places where no site is. [`compensating_translation`](@ref) finds the ``\tau`` that fixes this,
and ``\{R \mid \tau\}`` is the actual symmetry:

```@example threed
diamond = build(Diamond([2, 2, 2], 1.0; periodic=true))
inversion = -1.0 * [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]

is_symmetry(diamond, PointOperation(inversion))
```

```@example threed
τ = compensating_translation(diamond, inversion)
```

```@example threed
is_symmetry(diamond, SpaceOperation(PointOperation(inversion), τ))
```

Without that mechanism, diamond would report a proper subgroup of ``O_h`` and every
symmetry-reduced calculation built on it would silently use too few symmetries.

## Half the elements are proper rotations

For each of these, exactly half the point-group elements have determinant `+1`:

```@example threed
using LinearAlgebra: det

for (name, spec) in specs
    pg = point_group(build(spec))
    println(rpad(name, 12), " rotations: ", count(o -> det(o.matrix) ≈ 1, pg),
            " of ", length(pg))
end
```

## Splitting bonds by shell

Build with `max_order` above 1 to keep further neighbour shells, then select one with the
`order` keyword of [`bonds`](@ref):

```@example threed
bcc = build(BCC([3, 3, 3], 1.0; periodic=true), max_order=2)

length(bonds(bcc; order=1)), length(bonds(bcc; order=2))
```

BCC has 8 nearest neighbours and 6 next-nearest, so on 27 sites that is `27 * 8 / 2` and
`27 * 6 / 2` bonds respectively.
