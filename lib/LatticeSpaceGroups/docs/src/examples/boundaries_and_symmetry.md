```@meta
CurrentModule = LatticeSpaceGroups
```

# Boundary conditions and symmetry

Boundary conditions are not a detail you set and forget: they decide which symmetries a finite
lattice actually has. A periodic direction carries translation symmetry, an open one does not,
and a mirror that sits at the origin on one lattice sits at the midpoint on another. All of that
is worked out from the geometry, but it is worth seeing explicitly.

```@example bc
using LatticeSpaceGroups
```

## Translations need periodicity

An open chain has no translation symmetry at all — shifting it moves sites off the end:

```@example bc
open_chain = build(Hypercube([6]; periodic=false))
translation_generators(open_chain)
```

Close it into a ring and the generator appears:

```@example bc
ring = build(Hypercube([6]; periodic=true))
translation_generators(ring)
```

```@example bc
translation_permutation(ring, 1)
```

Asking for a translation along an open direction is an error rather than a wrong answer:

```@example bc
try
    translation_permutation(open_chain, 1)
catch err
    err
end
```

## Mixed boundaries: the cylinder

Periodicity is per dimension, so a 4 × 4 lattice periodic in one direction only is a cylinder.
It has translation symmetry around the circumference and none along the axis:

```@example bc
cylinder = build(Hypercube([4, 4]; periodic=[true, false]))

length(translation_generators(cylinder)), length(translation_group(cylinder))
```

```@example bc
is_symmetry(cylinder, Translation([1, 0])), is_symmetry(cylinder, Translation([0, 1]))
```

Compare the torus, where the translation group is the product of both cyclic groups:

```@example bc
torus = build(Hypercube([4, 4]; periodic=true))
length(translation_group(torus))
```

## The mirror is not always at the origin

Reflections survive open boundaries, but the mirror plane moves. On an open chain it sits at the
**midpoint**, so the permutation is the reversal:

```@example bc
reflection_permutation(open_chain, 1)
```

Nothing about that was hard-coded. The reflection about the origin alone is *not* a symmetry —
it sends site 1 to a position where no site exists — and
[`compensating_translation`](@ref) finds the shift that repairs it:

```@example bc
mirror = [-1.0 ;;]

is_symmetry(open_chain, PointOperation(mirror))
```

```@example bc
compensating_translation(open_chain, mirror)
```

That is the same mechanism that gives honeycomb and diamond their full point groups, where the
symmetry centre is fixed by the site offsets rather than by the boundary.

## Point groups survive opening the boundary

Opening the boundary costs translations, but not necessarily the point group. Both the periodic
and the open square lattice have ``D_4``, of order 8:

```@example bc
length(point_group(torus)), length(point_group(build(Hypercube([4, 4]; periodic=false))))
```

The reflection is its own inverse in both cases, as a reflection should be:

```@example bc
compose(b, a) = [b[a[i]] for i in eachindex(a)]

p = reflection_permutation(torus, 1)
compose(p, p) == collect(1:n_sites(torus))
```

## Why `preserves_edges` exists

A permutation can move the sites correctly and still not be a symmetry of anything you would
build on the lattice, because it maps a bond onto a non-bond. Every operation this package
returns is checked against the edge set, not just the site set:

```@example bc
all(preserves_edges(cylinder, site_permutation(cylinder, op)) for op in space_group(cylinder))
```

That check is on by default; pass `check_edges=false` to [`point_group`](@ref) or
[`space_group`](@ref) if you genuinely want site-set symmetries only.
