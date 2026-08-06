"""
    LatticeSpaceGroups

Lattice geometry and the space-group machinery that turns it into **site permutations**.

This is the Julia counterpart of NetKet's `netket/graph`. Graphs.jl supplies the graph half of
that; what this package adds is everything else a quantum lattice needs:

- Bravais lattice bases and site positions ([`LatticeBasis`](@ref), [`Lattice`](@ref))
- A periodic-boundary-aware distance metric, so neighbour *orders* (nearest, next-nearest, ...)
  are well defined on a torus
- The predefined lattice zoo: [`build`](@ref) for hypercubic, triclinic, triangular, honeycomb
  and kagome lattices
- Translation and point groups, and the site permutations they induce
  ([`site_permutation`](@ref), [`translation_group`](@ref), [`point_group`](@ref))

That last item is the reason this package exists as its own package rather than as part of the
NQS stack. `SymBasis.Translational`, `SymBasis.SpatialReflection`, and `SymBasis.Rotational` all
take a raw permutation vector; writing those by hand only ever works for a 1-D chain. Generating
them from lattice geometry is useful to anyone doing symmetry-reduced exact diagonalization,
with or without a neural network anywhere in sight — so this package must not depend on the
neural-network stack, and does not.

Loading SymBasis alongside this package activates `LatticeSpaceGroupsSymBasisExt`, which lets
the symmetry constructors take a lattice directly:

```julia
using LatticeSpaceGroups, SymBasis

lat = build(:Hypercube, [8], 1.0; periodic=[true])
dofo = dof_object(Spin(1 // 2))
sg = sym(Translational(0, lat), dofo)      # momentum-zero sector
basis(dofo, 8, sg)
```

!!! note "Scope"
    Irreducible representations and character tables at given wave vectors are deliberately out
    of scope for this first cut.
"""
module LatticeSpaceGroups

using LinearAlgebra: Diagonal

# Geometry: lattice bases, the periodic metric, neighbour-order search.
include("lattice.jl")
# The predefined lattice zoo, reached through `build`.
include("predefined_lattices.jl")
# Translations, point group, and the site permutations they induce.
include("space_group.jl")

# geometry
export AbstractLatticeBasis, LatticeBasis
export AbstractLattice, Lattice
export vertices, nv

# predefined lattice specifications, and the builder that turns one into a `Lattice`
export AbstractLatticeSpec
export Hypercube, Triclinic, Triangular, Honeycomb, Kagome
export build

# symmetry operations
export AbstractSymmetryOperation, Translation, PointOperation, SpaceOperation, point_part

# site indexing and permutations
export site_positions, site_labels, site_permutation, bonds
export is_symmetry, preserves_edges

# groups
export translation_generators, translation_group, point_group, space_group
export compensating_translation

# ready-made generators for SymBasis
export translation_permutation, reflection_permutation, rotation_permutations

end # module LatticeSpaceGroups
