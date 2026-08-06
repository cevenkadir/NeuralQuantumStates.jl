"""
    LatticeSpaceGroups

Lattice geometry and the space-group machinery that turns it into **site permutations**.

A quantum lattice needs more than a graph. This package supplies:

- Bravais lattice bases and site positions ([`LatticeBasis`](@ref), [`Lattice`](@ref))
- A periodic-boundary-aware minimum-image metric, so neighbour *orders* (nearest, next-nearest,
  …) are well defined on a torus
- The predefined lattice zoo, reached through [`build`](@ref): [`Hypercube`](@ref),
  [`Square`](@ref), [`Cube`](@ref), [`Triclinic`](@ref),
  [`Triangular`](@ref), [`Honeycomb`](@ref), [`Kagome`](@ref), [`BCC`](@ref), [`FCC`](@ref),
  [`Diamond`](@ref), and [`Pyrochlore`](@ref)
- Translation and point groups, and the site permutations they induce
  ([`site_permutation`](@ref), [`translation_group`](@ref), [`point_group`](@ref))

That last item is the reason this package exists on its own. `SymBasis.Translational`,
`SymBasis.SpatialReflection`, and `SymBasis.Rotational` all take a raw permutation vector;
writing those by hand only ever works for a 1-D chain. Generating them from lattice geometry is
useful to anyone doing symmetry-reduced exact diagonalization, with or without a neural network
anywhere in sight — so this package depends on **StaticArrays and nothing else**.

# Example
```julia
using LatticeSpaceGroups

lat = build(Honeycomb([3, 3], 1.0; periodic=true))
n_sites(lat)          # 18
bonds(lat)            # nearest-neighbour pairs, as site indices
length(point_group(lat))   # 12, the D₆ point group
```

# Optional integrations

Loading any of these packages alongside this one activates an extension; none is a dependency.

- **SymBasis.jl** — the symmetry constructors gain methods taking a lattice directly:
  ```julia
  using LatticeSpaceGroups, SymBasis
  lat = build(Hypercube([8]; periodic=true))
  dofo = dof_object(Spin(1 // 2))
  basis(dofo, 8, sym(Translational(0, lat), dofo))     # momentum-zero sector
  ```
- **Graphs.jl** — `SimpleGraph(lattice)`, plus `Graphs.nv`/`Graphs.ne` methods.
- **MetaGraphsNext.jl** — `MetaGraph(lattice)`, carrying site labels, positions, and
  neighbour-shell orders.

!!! note "Scope"
    Irreducible representations and character tables at given wave vectors are deliberately out
    of scope.
"""
module LatticeSpaceGroups

# Geometry: lattice bases, the minimum-image metric, neighbour-shell search.
include("lattice.jl")
# The predefined lattice zoo, reached through `build`.
include("predefined_lattices.jl")
# Translations, point group, and the site permutations they induce.
include("space_group.jl")

# geometry
export AbstractLatticeBasis, LatticeBasis
export AbstractLattice, Lattice
export n_sites, site_positions, site_labels, bonds

# predefined lattice specifications, and the builder that turns one into a `Lattice`
export AbstractLatticeSpec
export Hypercube, Square, Cube, Triclinic
export Triangular, Honeycomb, Kagome
export BCC, FCC, Diamond, Pyrochlore
export build

# symmetry operations
export AbstractSymmetryOperation, Translation, PointOperation, SpaceOperation, point_part

# site permutations
export site_permutation, is_symmetry, preserves_edges

# groups
export translation_generators, translation_group, point_group, space_group
export compensating_translation

# ready-made generators for SymBasis
export translation_permutation, reflection_permutation, rotation_permutations

end # module LatticeSpaceGroups
