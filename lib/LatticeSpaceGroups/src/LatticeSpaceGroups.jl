"""
    LatticeSpaceGroups

Lattice geometry and the space-group machinery that turns it into **site permutations**.

Symmetry-reduced exact diagonalization needs a **site permutation**: which site does site `i`
become under a translation, a reflection, a rotation? Writing one by hand only ever works for a
one-dimensional chain. This package derives them from lattice geometry instead, and depends on
**StaticArrays and nothing else**.

Describe a lattice with a spec — [`Hypercube`](@ref), [`Honeycomb`](@ref), [`Kagome`](@ref),
[`FCC`](@ref), [`Pyrochlore`](@ref) and the rest — turn it into a [`Lattice`](@ref) with
[`build`](@ref), then ask for its symmetries:

```julia
using LatticeSpaceGroups

lat = build(Honeycomb([3, 3], 1.0; periodic=true))
n_sites(lat)                  # 18
bonds(lat)                    # nearest-neighbour pairs, as site indices
length(point_group(lat))      # 12, the D₆ point group
translation_permutation(lat)  # feed this to SymBasis.Translational
```

Loading SymBasis.jl, Graphs.jl or MetaGraphsNext.jl activates an extension for each; none is a
dependency.

Full documentation, including worked examples:
<https://cevenkadir.github.io/NeuralQuantumStates.jl/LatticeSpaceGroups/>

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
