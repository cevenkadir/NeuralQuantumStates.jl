"""
    NQSAnsatzeLatticeSpaceGroupsExt

Builds symmetric ansätze straight from a lattice, deriving the site permutations from its
geometry rather than requiring them to be written by hand.

Hand-written permutations only ever work for a one-dimensional chain. This is what makes a
translation-invariant ansatz on a kagome torus no harder to write than one on a chain.
"""
module NQSAnsatzeLatticeSpaceGroupsExt

using NQSAnsatze
# `import`, not `using`: this adds a method to SymmetricRBM's constructor rather than
# defining a shadowing function in this module.
import NQSAnsatze: SymmetricRBM

using LatticeSpaceGroups
using LatticeSpaceGroups: Lattice, Translation, site_permutation, translation_group

"""
    SymmetricRBM(lattice, alpha; group=:translation, T=ComplexF64)

An RBM invariant under a symmetry group of `lattice`.

`group` selects which symmetries to impose:
- `:translation` — the full translation group, giving a translation-invariant wavefunction.
  This is the usual choice for a periodic lattice and reduces the parameter count by the number
  of cells.
- `:space` — the full space group, adding point-group operations. A stronger constraint, and
  correspondingly fewer parameters.

# Example
```julia
lat = build(Square(4; periodic=true))
model = SymmetricRBM(lat, 2)
```
"""
function SymmetricRBM(
    lattice::Lattice, alpha::Real=1; group::Symbol=:translation, T::Type=ComplexF64
)
    perms = if group === :translation
        [site_permutation(lattice, t) for t in translation_group(lattice)]
    elseif group === :space
        unique(site_permutation(lattice, op) for op in space_group(lattice))
    else
        throw(ArgumentError("group must be :translation or :space, got :$group"))
    end
    return SymmetricRBM(perms, alpha; T=T)
end

end # module NQSAnsatzeLatticeSpaceGroupsExt
