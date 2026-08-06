using Graphs: edges as _edges_of, src as _src, dst as _dst
using LinearAlgebra: det
using MetaGraphsNext: code_for, label_for, labels
using StaticArrays

"""
    AbstractSymmetryOperation{D}

A symmetry operation acting on the sites of a `D`-dimensional lattice.

Every concrete operation can be turned into a site permutation with
[`site_permutation`](@ref); that permutation is the object the rest of the ecosystem — and
SymBasis.jl in particular — actually consumes.
"""
abstract type AbstractSymmetryOperation{D} end

"""
    Translation{D,Tᵢ<:Integer} <: LatticeSpaceGroups.AbstractSymmetryOperation{D}

Translation by an integer number of primitive cells.

# Fields
- `displacement::SVector{D,Tᵢ}`: The displacement in units of the primitive vectors, so a
    displacement of `[1, 0]` shifts by one cell along the first primitive vector.
"""
struct Translation{D,Tᵢ<:Integer} <: AbstractSymmetryOperation{D}
    displacement::SVector{D,Tᵢ}
end
Translation(d::AbstractVector{Tᵢ}) where {Tᵢ<:Integer} =
    Translation{length(d),Tᵢ}(SVector{length(d),Tᵢ}(d))

"""
    PointOperation{D,T<:Real} <: LatticeSpaceGroups.AbstractSymmetryOperation{D}

An orthogonal transformation about the Cartesian origin — a rotation, reflection, inversion,
or a product of those.

# Fields
- `matrix::SMatrix{D,D,T}`: The orthogonal matrix acting on Cartesian site positions.
"""
struct PointOperation{D,T<:Real} <: AbstractSymmetryOperation{D}
    matrix::SMatrix{D,D,T}
end
PointOperation(m::AbstractMatrix{T}) where {T<:Real} =
    PointOperation{size(m, 1),T}(SMatrix{size(m, 1),size(m, 1),T}(m))

"""
    SpaceOperation{D,T} <: LatticeSpaceGroups.AbstractSymmetryOperation{D}

The general element of a space group, in Seitz notation ``\\{R \\mid \\tau\\}``: apply the
orthogonal part `R`, then translate by `τ`.

The translation is a **Cartesian** vector, not an integer number of cells, and that is
essential rather than incidental. In a lattice with more than one site per unit cell the
symmetry centre generally does not sit at the Cartesian origin — a honeycomb lattice's
six-fold axis passes through a hexagon centre, not through the coordinate origin — so its
rotations are only symmetries when paired with a compensating fractional translation. The same
mechanism handles open boundaries, where a reflection about the origin composed with a
translation is a reflection about the lattice's midpoint. Restricting `τ` to whole cells would
lose the point group of every such lattice.

# Fields
- `matrix::SMatrix{D,D,T}`: The orthogonal part ``R``.
- `offset::SVector{D,T}`: The Cartesian translation ``τ``.
"""
struct SpaceOperation{D,T<:Real} <: AbstractSymmetryOperation{D}
    matrix::SMatrix{D,D,T}
    offset::SVector{D,T}
end
SpaceOperation(p::PointOperation{D,T}, offset::AbstractVector) where {D,T} =
    SpaceOperation{D,T}(p.matrix, SVector{D,T}(offset))

"""
    point_part(op) -> PointOperation

The orthogonal part of a space-group operation, discarding its translation.
"""
point_part(op::SpaceOperation{D,T}) where {D,T} = PointOperation{D,T}(op.matrix)

# ---------------------------------------------------------------------------- site geometry

"""
    site_positions(lattice) -> Vector{SVector{D,T}}

Cartesian positions of the lattice sites, indexed by graph vertex number.

The ordering matters: it is the indexing that [`site_permutation`](@ref) permutes, and the
same one a Hilbert space built on this lattice will use for its degrees of freedom.
"""
function site_positions(lattice::Lattice{Tᵢ,T,D,O}) where {Tᵢ<:Integer,T<:Real,D,O}
    mg = lattice.metagraph
    positions = Vector{SVector{D,T}}(undef, _nv(mg))
    for l in labels(mg)
        positions[code_for(mg, l)] = mg[l]
    end
    return positions
end

"""
    site_labels(lattice) -> Vector

Lattice site labels, indexed by graph vertex number, in the same order as
[`site_positions`](@ref).
"""
function site_labels(lattice::Lattice{Tᵢ,T,D,O}) where {Tᵢ<:Integer,T<:Real,D,O}
    mg = lattice.metagraph
    return [label_for(mg, v) for v in 1:_nv(mg)]
end

"""
    bonds(lattice; order=nothing) -> Vector{Tuple{Int,Int}}

The lattice's edges as pairs of site indices, each listed once with `i < j`.

Site indices are graph vertex numbers, the same indexing as [`site_positions`](@ref) — so a
bond can be used directly as the site identifiers of an operator term.

`order` selects a neighbour shell: `1` for nearest neighbours, `2` for next-nearest, and so on,
matching the `max_order` the lattice was built with. The default, `nothing`, returns every bond
regardless of shell.

# Example
```julia
lat = build(Hypercube([4], 1.0; periodic=[true]))
bonds(lat)    # [(1,2), (2,3), (3,4), (1,4)]
```
"""
function bonds(
    lattice::Lattice{Tᵢ,T,D,O}; order::Union{Nothing,Integer}=nothing
) where {Tᵢ<:Integer,T<:Real,D,O}
    mg = lattice.metagraph
    out = Tuple{Int,Int}[]
    for e in _edges_of(mg)
        i, j = Int(_src(e)), Int(_dst(e))
        if order !== nothing
            # Edge data records which neighbour shell the bond belongs to.
            mg[label_for(mg, i), label_for(mg, j)] == order || continue
        end
        push!(out, (min(i, j), max(i, j)))
    end
    return out
end

"""
    _site_key(lattice, position; tol_digits) -> NTuple{D,Int}

A hashable key identifying a site, invariant under the lattice's periodic boundary conditions.

The position is expressed in fractional (primitive-cell) coordinates and, along each periodic
direction, folded back into `[0, shape)`. Two positions differing by a supercell translation
therefore produce the same key, which makes matching an operation's image back onto a site a
dictionary lookup rather than a search.

The key is built from **integers**, not rounded floats, and that matters. Folding a coordinate
in floating point computes the same lattice site two different ways — a site stored directly at
``8/3`` versus one reached as ``-1/3`` and wrapped by ``+3`` — and those can differ in the last
bit, so the two never compare equal and a genuine symmetry gets rejected. Scaling to integers
first and folding with integer `mod` makes the two paths agree exactly.
"""
function _site_key(
    lattice::Lattice{Tᵢ,T,D,O}, position::AbstractVector; tol_digits::Integer=TOL_DIGITS
) where {Tᵢ<:Integer,T<:Real,D,O}
    fractional = lattice.basis.vectors \ SVector{D,T}(position)
    scale = 10^tol_digits
    key = ntuple(D) do d
        n = round(Int, Float64(fractional[d]) * scale)
        lattice.periodic[d] ? mod(n, Int(lattice.shape[d]) * scale) : n
    end
    return key
end

"""
    apply(operation, position, lattice) -> SVector

Image of a Cartesian `position` under `operation`. Point operations act about the Cartesian
origin.
"""
function apply(
    op::Translation{D}, position::AbstractVector, lattice::Lattice{Tᵢ,T,D,O}
) where {Tᵢ<:Integer,T<:Real,D,O}
    return SVector{D,T}(position) + lattice.basis.vectors * SVector{D,T}(op.displacement)
end
function apply(
    op::PointOperation{D}, position::AbstractVector, lattice::Lattice{Tᵢ,T,D,O}
) where {Tᵢ<:Integer,T<:Real,D,O}
    return SMatrix{D,D,T}(op.matrix) * SVector{D,T}(position)
end
function apply(
    op::SpaceOperation{D}, position::AbstractVector, lattice::Lattice{Tᵢ,T,D,O}
) where {Tᵢ<:Integer,T<:Real,D,O}
    return SMatrix{D,D,T}(op.matrix) * SVector{D,T}(position) + SVector{D,T}(op.offset)
end

# ------------------------------------------------------------------------ site permutations

"""
    site_permutation(lattice, operation; tol_digits=12) -> Vector{Int}

The site permutation induced by `operation`, with `perm[i]` the site that site `i` is mapped
*to*.

This is the bridge to SymBasis.jl: `SymBasis.Translational`, `SymBasis.SpatialReflection`, and
`SymBasis.Rotational` each take exactly such a vector as their generator.

Throws an `ArgumentError` if `operation` is not a symmetry of the lattice — that is, if it
moves some site to a position where no site exists (after folding through the periodic
boundary conditions), or if it maps two sites onto the same one. Use
[`is_symmetry`](@ref) to test without throwing.

# Example
For a periodic chain, translating by one cell is the cyclic shift:
```julia
lat = build(:Hypercube, [4], 1.0; periodic=[true])
site_permutation(lat, Translation([1])) == [2, 3, 4, 1]
```
"""
function site_permutation(
    lattice::Lattice{Tᵢ,T,D,O},
    operation::AbstractSymmetryOperation{D};
    tol_digits::Integer=TOL_DIGITS
) where {Tᵢ<:Integer,T<:Real,D,O}
    positions = site_positions(lattice)
    n = length(positions)

    lookup = Dict{NTuple{D,Int},Int}()
    for (i, p) in enumerate(positions)
        lookup[_site_key(lattice, p; tol_digits=tol_digits)] = i
    end

    perm = Vector{Int}(undef, n)
    for (i, p) in enumerate(positions)
        image = apply(operation, p, lattice)
        key = _site_key(lattice, image; tol_digits=tol_digits)
        j = get(lookup, key, 0)
        if j == 0
            throw(ArgumentError(
                "$operation is not a symmetry of this lattice: site $i is mapped to " *
                "$(Vector(image)), where there is no lattice site"
            ))
        end
        perm[i] = j
    end

    if !isperm(perm)
        throw(ArgumentError(
            "$operation is not a symmetry of this lattice: it does not map the sites " *
            "one-to-one"
        ))
    end
    return perm
end

"""
    is_symmetry(lattice, operation; check_edges=true, tol_digits=12) -> Bool

Whether `operation` is a symmetry of `lattice`.

With `check_edges=true` (the default) the operation must preserve the lattice's edge set as
well as its site set. That distinction matters: a lattice built with `custom_edges`, or one
with open boundaries in some directions, can admit an operation that permutes the sites
correctly while mapping a bond onto a non-bond — such an operation is not a symmetry of any
Hamiltonian defined on those bonds.
"""
function is_symmetry(
    lattice::Lattice{Tᵢ,T,D,O},
    operation::AbstractSymmetryOperation{D};
    check_edges::Bool=true,
    tol_digits::Integer=TOL_DIGITS
) where {Tᵢ<:Integer,T<:Real,D,O}
    perm = try
        site_permutation(lattice, operation; tol_digits=tol_digits)
    catch err
        err isa ArgumentError && return false
        rethrow()
    end
    return check_edges ? preserves_edges(lattice, perm) : true
end

"""
    preserves_edges(lattice, perm) -> Bool

Whether the site permutation `perm` maps the lattice's edge set onto itself.
"""
function preserves_edges(
    lattice::Lattice{Tᵢ,T,D,O}, perm::AbstractVector{<:Integer}
) where {Tᵢ<:Integer,T<:Real,D,O}
    mg = lattice.metagraph
    # Edges are undirected, so store each as a sorted pair.
    edge_set = Set((min(_src(e), _dst(e)), max(_src(e), _dst(e))) for e in _edges_of(mg))
    for e in _edges_of(mg)
        i, j = perm[_src(e)], perm[_dst(e)]
        (min(i, j), max(i, j)) in edge_set || return false
    end
    return true
end

# ------------------------------------------------------------------------------- the groups

"""
    translation_generators(lattice) -> Vector{Translation}

One unit translation per periodic direction — the generators of the translation group.

These are what you hand to SymBasis: each generates a cyclic group of order `shape[d]`, and
`SymBasis.Translational` derives that cyclic group from the single generator permutation.
Non-periodic directions admit no translation symmetry and are skipped, so an open lattice
gives an empty vector.
"""
function translation_generators(
    lattice::Lattice{Tᵢ,T,D,O}
) where {Tᵢ<:Integer,T<:Real,D,O}
    generators = Translation{D,Int}[]
    for d in 1:D
        # A direction of extent 1 wraps onto itself: the "translation" is the identity, which
        # is not a usable generator.
        (lattice.periodic[d] && lattice.shape[d] > 1) || continue
        push!(generators, Translation(SVector{D,Int}(ntuple(i -> i == d ? 1 : 0, D))))
    end
    return generators
end

"""
    translation_group(lattice) -> Vector{Translation}

Every translation that is a symmetry of `lattice`, including the identity.

The group is the direct product of the cyclic groups along each periodic direction, so its
order is `prod(shape[d] for d in periodic directions)`.
"""
function translation_group(lattice::Lattice{Tᵢ,T,D,O}) where {Tᵢ<:Integer,T<:Real,D,O}
    ranges = ntuple(d -> lattice.periodic[d] ? (0:(lattice.shape[d]-1)) : (0:0), D)
    return [
        Translation(SVector{D,Int}(δ)) for δ in Iterators.product(ranges...)
    ] |> vec
end

"""
    _integer_point_candidates(basis) -> Vector{SMatrix}

Candidate point operations, as Cartesian matrices.

A point operation maps the Bravais lattice onto itself exactly when, expressed in the basis
of primitive vectors, it is an **integer** matrix `M` — and it is orthogonal exactly when `M`
preserves the metric tensor `G = AᵀA`, that is `MᵀGM = G`. So instead of guessing rotation
angles, enumerate integer matrices and keep those satisfying that condition; the Cartesian
operation is then `R = A M A⁻¹`.

Entries are restricted to `{-1, 0, 1}`, which is sufficient for every crystallographic point
group in the primitive basis (the 60° rotation of a triangular lattice, for instance, is
`[1 -1; 1 0]`).
"""
function _integer_point_candidates(basis::AbstractLatticeBasis{T,D,O}) where {T<:Real,D,O}
    A = basis.vectors
    G = A' * A

    candidates = SMatrix{D,D,Float64}[]
    for entries in Iterators.product(ntuple(_ -> (-1, 0, 1), D * D)...)
        M = SMatrix{D,D,Int}(entries...)
        abs(det(M)) ≈ 1 || continue
        isapprox(M' * G * M, G; atol=1e-10, rtol=1e-10) || continue
        push!(candidates, SMatrix{D,D,Float64}(A * M / A))
    end
    return candidates
end

"""
    compensating_translation(lattice, R; check_edges=true) -> Union{SVector,Nothing}

The Cartesian translation ``τ`` that makes ``\\{R \\mid τ\\}`` a symmetry of `lattice`, or
`nothing` if no such translation exists.

`R` alone is rarely a symmetry of a lattice with a multi-site unit cell, because the symmetry
centre need not be the coordinate origin. Any valid ``τ`` must carry the image of site 1 onto
*some* site, so it suffices to try the `nv` candidates ``τ = r_j - R r_1`` and keep the first
that works — turning an unbounded search into a linear one.
"""
function compensating_translation(
    lattice::Lattice{Tᵢ,T,D,O},
    R::AbstractMatrix;
    check_edges::Bool=true,
    tol_digits::Integer=TOL_DIGITS
) where {Tᵢ<:Integer,T<:Real,D,O}
    positions = site_positions(lattice)
    Rmat = SMatrix{D,D,Float64}(R)
    image₁ = Rmat * SVector{D,Float64}(positions[1])

    for target in positions
        τ = SVector{D,Float64}(target) - image₁
        op = SpaceOperation{D,Float64}(Rmat, τ)
        if is_symmetry(lattice, op; check_edges=check_edges, tol_digits=tol_digits)
            return τ
        end
    end
    return nothing
end

"""
    point_group(lattice; check_edges=true) -> Vector{PointOperation}

The point group of `lattice`: every orthogonal operation that is a symmetry of the lattice,
possibly once paired with a compensating translation. Includes the identity.

Candidates come from the integer-matrix condition described in `_integer_point_candidates` and
are then filtered by whether some ``\\{R \\mid τ\\}`` actually permutes this lattice's sites,
which accounts for the site offsets of a multi-site unit cell, the supercell shape, and the
boundary conditions. With `check_edges=true` the edge set must be preserved too.

The returned operations carry only the orthogonal part `R`; use [`space_group`](@ref) to get
the operations complete with their translations, or [`compensating_translation`](@ref) to
recover the translation for a particular `R`.
"""
function point_group(
    lattice::Lattice{Tᵢ,T,D,O}; check_edges::Bool=true, tol_digits::Integer=TOL_DIGITS
) where {Tᵢ<:Integer,T<:Real,D,O}
    operations = PointOperation{D,Float64}[]
    for R in _integer_point_candidates(lattice.basis)
        τ = compensating_translation(
            lattice, R; check_edges=check_edges, tol_digits=tol_digits
        )
        τ === nothing && continue
        push!(operations, PointOperation{D,Float64}(R))
    end
    return operations
end

"""
    space_group(lattice; check_edges=true) -> Vector{SpaceOperation}

The space group of `lattice`: every ``\\{R \\mid τ\\}`` that is a symmetry, built as the
point group combined with the translation group.

For a symmorphic lattice this has `length(point_group(lattice)) * length(translation_group(lattice))`
elements. Operations that induce the same site permutation are returned once each — distinct
group elements can coincide on a finite lattice, and counting them twice would misreport the
order of the group actually acting on the sites.
"""
function space_group(
    lattice::Lattice{Tᵢ,T,D,O}; check_edges::Bool=true, tol_digits::Integer=TOL_DIGITS
) where {Tᵢ<:Integer,T<:Real,D,O}
    translations = translation_group(lattice)

    operations = SpaceOperation{D,Float64}[]
    seen = Set{Vector{Int}}()
    for R in _integer_point_candidates(lattice.basis)
        τ₀ = compensating_translation(
            lattice, R; check_edges=check_edges, tol_digits=tol_digits
        )
        τ₀ === nothing && continue

        for t in translations
            τ = τ₀ + lattice.basis.vectors * SVector{D,Float64}(t.displacement)
            op = SpaceOperation{D,Float64}(SMatrix{D,D,Float64}(R), τ)
            perm = try
                site_permutation(lattice, op; tol_digits=tol_digits)
            catch err
                err isa ArgumentError && continue
                rethrow()
            end
            check_edges && !preserves_edges(lattice, perm) && continue
            perm in seen && continue
            push!(seen, perm)
            push!(operations, op)
        end
    end
    return operations
end

# ------------------------------------------------------- convenience permutation generators

"""
    translation_permutation(lattice, axis=1; cells=1) -> Vector{Int}

Site permutation for translating by `cells` primitive cells along `axis`.

This is the vector to feed `SymBasis.Translational`.
"""
function translation_permutation(
    lattice::Lattice{Tᵢ,T,D,O}, axis::Integer=1; cells::Integer=1, tol_digits::Integer=TOL_DIGITS
) where {Tᵢ<:Integer,T<:Real,D,O}
    1 <= axis <= D || throw(ArgumentError("axis $axis is out of range for a $D-D lattice"))
    lattice.periodic[axis] || throw(ArgumentError(
        "axis $axis is not periodic, so it admits no translation symmetry"
    ))
    δ = SVector{D,Int}(ntuple(i -> i == axis ? Int(cells) : 0, D))
    return site_permutation(lattice, Translation(δ); tol_digits=tol_digits)
end

"""
    reflection_permutation(lattice, axis=1) -> Vector{Int}

Site permutation for the reflection that reverses `axis` about the Cartesian origin.

This is the vector to feed `SymBasis.SpatialReflection`. Throws if the reflection is not a
symmetry of the lattice.
"""
function reflection_permutation(
    lattice::Lattice{Tᵢ,T,D,O}, axis::Integer=1; tol_digits::Integer=TOL_DIGITS
) where {Tᵢ<:Integer,T<:Real,D,O}
    1 <= axis <= D || throw(ArgumentError("axis $axis is out of range for a $D-D lattice"))
    R = SMatrix{D,D,Float64}(Diagonal([i == axis ? -1.0 : 1.0 for i in 1:D]))

    # The reflection plane need not pass through the origin: on an open chain the mirror sits
    # at the midpoint, and on a lattice with a basis it is fixed by the site offsets.
    τ = compensating_translation(lattice, R; tol_digits=tol_digits)
    τ === nothing && throw(ArgumentError(
        "reflection about axis $axis is not a symmetry of this lattice"
    ))
    return site_permutation(
        lattice, SpaceOperation{D,Float64}(R, τ); tol_digits=tol_digits
    )
end

"""
    rotation_permutations(lattice) -> Vector{Vector{Int}}

Site permutations for the proper rotations of `lattice` — the space-group elements whose
orthogonal part has determinant `+1` — excluding the identity permutation.

These are the vectors to feed `SymBasis.Rotational`, which requires a non-identity generator.
"""
function rotation_permutations(
    lattice::Lattice{Tᵢ,T,D,O}; tol_digits::Integer=TOL_DIGITS
) where {Tᵢ<:Integer,T<:Real,D,O}
    identity_perm = collect(1:_nv(lattice.metagraph))
    perms = Vector{Int}[]
    for R in _integer_point_candidates(lattice.basis)
        det(R) ≈ 1 || continue
        τ = compensating_translation(lattice, R; tol_digits=tol_digits)
        τ === nothing && continue
        perm = site_permutation(
            lattice, SpaceOperation{D,Float64}(SMatrix{D,D,Float64}(R), τ);
            tol_digits=tol_digits
        )
        perm == identity_perm && continue
        perm in perms || push!(perms, perm)
    end
    return perms
end
