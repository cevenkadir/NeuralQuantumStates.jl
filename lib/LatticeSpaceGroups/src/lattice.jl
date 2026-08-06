using LinearAlgebra: det, norm
using StaticArrays

# Default neighbour order: nearest neighbours only.
const ORDER = 1
# Digits to round distances and fractional coordinates to when deciding whether two are equal.
const TOL_DIGITS = 12
# Absolute tolerance below which a distance counts as zero.
const DIST_TOL = 1.0e-12

"""
    AbstractLatticeBasis{T<:Real,D,O}

Supertype of unit cells: `O` sites in `D` dimensions with coordinates of type `T`.
The concrete implementation is [`LatticeBasis`](@ref).
"""
abstract type AbstractLatticeBasis{T<:Real,D,O} end

"""
    AbstractLattice{T<:Real,D,O}

Supertype of lattices: a finite `D`-dimensional arrangement of cells carrying `O` sites each,
with coordinates of type `T`. The concrete implementation is [`Lattice`](@ref).
"""
abstract type AbstractLattice{T<:Real,D,O} end

# ============================================================================ lattice basis

"""
    LatticeBasis{T<:Real,D,O} <: LatticeSpaceGroups.AbstractLatticeBasis{T,D,O}

The unit cell of a `D`-dimensional lattice carrying `O` sites.

# Fields
- `vectors::SMatrix{D,D,T}`: The primitive vectors, one per **column**. The cell at integer
    coordinates ``(n_1, \\dots, n_D)`` sits at ``\\sum_i \\mathrm{vectors}[:, i] \\, (n_i - 1)``.
- `site_offsets::SMatrix{D,O,T}`: Positions of the `O` sites within the cell, one per column.

# Constructors
    LatticeBasis(vectors)                        # one site at the cell origin
    LatticeBasis(vectors, site_offset::Vector)   # one site at a given offset
    LatticeBasis(vectors, site_offsets::Matrix)  # O sites, one per column

`vectors` may be a `D`×`D` matrix, a vector of `D` primitive vectors, or — in one dimension —
a single real number. The primitive vectors must be linearly independent; a singular set does
not define a lattice and is rejected.

# Example
```julia
LatticeBasis([1.0 0.5; 0.0 sqrt(0.75)])           # triangular Bravais cell
LatticeBasis(1.0, [0.0, 0.5])                     # 1-D cell with a two-site basis
```
"""
struct LatticeBasis{T<:Real,D,O} <: AbstractLatticeBasis{T,D,O}
    vectors::SMatrix{D,D,T}
    site_offsets::SMatrix{D,O,T}

    function LatticeBasis{T,D,O}(
        vectors::AbstractMatrix, site_offsets::AbstractMatrix
    ) where {T<:Real,D,O}
        D > 0 || throw(ArgumentError("dimension must be positive"))
        O > 0 || throw(ArgumentError("a unit cell must contain at least one site"))
        v = SMatrix{D,D,T}(vectors)
        # A singular set of primitive vectors spans fewer than D dimensions, so it describes no
        # D-dimensional lattice -- and `_site_key` could not invert it to get fractional
        # coordinates. Rejecting it here beats a confusing failure later.
        iszero(det(v)) && throw(ArgumentError("the primitive vectors must be independent"))
        return new{T,D,O}(v, SMatrix{D,O,T}(site_offsets))
    end
end

function LatticeBasis(vectors::AbstractMatrix{T}, site_offsets::AbstractMatrix{T}) where {T<:Real}
    size(vectors, 1) == size(vectors, 2) ||
        throw(ArgumentError("expected a square matrix of primitive vectors"))
    D, O = size(vectors, 1), size(site_offsets, 2)
    size(site_offsets, 1) == D || throw(ArgumentError(
        "site offsets are $(size(site_offsets, 1))-dimensional but the cell is $D-dimensional"
    ))
    return LatticeBasis{T,D,O}(vectors, site_offsets)
end
function LatticeBasis(vectors::AbstractMatrix{T}, site_offset::AbstractVector{T}) where {T<:Real}
    return LatticeBasis(vectors, reshape(collect(site_offset), :, 1))
end
function LatticeBasis(vectors::AbstractMatrix{T}) where {T<:Real}
    return LatticeBasis(vectors, zeros(T, size(vectors, 1)))
end

function LatticeBasis(
    vectors::AbstractVector{<:AbstractVector{T}},
    site_offsets::AbstractVector{<:AbstractVector{T}}
) where {T<:Real}
    return LatticeBasis(reduce(hcat, vectors), reduce(hcat, site_offsets))
end
function LatticeBasis(vectors::AbstractVector{<:AbstractVector{T}}) where {T<:Real}
    return LatticeBasis(reduce(hcat, vectors))
end

function LatticeBasis(vector::T, site_offsets::AbstractVector{T}) where {T<:Real}
    return LatticeBasis(fill(vector, 1, 1), reshape(collect(site_offsets), 1, :))
end
LatticeBasis(vector::T, site_offset::T) where {T<:Real} = LatticeBasis(vector, [site_offset])
LatticeBasis(vector::T) where {T<:Real} = LatticeBasis(vector, [zero(T)])

# ============================================================================ site geometry

"""
    _cell_ranges(shape) -> NTuple{D,UnitRange{Int}}

The cell coordinate ranges of a lattice, for iteration in site order.
"""
_cell_ranges(shape::SVector{D,Int}) where {D} = ntuple(d -> 1:shape[d], D)

"""
    _site_positions(shape, basis) -> Vector{SVector{D,T}}

Cartesian positions of every site, in **site order**: the sublattice index varies fastest,
then the first cell coordinate, then the second, and so on.

That ordering is load-bearing. It is the indexing that [`site_permutation`](@ref) permutes and
that a Hilbert space built on this lattice assigns its degrees of freedom to, so it must not
drift.
"""
function _site_positions(
    shape::SVector{D,Int}, basis::LatticeBasis{T,D,O}
) where {T<:Real,D,O}
    positions = Vector{SVector{D,T}}(undef, prod(shape) * O)
    i = 0
    for cell in Iterators.product(_cell_ranges(shape)...)
        origin = basis.vectors * (SVector{D,T}(cell) .- one(T))
        for o in 1:O
            positions[i+=1] = origin + basis.site_offsets[:, o]
        end
    end
    return positions
end

"""
    _supercell_images(shape, basis, periodic) -> Vector{SVector{D,T}}

Translations by whole supercells along the periodic directions, one per combination of
``\\{-1, 0, +1\\}``. Adding each to a separation vector and taking the shortest result is the
minimum-image convention, which is what makes "nearest neighbour" mean the right thing on a
torus.
"""
function _supercell_images(
    shape::SVector{D,Int}, basis::LatticeBasis{T,D,O}, periodic::SVector{D,Bool}
) where {T<:Real,D,O}
    ranges = ntuple(d -> periodic[d] ? (-1:1) : (0:0), D)
    return [
        basis.vectors * (SVector{D,T}(n) .* SVector{D,T}(shape))
        for n in Iterators.product(ranges...)
    ] |> vec
end

"""
    _min_image_distance(Δ, images) -> T

Shortest distance between two sites separated by `Δ`, minimized over supercell `images`.

A single pass over the ``3^p`` images is not always enough for a strongly sheared supercell —
shifting to the best image can expose a better one still — so the pass repeats until no image
improves on the current best.
"""
function _min_image_distance(Δ::SVector{D,T}, images::Vector{SVector{D,T}}) where {T<:Real,D}
    current = Δ
    best = norm(current)
    while true
        improved = false
        for t in images
            candidate = current + t
            d = norm(candidate)
            # A plain `<` would let two tied images swap forever.
            if d < best - eps(T) * max(one(T), best)
                best, current, improved = d, candidate, true
            end
        end
        improved || return best
    end
end

"""
    _neighbour_edges(positions, images, max_order; tol_digits, dist_tol)
        -> (edges, orders)

Edges of the lattice grouped into neighbour shells, out to `max_order`.

Every pair of sites is measured under the minimum-image convention; the distinct distances are
sorted, and the `max_order` smallest define the shells. `orders[k]` is the shell `edges[k]`
belongs to, with `1` the nearest neighbours.

This is `O(N²)` in the number of sites. For the lattices this package is for — exact
diagonalization and variational Monte Carlo, so hundreds of sites, not millions — that is far
cheaper than the spatial index it replaced, whose custom periodic metric had to be evaluated
at every tree comparison anyway.
"""
function _neighbour_edges(
    positions::Vector{SVector{D,T}}, images::Vector{SVector{D,T}}, max_order::Int;
    tol_digits::Int=TOL_DIGITS, dist_tol::T=T(DIST_TOL)
) where {T<:Real,D}
    n = length(positions)
    pairs = Tuple{Int,Int}[]
    distances = T[]
    for i in 1:(n-1), j in (i+1):n
        push!(pairs, (i, j))
        push!(distances, round(
            _min_image_distance(positions[j] - positions[i], images); digits=tol_digits
        ))
    end

    shells = sort!(unique(d for d in distances if d > dist_tol))
    length(shells) >= max_order || throw(ArgumentError(
        "this lattice has only $(length(shells)) distinct neighbour shell(s), " *
        "but max_order = $max_order was requested"
    ))
    kept = shells[1:max_order]

    edges = Tuple{Int,Int}[]
    orders = Int[]
    for (pair, d) in zip(pairs, distances)
        shell = findfirst(==(d), kept)
        shell === nothing && continue
        push!(edges, pair)
        push!(orders, shell)
    end
    return edges, orders
end

# ================================================================================== lattice

"""
    Lattice{T<:Real,D,O} <: LatticeSpaceGroups.AbstractLattice{T,D,O}

A finite `D`-dimensional lattice: `shape` unit cells of `basis`, under `periodic` boundary
conditions, with its sites and bonds enumerated.

Sites are numbered in **site order** — sublattice fastest, then cell coordinate 1, 2, … — and
that numbering is the lattice's interface to the rest of the ecosystem: it is what
[`site_permutation`](@ref) permutes and what [`bonds`](@ref) reports pairs of.

# Fields
- `basis::LatticeBasis{T,D,O}`: The unit cell.
- `shape::SVector{D,Int}`: Number of cells along each dimension.
- `periodic::SVector{D,Bool}`: Boundary condition per dimension.
- `positions::Vector{SVector{D,T}}`: Cartesian position of each site.
- `edges::Vector{Tuple{Int,Int}}`: Bonds as site-index pairs, each listed once with `i < j`,
    in lexicographic order.
- `edge_orders::Vector{Int}`: Neighbour shell of each bond — `1` nearest, `2` next-nearest, …

# Constructors
    Lattice(shape, basis, periodic=false; max_order=1, tol_digits=12, dist_tol=1e-12)
    Lattice(shape, basis, edges, periodic=false; orders=ones(Int, length(edges)))

The first form derives bonds from distances, keeping every shell out to `max_order`. The
second takes them literally as site-index pairs, for connectivity that is not distance-derived.

`periodic` may be a single `Bool` applying to every dimension or one flag per dimension.

Prefer [`build`](@ref) with one of the predefined specs — [`Hypercube`](@ref),
[`Square`](@ref), [`Honeycomb`](@ref), … — over calling these directly.
"""
struct Lattice{T<:Real,D,O} <: AbstractLattice{T,D,O}
    basis::LatticeBasis{T,D,O}
    shape::SVector{D,Int}
    periodic::SVector{D,Bool}
    positions::Vector{SVector{D,T}}
    edges::Vector{Tuple{Int,Int}}
    edge_orders::Vector{Int}
end

"""
    _as_boundary(periodic, shape) -> SVector{D,Bool}

Normalize a boundary-condition argument to one flag per dimension, rejecting a periodic
dimension of extent 1 — wrapping it would make a site its own neighbour.
"""
function _as_boundary(periodic::AbstractVector, shape::SVector{D,Int}) where {D}
    length(periodic) == D || throw(ArgumentError(
        "periodic has $(length(periodic)) entries but the shape is $D-dimensional"
    ))
    p = SVector{D,Bool}(periodic)
    bad = findall(p .& (shape .== 1))
    isempty(bad) || throw(ArgumentError(
        "periodic must be false where shape == 1 (offending dimension(s): $bad)"
    ))
    return p
end
_as_boundary(periodic::Bool, shape::SVector{D,Int}) where {D} =
    _as_boundary(fill(periodic, D), shape)

"""
    _as_shape(shape) -> SVector{D,Int}

Normalize and validate a lattice shape.
"""
function _as_shape(shape::AbstractVector{<:Integer})
    D = length(shape)
    D > 0 || throw(ArgumentError("shape must have at least one entry"))
    all(shape .> 0) || throw(ArgumentError("shape must contain positive integers"))
    return SVector{D,Int}(shape)
end

function Lattice(
    shape::AbstractVector{<:Integer},
    basis::LatticeBasis{T,D,O},
    periodic::Union{Bool,AbstractVector}=false;
    max_order::Integer=ORDER, tol_digits::Integer=TOL_DIGITS, dist_tol::Real=DIST_TOL
) where {T<:Real,D,O}
    s = _as_shape(shape)
    length(s) == D || throw(ArgumentError(
        "shape is $(length(s))-dimensional but the basis is $D-dimensional"
    ))
    max_order > 0 || throw(ArgumentError("max_order must be positive"))

    p = _as_boundary(periodic, s)
    positions = _site_positions(s, basis)
    edges, orders = _neighbour_edges(
        positions, _supercell_images(s, basis, p), Int(max_order);
        tol_digits=Int(tol_digits), dist_tol=T(dist_tol)
    )
    return Lattice{T,D,O}(basis, s, p, positions, edges, orders)
end

function Lattice(
    shape::AbstractVector{<:Integer},
    basis::LatticeBasis{T,D,O},
    edges::AbstractVector{<:Tuple{Integer,Integer}},
    periodic::Union{Bool,AbstractVector}=false;
    orders::AbstractVector{<:Integer}=fill(1, length(edges))
) where {T<:Real,D,O}
    s = _as_shape(shape)
    length(s) == D || throw(ArgumentError(
        "shape is $(length(s))-dimensional but the basis is $D-dimensional"
    ))
    length(orders) == length(edges) ||
        throw(ArgumentError("got $(length(orders)) orders for $(length(edges)) edges"))

    positions = _site_positions(s, basis)
    n = length(positions)
    normalized = map(edges) do (i, j)
        (1 <= i <= n && 1 <= j <= n) ||
            throw(ArgumentError("edge ($i, $j) refers to a site outside 1:$n"))
        i == j && throw(ArgumentError("edge ($i, $j) connects a site to itself"))
        (min(Int(i), Int(j)), max(Int(i), Int(j)))
    end

    # Sorting makes `bonds` independent of how the bonds were arrived at: the grid path and
    # the distance search produce the same set in different sequences, and a lattice's bond
    # list should not depend on which one built it.
    permutation = sortperm(normalized)
    return Lattice{T,D,O}(
        basis, s, _as_boundary(periodic, s), positions,
        collect(normalized)[permutation], collect(Int, orders)[permutation]
    )
end

function Base.show(io::IO, lattice::Lattice{T,D,O}) where {T<:Real,D,O}
    print(io, "Lattice{", T, ",", D, ",", O, "}(", Vector(lattice.shape),
        "; periodic=", Vector(lattice.periodic), ", ",
        n_sites(lattice), " sites, ", length(lattice.edges), " bonds)")
    return nothing
end

# ================================================================================= accessors

"""
    n_sites(lattice) -> Int

Number of sites in the lattice, `prod(shape) * O`.
"""
n_sites(lattice::Lattice) = length(lattice.positions)

"""
    site_positions(lattice) -> Vector{SVector{D,T}}

Cartesian positions of the lattice sites, indexed by site number.
"""
site_positions(lattice::Lattice) = lattice.positions

"""
    site_labels(lattice) -> Vector{NTuple{D+1,Int}}

Lattice site labels `(sublattice, n₁, …, n_D)`, indexed by site number and in the same order
as [`site_positions`](@ref).

The sublattice index runs `1:O` and the cell coordinates run `1:shape[d]`, so the label says
where a site sits without reference to its Cartesian position.
"""
function site_labels(lattice::Lattice{T,D,O}) where {T<:Real,D,O}
    labels = Vector{NTuple{D + 1,Int}}(undef, n_sites(lattice))
    i = 0
    for cell in Iterators.product(_cell_ranges(lattice.shape)...)
        for o in 1:O
            labels[i+=1] = (o, cell...)
        end
    end
    return labels
end

"""
    bonds(lattice; order=nothing) -> Vector{Tuple{Int,Int}}

The lattice's bonds as pairs of site indices, each listed once with `i < j`, in lexicographic
order.

Site indices are the numbering of [`site_positions`](@ref), so a bond can be used directly as
the site identifiers of an operator term.

`order` selects a neighbour shell: `1` for nearest neighbours, `2` for next-nearest, and so on,
matching the `max_order` the lattice was built with. The default, `nothing`, returns every bond.

# Example
```julia
lat = build(Hypercube([4]; periodic=true))
bonds(lat)    # [(1,2), (1,4), (2,3), (3,4)]
```
"""
function bonds(lattice::Lattice; order::Union{Nothing,Integer}=nothing)
    order === nothing && return copy(lattice.edges)
    return [e for (e, o) in zip(lattice.edges, lattice.edge_orders) if o == order]
end
