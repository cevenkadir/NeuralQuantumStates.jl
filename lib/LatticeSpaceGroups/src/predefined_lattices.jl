using LinearAlgebra: Diagonal, I
using StaticArrays

"""
    AbstractLatticeSpec{D}

A specification of a `D`-dimensional predefined lattice.

A spec is a *description* of a lattice — its shape, size, and boundary conditions — not the
lattice itself. Pass one to [`build`](@ref) to get a [`Lattice`](@ref):

```julia
build(Square(4; periodic=true))
```

Each spec validates its arguments on construction, so an invalid lattice is rejected at the
point where it is described rather than deep inside the build.

# Available specs

| Spec | Dim | Sites/cell | Coordination |
|---|---|---|---|
| [`Hypercube`](@ref) (and [`Square`](@ref), [`Cube`](@ref)) | any | 1 | `2D` |
| [`Triclinic`](@ref) | 3 | 1 | 6 |
| [`Triangular`](@ref) | 2 | 1 | 6 |
| [`Honeycomb`](@ref) | 2 | 2 | 3 |
| [`Kagome`](@ref) | 2 | 3 | 4 |
| [`BCC`](@ref) | 3 | 1 | 8 |
| [`FCC`](@ref) | 3 | 1 | 12 |
| [`Diamond`](@ref) | 3 | 2 | 4 |
| [`Pyrochlore`](@ref) | 3 | 4 | 6 |

The type parameter `D` is the spatial dimension, which is fixed for most specs — a honeycomb is
always two-dimensional — and free only for [`Hypercube`](@ref).
"""
abstract type AbstractLatticeSpec{D} end

function Base.show(io::IO, spec::AbstractLatticeSpec)
    print(io, nameof(typeof(spec)), "(", Vector(spec.shape), "; periodic=",
        Vector(spec.periodic), ")")
    return nothing
end

# ------------------------------------------------------------------------- generic 1-cell specs

"""
    @simple_spec name dim doc

Define a lattice spec that is fully described by a shape, one edge length, and boundary
conditions — which is all of them except [`Triclinic`](@ref).
"""
macro simple_spec(name, dim, doc)
    quote
        @doc $doc struct $(esc(name)){T<:Real} <: AbstractLatticeSpec{$dim}
            shape::SVector{$dim,Int}
            edge_length::T
            periodic::SVector{$dim,Bool}

            function $(esc(name))(
                shape::AbstractVector{<:Integer},
                edge_length::T=1.0;
                periodic::Union{Bool,AbstractVector}=false
            ) where {T<:Real}
                s = _as_shape(shape)
                length(s) == $dim || throw(ArgumentError(
                    string($(QuoteNode(name)), " lattices are ", $dim, "-dimensional, but ",
                        "the given shape is ", length(s), "-dimensional")
                ))
                edge_length > 0 || throw(ArgumentError("edge length must be positive"))
                new{T}(s, edge_length, _as_boundary(periodic, s))
            end
        end
    end
end

"""
    Hypercube{D,T<:Real} <: LatticeSpaceGroups.AbstractLatticeSpec{D}

A `D`-dimensional hypercubic lattice: a chain in 1-D, a square lattice in 2-D, a cubic lattice
in 3-D, and so on. One site per cell, coordination `2D` under periodic boundaries.

# Fields
- `shape::SVector{D,Int}`: Number of cells along each dimension.
- `edge_length::T`: The lattice spacing, which is also the nearest-neighbour distance.
- `periodic::SVector{D,Bool}`: Boundary condition per dimension.

# Constructor
    Hypercube(shape, edge_length=1.0; periodic=false)

`periodic` may be a single `Bool` applying to every dimension, or one flag per dimension.
See [`Square`](@ref) and [`Cube`](@ref) for those two cases.

!!! note "Why there is no `Chain` shorthand"
    A one-dimensional chain is `Hypercube([n])`. There is deliberately no exported `Chain`,
    because the name collides with `Lux.Chain` and `Flux.Chain` — and this package is meant to
    be loaded alongside them.

# Example
```julia
build(Hypercube([8], 1.0; periodic=true))              # periodic chain
build(Hypercube([4, 4], 1.0; periodic=[true, false]))  # cylinder
```
"""
struct Hypercube{D,T<:Real} <: AbstractLatticeSpec{D}
    shape::SVector{D,Int}
    edge_length::T
    periodic::SVector{D,Bool}

    function Hypercube(
        shape::AbstractVector{<:Integer},
        edge_length::T=1.0;
        periodic::Union{Bool,AbstractVector}=false
    ) where {T<:Real}
        s = _as_shape(shape)
        edge_length > 0 || throw(ArgumentError("edge length must be positive"))
        return new{length(s),T}(s, edge_length, _as_boundary(periodic, s))
    end
end

"""
    Square(n, edge_length=1.0; periodic=false) -> Hypercube{2}
    Square((n₁, n₂), edge_length=1.0; periodic=false) -> Hypercube{2}

A two-dimensional square lattice, `n`×`n` cells or `n₁`×`n₂` if given a pair.
"""
Square(n::Integer, edge_length::Real=1.0; kwargs...) =
    Hypercube([n, n], edge_length; kwargs...)
Square(shape::AbstractVector{<:Integer}, edge_length::Real=1.0; kwargs...) =
    Hypercube(_require_dim(shape, 2, "Square"), edge_length; kwargs...)

"""
    Cube(n, edge_length=1.0; periodic=false) -> Hypercube{3}
    Cube((n₁, n₂, n₃), edge_length=1.0; periodic=false) -> Hypercube{3}

A three-dimensional cubic lattice, `n`×`n`×`n` cells or `n₁`×`n₂`×`n₃` if given a triple.
"""
Cube(n::Integer, edge_length::Real=1.0; kwargs...) =
    Hypercube([n, n, n], edge_length; kwargs...)
Cube(shape::AbstractVector{<:Integer}, edge_length::Real=1.0; kwargs...) =
    Hypercube(_require_dim(shape, 3, "Cube"), edge_length; kwargs...)

function _require_dim(shape::AbstractVector{<:Integer}, D::Int, name::AbstractString)
    length(shape) == D || throw(ArgumentError(
        "$name lattices are $D-dimensional, but the given shape is " *
        "$(length(shape))-dimensional"
    ))
    return shape
end

@simple_spec Triangular 2 """
    Triangular{T<:Real} <: LatticeSpaceGroups.AbstractLatticeSpec{2}

A two-dimensional triangular lattice: one site per cell, six nearest neighbours.

# Fields
- `shape::SVector{2,Int}`: Number of unit cells along each dimension.
- `edge_length::T`: Length of the primitive vectors, which is also the bond length here.
- `periodic::SVector{2,Bool}`: Boundary condition per dimension.

# Constructor
    Triangular(shape, edge_length=1.0; periodic=false)
"""

@simple_spec Honeycomb 2 """
    Honeycomb{T<:Real} <: LatticeSpaceGroups.AbstractLatticeSpec{2}

A two-dimensional honeycomb lattice: a triangular Bravais lattice with a two-site basis,
giving three nearest neighbours per site — the graphene lattice.

# Fields
- `shape::SVector{2,Int}`: Number of unit cells along each dimension.
- `edge_length::T`: Length of the primitive vectors. The **bond** length is
    `edge_length / √3`, since the two sublattices sit at thirds of the cell.
- `periodic::SVector{2,Bool}`: Boundary condition per dimension.

# Constructor
    Honeycomb(shape, edge_length=1.0; periodic=false)
"""

@simple_spec Kagome 2 """
    Kagome{T<:Real} <: LatticeSpaceGroups.AbstractLatticeSpec{2}

A two-dimensional kagome lattice: a triangular Bravais lattice with a three-site basis forming
corner-sharing triangles, giving four nearest neighbours per site.

# Fields
- `shape::SVector{2,Int}`: Number of unit cells along each dimension.
- `edge_length::T`: Length of the primitive vectors. The **bond** length is `edge_length / 2`.
- `periodic::SVector{2,Bool}`: Boundary condition per dimension.

# Constructor
    Kagome(shape, edge_length=1.0; periodic=false)
"""

@simple_spec BCC 3 """
    BCC{T<:Real} <: LatticeSpaceGroups.AbstractLatticeSpec{3}

A body-centred cubic lattice: one site per primitive cell, eight nearest neighbours.

Built in the **primitive** basis, so `shape` counts primitive cells and the lattice has
`prod(shape)` sites — not the two-site conventional cubic cell.

# Fields
- `shape::SVector{3,Int}`: Number of primitive cells along each dimension.
- `edge_length::T`: The **conventional** cubic cell parameter `a`. The nearest-neighbour
    distance is `√3 a / 2`.
- `periodic::SVector{3,Bool}`: Boundary condition per dimension.

# Constructor
    BCC(shape, edge_length=1.0; periodic=false)
"""

@simple_spec FCC 3 """
    FCC{T<:Real} <: LatticeSpaceGroups.AbstractLatticeSpec{3}

A face-centred cubic lattice: one site per primitive cell, twelve nearest neighbours — the
densest packing, and the Bravais lattice underlying both [`Diamond`](@ref) and
[`Pyrochlore`](@ref).

Built in the **primitive** basis, so `shape` counts primitive cells and the lattice has
`prod(shape)` sites — not the four-site conventional cubic cell.

# Fields
- `shape::SVector{3,Int}`: Number of primitive cells along each dimension.
- `edge_length::T`: The **conventional** cubic cell parameter `a`. The nearest-neighbour
    distance is `a / √2`.
- `periodic::SVector{3,Bool}`: Boundary condition per dimension.

# Constructor
    FCC(shape, edge_length=1.0; periodic=false)
"""

@simple_spec Diamond 3 """
    Diamond{T<:Real} <: LatticeSpaceGroups.AbstractLatticeSpec{3}

The diamond lattice: an FCC Bravais lattice with a two-site basis, giving four nearest
neighbours in a tetrahedral arrangement. Silicon, germanium, and diamond itself.

# Fields
- `shape::SVector{3,Int}`: Number of primitive cells along each dimension; the lattice has
    `2 * prod(shape)` sites.
- `edge_length::T`: The **conventional** cubic cell parameter `a`. The bond length is
    `√3 a / 4`.
- `periodic::SVector{3,Bool}`: Boundary condition per dimension.

# Constructor
    Diamond(shape, edge_length=1.0; periodic=false)
"""

@simple_spec Pyrochlore 3 """
    Pyrochlore{T<:Real} <: LatticeSpaceGroups.AbstractLatticeSpec{3}

The pyrochlore lattice: an FCC Bravais lattice with a four-site basis forming corner-sharing
tetrahedra, giving six nearest neighbours. The canonical three-dimensional geometrically
frustrated lattice, and the 3-D analogue of [`Kagome`](@ref).

# Fields
- `shape::SVector{3,Int}`: Number of primitive cells along each dimension; the lattice has
    `4 * prod(shape)` sites.
- `edge_length::T`: The **conventional** cubic cell parameter `a`. The bond length is
    `a / (2√2)`.
- `periodic::SVector{3,Bool}`: Boundary condition per dimension.

# Constructor
    Pyrochlore(shape, edge_length=1.0; periodic=false)
"""

"""
    Triclinic{T<:Real} <: LatticeSpaceGroups.AbstractLatticeSpec{3}

A three-dimensional triclinic lattice, the least symmetric Bravais lattice: three independent
edge lengths at three independent angles.

# Fields
- `shape::SVector{3,Int}`: Number of cells along each dimension.
- `edge_lengths::SVector{3,T}`: The three edge lengths.
- `angles::SVector{3,T}`: The three angles in **degrees**, where `angles[i]` is the angle
    between `edge_lengths[j]` and `edge_lengths[k]` for `(i,j,k)` a permutation of `(1,2,3)`.
- `periodic::SVector{3,Bool}`: Boundary condition per dimension.

# Constructor
    Triclinic(shape, edge_lengths, angles; periodic=false)
"""
struct Triclinic{T<:Real} <: AbstractLatticeSpec{3}
    shape::SVector{3,Int}
    edge_lengths::SVector{3,T}
    angles::SVector{3,T}
    periodic::SVector{3,Bool}

    function Triclinic(
        shape::AbstractVector{<:Integer},
        edge_lengths::AbstractVector{T},
        angles::AbstractVector{T};
        periodic::Union{Bool,AbstractVector}=false
    ) where {T<:Real}
        s = _as_shape(shape)
        length(s) == 3 || throw(ArgumentError("a triclinic lattice must be 3-dimensional"))
        length(edge_lengths) == 3 || throw(ArgumentError("expected three edge lengths"))
        length(angles) == 3 || throw(ArgumentError("expected three angles"))
        all(edge_lengths .> 0) || throw(ArgumentError("edge lengths must be positive"))
        return new{T}(
            s, SVector{3,T}(edge_lengths), SVector{3,T}(angles), _as_boundary(periodic, s)
        )
    end
end

# --------------------------------------------------------------------------------- building

"""
    _cartesian_edges(shape, periodic) -> Vector{Tuple{Int,Int}}

Nearest-neighbour bonds of a `shape`-sized grid of one-site cells: each cell joined to its
successor along every dimension, wrapping where periodic.

Used by the lattices whose connectivity *is* the grid, so their bond set does not depend on
the distance search at all — exact by construction, and `O(N)` rather than `O(N²)`.
"""
function _cartesian_edges(shape::SVector{D,Int}, periodic::SVector{D,Bool}) where {D}
    strides = ntuple(d -> d == 1 ? 1 : prod(shape[k] for k in 1:(d-1)), D)
    linear(cell) = 1 + sum(ntuple(d -> (cell[d] - 1) * strides[d], D))

    edges = Tuple{Int,Int}[]
    for cell in Iterators.product(_cell_ranges(shape)...)
        i = linear(cell)
        for d in 1:D
            j = if cell[d] < shape[d]
                linear(ntuple(k -> k == d ? cell[k] + 1 : cell[k], D))
            elseif periodic[d] && shape[d] > 2
                # At extent 2 the forward bond and the wrap-around bond are the same pair;
                # adding both would double-count it.
                linear(ntuple(k -> k == d ? 1 : cell[k], D))
            else
                continue
            end
            push!(edges, (min(i, j), max(i, j)))
        end
    end
    return edges
end

"""
    build(spec::AbstractLatticeSpec; kwargs...) -> Lattice

Build the lattice described by `spec`.

# Keywords
- `max_order::Integer=1`: How many neighbour shells to include as bonds. `1` means nearest
    neighbours only, `2` adds next-nearest, and so on.
- `tol_digits::Integer=12`: Digits to round computed distances to when grouping neighbour
    shells.
- `dist_tol::Real=1.0e-12`: Tolerance below which a distance counts as zero.

# Example
```julia
lat = build(Square(4; periodic=true))
lat = build(Kagome([3, 3], 1.0; periodic=true))
lat = build(Pyrochlore([2, 2, 2], 1.0; periodic=true))
```
"""
function build end

"""
    _build_grid(spec, basis; max_order=1, kwargs...) -> Lattice

Build a lattice whose nearest-neighbour bonds *are* the cell grid.

At `max_order == 1` the bonds follow from the grid directly, which is both exact and `O(N)` —
no distances involved. Beyond that the shells stop being grid steps (the second shell of a
square lattice is its diagonals), so the general distance search takes over.
"""
function _build_grid(
    spec::AbstractLatticeSpec, basis::LatticeBasis{T};
    max_order::Integer=ORDER, kwargs...
) where {T<:Real}
    max_order == 1 || return _build_from_basis(spec, basis; max_order=max_order, kwargs...)
    return Lattice(
        spec.shape, basis, _cartesian_edges(spec.shape, spec.periodic), spec.periodic
    )
end

build(spec::Hypercube{D,T}; kwargs...) where {D,T<:Real} =
    _build_grid(spec, LatticeBasis(SMatrix{D,D,T}(spec.edge_length * I(D))); kwargs...)

function build(spec::Triclinic{T}; kwargs...) where {T<:Real}
    α, β, γ = spec.angles
    vectors = Diagonal(spec.edge_lengths) * SMatrix{3,3,T}(
        [
            1.0 0.0 0.0
            cosd(γ) sind(γ) 0.0
            cosd(β) (cosd(α)-cosd(β)*cosd(γ))/sind(γ) sqrt(1 - (cosd(α)^2 + cosd(β)^2) / sind(γ)^2)
        ]
    )
    return _build_grid(spec, LatticeBasis(SMatrix{3,3,T}(transpose(vectors))); kwargs...)
end

"""Primitive vectors shared by the triangular, honeycomb, and kagome lattices."""
_triangular_vectors(a::T) where {T<:Real} = a * SMatrix{2,2,T}([1.0 0.5; 0.0 sqrt(0.75)])

"""
Primitive vectors of the face-centred cubic lattice, in terms of the conventional cubic cell
parameter `a`. Shared by [`FCC`](@ref), [`Diamond`](@ref), and [`Pyrochlore`](@ref), which
differ only in what sits inside the cell.
"""
_fcc_vectors(a::T) where {T<:Real} =
    (a / 2) * SMatrix{3,3,T}([0.0 1.0 1.0; 1.0 0.0 1.0; 1.0 1.0 0.0])

"""Primitive vectors of the body-centred cubic lattice."""
_bcc_vectors(a::T) where {T<:Real} =
    (a / 2) * SMatrix{3,3,T}([-1.0 1.0 1.0; 1.0 -1.0 1.0; 1.0 1.0 -1.0])

"""
    _build_from_basis(spec, basis; max_order, tol_digits, dist_tol) -> Lattice

Shared tail of every distance-derived `build` method.
"""
function _build_from_basis(
    spec::AbstractLatticeSpec, basis::LatticeBasis{T};
    max_order::Integer=ORDER, tol_digits::Integer=TOL_DIGITS, dist_tol::Real=DIST_TOL
) where {T<:Real}
    return Lattice(
        spec.shape, basis, spec.periodic;
        max_order=max_order, tol_digits=tol_digits, dist_tol=T(dist_tol)
    )
end

build(spec::Triangular{T}; kwargs...) where {T<:Real} =
    _build_from_basis(spec, LatticeBasis(_triangular_vectors(spec.edge_length)); kwargs...)

function build(spec::Honeycomb{T}; kwargs...) where {T<:Real}
    vectors = _triangular_vectors(spec.edge_length)
    # The two sublattices sit at fractional (1/3, 1/3) and (2/3, 2/3) of the cell, which puts
    # them a distance `edge_length / sqrt(3)` apart -- the honeycomb bond length. Deriving the
    # offsets from the primitive vectors is what makes them scale with `edge_length`.
    offsets = hcat(
        vectors * SVector{2,T}(1 / 3, 1 / 3), vectors * SVector{2,T}(2 / 3, 2 / 3)
    )
    return _build_from_basis(spec, LatticeBasis(vectors, SMatrix{2,2,T}(offsets)); kwargs...)
end

function build(spec::Kagome{T}; kwargs...) where {T<:Real}
    vectors = _triangular_vectors(spec.edge_length)
    # A cell corner plus the midpoints of the two primitive vectors, which is what makes the
    # corner-sharing triangles.
    offsets = SMatrix{2,3,T}([[0.0; 0.0] vectors ./ 2])
    return _build_from_basis(spec, LatticeBasis(vectors, offsets); kwargs...)
end

build(spec::BCC{T}; kwargs...) where {T<:Real} =
    _build_from_basis(spec, LatticeBasis(_bcc_vectors(spec.edge_length)); kwargs...)

build(spec::FCC{T}; kwargs...) where {T<:Real} =
    _build_from_basis(spec, LatticeBasis(_fcc_vectors(spec.edge_length)); kwargs...)

function build(spec::Diamond{T}; kwargs...) where {T<:Real}
    vectors = _fcc_vectors(spec.edge_length)
    # The second sublattice sits at a quarter of the conventional cube's body diagonal, which
    # is what gives each site four tetrahedrally arranged neighbours.
    offsets = SMatrix{3,2,T}([
        [0.0; 0.0; 0.0] fill(spec.edge_length / 4, 3)
    ])
    return _build_from_basis(spec, LatticeBasis(vectors, offsets); kwargs...)
end

function build(spec::Pyrochlore{T}; kwargs...) where {T<:Real}
    vectors = _fcc_vectors(spec.edge_length)
    # A cell corner plus the midpoints of the three primitive vectors. Those four sites form a
    # tetrahedron, and each shares its corners with the neighbouring cells' tetrahedra.
    offsets = SMatrix{3,4,T}([[0.0; 0.0; 0.0] vectors ./ 2])
    return _build_from_basis(spec, LatticeBasis(vectors, offsets); kwargs...)
end
