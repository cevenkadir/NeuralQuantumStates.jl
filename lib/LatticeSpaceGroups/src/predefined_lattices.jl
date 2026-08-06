using Graphs: adjacency_matrix, cartesian_product, cycle_graph, path_graph
using LinearAlgebra: Diagonal, I, triu
using SparseArrays: findnz
using StaticArrays

"""
    AbstractLatticeSpec{D}

A specification of a `D`-dimensional predefined lattice.

A spec is a *description* of a lattice — its shape, size, and boundary conditions — not the
lattice itself. Pass one to [`build`](@ref) to get a [`Lattice`](@ref):

```julia
build(Hypercube([4, 4], 1.0; periodic=true))
```

Concrete specs are [`Hypercube`](@ref), [`Triclinic`](@ref), [`Triangular`](@ref),
[`Honeycomb`](@ref), and [`Kagome`](@ref). Each validates its arguments on construction, so an
invalid lattice is rejected at the point where it is described rather than deep inside the
build.

The type parameter `D` is the spatial dimension, which is fixed for most specs — a honeycomb is
always two-dimensional — and free only for [`Hypercube`](@ref).
"""
abstract type AbstractLatticeSpec{D} end

# ------------------------------------------------------------------------ argument handling

function _as_shape(shape::AbstractVector{<:Integer})
    D = length(shape)
    D > 0 || throw(ArgumentError("shape must have at least one entry"))
    all(shape .> 0) || throw(ArgumentError("shape must contain positive integers"))
    return SVector{D,Int}(shape)
end

"""
    _as_periodic(periodic, shape) -> SVector{D,Bool}

Normalize a boundary-condition argument to one flag per dimension.

A dimension of extent 1 cannot be periodic: wrapping it would make a site its own neighbour.
"""
function _as_periodic(periodic::AbstractVector{Bool}, shape::SVector{D,Int}) where {D}
    length(periodic) == D || throw(ArgumentError(
        "periodic has $(length(periodic)) entries but the shape is $D-dimensional"
    ))
    p = SVector{D,Bool}(periodic)
    if any(p .& (shape .== 1))
        bad = findall(p .& (shape .== 1))
        throw(ArgumentError(
            "periodic must be false where shape == 1 (offending dimension(s): $bad)"
        ))
    end
    return p
end
_as_periodic(periodic::Bool, shape::SVector{D,Int}) where {D} =
    _as_periodic(fill(periodic, D), shape)

# --------------------------------------------------------------------------- the spec types

"""
    Hypercube{D,T<:Real} <: LatticeSpaceGroups.AbstractLatticeSpec{D}

A `D`-dimensional hypercubic lattice: a chain in 1-D, a square lattice in 2-D, a cubic lattice
in 3-D, and so on.

# Fields
- `shape::SVector{D,Int}`: Number of cells along each dimension.
- `edge_length::T`: The lattice spacing.
- `periodic::SVector{D,Bool}`: Boundary condition per dimension.

# Constructor
    Hypercube(shape, edge_length; periodic=false)

`periodic` may be a single `Bool` applying to every dimension, or one flag per dimension.

# Example
```julia
build(Hypercube([8], 1.0; periodic=true))         # periodic chain
build(Hypercube([4, 4], 1.0; periodic=[true, false]))  # cylinder
```
"""
struct Hypercube{D,T<:Real} <: AbstractLatticeSpec{D}
    shape::SVector{D,Int}
    edge_length::T
    periodic::SVector{D,Bool}

    function Hypercube(
        shape::AbstractVector{<:Integer},
        edge_length::T;
        periodic::Union{Bool,AbstractVector{Bool}}=false
    ) where {T<:Real}
        s = _as_shape(shape)
        edge_length > 0 || throw(ArgumentError("edge length must be positive"))
        return new{length(s),T}(s, edge_length, _as_periodic(periodic, s))
    end
end

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
        periodic::Union{Bool,AbstractVector{Bool}}=false
    ) where {T<:Real}
        s = _as_shape(shape)
        length(s) == 3 || throw(ArgumentError("a triclinic lattice must be 3-dimensional"))
        length(edge_lengths) == 3 || throw(ArgumentError("expected three edge lengths"))
        length(angles) == 3 || throw(ArgumentError("expected three angles"))
        all(edge_lengths .> 0) || throw(ArgumentError("edge lengths must be positive"))
        return new{T}(
            s, SVector{3,T}(edge_lengths), SVector{3,T}(angles), _as_periodic(periodic, s)
        )
    end
end

for (name, doc) in (
    (:Triangular, "A two-dimensional triangular lattice, with one site per unit cell and six nearest neighbours."),
    (:Honeycomb, "A two-dimensional honeycomb lattice: a triangular Bravais lattice with a two-site basis, giving three nearest neighbours per site."),
    (:Kagome, "A two-dimensional kagome lattice: a triangular Bravais lattice with a three-site basis, forming corner-sharing triangles."),
)
    @eval begin
        """
            $($name){T<:Real} <: LatticeSpaceGroups.AbstractLatticeSpec{2}

        $($doc)

        # Fields
        - `shape::SVector{2,Int}`: Number of unit cells along each dimension.
        - `edge_length::T`: The lattice spacing.
        - `periodic::SVector{2,Bool}`: Boundary condition per dimension.

        # Constructor
            $($name)(shape, edge_length; periodic=false)
        """
        struct $name{T<:Real} <: AbstractLatticeSpec{2}
            shape::SVector{2,Int}
            edge_length::T
            periodic::SVector{2,Bool}

            function $name(
                shape::AbstractVector{<:Integer},
                edge_length::T;
                periodic::Union{Bool,AbstractVector{Bool}}=false
            ) where {T<:Real}
                s = _as_shape(shape)
                length(s) == 2 ||
                    throw(ArgumentError("a $($name) lattice must be 2-dimensional"))
                edge_length > 0 || throw(ArgumentError("edge length must be positive"))
                return new{T}(s, edge_length, _as_periodic(periodic, s))
            end
        end
    end
end

function Base.show(io::IO, spec::AbstractLatticeSpec{D}) where {D}
    print(io, nameof(typeof(spec)), "(", Vector(spec.shape), "; periodic=",
        Vector(spec.periodic), ")")
    return nothing
end

# --------------------------------------------------------------------------------- building

"""
    _cartesian_edges(shape, periodic) -> (edge_labels, edge_orders)

Nearest-neighbour edges of a `shape`-sized grid, as a Cartesian product of chains — cyclic
along periodic dimensions and open along the rest.

Used by the lattices whose connectivity is a grid rather than something distance-derived, so
that their edge set does not depend on the neighbour-order search at all.
"""
function _cartesian_edges(shape::SVector{D,Int}, periodic::SVector{D,Bool}) where {D}
    g = periodic[1] ? cycle_graph(shape[1]) : path_graph(shape[1])
    for (index, s) in enumerate(shape[2:end])
        g = cartesian_product(periodic[index+1] ? cycle_graph(s) : path_graph(s), g)
    end

    pos_labels = collect(Iterators.product([1:i for i in shape]...))[:]
    i_s, j_s, = findnz(triu(adjacency_matrix(g)))
    edge_labels = map((i, j) -> (pos_labels[i], pos_labels[j]), i_s, j_s)
    return (edge_labels, fill(1, length(edge_labels)))
end

"""
    build(spec::AbstractLatticeSpec; kwargs...) -> Lattice

Build the lattice described by `spec`.

# Keywords
- `max_order::Integer=1`: How many neighbour shells to include as edges. `1` means nearest
    neighbours only, `2` adds next-nearest, and so on. Ignored by [`Hypercube`](@ref) and
    [`Triclinic`](@ref), whose connectivity is a grid product rather than distance-derived.
- `tol_digits::Integer=12`: Digits to round computed distances to when grouping neighbour
    shells.
- `dist_tol::Real=1.0e-12`: Tolerance for treating two distances as equal.

# Example
```julia
lat = build(Hypercube([4, 4], 1.0; periodic=true))
lat = build(Kagome([3, 3], 1.0; periodic=true))
```
"""
function build end

function build(spec::Hypercube{D,T}; kwargs...) where {D,T<:Real}
    basis = LatticeBasis(SMatrix{D,D,T}(spec.edge_length * I(D)))
    edges = _cartesian_edges(spec.shape, spec.periodic)
    return Lattice(spec.shape, basis, edges, spec.periodic)
end

function build(spec::Triclinic{T}; kwargs...) where {T<:Real}
    α, β, γ = spec.angles
    vectors = Diagonal(spec.edge_lengths) * SMatrix{3,3,T}(
        [
            1.0 0.0 0.0
            cosd(γ) sind(γ) 0
            cosd(β) (cosd(α)-cosd(β)*cosd(γ))/sind(γ) sqrt(1 - (cosd(α)^2 + cosd(β)^2) / sind(γ)^2)
        ]
    )
    basis = LatticeBasis(SMatrix{3,3,T}(transpose(vectors)))
    edges = _cartesian_edges(spec.shape, spec.periodic)
    return Lattice(spec.shape, basis, edges, spec.periodic)
end

"""Primitive vectors shared by the triangular, honeycomb, and kagome lattices."""
_triangular_vectors(edge_length::T) where {T<:Real} =
    edge_length * SMatrix{2,2,T}([1.0 0.5; 0.0 sqrt(0.75)])

function build(
    spec::Triangular{T};
    max_order::Integer=ORDER, tol_digits::Integer=TOL_DIGITS, dist_tol::Real=DIST_TOL
) where {T<:Real}
    basis = LatticeBasis(_triangular_vectors(spec.edge_length))
    return Lattice(
        spec.shape, basis, spec.periodic;
        max_order=Int(max_order), tol_digits=Int(tol_digits), dist_tol=T(dist_tol)
    )
end

function build(
    spec::Honeycomb{T};
    max_order::Integer=ORDER, tol_digits::Integer=TOL_DIGITS, dist_tol::Real=DIST_TOL
) where {T<:Real}
    vectors = _triangular_vectors(spec.edge_length)
    # The two sublattices sit at fractional (1/3, 1/3) and (2/3, 2/3) of the cell, which puts
    # them a distance `edge_length / sqrt(3)` apart -- the honeycomb bond length.
    site_offsets = hcat(vectors * SVector{2,T}(1 / 3, 1 / 3), vectors * SVector{2,T}(2 / 3, 2 / 3))
    basis = LatticeBasis(vectors, SMatrix{2,2,T}(site_offsets))
    return Lattice(
        spec.shape, basis, spec.periodic;
        max_order=Int(max_order), tol_digits=Int(tol_digits), dist_tol=T(dist_tol)
    )
end

function build(
    spec::Kagome{T};
    max_order::Integer=ORDER, tol_digits::Integer=TOL_DIGITS, dist_tol::Real=DIST_TOL
) where {T<:Real}
    vectors = _triangular_vectors(spec.edge_length)
    # Sites at a cell corner and at the midpoints of the two primitive vectors, which is what
    # makes the corner-sharing triangles.
    site_offsets = SMatrix{2,3,T}([[0.0; 0.0] vectors ./ 2.0])
    basis = LatticeBasis(vectors, site_offsets)
    return Lattice(
        spec.shape, basis, spec.periodic;
        max_order=Int(max_order), tol_digits=Int(tol_digits), dist_tol=T(dist_tol)
    )
end
