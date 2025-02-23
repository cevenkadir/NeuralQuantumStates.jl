import Distances: evaluate

using Distances: Metric
using Graphs: Graph, nv as _nv
using MetaGraphsNext: MetaGraph
using NearestNeighbors: BallTree, inrange, knn
using Parameters
using StaticArrays

# define constants
const ORDER = 1
const TOL_DIGITS = 12
const DIST_TOL = 1.0e-12

# define abstract types
abstract type AbstractLatticeBasis{T<:Real,D,O} end
abstract type AbstractLattice{Tᵢ<:Integer,T<:Real,D,O} end

"""
    LatticeBasis{T<:Real,D,O} <: NeuralQuantumStates.Lattices.AbstractLatticeBasis{T,D,O}

A lattice basis for representing the unit cell of a `D`-dimensional lattice with `O` site
    offsets.

# Fields
- `vectors::SMatrix{D,D,T}`: A `D```\\times```D` square matrix for the primitive vectors
    defining the unit cell.
- `site_offsets::SMatrix{D,O,T}`: A `D```\\times```O` matrix for the site offsets of the
    lattice basis in the unit cell.
"""
@with_kw struct LatticeBasis{T<:Real,D,O} <: AbstractLatticeBasis{T,D,O}
    vectors::SMatrix{D,D,T}
    site_offsets::SMatrix{D,O,T}
    @assert D isa Integer "dimension must be an integer"
    @assert O isa Integer "number of site offsets must be an integer"
    @assert D > 0 "dimension must be positive"
    @assert O > 0 "number of site offsets must be positive"
end
"""
    LatticeBasis(
        vector::T, site_offsets::AbstractVector{T}
    ) where {T<:Real} -> NeuralQuantumStates.Lattices.LatticeBasis{T,1,length(site_offsets)}

Define a lattice basis for representing the unit cell of a 1D lattice with
    `length(site_offsets)` site offsets.

# Arguments
- `vector::T`: A real number for a primitive vector defining the unit cell.
- `site_offsets::AbstractVector{T}`: A vector for the site offsets of the lattice basis in
    the unit cell.

# Returns
- `NeuralQuantumStates.Lattices.LatticeBasis{T,1,length(site_offsets)}`: The defined lattice
    basis of a 1D lattice with `length(site_offsets)` site offsets.
"""
function LatticeBasis(vector::T, site_offsets::AbstractVector{T}) where {T<:Real}
    site_offsets = unique(collect(site_offsets))

    n_offsets = length(site_offsets)

    vector = SMatrix{1,1,T}(vector)
    site_offsets = SMatrix{1,n_offsets,T}(site_offsets)

    return LatticeBasis{T,1,n_offsets}(vector, site_offsets)
end
"""
    LatticeBasis(vector::T, site_offset::T=T(0.0)) where {T<:Real}
        -> NeuralQuantumStates.Lattices.LatticeBasis{T,1,1}

Define a `NeuralQuantumStates.Lattices.LatticeBasis` for representing the unit cell of a 1D
    lattice with the given primitive `vector` and `site_offset`.

# Arguments
- `vector::T`: A real number for a primitive vector defining the unit cell.
- `site_offset::T`: A real number for the site offset of the lattice basis in the unit cell.
    Defaults to `T(0.0)`.

# Returns
- `NeuralQuantumStates.Lattices.LatticeBasis{T,1,1}`: The defined lattice basis of a 1D
    lattice with the given primitive `vector` and `site_offset`.
"""
function LatticeBasis(vector::T, site_offset::T) where {T<:Real}
    return LatticeBasis(vector, [site_offset])
end
LatticeBasis(vector::T) where {T<:Real} = LatticeBasis(vector, T(0.0))
"""
    LatticeBasis(
        vectors::AbstractMatrix{T},
        site_offset::AbstractVector{T}=zeros(T, size(vectors)[1])
    ) where {T<:Real} -> NeuralQuantumStates.Lattices.LatticeBasis{T,size(vectors)[1],1}

Define a `NeuralQuantumStates.Lattices.LatticeBasis` for representing the unit cell of a
    `size(vectors)[1]`-dimensional lattice with the given primitive `vectors` and
    `site_offset`.

# Arguments
- `vectors::AbstractMatrix{T}`: A square matrix for the primitive vectors defining the unit
    cell. Each column of the matrix is a primitive vector.
- `site_offset::AbstractVector{T}`: A vector for one site offset of the lattice basis in the
    unit cell. Defaults to `zeros(T, size(vectors)[1])`.

# Returns
- `NeuralQuantumStates.Lattices.LatticeBasis{T,size(vectors)[1],1}`: The defined
    `NeuralQuantumStates.Lattices.LatticeBasis` for representing the unit cell of a
    `size(vectors)[1]`-dimensional lattice with the given primitive `vectors` and
    `site_offset`.
"""
function LatticeBasis(
    vectors::AbstractMatrix{T}, site_offset::AbstractVector{T}
) where {T<:Real}
    vectors = unique(collect(vectors), dims=2)

    n_dims, = size(vectors)

    site_offset = SMatrix{n_dims,1,T}(site_offset)

    return LatticeBasis{T,n_dims,1}(vectors, site_offset)
end
"""
    LatticeBasis(
        vectors::AbstractMatrix{T},
        site_offsets::AbstractMatrix{T}=zeros(T, size(vectors)[1], 1)
    ) where {T<:Real}
        -> NeuralQuantumStates.Lattices.LatticeBasis{T,size(vectors)[1],size(site_offsets)[2]}

Define a `NeuralQuantumStates.Lattices.LatticeBasis` for representing the unit cell of a
    `size(vectors)[1]`-dimensional lattice with the given primitive `vectors` and
    `site_offsets`.

# Arguments
- `vectors::AbstractMatrix{T}`: A square matrix for the primitive vectors defining the unit
    cell. Each column of the matrix is a primitive vector.
- `site_offsets::AbstractMatrix{T}`: A matrix for the site offsets of the lattice basis in
    the unit cell. Each column of the matrix is a site offset vector. Defaults to
    `zeros(T, size(vectors)[1], 1)`.

# Returns
- `NeuralQuantumStates.Lattices.LatticeBasis{T,size(vectors)[1],size(site_offsets)[2]}`: The
    defined `NeuralQuantumStates.Lattices.LatticeBasis` for representing the unit cell of a
    `size(vectors)[1]`-dimensional lattice with the given primitive `vectors` and
    `site_offsets`.
"""
function LatticeBasis(
    vectors::AbstractMatrix{T}, site_offsets::AbstractMatrix{T}
) where {T<:Real}
    vectors = unique(collect(vectors), dims=2)
    site_offsets = unique(collect(site_offsets), dims=2)

    n_dims, = size(vectors)

    n_offsets = size(site_offsets)[2]

    return LatticeBasis{T,n_dims,n_offsets}(vectors, site_offsets)
end
function LatticeBasis(vectors::AbstractMatrix{T}) where {T<:Real}
    return LatticeBasis(vectors, zeros(T, size(vectors)[1]))
end
"""
    LatticeBasis(
        vectors::AbstractVector{<:AbstractVector{T}},
        site_offsets::AbstractVector{<:AbstractVector{T}}=[zeros(T, length(vectors[1]))]
    ) where {T<:Real}
        -> NeuralQuantumStates.Lattices.LatticeBasis{T,length(vectors[1]),length(site_offsets)}

Define a `NeuralQuantumStates.Lattices.LatticeBasis` for representing the unit cell of a
    `length(vectors[1])`-dimensional lattice with the given primitive `vectors` and
    `site_offsets`.

# Arguments
- `vectors::AbstractVector{<:AbstractVector{T}}`: A vector of primitive vectors defining the
    unit cell. Each element of the vector should be a primitive vector with the same
    dimension.
- `site_offsets::AbstractVector{<:AbstractVector{T}}`: A vector of site offsets of the
    lattice basis in the unit cell. Each element of the vector should be a site offset
    vector with the same dimension. Defaults to `[zeros(T, length(vectors[1]))]`.

# Returns
- `NeuralQuantumStates.Lattices.LatticeBasis{T,length(vectors[1]),length(site_offsets)}`:
    The defined `NeuralQuantumStates.Lattices.LatticeBasis` for representing the unit cell
    of a `length(vectors[1])`-dimensional lattice with the given primitive `vectors` and
    `site_offsets`.
"""
function LatticeBasis(
    vectors::AbstractVector{<:AbstractVector{T}},
    site_offsets::AbstractVector{<:AbstractVector{T}}
) where {T<:Real}
    return LatticeBasis(reduce(hcat, vectors), reduce(hcat, site_offsets))
end
function LatticeBasis(vectors::AbstractVector{<:AbstractVector{T}}) where {T<:Real}
    return LatticeBasis(reduce(hcat, vectors))
end

function vertices(
    shape::SVector{D,Tᵢ}, basis::AbstractLatticeBasis{T,D,O}
) where {Tᵢ<:Integer,T<:Real,D,O}
    @assert D > 0 "dimension must be positive"
    @assert O > 0 "number of site offsets must be positive"

    n_vertices = Tᵢ(prod(shape) * O)

    site_symbols = Symbol.('A':('A'+O-1))
    site_offsets_dict = Dict(site_symbols[i] => basis.site_offsets[:, i] for i in 1:O)

    labels_iter = Iterators.product(site_symbols, [1:i for i in shape]...)

    _position = label -> site_offsets_dict[label[1]] .+
                         sum(basis.vectors' .* (label[2:end] .- T(1.0)), dims=1)[1, :]

    data_iter = Iterators.map(_position, labels_iter)

    labels = SVector{n_vertices,eltype(labels_iter)}(collect(labels_iter)[:])
    label_data = SVector{n_vertices,SVector{D,T}}(collect(data_iter)[:])

    return labels, label_data
end
function vertices(
    shape::SVector{D,Tᵢ}, basis::AbstractLatticeBasis{T,D,1}
) where {Tᵢ<:Integer,T<:Real,D}
    @assert D > 0 "dimension must be positive"

    n_vertices = prod(shape)

    labels_iter = Iterators.product([1:i for i in shape]...)

    _position = label -> basis.site_offsets[:] .+
                         sum(basis.vectors' .* (label .- T(1.0)), dims=1)[1, :]

    data_iter = Iterators.map(_position, labels_iter)

    labels = SVector{n_vertices,eltype(labels_iter)}(collect(labels_iter)[:])
    label_data = SVector{n_vertices,SVector{D,T}}(collect(data_iter)[:])

    return labels, label_data
end
# function vertices(
#     shape::SVector{1,Tᵢ}, basis::AbstractLatticeBasis{T,1,1}
# ) where {Tᵢ<:Integer,T<:Real}
#     n_vertices = prod(shape)

#     labels_iter = 1:shape[1]

#     _position = label -> basis.site_offsets[:] .+
#                          sum(basis.vectors' .* (label .- T(1.0)), dims=1)[1, :]

#     data_iter = Iterators.map(_position, labels_iter)

#     labels = SVector{n_vertices,eltype(labels_iter)}(collect(labels_iter)[:])
#     label_data = SVector{n_vertices,SVector{1,T}}(collect(data_iter)[:])

#     return labels, label_data
# end

@with_kw struct _LatticeMetric{Tᵢ<:Integer,T<:Real,D,O} <: Metric
    shape::SVector{D,Tᵢ}
    basis::AbstractLatticeBasis{T,D,O}
    periodic::SVector{D,Bool}
    @assert D isa Integer "dimension must be an integer"
    @assert O isa Integer "number of site offsets must be an integer"
    @assert D > 0 "dimension must be positive"
    @assert O > 0 "number of site offsets must be positive"
    @assert all(shape .> 0) "shape must contain positive integers"
end
function _LatticeMetric(
    shape::AbstractVector{Tᵢ},
    basis::AbstractLatticeBasis{T,D,O},
    periodic::AbstractVector{Bool}
) where {Tᵢ<:Integer,T<:Real,D,O}
    _LatticeMetric{Tᵢ,T,D,O}(shape, basis, periodic)
end
function evaluate(
    dist::_LatticeMetric{Tᵢ,T,D,O}, a::AbstractVector{T}, b::AbstractVector{T}
) where {Tᵢ<:Integer,T<:Real,D,O}
    index_iter = Iterators.product([p ? (-1:1) : 0 for p in dist.periodic]...)

    T_mat = Iterators.map(
        x -> sum(dist.basis.vectors' .* dist.shape .* x, dims=1)[1, :],
        index_iter
    ) |> collect

    rel_pos = a .- b

    min_dist = sqrt(sum(rel_pos .^ 2))

    while true
        new_rel_pos = map(x -> x .+ rel_pos, T_mat)
        new_dist = map(p -> sqrt(sum(p .^ 2)), new_rel_pos)
        new_min_dist, min_index = findmin(new_dist)

        if CartesianIndex(min_index) == CartesianIndex(fill(Tᵢ(2), count(dist.periodic))...)
            min_dist = new_min_dist
            break
        else
            rel_pos = new_rel_pos[min_index]
        end
    end

    return min_dist
end

function _all_neighbor_dists(
    tree::BallTree{SVector{D,T},D,T,_LatticeMetric{Tᵢ,T,D,O}}, max_order::Tᵢ=ORDER;
    tol_digits::Tᵢ=TOL_DIGITS
) where {Tᵢ<:Integer,T<:Real,D,O}
    i = 1
    c = max_order

    dists = SVector{max_order,T}[]

    while true
        dists = unique(
            x -> round(x, digits=tol_digits),
            knn(tree, tree.data[i], c, true, j -> j == i)[2]
        )
        dists = filter(x -> x != T(0.0), dists)

        if length(dists) == max_order
            break
        else
            c += 1
        end
    end

    return SVector{max_order,T}(dists)
end
function _neighbor_dist(
    tree::BallTree{SVector{D,T},D,T,_LatticeMetric{Tᵢ,T,D,O}}, order::Tᵢ=ORDER;
    tol_digits::Tᵢ=TOL_DIGITS
) where {Tᵢ<:Integer,T<:Real,D,O}
    dists = _all_neighbor_dists(tree, order; tol_digits=tol_digits)
    dist = dists[order]

    return dist
end

function _edges(
    tree::BallTree{SVector{D,T},D,T,_LatticeMetric{Tᵢ,T,D,O}},
    superindeces::SVector{Nᵥ,Tᵢ},
    max_order::Tᵢ=ORDER;
    tol_digits::Tᵢ=TOL_DIGITS,
    dist_tol::T=DIST_TOL
) where {Tᵢ<:Integer,T<:Real,D,O,Nᵥ}
    @assert D > 0 "dimension must be positive"
    @assert O > 0 "number of site offsets must be positive"
    @assert max_order > 0 "`max_order` must be a positive integer"

    S = length(superindeces)
    @assert S > 0 "`superindeces` must have at least one element"

    reordered_superindices = Vector{Tᵢ}(indexin(superindeces, tree.indices))

    last_indeces = inrange(
        tree,
        tree.data[reordered_superindices],
        _neighbor_dist(tree, 1; tol_digits=tol_digits) + dist_tol,
        true
    )
    indeces = Tᵢ(1) .=> map((x, i) -> setdiff(x, i), last_indeces, superindeces)

    if max_order > 1
        for k in 1:(max_order-1)
            new_indices = inrange(
                tree,
                tree.data[reordered_superindices],
                _neighbor_dist(tree, k + 1; tol_digits=tol_digits) + dist_tol,
                true
            )

            indeces = [indeces Tᵢ(k + 1) .=> map(
                (x, y) -> setdiff(y, x),
                last_indeces,
                new_indices,
            )]

            last_indeces = new_indices
        end
    end

    return SVector{S,Vector{Pair{Tᵢ,Vector{Tᵢ}}}}(eachrow(indeces))
end
function _edges(
    tree::BallTree{SVector{D,T},D,T,_LatticeMetric{Tᵢ,T,D,O}},
    superindex::Tᵢ,
    max_order::Tᵢ=ORDER;
    tol_digits::Tᵢ=TOL_DIGITS,
    dist_tol::T=DIST_TOL
) where {Tᵢ<:Integer,T<:Real,D,O}
    return _edges(
        tree,
        SVector{1,Tᵢ}([superindex]),
        max_order;
        tol_digits=tol_digits,
        dist_tol=dist_tol
    )
end

"""
    Lattice{Tᵢ<:Integer,T<:Real,D,O}
        <: NeuralQuantumStates.Lattices.AbstractLattice{Tᵢ,T,D,O}

A `D`-dimensional `NeuralQuantumStates.Lattices.Lattice` of `shape` with given lattice
    `basis` and `periodic` boundary conditions.

# Fields
- `metagraph::MetaGraph{Tᵢ}`: A `MetaGraphsNext.MetaGraph` for the lattice to store its
    vertices and edges.
- `shape::SVector{D,Tᵢ}`: A vector for the shape of the lattice. It must contain `D`
    positive integers.
- `basis::NeuralQuantumStates.Lattices.AbstractLatticeBasis{T,D,O}`: A lattice basis for
    representing the unit cell of the lattice.
- `periodic::SVector{D,Bool}`: A vector for the periodic boundary condition of the lattice
    in each dimension.
"""
@with_kw struct Lattice{Tᵢ<:Integer,T<:Real,D,O} <: AbstractLattice{Tᵢ,T,D,O}
    metagraph::MetaGraph{Tᵢ}
    shape::SVector{D,Tᵢ}
    basis::AbstractLatticeBasis{T,D,O}
    periodic::SVector{D,Bool}
    @assert all(shape .> 0) "shape must contain positive integers"
    @assert D isa Integer "dimension must be an integer"
    @assert O isa Integer "number of site offsets must be an integer"
    @assert D > 0 "dimension must be positive"
    @assert O > 0 "number of site offsets must be positive"
end
"""
    Lattice(
        shape::SVector{D,Tᵢ},
        basis::NeuralQuantumStates.Lattices.AbstractLatticeBasis{T,D,O},
        periodic::SVector{D,Bool}=SVector{D,Bool}(fill(false, D));
        max_order::Tᵢ=1,
        tol_digits::Tᵢ=12,
        dist_tol::T=1.0e-12
    ) where {Tᵢ<:Integer,T<:Real,D,O} -> NeuralQuantumStates.Lattices.Lattice{Tᵢ,T,D,O}

Build a `D`-dimensional `NeuralQuantumStates.Lattices.Lattice` of `shape` by using the given
    lattice `basis` and `periodic` boundary conditions.

# Arguments
- `shape::SVector{D,Tᵢ}`: A vector for the shape of the lattice. It must contain `D`
    positive integers.
- `basis::NeuralQuantumStates.Lattices.AbstractLatticeBasis{T,D,O}`: A lattice basis for
    representing the unit cell of the lattice.
- `periodic::SVector{D,Bool}`: A vector for the periodic boundary condition of the lattice
    in each dimension. Defaults to `SVector{D,Bool}(fill(false, D))`.

# Keywords
- `max_order::Tᵢ`: An integer for the maximum order of the edges to be included in the
    lattice as `max_order`-nearest neighbors. Defaults to `1`, which means only nearest
    neighbors are included. For example, if it is set to `2`, then nearest and next-nearest
    neighbors are included.
- `tol_digits::Tᵢ`: An integer for the number of digits to round the calculated distances
    to. Defaults to `12`.
- `dist_tol::T`: A positive number for the tolerance of the distance between two lattice
    sites to be considered as the same site. Defaults to `1.0e-12`.

# Returns
- `NeuralQuantumStates.Lattices.Lattice{Tᵢ,T,D,O}`: The built `D`-dimensional
    `NeuralQuantumStates.Lattices.Lattice` of `shape` by using the given lattice `basis` and
    `periodic` boundary conditions.
"""
function Lattice(
    shape::SVector{D,Tᵢ}, basis::AbstractLatticeBasis{T,D,O}, periodic::SVector{D,Bool};
    max_order::Tᵢ=ORDER, tol_digits::Tᵢ=TOL_DIGITS, dist_tol::T=DIST_TOL
) where {Tᵢ<:Integer,T<:Real,D,O}
    @assert max_order > 0 "order must be a positive integer"

    ver_labels, ver_data = vertices(shape, basis)
    n_vertices = length(ver_labels)

    metric = _LatticeMetric(shape, basis, periodic)

    tree = BallTree(Vector(ver_data), metric)

    g = Graph(0)
    mg = MetaGraph(
        g;
        label_type=eltype(ver_labels),
        vertex_data_type=SVector{D,T},
        edge_data_type=Tᵢ
    )

    foreach((l, d) -> setindex!(mg, d, l), ver_labels, ver_data)

    foreach(
        (lᵢ, e) -> foreach(
            p -> foreach(jₛ -> setindex!(mg, p[1], lᵢ, ver_labels[jₛ]), p[2]),
            e
        ),
        ver_labels,
        _edges(
            tree, SVector{n_vertices,Tᵢ}(1:n_vertices), max_order;
            tol_digits=tol_digits, dist_tol=dist_tol
        )
    )

    return Lattice(mg, shape, basis, periodic)
end
function Lattice(
    shape::SVector{D,Tᵢ}, basis::AbstractLatticeBasis{T,D,O};
    max_order::Tᵢ=ORDER, tol_digits::Tᵢ=TOL_DIGITS, dist_tol::T=DIST_TOL
) where {Tᵢ<:Integer,T<:Real,D,O}
    return Lattice(
        shape, basis, SVector{D,Bool}(fill(false, D));
        max_order=max_order, tol_digits=tol_digits, dist_tol=dist_tol
    )
end
"""
    Lattice(
        shape::AbstractVector{Tᵢ},
        basis::NeuralQuantumStates.Lattices.AbstractLatticeBasis{T,D,O},
        periodic::AbstractVector{Bool}=fill(false, D);
        max_order::Tᵢ=1,
        tol_digits::Tᵢ=12,
        dist_tol::T=1.0e-12
    ) where {Tᵢ<:Integer,T<:Real,D,O} -> NeuralQuantumStates.Lattices.Lattice{Tᵢ,T,D,O}

Build a `D`-dimensional `NeuralQuantumStates.Lattices.Lattice` of `shape` by using the given
    lattice `basis` and `periodic` boundary conditions.

# Arguments
- `shape::AbstractVector{Tᵢ}`: A vector for the shape of the lattice. It must contain `D`
    positive integers.
- `basis::NeuralQuantumStates.Lattices.AbstractLatticeBasis{T,D,O}`: A lattice basis for
    representing the unit cell of the lattice.
- `periodic::AbstractVector{Bool}`: A vector for the periodic boundary condition of the
    lattice in each dimension. Defaults to `fill(false, D)`.

# Keywords
- `max_order::Tᵢ`: An integer for the maximum order of the edges to be included in the
    lattice as `max_order`-nearest neighbors. Defaults to `1`, which means only nearest
    neighbors are included. For example, if it is set to `2`, then nearest and next-nearest
    neighbors are included.
- `tol_digits::Tᵢ`: An integer for the number of digits to round the calculated distances
    to. Defaults to `12`.
- `dist_tol::T`: A positive number for the tolerance of the distance between two lattice
    sites to be considered as the same site. Defaults to `1.0e-12`.

# Returns
- `NeuralQuantumStates.Lattices.Lattice{Tᵢ,T,D,O}`: The built `D`-dimensional
    `NeuralQuantumStates.Lattices.Lattice` of `shape` by using the given lattice `basis` and
    `periodic` boundary conditions.
"""
function Lattice(
    shape::AbstractVector{Tᵢ},
    basis::AbstractLatticeBasis{T,D,O},
    periodic::AbstractVector{Bool};
    max_order::Tᵢ=ORDER,
    tol_digits::Tᵢ=TOL_DIGITS,
    dist_tol::T=DIST_TOL
) where {Tᵢ<:Integer,T<:Real,D,O}
    return Lattice(
        SVector{length(shape),Tᵢ}(shape), basis, SVector{length(periodic),Bool}(periodic);
        max_order=max_order, tol_digits=tol_digits, dist_tol=dist_tol
    )
end
function Lattice(
    shape::AbstractVector{Tᵢ}, basis::AbstractLatticeBasis{T,D,O};
    max_order::Tᵢ=ORDER, tol_digits::Tᵢ=TOL_DIGITS, dist_tol::T=DIST_TOL
) where {Tᵢ<:Integer,T<:Real,D,O}
    return Lattice(
        SVector{length(shape),Tᵢ}(shape), basis, SVector{D,Bool}(fill(false, D));
        max_order=max_order, tol_digits=tol_digits, dist_tol=dist_tol
    )
end
"""
    Lattice(
        shape::SVector{D,Tᵢ},
        basis::NeuralQuantumStates.Lattices.AbstractLatticeBasis{T,D,O},
        custom_edges::Tuple{AbstractVector{NTuple{2,L}},AbstractVector{Tᵢ}},
        periodic::SVector{D,Bool}=SVector{D,Bool}(fill(false, D))
    ) where {Tᵢ<:Integer,T<:Real,D,O,L} -> NeuralQuantumStates.Lattices.Lattice{Tᵢ,T,D,O}

Build a `D`-dimensional `NeuralQuantumStates.Lattices.Lattice` of `shape` by using the given
    lattice `basis`, `periodic` boundary conditions and `custom_edges`.

# Arguments
- `shape::SVector{D,Tᵢ}`: A vector for the shape of the lattice. It must contain `D`
    positive integers.
- `basis::NeuralQuantumStates.Lattices.AbstractLatticeBasis{T,D,O}`: A lattice basis for
    representing the unit cell of the lattice.
- `custom_edges::Tuple{AbstractVector{NTuple{2,L}},AbstractVector{Tᵢ}}`: A tuple of two
    vectors for the custom edges to be added to the lattice. The first vector contains the
    lattice site labels of the edges to be added in the form of `NTuple{2,L}` where `L` is
    the type of the lattice site labels. The second vector is of positive integers for
    distingushing the edges to be added. These two vectors must have the same length.
- `periodic::SVector{D,Bool}`: A vector for the periodic boundary condition of the lattice
    in each dimension. Defaults to `SVector{D,Bool}(fill(false, D))`.

# Returns
- `NeuralQuantumStates.Lattices.Lattice{Tᵢ,T,D,O}`: The built `D`-dimensional
    `NeuralQuantumStates.Lattices.Lattice` of `shape` by using the given lattice `basis`,
    `periodic` boundary conditions and `custom_edges`.
"""
function Lattice(
    shape::SVector{D,Tᵢ},
    basis::AbstractLatticeBasis{T,D,O},
    custom_edges::Tuple{AbstractVector{NTuple{2,L}},AbstractVector{Tᵢ}},
    periodic::SVector{D,Bool}
) where {Tᵢ<:Integer,T<:Real,D,O,L}
    @assert length(custom_edges[1]) ==
            length(custom_edges[2]) "vectors in `custom_edges` must have the same length"

    ver_labels, ver_data = vertices(shape, basis)

    @assert eltype(ver_labels) == L "vertex label type must match"

    g = Graph(0)
    mg = MetaGraph(
        g;
        label_type=eltype(ver_labels),
        vertex_data_type=SVector{D,T},
        edge_data_type=Tᵢ
    )

    foreach((l, d) -> setindex!(mg, d, l), ver_labels, ver_data)

    foreach((ll, c) -> setindex!(mg, c, ll...), custom_edges...)

    return Lattice(mg, shape, basis, periodic)
end
function Lattice(
    shape::SVector{D,Tᵢ},
    basis::AbstractLatticeBasis{T,D,O},
    custom_edges::Tuple{AbstractVector{NTuple{2,L}},AbstractVector{Tᵢ}}
) where {Tᵢ<:Integer,T<:Real,D,O,L}
    return Lattice(shape, basis, custom_edges, SVector{D,Bool}(fill(false, D)))
end
"""
    Lattice(
        shape::AbstractVector{Tᵢ},
        basis::NeuralQuantumStates.Lattices.AbstractLatticeBasis{T,D,O},
        custom_edges::Tuple{AbstractVector{NTuple{2,L}},AbstractVector{Tᵢ}},
        periodic::SVector{Bool}=fill(false, D)
    ) where {Tᵢ<:Integer,T<:Real,D,O,L} -> NeuralQuantumStates.Lattices.Lattice{Tᵢ,T,D,O}

Build a `D`-dimensional `NeuralQuantumStates.Lattices.Lattice` of `shape` by using the given
    lattice `basis`, `periodic` boundary conditions and `custom_edges`.

# Arguments
- `shape::AbstractVector{Tᵢ}`: A vector for the shape of the lattice. It must contain `D`
    positive integers.
- `basis::NeuralQuantumStates.Lattices.AbstractLatticeBasis{T,D,O}`: A lattice basis for
    representing the unit cell of the lattice.
- `custom_edges::Tuple{AbstractVector{NTuple{2,L}},AbstractVector{Tᵢ}}`: A tuple of two
    vectors for the custom edges to be added to the lattice. The first vector contains the
    lattice site labels of the edges to be added in the form of `NTuple{2,L}` where `L` is
    the type of the lattice site labels. The second vector is of positive integers for
    distingushing the edges to be added. These two vectors must have the same length.
- `periodic::AbstractVector{Bool}`: A vector for the periodic boundary condition of the
    lattice in each dimension. Defaults to `fill(false, D)`.

# Returns
- `NeuralQuantumStates.Lattices.Lattice{Tᵢ,T,D,O}`: The built `D`-dimensional
    `NeuralQuantumStates.Lattices.Lattice` of `shape` by using the given lattice `basis`,
    `periodic` boundary conditions and `custom_edges`.
"""
function Lattice(
    shape::AbstractVector{Tᵢ},
    basis::AbstractLatticeBasis{T,D,O},
    custom_edges::Tuple{AbstractVector{NTuple{2,L}},AbstractVector{Tᵢ}},
    periodic::AbstractVector{Bool}
) where {Tᵢ<:Integer,T<:Real,D,O,L}
    return Lattice(
        SVector{length(shape),Tᵢ}(shape),
        basis,
        custom_edges,
        SVector{length(periodic),Bool}(periodic)
    )
end
function Lattice(
    shape::AbstractVector{Tᵢ},
    basis::AbstractLatticeBasis{T,D,O},
    custom_edges::Tuple{AbstractVector{NTuple{2,L}},AbstractVector{Tᵢ}}
) where {Tᵢ<:Integer,T<:Real,D,O,L}
    return Lattice(shape, basis, custom_edges, fill(false, D))
end

function nv(lattice::Lattice{Tᵢ,T,D,O}) where {Tᵢ<:Integer,T<:Real,D,O}
    return _nv(lattice.metagraph)
end
