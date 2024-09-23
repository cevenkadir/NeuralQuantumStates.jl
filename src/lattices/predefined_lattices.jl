using Graphs: adjacency_matrix, cartesian_product, cycle_graph, path_graph
using LinearAlgebra: Diagonal, I, triu
using NeuralQuantumStates.Lattices
using SparseArrays: findnz
using StaticArrays

function _vectorize(periodic::Bool, shape::SVector{D,Tᵢ}) where {Tᵢ<:Integer,D}
    periodic = SVector{D,Bool}(fill(periodic, D))

    # check whether `periodic[i]` is false if `shape[i] == 1`
    @assert !any(periodic .== true * shape .== 1) "periodic must be false for shape == 1"

    return periodic
end
function _vectorize(periodic::Bool, shape::AbstractVector{Tᵢ}) where {Tᵢ<:Integer}
    return _vectorize(periodic, SVector{length(shape),Tᵢ}(shape))
end
function _vectorize(periodic::SVector{D,Bool}, shape::SVector{D,Tᵢ}) where {Tᵢ<:Integer,D}
    # check whether `periodic[i]` is false if `shape[i] == 1`
    @assert !any(periodic .== true * shape .== 1) "periodic must be false for shape == 1"

    return periodic
end
function _vectorize(
    periodic::AbstractVector{Bool}, shape::AbstractVector{T}
) where {T<:Integer}
    return _vectorize(
        SVector{length(periodic),Bool}(periodic),
        SVector{length(shape),T}(shape),
    )
end

build(lattice_sym::Symbol, args...; kwargs...) = build(Val(lattice_sym), args...; kwargs...)

"""
    build(
        ::Val{:Triclinic},
        shape::AbstractVector{Tᵢ},
        edge_lengths::AbstractVector{T},
        angles::AbstractVector{T};
        periodic::Union{Bool,AbstractVector{Bool}}=false,
    ) where {Tᵢ<:Integer,T<:Real} -> NeuralQuantumStates.Lattices.Lattice{Tᵢ,T,3,1}

Build a triclinic lattice from the given parameters.

# Arguments
- `::Val{:Triclinic}`: A value to dispatch to this function.
- `shape::AbstractVector{Tᵢ}`: A vector for the shape of the lattice. It must contain three
    positive integers.
- `edge_lengths::AbstractVector{T}`: A vector for the edge lengths of the lattice. It must
    contain three positive numbers.
- `angles::AbstractVector{T}`: A vector for the angles in degrees between the edge lengths
    of the lattice. It must contain three numbers. `angles[i]` is the angle between
    `edge_lengths[j]` and `edge_lengths[k]` where `(i,j,k)` is a permutation of `(1,2,3)`.

# Keywords
- `periodic::Union{Bool,AbstractVector{Bool}}`: A boolean or a vector for the periodic
    boundary condition of the lattice in each dimension. If it is a boolean, then it is
    applied to all dimensions. If it is a vector, then it must contain three booleans.
    Defaults to `false`.

# Returns
- `NeuralQuantumStates.Lattices.Lattice{Tᵢ,T,3,1}`: The built triclinic lattice from the
    given parameters.
"""
function build(
    ::Val{:Triclinic},
    shape::AbstractVector{Tᵢ},
    edge_lengths::AbstractVector{T},
    angles::AbstractVector{T};
    periodic::Union{Bool,AbstractVector{Bool}}=false
) where {Tᵢ<:Integer,T<:Real}
    @assert all(shape .> 0) "shape must contain positive integers"
    @assert length(shape) == 3 "shape must contain three integers"
    @assert all(edge_lengths .> 0) "edge lengths must be positive"

    # vectorize the periodic boundary condition
    periodic = _vectorize(periodic, shape)

    vectors = Diagonal(edge_lengths)
    α, β, γ = angles
    vectors = vectors * SMatrix{3,3,T}(
        [
            1.0 0.0 0.0
            cosd(γ) sind(γ) 0
            cosd(β) (cosd(α)-cosd(β)*cosd(γ))/sind(γ) sqrt(1 - (cosd(α)^2 + cosd(β)^2) / sind(γ)^2)
        ]
    )
    vectors = vectors |> transpose

    basis = LatticeBasis(vectors)

    # create the graph depending on periodicity in the first dimension
    if periodic[1]
        g = cycle_graph(shape[1])
    else
        g = path_graph(shape[1])
    end

    # continue to create the graph depending on periodicity in the remaining dimensions
    for (index, s) in enumerate(shape[2:end])
        if periodic[index+1]
            g = cartesian_product(cycle_graph(s), g)
        else
            g = cartesian_product(path_graph(s), g)
        end
    end

    pos_labels = collect(Iterators.product([1:i for i in shape]...))[:]
    i_s, j_s, = findnz(triu(adjacency_matrix(g)))
    edge_labels = map((i, j) -> (pos_labels[i], pos_labels[j]), i_s, j_s)
    custom_edges = (edge_labels, fill(Tᵢ(1), length(edge_labels)))

    return Lattice(shape, basis, custom_edges, periodic)
end

"""
    build(
        ::Val{:Hypercube},
        shape::AbstractVector{Tᵢ},
        edge_length::T;
        periodic::Union{Bool,AbstractVector{Bool}}=false,
    ) where {Tᵢ<:Integer,T<:Real}
        -> NeuralQuantumStates.Lattices.Lattice{Tᵢ,T,length(shape),1}

Build a hypercubic lattice from the given parameters.

# Arguments
- `::Val{:Hypercube}`: A value to dispatch to this function.
- `shape::AbstractVector{Tᵢ}`: A vector for the shape of the lattice.
- `edge_length::T`: A positive number for the edge length of the lattice.

# Keywords
- `periodic::Union{Bool,AbstractVector{Bool}}`: A boolean or a vector for the periodic
    boundary condition of the lattice in each dimension. If it is a boolean, then it is
    applied to all dimensions. If it is a vector, then it must contain `length(shape)`
    booleans. Defaults to `false`.

# Returns
- `NeuralQuantumStates.Lattices.Lattice{Tᵢ,T,length(shape),1}`: The built hypercubic lattice
    from the given parameters.
"""
function build(
    ::Val{:Hypercube},
    shape::AbstractVector{Tᵢ},
    edge_length::T;
    periodic::Union{Bool,AbstractVector{Bool}}=false
) where {Tᵢ<:Integer,T<:Real}
    @assert edge_length > 0 "edge length must be positive"
    @assert all(shape .> 0) "shape must contain positive integers"
    @assert length(shape) > 0 "shape must contain at least one integer"

    periodic = _vectorize(periodic, shape)

    vectors = hypercube.edge_length * I(D)

    basis = LatticeBasis(vectors)

    # create the graph depending on periodicity in the first dimension
    if hypercube.periodic[1]
        g = cycle_graph(shape[1])
    else
        g = path_graph(shape[1])
    end

    # continue to create the graph depending on periodicity in the remaining dimensions
    for (index, s) in enumerate(shape[2:end])
        if periodic[index+1]
            g = cartesian_product(cycle_graph(s), g)
        else
            g = cartesian_product(path_graph(s), g)
        end
    end

    pos_labels = collect(Iterators.product([1:i for i in shape]...))[:]
    i_s, j_s, = findnz(triu(adjacency_matrix(g)))
    edge_labels = map((i, j) -> (pos_labels[i], pos_labels[j]), i_s, j_s)
    custom_edges = (edge_labels, fill(Tᵢ(1), length(edge_labels)))

    return Lattice(shape, basis, custom_edges, periodic)
end

"""
    build(
        ::Val{:Triangular},
        shape::AbstractVector{Tᵢ},
        edge_length::T;
        periodic::Union{Bool,AbstractVector{Bool}}=false,
        tol_digits::Tᵢ=TOL_DIGITS,
        dist_tol::T=DIST_TOL
    ) where {Tᵢ<:Integer,T<:Real} -> NeuralQuantumStates.Lattices.Lattice{Tᵢ,T,2,1}

Build a 2D triangular lattice from the given parameters.

# Arguments
- `::Val{:Triangular}`: A value to dispatch to this function.
- `shape::AbstractVector{Tᵢ}`: A vector for the shape of the lattice. It must contain two
    positive integers.
- `edge_length::T`: A positive number for the edge length of the lattice.

# Keywords
- `periodic::Union{Bool,AbstractVector{Bool}}`: A boolean or a vector for the periodic
    boundary condition of the lattice in each dimension. If it is a boolean, then it is
    applied to all dimensions. If it is a vector, then it must contain two booleans.
    Defaults to `false`.
- `tol_digits::Tᵢ`: An integer for the number of digits to round the calculated distances
    to. Defaults to `12`.
- `dist_tol::T`: A positive number for the tolerance of the distance between two lattice
    sites to be considered as the same site. Defaults to `1.0e-12`.

# Returns
- `NeuralQuantumStates.Lattices.Lattice{Tᵢ,T,2,1}`: The built 2D triangular lattice from the
    given parameters.
"""
function build(
    ::Val{:Triangular},
    shape::AbstractVector{Tᵢ},
    edge_length::T;
    periodic::Union{Bool,AbstractVector{Bool}}=false,
    tol_digits::Tᵢ=TOL_DIGITS,
    dist_tol::T=DIST_TOL
) where {Tᵢ<:Integer,T<:Real}
    @assert all(shape .> 0) "shape must contain positive integers"
    @assert length(shape) == 2 "shape must contain two integers"
    @assert edge_length > 0 "edge length must be positive"

    periodic = _vectorize(periodic, shape)

    vectors = edge_length * SMatrix{2,2,T}([1.0 0.5; 0.0 sqrt(0.75)])

    basis = LatticeBasis(vectors)

    return Lattice(
        shape, basis, periodic;
        max_order=Tᵢ(1), tol_digits=tol_digits, dist_tol=dist_tol
    )
end

"""
    build(
        ::Val{:Honeycomb},
        shape::AbstractVector{Tᵢ},
        edge_length::T;
        periodic::Union{Bool,AbstractVector{Bool}}=false,
        tol_digits::Tᵢ=TOL_DIGITS,
        dist_tol::T=DIST_TOL
    ) where {Tᵢ<:Integer,T<:Real} -> NeuralQuantumStates.Lattices.Lattice{Tᵢ,T,2,1}

Build a 2D honeycomb lattice from the given parameters.

# Arguments
- `shape::AbstractVector{Tᵢ}`: A vector for the shape of the lattice. It must contain two
    positive integers.
- `edge_length::T`: A positive number for the edge length of the lattice.

# Keywords
- `periodic::Union{Bool,AbstractVector{Bool}}`: A boolean or a vector for the periodic
    boundary condition of the lattice in each dimension. If it is a boolean, then it is
    applied to all dimensions. If it is a vector, then it must contain two booleans.
    Defaults to `false`.
- `tol_digits::Tᵢ`: An integer for the number of digits to round the calculated distances
    to. Defaults to `12`.
- `dist_tol::T`: A positive number for the tolerance of the distance between two lattice
    sites to be considered as the same site. Defaults to `1.0e-12`.

# Returns
- `NeuralQuantumStates.Lattices.Lattice{Tᵢ,T,2,2}`: The built 2D honeycomb lattice from the
    given parameters.
"""
function build(
    ::Val{:Honeycomb},
    shape::AbstractVector{Tᵢ},
    edge_length::T;
    periodic::Union{Bool,AbstractVector{Bool}}=false,
    tol_digits::Tᵢ=TOL_DIGITS,
    dist_tol::T=DIST_TOL
) where {Tᵢ<:Integer,T<:Real}
    @assert all(shape .> 0) "shape must contain positive integers"
    @assert length(shape) == 2 "shape must contain two integers"
    @assert edge_length > 0 "edge length must be positive"

    periodic = _vectorize(periodic, shape)

    vectors = edge_length * SMatrix{2,2,T}([1.0 0.5; 0.0 sqrt(0.75)])
    site_offsets = SMatrix{2,2,T}([0.5 1; 0.5/3^0.5 1/3^0.5])

    basis = LatticeBasis(vectors, site_offsets)

    return Lattice(
        shape, basis, periodic;
        max_order=Tᵢ(1), tol_digits=tol_digits, dist_tol=dist_tol
    )
end

"""
    build(
        ::Val{:Kagome},
        shape::AbstractVector{Tᵢ},
        edge_length::T;
        periodic::Union{Bool,AbstractVector{Bool}}=false,
        tol_digits::Tᵢ=TOL_DIGITS,
        dist_tol::T=DIST_TOL
    ) where {Tᵢ<:Integer,T<:Real} -> NeuralQuantumStates.Lattices.Lattice{Tᵢ,T,2,3}

Build a 2D kagome lattice from the given parameters.

# Arguments
- `shape::AbstractVector{Tᵢ}`: A vector for the shape of the lattice. It must contain two
    positive integers.
- `edge_length::T`: A positive number for the edge length of the lattice.

# Keywords
- `periodic::Union{Bool,AbstractVector{Bool}}`: A boolean or a vector for the periodic
    boundary condition of the lattice in each dimension. If it is a boolean, then it is
    applied to all dimensions. If it is a vector, then it must contain two booleans.
    Defaults to `false`.
- `tol_digits::Tᵢ`: An integer for the number of digits to round the calculated distances
    to. Defaults to `12`.
- `dist_tol::T`: A positive number for the tolerance of the distance between two lattice
    sites to be considered as the same site. Defaults to `1.0e-12`.

# Returns
- `NeuralQuantumStates.Lattices.Lattice{Tᵢ,T,2,3}`: The built 2D kagome lattice from the
    given parameters.
"""
function build(
    ::Val{:Kagome},
    shape::AbstractVector{Tᵢ},
    edge_length::T;
    periodic::Union{Bool,AbstractVector{Bool}}=false,
    tol_digits::Tᵢ=TOL_DIGITS,
    dist_tol::T=DIST_TOL
) where {Tᵢ<:Integer,T<:Real}
    @assert all(shape .> 0) "shape must contain positive integers"
    @assert edge_length > 0 "edge length must be positive"
    @assert length(shape) == 2 "shape must contain two integers"

    periodic = _vectorize(periodic, shape)

    vectors = edge_length * SMatrix{2,2,T}([1.0 0.5; 0.0 sqrt(3.0)/2.0])
    site_offsets = SMatrix{2,3,T}([[0.0; 0.0] vectors ./ 2.0])

    basis = LatticeBasis(vectors, site_offsets)

    return Lattice(
        shape, basis, periodic;
        max_order=Tᵢ(1), tol_digits=tol_digits, dist_tol=dist_tol
    )
end
