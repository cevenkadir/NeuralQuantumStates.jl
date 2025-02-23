using Parameters
using StaticArrays

using LinearAlgebra

using MetaGraphsNext: edges, src, dst, weights

using NeuralQuantumStates.Hilberts
using NeuralQuantumStates.Lattices
using NeuralQuantumStates.Operators
using NeuralQuantumStates.Extras

const DEFAULT_bosehubbard_T = Union{Complex,Real}
const DEFAULT_bosehubbard_T_h = Union{Integer}
const DEFAULT_bosehubbard_T_h_Array = AbstractArray
const DEFAULT_bosehubbard_Tᵢ = Integer
const DEFAULT_bosehubbard_T_l = Real
const DEFAULT_bosehubbard_T_p = Number

extend_cartsian_index(x, y) = CartesianIndex.(x, (Ref.(y))...)

"""
    BoseHubbardOperator{
        T<:DEFAULT_bosehubbard_T,
        T_h<:DEFAULT_bosehubbard_T_h,
        T_h_Array<:DEFAULT_bosehubbard_T_h_Array,
        Tᵢ<:DEFAULT_bosehubbard_Tᵢ,
        T_l<:DEFAULT_bosehubbard_T_l,
        T_p<:DEFAULT_bosehubbard_T_p,
        N,n_lDoF,D,O
    }(
        hilbert::NeuralQuantumStates.Hilberts.FiniteUniformHilbert{T_h,T_h_Array,N,n_lDoF},
        lattice::NeuralQuantumStates.Lattices.Lattice{Tᵢ,T_l,D,O},
        J::T_p,
        U::T_p,
        V::T_p,
        μ::T_p
    ) <: NeuralQuantumStates.Operators.DiscreteOperator{T}

Extended Bose-Hubbard operator instance according to the following formula:
``\\hat{H} = -J \\sum_{\\langle i, j \\rangle} \\left( \\hat{b}^\\dagger_i \\hat{b}_j +
\\mathrm{h.c.} \\right) + \\frac{U}{2}\\sum_{i} \\hat{n}_i \\left( \\hat{n}_i - 1 \\right) +
V\\sum_{\\langle i, j \\rangle} \\hat{n}_i \\hat{n}_j + \\mu \\sum_i \\hat{n}_i``

# Fields
- `hilbert::NeuralQuantumStates.Hilberts.FiniteUniformHilbert{T_h,T_h_Array,N,n_lDoF}`: An
    instance of `NeuralQuantumStates.Hilberts.FiniteUniformHilbert`.
- `lattice::NeuralQuantumStates.Lattices.Lattice{Tᵢ,T_l,D,O}`: An instance of
    `NeuralQuantumStates.Lattices.Lattice`.
- `J::T_p`: Hopping amplitude.
- `U::T_p`: On-site interaction strength.
- `V::T_p`: Extended-range interaction strength.
- `μ::T_p`: Chemical potential.
"""
@with_kw struct BoseHubbardOperator{
    T<:DEFAULT_bosehubbard_T,
    T_h<:DEFAULT_bosehubbard_T_h,
    T_h_Array<:DEFAULT_bosehubbard_T_h_Array,
    Tᵢ<:DEFAULT_bosehubbard_Tᵢ,
    T_l<:DEFAULT_bosehubbard_T_l,
    T_p<:DEFAULT_bosehubbard_T_p,
    N,n_lDoF,D,O
} <: DiscreteOperator{T}
    hilbert::Hilberts.FiniteUniformHilbert{T_h,T_h_Array,N,n_lDoF}
    lattice::Lattices.Lattice{Tᵢ,T_l,D,O}
    J::T_p
    U::T_p
    V::T_p
    μ::T_p
end

"""
    build(
        ::Val{:ExtendedBoseHubbard},
        hilbert::NeuralQuantumStates.Hilberts.FiniteUniformHilbert{T_h,T_h_Array,N,n_lDoF},
        lattice::NeuralQuantumStates.Lattices.Lattice{Tᵢ,T_l,D,O};
        J::T_p=1.0, U::T_p=1.0, V::T_p=1.0, μ::T_p=0.0, T=Float64
    ) where {
        T_h<:DEFAULT_bosehubbard_T_h,
        T_h_Array<:DEFAULT_bosehubbard_T_h_Array,
        Tᵢ<:DEFAULT_bosehubbard_Tᵢ,
        T_l<:DEFAULT_bosehubbard_T_l,
        T_p<:DEFAULT_bosehubbard_T_p,
        N,n_lDoF,D,O
    }

Build an extended Bose-Hubbard operator from the given parameters according to the following
    formula:
``\\hat{H} = -J \\sum_{\\langle i, j \\rangle} \\left( \\hat{b}^\\dagger_i \\hat{b}_j +
\\mathrm{h.c.} \\right) + \\frac{U}{2}\\sum_{i} \\hat{n}_i \\left( \\hat{n}_i - 1 \\right) +
V\\sum_{\\langle i, j \\rangle} \\hat{n}_i \\hat{n}_j + \\mu \\sum_i \\hat{n}_i``

# Arguments
- `::Val{:ExtendedBoseHubbard}`: A value to dispatch to this function.
- `hilbert::NeuralQuantumStates.Hilberts.FiniteUniformHilbert{T_h,T_h_Array,N,n_lDoF}`: An
    instance of `NeuralQuantumStates.Hilberts.FiniteUniformHilbert`.
- `lattice::NeuralQuantumStates.Lattices.Lattice{Tᵢ,T_l,D,O}`: An instance of
    `NeuralQuantumStates.Lattices.Lattice`.

# Keywords
- `J::T_p`: Hopping amplitude. Defaults to `1.0`.
- `U::T_p`: On-site interaction strength. Defaults to `1.0`.
- `V::T_p`: Extended-range interaction strength. Defaults to `1.0`.
- `μ::T_p`: Chemical potential. Defaults to `0.0`.
- `T::DataType`: Type of numbers in the operator. Defaults to `Float64`.
"""
function build(
    ::Val{:ExtendedBoseHubbard},
    hilbert::Hilberts.FiniteUniformHilbert{T_h,T_h_Array,N,n_lDoF},
    lattice::Lattices.Lattice{Tᵢ,T_l,D,O};
    J::T_p=1.0, U::T_p=1.0, V::T_p=1.0, μ::T_p=0.0, T::DataType=Float64
) where {
    T_h<:DEFAULT_bosehubbard_T_h,
    T_h_Array<:DEFAULT_bosehubbard_T_h_Array,
    Tᵢ<:DEFAULT_bosehubbard_Tᵢ,
    T_l<:DEFAULT_bosehubbard_T_l,
    T_p<:DEFAULT_bosehubbard_T_p,
    N,n_lDoF,D,O
}
    @assert hilbert.type == :Fock "The Hilbert space must be a Fock space."

    return BoseHubbardOperator{T,T_h,T_h_Array,Tᵢ,T_l,T_p,N,n_lDoF,D,O}(
        hilbert=hilbert,
        lattice=lattice,
        J=J,
        U=U,
        V=V,
        μ=μ
    )
end

function add_at!(
    operator::BoseHubbardOperator{T,T_h,T_h_Array,Tᵢ,T_l,T_p,N,n_lDoF,D,O},
    x::AbstractVector{Union{T_h,Missing}},
    idx::Tᵢ;
    addend::Integer=1
) where {T,T_h,T_h_Array,Tᵢ,T_l,T_p,N,n_lDoF,D,O}
    @views x[idx] += addend
end

"""
    connected_basis_configs(
        operator::NeuralQuantumStates.Operators.BoseHubbardOperator{
            T,T_h,T_h_Array,Tᵢ,T_l,T_p,N,n_lDoF,D,O
        },
        sample::AbstractVector{T_h}
    ) where {
        T<:DEFAULT_bosehubbard_T,
        T_h<:DEFAULT_bosehubbard_T_h,
        T_h_Array<:DEFAULT_bosehubbard_T_h_Array,
        Tᵢ<:DEFAULT_bosehubbard_Tᵢ,
        T_l<:DEFAULT_bosehubbard_T_l,
        T_p<:DEFAULT_bosehubbard_T_p,
        N,n_lDoF,D,O
    }

Return the connected basis configurations ``\\vert s^\\prime \\rangle`` to a given operator
    ``\\hat{O}``, *i.e.* `operator`, for a given configuration ``\\vert s \\rangle``, *i.e.*
    `sample`, with their matrix elements ``\\langle s^\\prime \\lvert \\hat{O} \\rvert s
    \\rangle``.
For this function, `sample` needs to be a vector that contains only one configuration.

# Arguments
- `operator::NeuralQuantumStates.Operators.BoseHubbardOperator{
    T,T_h,T_h_Array,Tᵢ,T_l,T_p,N,n_lDoF,D,O}`: An instance of
    `NeuralQuantumStates.Operators.BoseHubbardOperator`.
- `sample::AbstractVector{T_h}`: A representative vector for the given basis configuration.

# Returns
- `T_h_Array{Union{T_h,Missing}, 2}`: A ``N \\times M`` matrix of ``M`` connected basis
    configurations.
- `Vector{T}`: A vector of ``M`` matrix elements ``\\langle s^\\prime \\lvert \\hat{O}
    \\rvert s \\rangle``.
"""
function connected_basis_configs(
    operator::BoseHubbardOperator{T,T_h,T_h_Array,Tᵢ,T_l,T_p,N,n_lDoF,D,O},
    sample::AbstractVector{T_h}
) where {
    T<:DEFAULT_bosehubbard_T,
    T_h<:DEFAULT_bosehubbard_T_h,
    T_h_Array<:DEFAULT_bosehubbard_T_h_Array,
    Tᵢ<:DEFAULT_bosehubbard_Tᵢ,
    T_l<:DEFAULT_bosehubbard_T_l,
    T_p<:DEFAULT_bosehubbard_T_p,
    N,n_lDoF,D,O
}
    @assert length(sample) == N """
    The sample must have the same length as the number of degrees of freedom.
    """
    @assert Hilberts.check(operator.hilbert.constraint, sample)

    # get the maximum occupation number
    n_max = operator.hilbert.lDoF[end]

    array_type_params = Extras.get_array_type_params(T_h_Array, Union{T_h,Missing}, (N,))
    T_h_Array_w_params = T_h_Array{array_type_params...}

    if typeof(sample) != T_h_Array_w_params
        x = T_h_Array_w_params(sample)
    else
        x = sample
    end

    # get nearest neighbor pairs
    nn_idx_pairs = edges(operator.lattice.metagraph)

    i = [edge.src for edge in nn_idx_pairs]
    j = [edge.dst for edge in nn_idx_pairs]

    # get the occupation numbers
    nᵢ = stack(x[edge.src] for edge in nn_idx_pairs)
    nⱼ = stack(x[edge.dst] for edge in nn_idx_pairs)

    # get the weights
    ws = weights(operator.lattice.metagraph)

    # initialize the configurations and values
    configs₀ = T_h_Array_w_params[copy(x)]
    vals₀ = zeros(T, 1)

    # from diagonal terms
    vals₀ .+= 0.5 * operator.U * sum(x .* (x .- 1))
    vals₀ .+= operator.V * sum(nᵢ .* nⱼ .* getindex.(Ref(ws), i, j))
    vals₀ .-= operator.μ * sum(x)

    # from off-diagonal terms (bᵢ^† bⱼ)
    maskᵢⱼ = @. (nⱼ > 0) && (nᵢ < n_max)
    valsᵢⱼ = @. -operator.J * √(nⱼ[maskᵢⱼ]) * √(nᵢ[maskᵢⱼ] + 1)
    valsᵢⱼ .*= getindex.(Ref(ws), i, j)[maskᵢⱼ]
    configsᵢⱼ = T_h_Array_w_params[copy(x) for _ in 1:sum(maskᵢⱼ)]
    add_at!.(Ref(operator), configsᵢⱼ, i[maskᵢⱼ]; addend=1)
    add_at!.(Ref(operator), configsᵢⱼ, j[maskᵢⱼ]; addend=-1)

    # from off-diagonal terms (bⱼ^† bᵢ)
    maskⱼᵢ = @. (nᵢ > 0) && (nⱼ < n_max)
    valsⱼᵢ = @. -operator.J * √(nᵢ[maskⱼᵢ]) * √(nⱼ[maskⱼᵢ] + 1)
    valsⱼᵢ .*= getindex.(Ref(ws), j, i)[maskⱼᵢ]
    configsⱼᵢ = T_h_Array_w_params[copy(x) for _ in 1:sum(maskⱼᵢ)]
    add_at!.(Ref(operator), configsⱼᵢ, i[maskⱼᵢ]; addend=-1)
    add_at!.(Ref(operator), configsⱼᵢ, j[maskⱼᵢ]; addend=1)

    # concatenate the configurations and values
    configs = vcat(configs₀, configsᵢⱼ, configsⱼᵢ)
    vals = vcat(vals₀, valsᵢⱼ, valsⱼᵢ)

    return hcat(configs...), vals
end

"""
    connected_basis_configs(
        operator::NeuralQuantumStates.Operators.BoseHubbardOperator{
            T,T_h,T_h_Array,Tᵢ,T_l,T_p,N,n_lDoF,D,O
        },
        samples::AbstractVector{<:AbstractVector{T_h}};
        kwargs...
    ) where {
        T<:DEFAULT_bosehubbard_T,
        T_h<:DEFAULT_bosehubbard_T_h,
        T_h_Array<:DEFAULT_bosehubbard_T_h_Array,
        Tᵢ<:DEFAULT_bosehubbard_Tᵢ,
        T_l<:DEFAULT_bosehubbard_T_l,
        T_p<:DEFAULT_bosehubbard_T_p,
        N,n_lDoF,D,O
    }

Return the connected basis configurations ``\\vert s^\\prime \\rangle`` to a given operator
    ``\\hat{O}``, *i.e.* `operator`, for given configurations ``\\vert s \\rangle``, *i.e.*
    `samples`, with their matrix elements ``\\langle s^\\prime \\lvert \\hat{O} \\rvert s
    \\rangle``.
For this function, `samples` needs to be a vector of vectors that contains the given
    configurations.

# Arguments
- `operator::NeuralQuantumStates.Operators.BoseHubbardOperator{
    T,T_h,T_h_Array,Tᵢ,T_l,T_p,N,n_lDoF,D,O}`: An instance of
    `NeuralQuantumStates.Operators.BoseHubbardOperator`.
- `samples::AbstractVector{<:AbstractVector{T_h}}`: The representative vectors for the given
    basis configurations.

# Returns
- `Vector{Tuple{T_h_Array{Union{T_h,Missing}, 2}, Vector{T}}}`: The ``N \\times M_i``
    matrices of ``M_i`` connected basis configurations and the vectors of `M_i` matrix
    elements ``\\langle s^\\prime \\lvert \\hat{O} \\rvert s \\rangle`` for each the
    ``i``-th element of `samples`.
"""
function connected_basis_configs(
    operator::BoseHubbardOperator{T,T_h,T_h_Array,Tᵢ,T_l,T_p,N,n_lDoF,D,O},
    samples::AbstractVector{<:AbstractVector{T_h}};
    kwargs...
) where {
    T<:DEFAULT_bosehubbard_T,
    T_h<:DEFAULT_bosehubbard_T_h,
    T_h_Array<:DEFAULT_bosehubbard_T_h_Array,
    Tᵢ<:DEFAULT_bosehubbard_Tᵢ,
    T_l<:DEFAULT_bosehubbard_T_l,
    T_p<:DEFAULT_bosehubbard_T_p,
    N,n_lDoF,D,O
}
    return connected_basis_configs.(Ref(operator), samples; kwargs...)
end

"""
    connected_basis_configs(
        operator::NeuralQuantumStates.Operators.BoseHubbardOperator{
            T,T_h,T_h_Array,Tᵢ,T_l,T_p,N,n_lDoF,D,O
        },
        samples::AbstractArray{T_h};
        kwargs...
    ) where {
        T<:DEFAULT_bosehubbard_T,
        T_h<:DEFAULT_bosehubbard_T_h,
        T_h_Array<:DEFAULT_bosehubbard_T_h_Array,
        Tᵢ<:DEFAULT_bosehubbard_Tᵢ,
        T_l<:DEFAULT_bosehubbard_T_l,
        T_p<:DEFAULT_bosehubbard_T_p,
        N,n_lDoF,D,O
    }

Return the connected basis configurations ``\\vert s^\\prime \\rangle`` to a given operator
    ``\\hat{O}``, *i.e.* `operator`, for given configurations ``\\vert s \\rangle``, *i.e.*
    `samples`, with their matrix elements ``\\langle s^\\prime \\lvert \\hat{O} \\rvert s
    \\rangle``.
For this function, `samples` needs to be a ``\\dots \\times N`` array (or tensor) contains
    the given configurations.

# Arguments
- `operator::NeuralQuantumStates.Operators.BoseHubbardOperator{
    T,T_h,T_h_Array,Tᵢ,T_l,T_p,N,n_lDoF,D,O}`: An instance of
    `NeuralQuantumStates.Operators.BoseHubbardOperator`.
- `samples::AbstractArray{T_h}`: A ``\\dots \\times N`` array (or tensor) of the
    representative vectors for the given basis configurations.

# Returns
- `T_h_Array{Union{T_h,Missing}, ndims(samples) + 1}`: A ``M_\\mathrm{max} \\times \\dots
    \\times N`` array (or tensor) of connected basis configurations where
    ``M_\\mathrm{max}`` is the maximum of connections for each given configuration. Each
    element of the constructed configurations as extra will be assigned as `missing`.
- `Array{T, ndims(samples)}`: A ``M_\\mathrm{max} \\times \\dots`` array or (tensor) of
    matrix elements ``\\langle s^\\prime \\lvert \\hat{O} \\rvert s \\rangle``. Each matrix
    elements will be assigned as `missing`.
"""
function connected_basis_configs(
    operator::BoseHubbardOperator{T,T_h,T_h_Array,Tᵢ,T_l,T_p,N,n_lDoF,D,O},
    samples::AbstractArray{T_h};
    kwargs...
) where {
    T<:DEFAULT_bosehubbard_T,
    T_h<:DEFAULT_bosehubbard_T_h,
    T_h_Array<:DEFAULT_bosehubbard_T_h_Array,
    Tᵢ<:DEFAULT_bosehubbard_Tᵢ,
    T_l<:DEFAULT_bosehubbard_T_l,
    T_p<:DEFAULT_bosehubbard_T_p,
    N,n_lDoF,D,O
}
    batch_dims = ndims(samples) - 1
    total_dims = batch_dims + 1

    @assert size(samples, total_dims) == N """
      The sample must have the same length as the number of degrees of freedom.
      """

    size_samples = size(samples)

    # get the maximum occupation number
    n_max = operator.hilbert.lDoF[end]

    array_type_params = get_array_type_params(T_h_Array, T_h, size_samples)
    _array_type_params = get_array_type_params(T_h_Array, T_h, (1, size_samples...))
    T_h_Array_w_params = T_h_Array{array_type_params...}
    _T_h_Array_w_params = T_h_Array{_array_type_params...}

    if typeof(samples) != T_h_Array_w_params
        xs = T_h_Array_w_params(samples)
    else
        xs = samples
    end

    _xs = reshape(xs, :, size_samples...) |> _T_h_Array_w_params

    # get nearest neighbor pairs
    nn_idx_pairs = edges(operator.lattice.metagraph)
    i = [edge.src for edge in nn_idx_pairs]
    j = [edge.dst for edge in nn_idx_pairs]

    # get the occupation numbers
    nᵢ = @views xs[fill(:, batch_dims)..., src.(nn_idx_pairs)]
    nⱼ = @views xs[fill(:, batch_dims)..., dst.(nn_idx_pairs)]

    # get the weights
    ws = weights(operator.lattice.metagraph)

    # for off-diagonal terms (bᵢ^† bⱼ)
    maskᵢⱼ = @. (nⱼ > 0) && (nᵢ < n_max)
    wsᵢⱼ = reshape(getindex.(Ref(ws), i, j), (fill(1, batch_dims)..., N))
    Nᵢⱼ = sum(maskᵢⱼ, dims=total_dims)

    # for off-diagonal terms (bⱼ^† bᵢ)
    maskⱼᵢ = @. (nᵢ > 0) && (nⱼ < n_max)
    wsⱼᵢ = reshape(getindex.(Ref(ws), j, i), (fill(1, batch_dims)..., N))
    Nⱼᵢ = sum(maskⱼᵢ, dims=total_dims)

    N_configs = 1 .+ Nᵢⱼ .+ Nⱼᵢ
    N_max_configs = N_configs |> maximum

    # initialize matrices of "configs" and "vals" in which the output results are stored
    configs_array_type_params = get_array_type_params(
        T_h_Array, Union{T_h,Missing},
        (N_max_configs, size_samples...)
    )
    configs_T_h_Array_w_params = T_h_Array{configs_array_type_params...}
    configs = repeat(
        _xs;
        inner=(N_max_configs, fill(1, total_dims)...)
    ) |> configs_T_h_Array_w_params
    vals = zeros.(Ref(Union{T,Missing}), N_max_configs, size_samples[1:end-1]...)

    # set the irrelevant output results to "missing"
    batch_dims_ids = 1:batch_dims
    colon_prod = Iterators.product((:).(Ref(1), size_samples[batch_dims_ids])...)
    filter_missing = vcat(
        extend_cartsian_index.(
            dropdims(
                (:).(N_configs .+ 1, Ref(N_max_configs));
                dims=total_dims
            ),
            colon_prod)...
    )
    configs[filter_missing, :] .= missing
    vals[filter_missing] .= missing

    # from diagonal terms
    vals[1, fill(:, batch_dims)...] .+=
        0.5 * operator.U * sum(xs .* (xs .- 1); dims=total_dims) .+
        operator.V * sum(wsᵢⱼ .* (nᵢ .* nⱼ); dims=total_dims) .-
        operator.μ * sum(xs; dims=total_dims)

    # from off-diagonal terms (bᵢ^† bⱼ)
    cᵢⱼ_ids = sort(findall(maskᵢⱼ), by=x -> x.I[batch_dims_ids])
    filterᵢⱼ = vcat(
        extend_cartsian_index.(
            dropdims(
                (:).(Ref(2), 1 .+ Nᵢⱼ);
                dims=total_dims
            ),
            colon_prod)...
    )
    filtered_configsᵢⱼ = view.(Ref(configs), filterᵢⱼ, Ref(:))
    @views vals[filterᵢⱼ] .+=
        @. -operator.J * ((sqrt(nⱼ)*sqrt(nᵢ + 1)*wsᵢⱼ)*maskᵢⱼ)[cᵢⱼ_ids]
    eachslice_maskᵢⱼ = eachslice(maskᵢⱼ, dims=tuple(batch_dims_ids...))
    add_at!.(
        Ref(operator), filtered_configsᵢⱼ, vcat(getindex.(Ref(i), eachslice_maskᵢⱼ)...);
        addend=1
    )
    add_at!.(
        Ref(operator), filtered_configsᵢⱼ, vcat(getindex.(Ref(j), eachslice_maskᵢⱼ)...);
        addend=-1
    )

    # from off-diagonal terms (bⱼ^† bᵢ)
    cⱼᵢ_ids = sort(findall(maskⱼᵢ), by=x -> x.I[batch_dims_ids])
    filterⱼᵢ = vcat(
        extend_cartsian_index.(
            dropdims(
                (:).(2 .+ Nᵢⱼ, 1 .+ Nᵢⱼ .+ Nⱼᵢ);
                dims=batch_dims + 1
            ),
            colon_prod)...
    )
    filtered_configsⱼᵢ = view.(Ref(configs), filterⱼᵢ, Ref(:))
    @views vals[filterⱼᵢ] .+=
        @. -operator.J * ((sqrt(nᵢ)*sqrt(nⱼ + 1)*wsⱼᵢ)*maskⱼᵢ)[cⱼᵢ_ids]
    eachslice_maskⱼᵢ = eachslice(maskⱼᵢ, dims=tuple(batch_dims_ids...))
    add_at!.(
        Ref(operator), filtered_configsⱼᵢ, vcat(getindex.(Ref(i), eachslice_maskⱼᵢ)...);
        addend=-1
    )
    add_at!.(
        Ref(operator), filtered_configsⱼᵢ, vcat(getindex.(Ref(j), eachslice_maskⱼᵢ)...);
        addend=1
    )

    return configs, vals
end
