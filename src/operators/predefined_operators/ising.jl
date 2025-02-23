using Parameters, Accessors

using LinearAlgebra

using MetaGraphsNext: edges, src, dst, weights

using NeuralQuantumStates.Hilberts
using NeuralQuantumStates.Lattices
using NeuralQuantumStates.Operators
using NeuralQuantumStates.Extras

const DEFAULT_ising_T = Union{Complex,Real}
const DEFAULT_ising_T_h = Rational
const DEFAULT_ising_T_h_Array = AbstractArray
const DEFAULT_ising_Tᵢ = Integer
const DEFAULT_ising_T_l = Real
const DEFAULT_ising_T_p = Number

"""
    TransverseFieldIsingOperator{
        T<:DEFAULT_ising_T,
        T_h<:DEFAULT_ising_T_h,
        T_h_Array<:DEFAULT_ising_T_h_Array,
        Tᵢ<:DEFAULT_ising_Tᵢ,
        T_l<:DEFAULT_ising_T_l,
        T_p<:DEFAULT_ising_T_p,
        N,n_lDoF,D,O
    }(
        hilbert::NeuralQuantumStates.Hilberts.FiniteUniformHilbert{T_h,T_h_Array,N,n_lDoF}
        lattice::NeuralQuantumStates.Lattices.Lattice{Tᵢ,T_l,D,O}
        J::T_p
        h_x::T_p
        h_z::T_p
    ) <: NeuralQuantumStates.Operators.DiscreteOperator{T}

Transverse-field Ising operator instance according to the following formula:
``\\hat{H} = J \\sum_{\\langle i, j \\rangle} \\left( \\hat{\\sigma}^z_i \\hat{\\sigma}^z_j
\\right) + h_x \\sum_i \\hat{\\sigma}^x_i + h_z \\sum_i \\hat{\\sigma}^z_i``

# Fields
- `hilbert::NeuralQuantumStates.Hilberts.FiniteUniformHilbert{T_h,T_h_Array,N,n_lDoF}`: An
    instance of `NeuralQuantumStates.Hilberts.FiniteUniformHilbert`.
- `lattice::NeuralQuantumStates.Lattices.Lattice{Tᵢ,T_l,D,O}`: An instance of
    `NeuralQuantumStates.Lattices.Lattice`.
- `J::T_p`: Coupling constant.
- `h_x::T_p`: Transverse-field strength.
- `h_z::T_p`: Strength of external magnetic field.
"""
@with_kw struct TransverseFieldIsingOperator{
    T<:DEFAULT_ising_T,
    T_h<:DEFAULT_ising_T_h,
    T_h_Array<:DEFAULT_ising_T_h_Array,
    Tᵢ<:DEFAULT_ising_Tᵢ,
    T_l<:DEFAULT_ising_T_l,
    T_p<:DEFAULT_ising_T_p,
    N,n_lDoF,D,O
} <: DiscreteOperator{T}
    hilbert::Hilberts.FiniteUniformHilbert{T_h,T_h_Array,N,n_lDoF}
    lattice::Lattices.Lattice{Tᵢ,T_l,D,O}
    J::T_p
    h_x::T_p
    h_z::T_p
end

"""
    build(
        ::Val{:TransverseFieldIsing},
        hilbert::NeuralQuantumStates.Hilberts.FiniteUniformHilbert{T_h,T_h_Array,N,n_lDoF},
        lattice::NeuralQuantumStates.Lattices.Lattice{Tᵢ,T_l,D,O};
        J::T_p=1.0, h_x::T_p=1.0, h_z::T_p=1.0, T::DataType=Float64
    ) where {
        T_h<:DEFAULT_ising_T_h,
        T_h_Array<:DEFAULT_ising_T_h_Array,
        Tᵢ<:DEFAULT_ising_Tᵢ,
        T_l<:DEFAULT_ising_T_l,
        T_p<:DEFAULT_ising_T_p,
        N,n_lDoF,D,O
    }

Build a Transverse-field Ising operator from the given parameters according to the following
    formula:
``\\hat{H} = J \\sum_{\\langle i, j \\rangle} \\left( \\hat{\\sigma}^z_i \\hat{\\sigma}^z_j
\\right) + h_x \\sum_i \\hat{\\sigma}^x_i + h_z \\sum_i \\hat{\\sigma}^z_i``

# Arguments
- `::Val{:TransverseFieldIsing}`: A value to dispatch to this function.
- `hilbert::NeuralQuantumStates.Hilberts.FiniteUniformHilbert{T_h,T_h_Array,N,n_lDoF}`: An
    instance of `NeuralQuantumStates.Hilberts.FiniteUniformHilbert`.
- `lattice::NeuralQuantumStates.Lattices.Lattice{Tᵢ,T_l,D,O}`: An instance of
    `NeuralQuantumStates.Lattices.Lattice`.

# Keywords
- `J::T_p`: Coupling constant. Defaults to `1.0`.
- `h_x::T_p`: Transverse-field strength. Defaults to `1.0`.
- `h_z::T_p`: Strength of external magnetic field. Defaults to `1.0`.
"""
function build(
    ::Val{:TransverseFieldIsing},
    hilbert::Hilberts.FiniteUniformHilbert{T_h,T_h_Array,N,n_lDoF},
    lattice::Lattices.Lattice{Tᵢ,T_l,D,O};
    J::T_p=1.0, h_x::T_p=1.0, h_z::T_p=1.0, T::DataType=Float64
) where {
    T_h<:DEFAULT_ising_T_h,
    T_h_Array<:DEFAULT_ising_T_h_Array,
    Tᵢ<:DEFAULT_ising_Tᵢ,
    T_l<:DEFAULT_ising_T_l,
    T_p<:DEFAULT_ising_T_p,
    N,n_lDoF,D,O
}
    @assert hilbert.type == :Spin "The Hilbert space must be a spin space."
    @assert hilbert.lDoF[1] == -1 // 2
    @assert hilbert.lDoF[2] == 1 // 2

    return TransverseFieldIsingOperator{T,T_h,T_h_Array,Tᵢ,T_l,T_p,N,n_lDoF,D,O}(
        hilbert=hilbert,
        lattice=lattice,
        J=J,
        h_x=h_x,
        h_z=h_z
    )
end

"""
    connected_basis_configs(
        operator::NeuralQuantumStates.Operators.TransverseFieldIsingOperator{
            T,T_h,T_h_Array,Tᵢ,T_l,T_p,N,n_lDoF,D,O
        },
        sample::AbstractVector{T_h}
    ) where {
        T<:DEFAULT_ising_T,
        T_h<:DEFAULT_ising_T_h,
        T_h_Array<:DEFAULT_ising_T_h_Array,
        Tᵢ<:DEFAULT_ising_Tᵢ,
        T_l<:DEFAULT_ising_T_l,
        T_p<:DEFAULT_ising_T_p,
        N,n_lDoF,D,O
    }

Return the connected basis configurations ``\vert s^\\prime \\rangle`` to a given operator
    `O`, *i.e.* `operator`, for a given configuration ``\vert s \\rangle``, *i.e.* `sample`,
    with their matrix elements ``\\langle s^\\prime \\lvert \\hat{O} \\rvert s \\rangle``.
For this function, `sample` needs to be a vector that contains only one configuration.


# Arguments
- `operator::NeuralQuantumStates.Operators.TransverseFieldIsingOperator{
    T,T_h,T_h_Array,Tᵢ,T_l,T_p,N,n_lDoF,D,O}`: An instance of
    `NeuralQuantumStates.Operators.TransverseFieldIsingOperator`.
- `sample::AbstractVector{T_h}`: A representative vector for the given basis configuration.

# Returns
- `T_h_Array{Union{T_h,Missing}, 2}`: A ``N \\times M`` matrix of ``M`` connected basis
    configurations.
- `Vector{T}`: A vector of ``M`` matrix elements ``\\langle s^\\prime \\lvert \\hat{O}
    \\rvert s \\rangle``.
"""
function connected_basis_configs(
    operator::TransverseFieldIsingOperator{T,T_h,T_h_Array,Tᵢ,T_l,T_p,N,n_lDoF,D,O},
    sample::AbstractVector{T_h}
) where {
    T<:DEFAULT_ising_T,
    T_h<:DEFAULT_ising_T_h,
    T_h_Array<:DEFAULT_ising_T_h_Array,
    Tᵢ<:DEFAULT_ising_Tᵢ,
    T_l<:DEFAULT_ising_T_l,
    T_p<:DEFAULT_ising_T_p,
    N,n_lDoF,D,O
}
    @assert length(sample) == N """
    The sample must have the same length as the number of degrees of freedom.
    """
    @assert Hilberts.check(operator.hilbert.constraint, sample)

    array_type_params = get_array_type_params(T_h_Array, T_h, (N,))
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

    # get the spin quantum numbers
    zᵢ = stack(x[edge.src] for edge in nn_idx_pairs)
    zⱼ = stack(x[edge.dst] for edge in nn_idx_pairs)

    # get the weights
    ws = weights(operator.lattice.metagraph)

    # initialize the configurations and values
    configs₀ = T_h_Array_w_params[copy(x)]
    vals₀ = zeros(T, 1)

    # from diagonal terms
    vals₀ .+= 4 * operator.J * sum(zᵢ .* zⱼ .* getindex.(Ref(ws), i, j))
    vals₀ .+= 2 * operator.h_z * sum(x)

    # from off-diagonal terms
    if iszero(operator.h_x) ||
       operator.hilbert.constraint != Hilberts.NoDiscreteHilbertConstraint()
        configs₁ = T_h_Array_w_params[]
        vals₁ = T[]
    else
        configs₁ = T_h_Array_w_params[@set x[i] *= -1 for i in 1:N]
        vals₁ = fill(T(operator.h_x), N)
    end

    # concatenate the configurations and values
    configs = vcat(configs₀, configs₁)
    vals = vcat(vals₀, vals₁)

    return configs, vals
end

"""
    connected_basis_configs(
        operator::NeuralQuantumStates.Operators.TransverseFieldIsingOperator{
            T,T_h,T_h_Array,Tᵢ,T_l,T_p,N,n_lDoF,D,O
        },
        samples::AbstractVector{<:AbstractVector{T_h}};
        kwargs...
    ) where {
        T<:DEFAULT_ising_T,
        T_h<:DEFAULT_ising_T_h,
        T_h_Array<:DEFAULT_ising_T_h_Array,
        Tᵢ<:DEFAULT_ising_Tᵢ,
        T_l<:DEFAULT_ising_T_l,
        T_p<:DEFAULT_ising_T_p,
        N,n_lDoF,D,O
    }

Return the connected basis configurations ``\vert s^\\prime \\rangle`` to a given operator
    `O`, *i.e.* `operator`, for given configurations ``\vert s \\rangle``, *i.e.* `samples`,
    with their matrix elements ``\\langle s^\\prime \\lvert \\hat{O} \\rvert s \\rangle``.
For this function, `samples` needs to be a vector of vectors that contains the given
    configurations.

# Arguments
- `operator::NeuralQuantumStates.Operators.TransverseFieldIsingOperator{
    T,T_h,T_h_Array,Tᵢ,T_l,T_p,N,n_lDoF,D,O}`: An instance of
    `NeuralQuantumStates.Operators.TransverseFieldIsingOperator`.
- `samples::AbstractVector{<:AbstractVector{T_h}}`: The representative vectors for the given
    basis configurations.

# Returns
- `Vector{Tuple{T_h_Array{Union{T_h,Missing}, 2}, Vector{T}}}`: The ``N \\times M_i``
    matrices of ``M_i`` connected basis configurations and the vectors of `M_i` matrix
    elements ``\\langle s^\\prime \\lvert \\hat{O} \\rvert s \\rangle`` for each the `i`-th
    element of `samples`.
"""
function connected_basis_configs(
    operator::TransverseFieldIsingOperator{T,T_h,T_h_Array,Tᵢ,T_l,T_p,N,n_lDoF,D,O},
    samples::AbstractVector{<:AbstractVector{T_h}};
    kwargs...
) where {
    T<:DEFAULT_ising_T,
    T_h<:DEFAULT_ising_T_h,
    T_h_Array<:DEFAULT_ising_T_h_Array,
    Tᵢ<:DEFAULT_ising_Tᵢ,
    T_l<:DEFAULT_ising_T_l,
    T_p<:DEFAULT_ising_T_p,
    N,n_lDoF,D,O
}
    return connected_basis_configs.(Ref(operator), samples; kwargs...)
end

"""
    connected_basis_configs(
        operator::NeuralQuantumStates.Operators.TransverseFieldIsingOperator{
            T,T_h,T_h_Array,Tᵢ,T_l,T_p,N,n_lDoF,D,O
        },
        samples::AbstractArray{T_h};
    ) where {
        T<:DEFAULT_ising_T,
        T_h<:DEFAULT_ising_T_h,
        T_h_Array<:DEFAULT_ising_T_h_Array,
        Tᵢ<:DEFAULT_ising_Tᵢ,
        T_l<:DEFAULT_ising_T_l,
        T_p<:DEFAULT_ising_T_p,
        N,n_lDoF,D,O
    }

Return the connected basis configurations ``\vert s^\\prime \\rangle`` to a given operator
    `O`, *i.e.* `operator`, for given configurations ``\vert s \\rangle``, *i.e.* `samples`,
    with their matrix elements ``\\langle s^\\prime \\lvert \\hat{O} \\rvert s \\rangle``.
For this function, `samples` needs to be a ``\\dots \\times N`` array (or tensor) contains
    the given configurations.

# Arguments
- `operator::NeuralQuantumStates.Operators.TransverseFieldIsingOperator{
    T,T_h,T_h_Array,Tᵢ,T_l,T_p,N,n_lDoF,D,O}`: An instance of
    `NeuralQuantumStates.Operators.TransverseFieldIsingOperator`.
- `samples::AbstractArray{T_h}`: A `\\dots \\times N` array (or tensor) of the
    representative vectors for the given basis configurations.

# Returns
- `T_h_Array{Union{T_h,Missing}, ndims(samples) + 1}`: A `M_\\mathrm{max} \\times \\dots \\times N` array
    (or tensor) of connected basis configurations where `M_\\mathrm{max}` is the maximum
    of connections for each given configuration. Each element of the constructed
        configurations as extra will be assigned as `missing`.
- `Array{T, ndims(samples)}`: A `M_\\mathrm{max} \\times \\dots` array or (tensor) of matrix
    elements ``\\langle s^\\prime \\lvert \\hat{O} \\rvert s \\rangle``. Each matrix
    elements will be assigned as `missing`.
"""
function connected_basis_configs(
    operator::TransverseFieldIsingOperator{T,T_h,T_h_Array,Tᵢ,T_l,T_p,N,n_lDoF,D,O},
    samples::AbstractArray{T_h};
) where {
    T<:DEFAULT_ising_T,
    T_h<:DEFAULT_ising_T_h,
    T_h_Array<:DEFAULT_ising_T_h_Array,
    Tᵢ<:DEFAULT_ising_Tᵢ,
    T_l<:DEFAULT_ising_T_l,
    T_p<:DEFAULT_ising_T_p,
    N,n_lDoF,D,O
}
    batch_dims = ndims(samples) - 1
    total_dims = batch_dims + 1

    @assert size(samples, total_dims) == N """
      The sample must have the same length as the number of degrees of freedom.
      """
    size_samples = size(samples)

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
    zᵢ = @views xs[fill(:, batch_dims)..., src.(nn_idx_pairs)]
    zⱼ = @views xs[fill(:, batch_dims)..., dst.(nn_idx_pairs)]

    # get the weights
    ws = weights(operator.lattice.metagraph)
    wsᵢⱼ = reshape(getindex.(Ref(ws), i, j), (fill(1, batch_dims)..., N))

    are_there_any_off_diagonal_terms = !iszero(operator.h_x) &&
                                       (operator.hilbert.constraint ==
                                        Hilberts.NoDiscreteHilbertConstraint())
    N_configs = fill(N * are_there_any_off_diagonal_terms, size_samples[1:batch_dims])
    N_configs .+= 1
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

    # from diagonal terms
    vals[1, fill(:, batch_dims)...] .+=
        4 * operator.J * sum(wsᵢⱼ .* (zᵢ .* zⱼ); dims=total_dims) .+
        2 * operator.h_z * sum(xs; dims=total_dims)

    # from off-diagonal terms
    if are_there_any_off_diagonal_terms
        for i in 1:N
            configs[i+1, fill(:, batch_dims)..., i] .*= -1
        end
        vals[2:end, fill(:, batch_dims)...] .+= operator.h_x
    end

    return configs, vals
end
