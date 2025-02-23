using Parameters
using StaticArrays

abstract type AbstractHilbertConstraint end
abstract type AbstractDiscreteHilbertConstraint <: AbstractHilbertConstraint end
abstract type DiscreteHilbertConstraint <: AbstractDiscreteHilbertConstraint end

"""
    NoDiscreteHilbertConstraint <: NeuralQuantumStates.Hilberts.AbstractDiscreteHilbertConstraint

A constraint that does not introduce any constraint to the discrete Hilbert space.
"""
struct NoDiscreteHilbertConstraint <: AbstractDiscreteHilbertConstraint end

"""
    check(
        constraint::NeuralQuantumStates.Hilberts.NoDiscreteHilbertConstraint,
        x::AbstractVector{T_lDoF}
    ) where {T_lDoF<:Real} -> Bool

Always return `true` if there is no Hilbert space constraint.

# Arguments
- `constraint::NeuralQuantumStates.Hilberts.NoDiscreteHilbertConstraint`: A
    `NeuralQuantumStates.Hilberts.NoDiscreteHilbertConstraint` object.
- `x::AbstractVector{T_lDoF}`: A state to be checked.

# Returns
- `Bool`: Always `true`.
"""
function check(
    constraint::NoDiscreteHilbertConstraint, x::AbstractVector{T_lDoF}
) where {T_lDoF<:Real}
    return true
end

"""
    SumConstraint{T<:Real}(sum_value::T)
        <: NeuralQuantumStates.Hilberts.DiscreteHilbertConstraint
A constraint introducing that the sum of the elements of a given state is equal to a given
    value.

# Fields
- `sum_value::T`: The value that the sum of the elements of a given state must be equal to.
"""
@with_kw struct SumConstraint{T<:Real} <: DiscreteHilbertConstraint
    sum_value::T
end

"""
    check(
        constraint::NeuralQuantumStates.Hilberts.SumConstraint{T}, x::AbstractVector{T_lDoF}
    ) where {T<:Real,T_lDoF<:Real} -> Bool

Check if the sum of the elements of a given state is equal to the given sum value.

# Arguments
- `constraint::NeuralQuantumStates.Hilberts.SumConstraint{T}`: A
    `NeuralQuantumStates.Hilberts.SumConstraint` object.
- `x::AbstractVector{T_lDoF}`: A state to be checked.

# Returns
- `Bool`: `true` if the sum of the elements of the given state is equal to the given sum
    value.
"""
function check(
    constraint::SumConstraint{T}, x::AbstractVector{T_lDoF}
) where {T<:Real,T_lDoF<:Real}
    return sum(x) ≈ constraint.sum_value
end

abstract type AbstractCompositeHilbertConstraint end
abstract type AbstractCompositeDiscreteHilbertConstraint <: AbstractCompositeHilbertConstraint end
abstract type CompositeDiscreteHilbertConstraint <: AbstractCompositeDiscreteHilbertConstraint end

"""
    NoCompositeDiscreteHilbertConstraint
        <: NeuralQuantumStates.Hilberts.AbstractCompositeDiscreteHilbertConstraint

A constraint that does not introduce any constraint to the composite discrete Hilbert space.
"""
struct NoCompositeDiscreteHilbertConstraint <: AbstractCompositeDiscreteHilbertConstraint end

"""
    check(
        constraint::NeuralQuantumStates.Hilberts.NoCompositeDiscreteHilbertConstraint,
        x::AbstractVector{AbstractVector}
    ) -> Bool

Always return `true` if there is no composite discrete Hilbert space constraint.

# Arguments
- `constraint::NeuralQuantumStates.Hilberts.NoCompositeDiscreteHilbertConstraint`: A
    `NeuralQuantumStates.Hilberts.NoCompositeDiscreteHilbertConstraint` object.
- `x::AbstractVector{AbstractVector}`: A composite state to be checked.

# Returns
- `Bool`: Always `true`.
"""
function check(
    constraint::NoCompositeDiscreteHilbertConstraint, x::AbstractVector{<:AbstractVector}
)
    return true
end

"""
    CompositeSumConstraint{T<:Real,N_HS}
        <: NeuralQuantumStates.Hilberts.CompositeDiscreteHilbertConstraint
"""
@with_kw struct CompositeSumConstraint{T<:Real,N_HS} <: CompositeDiscreteHilbertConstraint
    sum_values::NTuple{N_HS,T}
end

"""
    check(
        constraint::NeuralQuantumStates.Hilberts.CompositeSumConstraint{T,N_HS},
        x::AbstractVector{AbstractVector}
    ) where {T<:Real,N_HS} -> Bool

Check if the sum of the elements of a given composite state is equal to the given sum
    values.

# Arguments
- `constraint::NeuralQuantumStates.Hilberts.CompositeSumConstraint{T,N_HS}`: A
    `NeuralQuantumStates.Hilberts.CompositeSumConstraint` object.
- `x::AbstractVector{AbstractVector}`: A composite state to be checked.

# Returns
- `Bool`: `true` if the sum of the elements of the given composite state is equal to the
    given sum values.
"""
function check(
    constraint::CompositeSumConstraint{T,N_HS}, x::AbstractVector{AbstractVector}
) where {T<:Real,N_HS}
    @assert N_HS == length(x)
    return all(sum(x[i]) ≈ constraint.sum_values[i] for i in 1:N_HS)
end

abstract type AbstractHilbert end
abstract type DiscreteHilbert{N_DoF} <: AbstractHilbert end
abstract type UniformHilbert{N_DoF} <: DiscreteHilbert{N_DoF} end

"""
    FiniteUniformHilbert{T<:Real,T_Array<:AbstractArray,N_DoF,N_lDoF}
        <: NeuralQuantumStates.Hilberts.UniformHilbert{N_DoF}

A finite uniform Hilbert space.

# Fields
- `lDoF::SVector{N_lDoF,T}`: Local degrees of freedom.
- `constraint::NeuralQuantumStates.Hilberts.AbstractDiscreteHilbertConstraint`: A constraint
    to the uniform Hilbert space.
- `type::Symbol`: The type of the uniform Hilbert space.
"""
@with_kw struct FiniteUniformHilbert{T<:Real,T_Array<:AbstractArray,N_DoF,N_lDoF} <: UniformHilbert{N_DoF}
    lDoF::SVector{N_lDoF,T}
    constraint::AbstractDiscreteHilbertConstraint
    type::Symbol
    @assert N_DoF isa Integer "number of degrees of freedom must be an integer"
    @assert N_lDoF isa Integer "number of local degrees of freedom must be an integer"
    @assert N_DoF > 0 "number of degrees of freedom must be positive"
    @assert N_lDoF > 0 "number of local degrees of freedom must be positive"
    @assert length(lDoF) == N_lDoF "number of local degrees of freedom must match the given number of local degrees of freedom"
end

"""
    InfiniteUniformHilbert{T<:Real,T_Array<:AbstractArray,N_DoF}
        <: NeuralQuantumStates.Hilberts.UniformHilbert{N_DoF}

An infinite uniform Hilbert space.

# Fields
- `type::Symbol`: The type of the uniform Hilbert space.
"""
@with_kw struct InfiniteUniformHilbert{T<:Real,T_Array<:AbstractArray,N_DoF} <: UniformHilbert{N_DoF}
    type::Symbol
    @assert N_DoF isa Integer "number of degrees of freedom must be an integer"
    @assert N_DoF > 0 "number of degrees of freedom must be positive"
end

"""
    n_DoF(hilbert::NeuralQuantumStates.Hilberts.UniformHilbert{N_DoF}) where {N_DoF}
        -> typeof(N_DoF)

Get the number of degrees of freedom of the given uniform Hilbert space.

# Arguments
- `hilbert::NeuralQuantumStates.Hilberts.UniformHilbert{N_DoF}`: A uniform Hilbert space.

# Returns
- `typeof(N_DoF)`: The number of degrees of freedom of the given uniform Hilbert space.
"""
n_DoF(hilbert::UniformHilbert{N_DoF}) where {N_DoF} = N_DoF

"""
    n_lDoF(hilbert::NeuralQuantumStates.Hilberts.FiniteUniformHilbert{T,T_Array,N_DoF,N_lDoF})
        where {T<:Real,T_Array<:AbstractArray,N_DoF,N_lDoF} -> typeof(N_lDoF)

Get the number of local degrees of freedom of the given finite uniform Hilbert space.

# Arguments
- `hilbert::NeuralQuantumStates.Hilberts.FiniteUniformHilbert{T,T_Array,N_DoF,N_lDoF}`: A finite
    uniform Hilbert space.

# Returns
- `typeof(N_lDoF)`: The number of local degrees of freedom of the given finite uniform
    Hilbert space.
"""
n_lDoF(hilbert::FiniteUniformHilbert{T,T_Array,N_DoF,N_lDoF}) where {T<:Real,T_Array<:AbstractArray,N_DoF,N_lDoF} = N_lDoF

"""
    n_lDoF(hilbert::NeuralQuantumStates.Hilberts.InfiniteUniformHilbert{T,T_Array,N_DoF})
        where {T<:Real,T_Array<:AbstractArray,N_DoF} -> typeof(N_DoF)

Throws an `OverflowError` since the number of local degrees of freedom of an infinite
    uniform Hilbert space is infinite.

# Arguments
- `hilbert::NeuralQuantumStates.Hilberts.InfiniteUniformHilbert{T,T_Array,N_DoF}`: An infinite
    uniform Hilbert space.
"""
function n_lDoF(hilbert::InfiniteUniformHilbert{T,T_Array,N_DoF}) where {T<:Real,T_Array<:AbstractArray,N_DoF}
    throw(OverflowError("Infinite local Hilbert space"))
end

"""
    n_states(hilbert::NeuralQuantumStates.Hilberts.FiniteUniformHilbert{T,N_DoF,N_lDoF})
        where {T<:Real,N_DoF,N_lDoF} -> Integer

Get the number of states of the given finite uniform Hilbert space.

# Arguments
- `hilbert::NeuralQuantumStates.Hilberts.FiniteUniformHilbert{T,N_DoF,N_lDoF}`: A finite
    uniform Hilbert space.

# Returns
- `Integer`: The number of states of the given finite uniform Hilbert space.
"""
function n_states(hilbert::FiniteUniformHilbert{T,T_Array,N_DoF,N_lDoF}) where
{T<:Real,T_Array<:AbstractArray,N_DoF,N_lDoF}
    return length(all_states(hilbert))
end

"""
    n_states(hilbert::NeuralQuantumStates.Hilberts.InfiniteUniformHilbert{T,T_Array,N_DoF})
        where {T<:Real,T_Array<:AbstractArray,N_DoF} -> Integer

Throws an `OverflowError` since the number of states of an infinite uniform Hilbert space
    is infinite.

# Arguments
- `hilbert::NeuralQuantumStates.Hilberts.InfiniteUniformHilbert{T,T_Array,N_DoF}`: An infinite
    uniform Hilbert space.
"""
function n_states(hilbert::InfiniteUniformHilbert{T,T_Array,N_DoF}) where {T<:Real,T_Array<:AbstractArray,N_DoF}
    throw(OverflowError("Infinite Hilbert space"))
end

"""
    CompositeUniformHilbert{T_C_Array<:AbstractArray,N,N_DoF} <:
    NeuralQuantumStates.Hilberts.UniformHilbert{N_DoF}

A composite uniform Hilbert space.

# Fields
- `hilberts::NTuple{N,NeuralQuantumStates.Hilberts.UniformHilbert}`: Uniform Hilbert
    spaces to be composed.
- `constraint::NeuralQuantumStates.Hilberts.AbstractCompositeDiscreteHilbertConstraint`:
    A constraint to the composite uniform Hilbert space.
"""
@with_kw struct CompositeUniformHilbert{T_C_Array<:AbstractArray,N,N_DoF} <:
                UniformHilbert{N_DoF}
    hilberts::NTuple{N,UniformHilbert}
    constraint::AbstractCompositeDiscreteHilbertConstraint
    @assert N_DoF isa Integer "number of degrees of freedom must be an integer"
    @assert N_DoF > 0 "number of degrees of freedom must be positive"
    @assert N_DoF == sum(n_DoF.(hilberts)) "number of degrees of freedom must match the sum of the number of degrees of freedom of the given uniform Hilbert spaces"
end

"""
    n_lDoF(
        hilbert::NeuralQuantumStates.Hilberts.CompositeUniformHilbert{T_C_Array,N,N_DoF},
        dof_index::integer
    ) where {T_C_Array<:AbstractArray,N,N_DoF} -> Integer

Get the number of local degrees of freedom of the given composite uniform Hilbert space at
    the given degree of freedom index.

# Arguments
- `hilbert::NeuralQuantumStates.Hilberts.CompositeUniformHilbert{T_C_Array,N,N_DoF}`: A
    composite uniform Hilbert space.
- `dof_index::Integer`: A degree of freedom index.

# Returns
- `Integer`: The number of local degrees of freedom of the given composite uniform Hilbert
    space at the given degree of freedom index.
"""
function n_lDoF(
    hilbert::CompositeUniformHilbert{T_C_Array,N,N_DoF}, dof_index::Integer
) where {T_C_Array<:AbstractArray,N,N_DoF}
    lengths = n_DoF.(hilbert.hilberts)
    cumsum_lengths = cumsum(lengths)

    in_which_uniform_hilbert = findfirst(x -> dof_index ≤ x, cumsum_lengths)

    return n_lDoF(hilbert.hilberts[in_which_uniform_hilbert])
end
