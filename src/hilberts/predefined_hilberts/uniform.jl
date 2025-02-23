using NeuralQuantumStates.Hilberts
using NeuralQuantumStates.Extras
using Random

function build(hilbert_type::Symbol, args...; kwargs...)
    return build(Val(hilbert_type), args...; kwargs...)
end

"""
    build(
        ::Val{:Spin}, s::T, N::Integer;
        ∑Sz::Union{T_Sz,Nothing}=nothing, array_type::Type=Array
    ) where {T<:Union{Rational,Integer},T_Sz<:Real}
        -> NeuralQuantumStates.Hilberts.FiniteUniformHilbert

Build a finite uniform Hilbert space for a spin-`s` system.

# Arguments
- `s::T`: Spin of the system.
- `N::Integer`: Number of degrees of freedom.
- `∑Sz::Union{T_Sz,Nothing}`: The total magnetization to be conserved. Default is `nothing`.
- `array_type::Type`: The Vector type to be used in each state. Default is `Array`.

# Returns
- `NeuralQuantumStates.Hilberts.FiniteUniformHilbert`: The finite uniform Hilbert space for
    a spin-`s` system.
"""
function build(
    ::Val{:Spin}, s::T, N::Integer;
    ∑Sz::Union{T_Sz,Nothing}=nothing, array_type::Type=Array
) where {T<:Union{Rational,Integer},T_Sz<:Real}
    n_lDoF = 2s + 1 |> Integer

    lDoF = range(-s, s; step=1)
    lDoF = SVector{length(lDoF)}(lDoF)


    if ∑Sz === nothing
        constraint = NoDiscreteHilbertConstraint()
    else
        constraint = SumConstraint{T_Sz}(∑Sz)
    end

    return FiniteUniformHilbert{T,array_type,N,n_lDoF}(
        lDoF=lDoF,
        constraint=constraint,
        type=:Spin
    )
end

"""
    build(
        ::Val{:Fock}, n_max::T, N::Integer;
        ∑n::Union{T_n,Nothing}=nothing, array_type::Type=Array
    ) where {T<:Integer,T_n<:Real} -> NeuralQuantumStates.Hilberts.FiniteUniformHilbert

Build a finite uniform Hilbert space for a Fock system.

# Arguments
- `n_max::T`: Allowed maximum number of particles in each lattice site.
- `N::Integer`: Number of degrees of freedom.
- `∑n::Union{T_n,Nothing}`: The total number of particles to be conserved. Default is
    `nothing`.
- `array_type::Type`: The Vector type to be used in each state. Default is `Array`.

# Returns
- `NeuralQuantumStates.Hilberts.FiniteUniformHilbert`: The finite uniform Hilbert space for
    a Fock system.
"""
function build(
    ::Val{:Fock}, n_max::T, N::Integer;
    ∑n::Union{T_n,Nothing}=nothing, array_type::Type=Array
) where {T<:Integer,T_n<:Real}
    lDoF = range(0, n_max; step=1)

    n_lDoF = length(lDoF)

    lDoF = SVector{n_lDoF}(lDoF)

    if ∑n === nothing
        constraint = NoDiscreteHilbertConstraint()
    else
        constraint = SumConstraint{T_n}(∑n)
    end

    return FiniteUniformHilbert{T,array_type,N,n_lDoF}(
        lDoF=lDoF,
        constraint=constraint,
        type=:Fock
    )
end

"""
    build(::Val{:Qubit}, N::Integer, array_type::Type=Array) where {T<:Integer}
        -> NeuralQuantumStates.Hilberts.FiniteUniformHilbert

Build a finite uniform Hilbert space for a qubit system.

# Arguments
- `N::Integer`: Number of degrees of freedom.
- `array_type::Type`: The Vector type to be used in each state. Default is `Array`.

# Returns
- `NeuralQuantumStates.Hilberts.FiniteUniformHilbert`: The finite uniform Hilbert space for
    a qubit system.
"""
function build(::Val{:Qubit}, N::Integer; array_type::Type=Array)
    lDoF = range(0, 1; step=1)
    lDoF = SVector{length(lDoF)}(lDoF)

    n_lDoF = 2

    constraint = NoDiscreteHilbertConstraint()

    return FiniteUniformHilbert{eltype(lDoF),array_type,N,n_lDoF}(
        lDoF=lDoF,
        constraint=constraint,
        type=:Qubit
    )
end

"""
    all_states(
        hilbert::NeuralQuantumStates.Hilberts.FiniteUniformHilbert{T,T_Array,N_DoF,N_lDoF}
    ) where {T<:Real,T_Array<:AbstractArray,N_DoF,N_lDoF} -> Vector{Vector}

Return all states in the given finite uniform Hilbert space.

# Arguments
- `hilbert::NeuralQuantumStates.Hilberts.FiniteUniformHilbert{T,T_Array,N_DoF,N_lDoF}`: The
    finite uniform Hilbert space.

# Returns
- `Vector{T_Array}`: All states in the given finite uniform Hilbert space.
"""
function all_states(
    hilbert::Hilberts.FiniteUniformHilbert{T,T_Array,N_DoF,N_lDoF}
) where {T<:Real,T_Array<:AbstractArray,N_DoF,N_lDoF}
    array_type_params = get_array_type_params(T_Array, T, (N_DoF,))
    T_Array_w_params = T_Array{array_type_params...}
    combinations = T_Array_w_params[]

    for state in Iterators.product([hilbert.lDoF for _ in 1:N_DoF]...)
        state = T_Array_w_params(state |> collect)
        if Hilberts.check(hilbert.constraint, state)
            push!(combinations, state)
        end
    end

    unique!(combinations)

    return combinations
end

"""
    all_states(hilbert::Hilberts.InfiniteUniformHilbert{T,T_Array,N_DoF}) where
        {T<:Real,T_Array<:AbstractArray,N_DoF}

Throws an `OverflowError` since the Hilbert space is infinite.

# Arguments
- `hilbert::NeuralQuantumStates.Hilberts.InfiniteUniformHilbert{T,T_Array,N_DoF}`: The
    infinite uniform Hilbert space.
"""
function all_states(hilbert::Hilberts.InfiniteUniformHilbert{T,T_Array,N_DoF}) where
{T<:Real,T_Array<:AbstractArray,N_DoF}
    throw(OverflowError("Infinite Hilbert space"))
end

"""
    all_states(hilbert::NeuralQuantumStates.Hilberts.CompositeUniformHilbert{T_C_Array,N,N_DoF}) where
        {T_C_Array<:AbstractArray,N,N_DoF} -> Vector{NTuple}

Return all states in the given composite uniform Hilbert space.

# Arguments
- `hilbert::NeuralQuantumStates.Hilberts.CompositeUniformHilbert{T_C_Array,N,N_DoF}`: The
    composite uniform Hilbert space.

# Returns
- `Vector{T_C_Array}`: All states in the given composite uniform Hilbert space.
"""
function all_states(composite_hilbert::CompositeUniformHilbert{T_C_Array,N,N_DoF}) where
{T_C_Array<:AbstractArray,N,N_DoF}
    c_array_type_params = get_array_type_params(
        T_C_Array,
        Union{[
            begin
                T_Arrayᵢ = typeof(hilbertᵢ).parameters[2]
                Tᵢ = typeof(hilbertᵢ).parameters[1]
                N_DoFᵢ = typeof(hilbertᵢ).parameters[3]
                array_type_paramsᵢ = get_array_type_params(T_Arrayᵢ, Tᵢ, (N_DoFᵢ,))
                T_Arrayᵢ{array_type_paramsᵢ...}
            end
            for hilbertᵢ in composite_hilbert.hilberts
        ]...},
        (length(composite_hilbert.hilberts),)
    )
    T_C_Array_w_params = T_C_Array{c_array_type_params...}
    combinations = T_C_Array_w_params[]

    for state in Iterators.product(
        [all_states(hilbert) for hilbert in composite_hilbert.hilberts]...
    )
        state = T_C_Array_w_params(state |> collect)

        if Hilberts.check(composite_hilbert.constraint, state) &&
           all(
            Hilberts.check(local_hilbert.constraint, state[id_hilbert])
            for (id_hilbert, local_hilbert) in enumerate(composite_hilbert.hilberts)
        )
            push!(combinations, state)
        end
    end

    unique!(combinations)

    return combinations
end

"""
    ⊗(
        (hilberts::NeuralQuantumStates.Hilberts.UniformHilbert)...;
        constraint::NeuralQuantumStates.Hilberts.AbstractCompositeDiscreteHilbertConstraint=NeuralQuantumStates.Hilberts.NoCompositeDiscreteHilbertConstraint(),
        array_type::Type=Array
    ) -> NeuralQuantumStates.Hilberts.CompositeUniformHilbert

Return a composite uniform Hilbert space by taking the tensor product of the given uniform
    Hilbert spaces.

# Arguments
- `hilberts::NeuralQuantumStates.Hilberts.UniformHilbert`: The uniform Hilbert spaces to be
    tensor producted.
- `constraint::NeuralQuantumStates.Hilberts.AbstractCompositeDiscreteHilbertConstraint`:
    The composite constraint to be applied. Default is
    `NeuralQuantumStates.Hilberts.NoCompositeDiscreteHilbertConstraint()`.
- `array_type::Type`: The Vector type to be used in each state. Default is `Array`.

# Returns
- `NeuralQuantumStates.Hilberts.CompositeUniformHilbert`: The generated composite uniform
    Hilbert space.
"""
function ⊗(
    (hilberts::Hilberts.UniformHilbert)...;
    constraint::Hilberts.AbstractCompositeDiscreteHilbertConstraint=Hilberts.NoCompositeDiscreteHilbertConstraint(),
    array_type::Type=Array
)
    N_hilberts = length(hilberts)
    total_N_DoF = sum(Hilberts.n_DoF.(hilberts))

    return Hilberts.CompositeUniformHilbert{array_type,N_hilberts,total_N_DoF}(
        hilberts=hilberts, constraint=constraint
    )
end

"""
    state_index_to_state(
        hilbert::NeuralQuantumStates.Hilberts.UniformHilbert{N_DoF},
        indices::AbstractVector{Int}
    ) where {N_DoF} -> Vector{NTuple}

Return the states corresponding to the given indices in the given uniform Hilbert space.

# Arguments
- `hilbert::NeuralQuantumStates.Hilberts.UniformHilbert{N_DoF}`: The uniform Hilbert space.
- `indices::AbstractVector{Int}`: The indices of the states to be returned in the state
    representation.

# Returns
- `Vector{Vector}`: A vector of states corresponding to the given indices in the state
    representation.
"""
function state_index_to_state(
    hilbert::Hilberts.UniformHilbert{N_DoF}, indices::AbstractVector{Int}
) where {N_DoF}
    return all_states(hilbert)[indices]
end

"""
    state_index_to_state(
        hilbert::NeuralQuantumStates.Hilberts.UniformHilbert{N_DoF}, index::Int
    ) where {N_DoF} -> NTuple

Return the state corresponding to the given index in the given uniform Hilbert space.

# Arguments
- `hilbert::NeuralQuantumStates.Hilberts.UniformHilbert{N_DoF}`: The uniform Hilbert space.
- `index::AbstractVector{Int}`: The index of the state to be returned in the state
    representation.

# Returns
- `Array`: A state corresponding to the given index in the state representation.
"""
function state_index_to_state(
    hilbert::Hilberts.UniformHilbert{N_DoF}, index::Int
) where {N_DoF}
    return all_states(hilbert)[index]
end

"""
    state_to_state_index(
        hilbert::NeuralQuantumStates.Hilberts.UniformHilbert{N_DoF},
        states::AbstractVector{<:AbstractVector{T}}
    ) where {N_DoF,T<:Real} -> Vector{Integer}

Return the indices of the given states in the given uniform Hilbert space.

# Arguments
- `hilbert::NeuralQuantumStates.Hilberts.UniformHilbert{N_DoF}`: The uniform Hilbert space.
- `states::AbstractVector{<:AbstractVector{T}}`: The states to be returned in the index
    representation.

# Returns
- `Vector{Integer}`: A vector of indices of the given states in the index representation.
"""
function state_to_state_index(
    hilbert::Hilberts.UniformHilbert{N_DoF}, states::AbstractVector{<:AbstractVector{T}}
) where {N_DoF,T<:Real}
    allstates = all_states(hilbert)
    return [findfirst(x -> x == state, allstates) for state in states]
end

"""
    state_to_state_index(
        hilbert::NeuralQuantumStates.Hilberts.UniformHilbert{N_DoF},
        state::AbstractVector{T}
    ) where {N_DoF,T<:Real} -> Integer

Return the index of the given state in the given uniform Hilbert space.

# Arguments
- `hilbert::NeuralQuantumStates.Hilberts.UniformHilbert{N_DoF}`: The uniform Hilbert space.
- `state::AbstractVector{T}`: The state to be returned in the index representation.

# Returns
- `Integer`: The index of the given state in the index representation.
"""
function state_to_state_index(
    hilbert::Hilberts.UniformHilbert{N_DoF}, state::AbstractVector{T}
) where {N_DoF,T<:Real}
    return state_to_state_index(hilbert, [state]).first
end

function random_state(
    hilbert::Hilberts.FiniteUniformHilbert{T,T_Array,N_DoF,N_lDoF}, seed::Int=0
) where {T<:Real,T_Array<:AbstractArray,N_DoF,N_lDoF}
    Random.seed!(seed)
    random_index = rand(1:n_states(hilbert))
    return all_states(hilbert)[random_index]
end

function random_state(
    hilbert::Hilberts.InfiniteUniformHilbert{T,T_Array,N_DoF}, seed::Int=0
) where {T<:Real,T_Array<:AbstractArray,N_DoF}
    throw("Not implemented")
end
