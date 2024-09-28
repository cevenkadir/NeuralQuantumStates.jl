using NeuralQuantumStates.Hilberts

# Uniform Hilbert space
## Fock
## Spin
## Qubit

build(hilbert_type::Symbol, args...; kwargs...) = build(Val(hilbert_type), args...; kwargs...)

"""
    build(::Val{:Spin}, s::T, N::Integer; ∑Sz::Union{T_Sz,Nothing}=nothing) where
        {T<:Union{Rational,Integer},T_Sz<:Real}
        -> NeuralQuantumStates.Hilberts.FiniteUniformHilbert

Build a finite uniform Hilbert space for a spin-`s` system.

# Arguments
- `s::T`: Spin of the system.
- `N::Integer`: Number of degrees of freedom.
- `∑Sz::Union{T_Sz,Nothing}`: The total magnetization to be conserved. Default is `nothing`.

# Returns
- `NeuralQuantumStates.Hilberts.FiniteUniformHilbert`: The finite uniform Hilbert space for
    a spin-`s` system.
"""
function build(
    ::Val{:Spin}, s::T, N::Integer;
    ∑Sz::Union{T_Sz,Nothing}=nothing
) where {T<:Union{Rational,Integer},T_Sz<:Real}
    n_lDoF = 2s + 1 |> Integer

    lDoF = range(-s, s; step=1)
    lDoF = SVector{length(lDoF)}(lDoF)


    if ∑Sz === nothing
        constraint = nothing
    else
        constraint = SumConstraint{typeof(∑Sz)}(∑Sz)
    end

    return FiniteUniformHilbert{T,N,n_lDoF}(
        lDoF=lDoF,
        constraint=constraint,
        type=:Spin
    )
end

"""
    build(::Val{:Fock}, n_max::T, N::Integer; ∑n::Union{T_n,Nothing}=nothing) where
        {T<:Integer,T_n<:Real}
        -> NeuralQuantumStates.Hilberts.FiniteUniformHilbert

Build a finite uniform Hilbert space for a Fock system.

# Arguments
- `n_max::T`: Allowed maximum number of particles in each lattice site.
- `N::Integer`: Number of degrees of freedom.
- `∑n::Union{T_n,Nothing}`: The total number of particles to be conserved. Default is
    `nothing`.

# Returns
- `NeuralQuantumStates.Hilberts.FiniteUniformHilbert`: The finite uniform Hilbert space for
    a Fock system.
"""
function build(
    ::Val{:Fock}, n_max::T, N::Integer;
    ∑n::Union{T_n,Nothing}=nothing
) where {T<:Integer,T_n<:Real}
    lDoF = range(0, n_max; step=1)

    n_lDoF = length(lDoF)

    lDoF = SVector{n_lDoF}(lDoF)

    if ∑n === nothing
        constraint = NoHilbertConstraint()
    else
        constraint = SumConstraint{typeof(∑n)}(∑n)
    end

    return FiniteUniformHilbert{T,N,n_lDoF}(
        lDoF=lDoF,
        constraint=constraint,
        type=:Fock
    )
end

"""
    build(::Val{:Qubit}, N::Integer) where {T<:Integer}
        -> NeuralQuantumStates.Hilberts.FiniteUniformHilbert

Build a finite uniform Hilbert space for a qubit system.

# Arguments
- `N::Integer`: Number of degrees of freedom.

# Returns
- `NeuralQuantumStates.Hilberts.FiniteUniformHilbert`: The finite uniform Hilbert space for
    a qubit system.
"""
function build(::Val{:Qubit}, N::Integer) where {T<:Integer}
    lDoF = range(0, 1; step=1)
    lDoF = SVector{length(lDoF)}(lDoF)

    n_lDoF = 2

    constraint = nothing

    return FiniteUniformHilbert{T,N,n_lDoF}(
        lDoF=lDoF,
        constraint=constraint,
        type=:Qubit
    )
end


Return all states in the given finite uniform Hilbert space.

# Arguments
- `hilbert::NeuralQuantumStates.Hilberts.FiniteUniformHilbert{T,N_DoF,N_lDoF}`: The finite uniform Hilbert
    space.

# Returns
- `Vector{NTuple{N_DoF,T}}`: All states in the given finite uniform Hilbert space.
"""
function all_states(
    hilbert::Hilberts.FiniteUniformHilbert{T,N_DoF,N_lDoF}
) where {T<:Real,N_DoF,N_lDoF}
    i = Iterators.product([hilbert.lDoF for _ in 1:N_DoF]...)
    return Iterators.filter(x -> Hilberts.check(hilbert.constraint, x), i) |> collect |> unique
end

"""
    all_states(hilbert::Hilberts.InfiniteUniformHilbert{T,N_DoF}) where
        {T<:Real,N_DoF}

Throws an `OverflowError` since the Hilbert space is infinite.

# Arguments
- `hilbert::NeuralQuantumStates.Hilberts.InfiniteUniformHilbert{T,N_DoF}`: The infinite uniform Hilbert space.
"""
function all_states(hilbert::Hilberts.InfiniteUniformHilbert{T,N_DoF}) where {T<:Real,N_DoF}
    throw(OverflowError("Infinite Hilbert space"))
end

"""
    all_states(hilbert::NeuralQuantumStates.Hilberts.CompositeUniformHilbert{N,N_DoF}) where
        {N,N_DoF} -> Vector{NTuple}

Return all states in the given composite uniform Hilbert space.

# Arguments
- `hilbert::NeuralQuantumStates.Hilberts.CompositeUniformHilbert{N,N_DoF}`: The composite
    uniform Hilbert space.

# Returns
- `Vector{NTuple}`: All states in the given composite uniform Hilbert space.
"""
function all_states(composite_hilbert::CompositeUniformHilbert{N,N_DoF}) where {N,N_DoF}
    i = Iterators.product([all_states(hilbert) for hilbert in composite_hilbert.hilberts]...)
    return Iterators.filter(x -> Hilberts.check(composite_hilbert.composite_constraint, x), i) |> collect |> unique
end

"""
    ⊗(
        (hilberts::NeuralQuantumStates.Hilberts.UniformHilbert)...;
        composite_constraint::NeuralQuantumStates.Hilberts.AbstractCompositeDiscreteHilbertConstraint=NeuralQuantumStates.Hilberts.NoCompositeDiscreteHilbertConstraint()
    ) -> NeuralQuantumStates.Hilberts.CompositeUniformHilbert

Return a composite uniform Hilbert space by taking the tensor product of the given uniform
    Hilbert spaces.

# Arguments
- `hilberts::NeuralQuantumStates.Hilberts.UniformHilbert`: The uniform Hilbert spaces to be
    tensor producted.
- `composite_constraint::NeuralQuantumStates.Hilberts.AbstractCompositeDiscreteHilbertConstraint`:
    The composite constraint to be applied. Default is
    `NeuralQuantumStates.Hilberts.NoCompositeDiscreteHilbertConstraint()`.

# Returns
- `NeuralQuantumStates.Hilberts.CompositeUniformHilbert`: The generated composite uniform
    Hilbert space.
"""
function ⊗((hilberts::Hilberts.UniformHilbert)...; composite_constraint::Hilberts.AbstractCompositeDiscreteHilbertConstraint=Hilberts.NoCompositeDiscreteHilbertConstraint())
    total_N_DoF = sum(Hilberts.n_DoF.(hilberts))

    return Hilberts.CompositeUniformHilbert{length(hilberts),total_N_DoF}(hilberts=hilberts, composite_constraint=composite_constraint)
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
- `Vector{NTuple}`: A vector of states corresponding to the given indices in the state
    representation.
"""
function state_index_to_state(hilbert::Hilberts.UniformHilbert{N_DoF}, indices::AbstractVector{Int}) where {N_DoF}
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
- `NTuple`: A state corresponding to the given index in the state representation.
"""
function state_index_to_state(hilbert::Hilberts.UniformHilbert{N_DoF}, index::Int) where {N_DoF}
    return all_states(hilbert)[index]
end

"""
    state_to_state_index(
        hilbert::NeuralQuantumStates.Hilberts.UniformHilbert{N_DoF},
        states::AbstractVector{NTuple{N_DoF,T}}
    ) where {N_DoF,T<:Real} -> Vector{Integer}

Return the indices of the given states in the given uniform Hilbert space.

# Arguments
- `hilbert::NeuralQuantumStates.Hilberts.UniformHilbert{N_DoF}`: The uniform Hilbert space.
- `states::AbstractVector{NTuple{N_DoF,T}}`: The states to be returned in the index
    representation.

# Returns
- `Vector{Integer}`: A vector of indices of the given states in the index representation.
"""
function state_to_state_index(hilbert::Hilberts.UniformHilbert{N_DoF}, states::AbstractVector{NTuple{N_DoF,T}}) where {N_DoF,T<:Real}
    allstates = all_states(hilbert)
    return [findfirst(x -> x == state, allstates) for state in states]
end

"""
    state_to_state_index(
        hilbert::NeuralQuantumStates.Hilberts.UniformHilbert{N_DoF},
        state::NTuple{N_DoF,T}
    ) where {N_DoF,T<:Real} -> Integer

Return the index of the given state in the given uniform Hilbert space.

# Arguments
- `hilbert::NeuralQuantumStates.Hilberts.UniformHilbert{N_DoF}`: The uniform Hilbert space.
- `state::NTuple{N_DoF,T}`: The state to be returned in the index representation.

# Returns
- `Integer`: The index of the given state in the index representation.
"""
function state_to_state_index(hilbert::Hilberts.UniformHilbert{N_DoF}, state::NTuple{N_DoF,T}) where {N_DoF,T<:Real}
    return state_to_state_index(hilbert, [state]).first
end

#! WRITE random_state function
