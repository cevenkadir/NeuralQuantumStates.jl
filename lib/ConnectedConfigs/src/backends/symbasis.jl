# Default state backend: SymBasis packed integers and symmetry-reduced bases.

@inline read_digit(state::BaseInt{T,Ti,B}, position::Integer) where {T,Ti,B} =
    Int(read(state, Ti(position)))

@inline write_digit(state::BaseInt{T,Ti,B}, position::Integer, digit::Integer) where {T,Ti,B} =
    write(state, Ti(position), digit)

fold_state(state, basis::SymBasis.Bases.Basis) = representative(state, basis)

state_lookup(basis::SymBasis.Bases.Basis) =
    Dict{eltype(basis.states),Int}(s => i for (i, s) in pairs(basis.states))

orbit_norms(basis::SymBasis.Bases.Basis) = basis.norms
