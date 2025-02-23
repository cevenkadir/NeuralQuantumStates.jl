module Operators

export AbstractOperator
export DiscreteOperator

export IsingOperator, LocalOperator

export BoseHubbardOperator

export connected_basis_configs

include("operator.jl")
include("predefined_operators/bosehubbard.jl")
include("predefined_operators/ising.jl")

end
