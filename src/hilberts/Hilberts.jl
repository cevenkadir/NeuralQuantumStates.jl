module Hilberts

export AbstractHilbert
export DiscreteHilbert, UniformHilbert
export FiniteUniformHilbert, InfiniteUniformHilbert, CompositeUniformHilbert

export AbstractHilbertConstraint, AbstractCompositeHilbertConstraint
export AbstractDiscreteHilbertConstraint, AbstractCompositeDiscreteHilbertConstraint
export NoDiscreteHilbertConstraint, NoCompositeDiscreteHilbertConstraint
export DiscreteHilbertConstraint, CompositeDiscreteHilbertConstraint

export SumConstraint

export n_DoF, n_lDoF
export all_states, n_states, ⊗, state_indices_to_states, states_to_state_indices
export random_state

include("hilbert.jl")
include("predefined_hilberts/uniform.jl")

end
