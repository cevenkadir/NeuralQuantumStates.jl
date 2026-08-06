"""
    NeuralQuantumStates

Umbrella package for the neural-quantum-states stack.

Installing this package gives you the whole ecosystem: it re-exports each component and adds
the predefined models and, in time, the optimization drivers, callbacks, and logging.

# The stack

Tier 1 carries no machine-learning dependencies and is usable on its own for exact
diagonalization:

| Package | Role |
|---|---|
| [SymBasis.jl](https://github.com/cevenkadir/SymBasis.jl) | States, degrees of freedom, symmetry groups, symmetry-reduced bases |
| [OperatorAlgebra.jl](https://github.com/h-mnzlr/OperatorAlgebra.jl) | `Op`/`OpChain`/`OpSum`, sparse and dense conversion, fermionic sites |
| `LatticeSpaceGroups` | Lattice geometry, bonds, and the site permutations symmetry groups need |
| `ConnectedConfigs` | The batched local-energy kernel |

| Package | Role |
|---|---|
| `NQSCore` | Interfaces, `MCState`/`FullSumState`, log-derivatives, statistics |
| `NQSAnsatze` | Ansätze on Lux: RBM, symmetric RBM, Jastrow |
| `NQSSamplers` | Metropolis samplers and transition rules |
| `NQSOptimisers` | Stochastic reconfiguration and MinSR |

This package adds the predefined models and the [`VMC`](@ref) driver on top.

# Predefined models

[`TransverseFieldIsing`](@ref) and [`ExtendedBoseHubbard`](@ref) are model *specifications*;
`build` turns one into a [`Model`](@ref), pairing an `OpSum` Hamiltonian with the space
it acts on.

```julia
using NeuralQuantumStates

lat = build(Hypercube([8]; periodic=true))
model = build(TransverseFieldIsing(lat; J=1.0, h_x=1.0))

states = basis(model).states
res = connected_padded(compile(model.hamiltonian), states)
```

# Migrating from the pre-split package

`Lattices`, `Hilberts`, and `Operators` are gone. Their replacements:

| Was | Now |
|---|---|
| `Lattices.build(:Hypercube, [8], 1.0; periodic=[true])` | `build(Hypercube([8]; periodic=true))` |
| `Hilberts.build(:Spin, 1//2, N)` | `basis(dof_object(Spin(1//2)), N)` |
| `Hilberts.build(:Fock, n_max, N; ∑n=n)` | `basis(dofo, N, sym(TotalBosonicNumber(n, N), dofo))` |
| `Operators.build(:TransverseFieldIsing, h, l; ...)` | `build(TransverseFieldIsing(l; ...))` |
| `Operators.connected_basis_configs(H, samples)` | `connected_padded(H, states)` |

Three behavioural differences are deliberate: connected configurations are padded with a
**zero** matrix element rather than `missing`; the degree-of-freedom axis is always **first**
(the old code put it last for Ising and first for Bose-Hubbard); and two terms reaching the
same configuration now produce two rows rather than one summed row, which every consumer sums
over anyway. See `ConnectedConfigs` for all three, and `compile` for how to stop re-flattening
a Hamiltonian on every call.
"""
module NeuralQuantumStates

using Reexport

# Tier 1: no machine-learning dependencies.
@reexport using LatticeSpaceGroups
@reexport using ConnectedConfigs
@reexport using OperatorAlgebra
@reexport using SymBasis

# Tier 2: the neural-network stack.
@reexport using NQSCore
@reexport using NQSAnsatze
@reexport using NQSSamplers
@reexport using NQSOptimisers

# The driver's `optimizer` argument is an Optimisers.jl rule, so `Descent`, `Adam` and friends
# have to be reachable from `using NeuralQuantumStates` alone.
@reexport using Optimisers

using Optimisers

using NQSCore: AbstractVariationalState, MCState, Stats
using NQSCore: parameters, setparameters!, resample!
using OperatorAlgebra: AbstractOp, Op, OpSum
using SymBasis: Boson, Spin, dof_object
using SymBasis.Bases: basis

# `import`, not `using`: `build` gains methods for model specs alongside the lattice-spec
# methods it already has, so `build` means the same verb throughout the stack. With `using` it
# would instead be shadowed by a new, separate generic here.
import LatticeSpaceGroups: build

include("models.jl")
include("driver.jl")

export AbstractModelSpec, Model
export TransverseFieldIsing, ExtendedBoseHubbard

export VMC, VMCLog, run!, final_energy
export EarlyStopping, InvalidLossStopping

end
