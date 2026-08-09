"""
    NQSOptimisers

Stochastic reconfiguration and natural-gradient preconditioning for variational quantum states.

This corresponds to NetKet's `netket/_src/ngd` together with `netket/_src/solvers`.

Plain gradient descent follows the steepest direction in *parameter* space, which is the wrong
geometry for a wavefunction: rescaling two parameters against each other leaves the state
untouched while changing the gradient. [`StochasticReconfiguration`](@ref) follows the steepest
direction in *state* space instead, with the quantum geometric tensor as the metric. In the
ordered regime of a transverse-field Ising chain that is the difference between converging and
stalling at the classical energy.

# What is here
- [`StochasticReconfiguration`](@ref) — in both the `P × P` and the kernel-trick `N × N` form,
  which give exactly the same update.
- [`Identity`](@ref) — no preconditioning, so plain descent and SR can be compared under
  identical conditions.
- Solvers: [`CholeskySolver`](@ref), [`PseudoInverseSolver`](@ref),
  [`ConjugateGradientSolver`](@ref) — a real choice, since the geometric tensor is routinely
  singular.
- [`optimize!`](@ref) — a minimal driver.

Linear solves go through [`AbstractLinearSolver`](@ref) rather than a hardcoded call, which is
where a GPU solver attaches as a package extension. `ConjugateGradientSolver` needs only
matrix-vector products, so paired with [`QuantumGeometricTensor`](@ref) it never asks for the
square matrix — which is what makes the geometric tensor of a real network tractable.

Time evolution (TDVP) uses the same geometric tensor with a time step, and is deliberately left
out rather than half-implemented.

# Example
```julia
using NQSCore, NQSOptimisers

sr = StochasticReconfiguration(; diag_shift=0.01)
history = optimize!(vs, H, sr; iterations=200, learning_rate=0.05)
```
"""
module NQSOptimisers

using LinearAlgebra
using Functors: fmap
using KrylovKit: linsolve

using NQSCore
using NQSCore: AbstractPreconditioner, AbstractVariationalState, MCState, Stats
using NQSCore: expect_and_grad, local_estimators, match_parameter_shape, parameters
using NQSCore: precondition, resample!, setparameters!, statistics, weighted_statistics

include("solvers.jl")
include("qgt.jl")
include("sr.jl")

export AbstractLinearSolver, solve
export CholeskySolver, PseudoInverseSolver, ConjugateGradientSolver
export QuantumGeometricTensor, to_dense
export StochasticReconfiguration, Identity, precondition, optimize!

end # module NQSOptimisers
