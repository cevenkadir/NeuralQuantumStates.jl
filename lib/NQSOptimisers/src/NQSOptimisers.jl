"""
    NQSOptimisers

Stochastic reconfiguration and natural-gradient preconditioning for variational quantum states.

This corresponds to NetKet's `netket/_src/ngd` together with `netket/_src/solvers`.

# Why preconditioning is not optional

Plain gradient descent follows the steepest direction in *parameter* space. That is the wrong
geometry for a wavefunction: rescaling two parameters against each other can leave the state
untouched while completely changing the gradient. [`StochasticReconfiguration`](@ref) instead
follows the steepest direction in *state* space, using the quantum geometric tensor as the
metric.

The practical consequence is stark. In the ordered regime of a transverse-field Ising chain,
plain descent stalls at the classical energy, because the gradient carries a factor of the Born
probability `p(s)` that vanishes for exactly the configurations whose amplitude must grow. The
geometric tensor carries the same factor, and dividing by it undoes the suppression.

# What is here
- [`StochasticReconfiguration`](@ref) — in both the `P × P` and the kernel-trick `N × N` form,
  which give exactly the same update.
- [`Identity`](@ref) — no preconditioning, so plain descent and SR can be compared under
  identical conditions.
- Solvers: [`CholeskySolver`](@ref), [`PseudoInverseSolver`](@ref),
  [`ConjugateGradientSolver`](@ref) — a real choice, because the geometric tensor is routinely
  singular.
- [`optimize!`](@ref) — a minimal driver.

# The solver seam

Linear solves go through [`AbstractLinearSolver`](@ref) rather than a hardcoded call, which is
where a GPU solver attaches as a package extension. `ConjugateGradientSolver` needs only
matrix-vector products, so paired with [`QuantumGeometricTensor`](@ref) it never asks for the
square matrix at all — which is what makes the geometric tensor of a real network tractable, on
a device or off one.

A note for anyone reading an earlier version of this file: the seam was once described as the
place a *batched* solver such as BatchSolve.jl would attach. That was a misreading. BatchSolve
solves many small independent problems at once; stochastic reconfiguration solves one large
Hermitian system per step, and has nothing to batch over.

# Not yet here

Time evolution (TDVP) uses the same geometric tensor with an imaginary or real time step. It is
deliberately left out rather than half-implemented.

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
