"""
    StochasticReconfiguration(; diag_shift=0.01, diag_scale=0.0, solver=CholeskySolver(),
                                mode=:auto, holomorphic=false) <: AbstractPreconditioner

Stochastic reconfiguration, also known as the natural gradient or imaginary-time projection.

Plain gradient descent follows the steepest direction in *parameter* space, which is the wrong
geometry: two parameters can be rescaled against each other without changing the wavefunction.
Stochastic reconfiguration follows the steepest direction in *state* space instead, by
preconditioning the gradient with the quantum geometric tensor

```math
S_{kk'} = \\mathrm{Re} \\left[ \\langle O_k^* O_{k'} \\rangle
          - \\langle O_k^* \\rangle \\langle O_{k'} \\rangle \\right]
```

and solving `S δ = ∇E` for the update direction.

!!! warning "Use a smaller step than for plain descent"
    The preconditioned update is routinely one to two orders of magnitude larger than the raw
    gradient, so a learning rate that is fine for plain gradient descent will overshoot here
    and can leave the run stalled at a non-eigenstate. Start around `0.05` with
    `diag_shift = 0.01` and reduce if the energy plateaus above the ground state while its
    variance stays large.

This is what rescues the ordered regime, where plain descent stalls: the gradient carries a
factor of the Born probability `p(s)`, which vanishes for exactly the configurations whose
amplitude needs to grow, and `S` carries the same factor.

# Fields
- `diag_shift`: absolute regularization `λ` added to the diagonal. The geometric tensor is
  routinely singular — redundant parameters and unexplored directions give exact zero modes —
  so some regularization is mandatory.
- `diag_scale`: regularization relative to each diagonal entry, tracking the tensor's own
  magnitude as it changes during a run.
- `solver`: an [`AbstractLinearSolver`](@ref).
- `mode`: `:sr` builds the `P × P` matrix `S`, `:minsr` builds the `N × N` matrix instead (see
  below), `:matrixfree` builds neither, and `:auto` picks whichever of the first two is smaller.
- `holomorphic`: passed through to `log_derivatives`.
- `chunk_size`: passed through to `local_estimators`, bounding the memory of the
  log-derivative differentiation pass. `nothing` disables it.

# `:matrixfree`

Both `:sr` and `:minsr` form a square matrix — `P × P` or `2N × 2N` — and for a real network
that matrix, not the sampling, exhausts memory first. `:matrixfree` wraps the design matrix in a
[`QuantumGeometricTensor`](@ref) whose only operation is multiplication and hands that to an
iterative solver, so nothing square is allocated. It requires
[`ConjugateGradientSolver`](@ref); a direct solver has nothing to factorize.

# The two forms are the same update

With `X` the centered, probability-weighted log-derivative matrix split into real and imaginary
parts and stacked, `S = XᵀX` and `∇E = 2Xᵀε`. The identity

```math
(X^T X + \\lambda I)^{-1} X^T = X^T (X X^T + \\lambda I)^{-1}
```

means the update can come either from the `P × P` matrix `XᵀX` or from the `N × N` matrix
`XXᵀ`. They agree exactly, so the choice is one of cost: with far more parameters than samples,
the usual case for a network, the second is dramatically cheaper. That is the kernel trick,
known here as MinSR or SRt.
"""
struct StochasticReconfiguration{S<:AbstractLinearSolver,C} <: AbstractPreconditioner
    diag_shift::Float64
    diag_scale::Float64
    solver::S
    mode::Symbol
    holomorphic::Bool
    chunk_size::C

    function StochasticReconfiguration(;
        diag_shift::Real=0.01, diag_scale::Real=0.0, solver::S=CholeskySolver(),
        mode::Symbol=:auto, holomorphic::Bool=false, chunk_size::C=nothing
    ) where {S<:AbstractLinearSolver,C}
        mode in (:sr, :minsr, :matrixfree, :auto) || throw(ArgumentError(
            "mode must be :sr, :minsr, :matrixfree or :auto, got :$mode"
        ))
        diag_shift >= 0 || throw(ArgumentError("diag_shift must be non-negative"))
        diag_scale >= 0 || throw(ArgumentError("diag_scale must be non-negative"))
        return new{S,C}(
            Float64(diag_shift), Float64(diag_scale), solver, mode, holomorphic, chunk_size
        )
    end
end

"""
    _regularize(A, scale) -> A + scale * Diagonal(diag(A))

Sorella's relative regularization, which tracks the tensor's own magnitude as it changes over a
run. The absolute `diag_shift` is applied by [`solve`](@ref), not here.
"""
function _regularize(A::AbstractMatrix, scale::Real)
    iszero(scale) && return A
    return A + scale * Diagonal(diag(A))
end

"""
    precondition(sr, state, operator) -> (Stats, update)

The stochastic-reconfiguration update direction, together with the energy it was computed from.

`update` has the same structure as the state's parameters, so a step is
`θ .- η .* update` (or `fmap` over a nested container). Both returned values come from a single
set of samples.
"""
function NQSCore.precondition(
    sr::StochasticReconfiguration, vs::AbstractVariationalState, operator
)
    est = local_estimators(vs, operator; holomorphic=sr.holomorphic, chunk_size=sr.chunk_size)

    # The design matrix and residual both forms of the update are built from. Real and
    # imaginary parts are stacked rather than kept complex, which makes `Re[X†X] = XᵀX`
    # literally true and lets the kernel-trick identity below apply with no special-casing.
    p = est.weights === nothing ? fill(1 / length(est.E), length(est.E)) :
        est.weights ./ sum(est.weights)
    sqrt_p = sqrt.(p)
    Xc = sqrt_p .* centered(est.O, est.weights)
    εc = sqrt_p .* (est.E .- sum(p .* est.E))
    X = vcat(real.(Xc), imag.(Xc))
    ε = vcat(real.(εc), imag.(εc))

    n_rows, n_params = size(X)
    mode = sr.mode === :auto ? (n_rows < n_params ? :minsr : :sr) : sr.mode

    gradient = 2 .* (transpose(X) * ε)

    δ = if mode === :matrixfree
        sr.solver isa ConjugateGradientSolver || throw(ArgumentError(
            "mode=:matrixfree never forms the geometric tensor, so it needs a solver that " *
            "only multiplies by it; got $(typeof(sr.solver)). Use ConjugateGradientSolver, " *
            "or :sr / :minsr for a direct solver."
        ))
        solve(sr.solver, QuantumGeometricTensor(X, sr.diag_scale), gradient, sr.diag_shift)
    elseif mode === :sr
        A = _regularize(transpose(X) * X, sr.diag_scale)
        solve(sr.solver, A, gradient, sr.diag_shift)
    else
        # (XᵀX + λI)⁻¹ Xᵀ b == Xᵀ (XXᵀ + λI)⁻¹ b, exactly.
        A = _regularize(X * transpose(X), sr.diag_scale)
        transpose(X) * solve(sr.solver, A, 2 .* ε, sr.diag_shift)
    end

    stats = est.weights === nothing ?
            statistics(est.E) : weighted_statistics(est.E, est.weights)
    return stats, match_parameter_shape(δ, parameters(vs))
end

"""
    Identity() <: AbstractPreconditioner

No preconditioning: the update is the plain energy gradient.

Present so a driver can treat plain gradient descent and stochastic reconfiguration uniformly,
and so the two can be compared under identical conditions.
"""
struct Identity <: AbstractPreconditioner end

function NQSCore.precondition(::Identity, vs::AbstractVariationalState, operator)
    return expect_and_grad(vs, operator)
end

"""
    optimize!(state, operator, preconditioner; iterations, learning_rate, callback=nothing)
        -> Vector{Stats}

Run a variational optimization, returning the energy at every step.

The smallest useful driver: no logging, no checkpointing, no early stopping — those belong with
the umbrella package. This is enough to test that a preconditioner descends.

`callback(iteration, stats, state)` runs after each step if given.
"""
function optimize!(
    vs::AbstractVariationalState, operator, preconditioner::AbstractPreconditioner;
    iterations::Integer=100, learning_rate::Real=0.05, callback=nothing
)
    history = Stats[]
    for it in 1:iterations
        stats, δ = precondition(preconditioner, vs, operator)
        push!(history, stats)
        setparameters!(vs, fmap((p, d) -> p .- learning_rate .* d, parameters(vs), δ))
        # A Monte Carlo state must redraw: its cached samples came from the old parameters.
        vs isa MCState && resample!(vs)
        callback === nothing || callback(it, stats, vs)
    end
    return history
end
