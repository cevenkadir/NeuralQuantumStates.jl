"""
    StochasticReconfiguration(; diag_shift=0.01, diag_scale=0.0, solver=CholeskySolver(),
                                mode=:auto, holomorphic=false) <: AbstractPreconditioner

Stochastic reconfiguration, also known as the natural gradient or imaginary-time projection.

Plain gradient descent follows the steepest direction in *parameter* space, which is the wrong
geometry: two parameters can be scaled arbitrarily relative to one another without changing the
wavefunction at all. Stochastic reconfiguration instead follows the steepest direction in
*state* space, by preconditioning the gradient with the quantum geometric tensor

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

This is what rescues the ordered regime. Plain descent stalls there because the gradient
carries a factor of the Born probability `p(s)`, which vanishes for exactly the configurations
whose amplitude needs to grow; `S` carries the same factor and dividing by it undoes the
suppression.

# Fields
- `diag_shift`: absolute regularization `λ` added to the diagonal. The geometric tensor is
  routinely singular — redundant parameters and unexplored directions both give exact zero
  modes — so some regularization is mandatory rather than optional.
- `diag_scale`: regularization relative to each diagonal entry, which tracks the tensor's own
  magnitude as it changes during a run. See `_regularize`.
- `solver`: an [`AbstractLinearSolver`](@ref).
- `mode`: `:sr` builds the `P × P` matrix `S`, `:minsr` builds the `N × N` matrix instead (see
  below), and `:auto` picks whichever is smaller.
- `holomorphic`: passed through to `log_derivatives`.

# The two forms are the same update

With `X` the centered, probability-weighted log-derivative matrix split into real and imaginary
parts and stacked, `S = XᵀX` and `∇E = 2Xᵀε`. The identity

```math
(X^T X + \\lambda I)^{-1} X^T = X^T (X X^T + \\lambda I)^{-1}
```

means the update can be computed either from the `P × P` matrix `XᵀX` or from the `N × N`
matrix `XXᵀ`. They agree exactly — not approximately — so the choice is purely one of cost.
For a neural network with far more parameters than samples, which is the usual case, the
second is dramatically cheaper. That form is the kernel trick, known in this context as MinSR
or SRt.
"""
struct StochasticReconfiguration{S<:AbstractLinearSolver} <: AbstractPreconditioner
    diag_shift::Float64
    diag_scale::Float64
    solver::S
    mode::Symbol
    holomorphic::Bool

    function StochasticReconfiguration(;
        diag_shift::Real=0.01, diag_scale::Real=0.0, solver::S=CholeskySolver(),
        mode::Symbol=:auto, holomorphic::Bool=false
    ) where {S<:AbstractLinearSolver}
        mode in (:sr, :minsr, :auto) ||
            throw(ArgumentError("mode must be :sr, :minsr or :auto, got :$mode"))
        diag_shift >= 0 || throw(ArgumentError("diag_shift must be non-negative"))
        diag_scale >= 0 || throw(ArgumentError("diag_scale must be non-negative"))
        return new{S}(Float64(diag_shift), Float64(diag_scale), solver, mode, holomorphic)
    end
end

"""
    _regularize(A, shift, scale) -> A + shift*I + scale*Diagonal(diag(A))

Apply both forms of regularization.

`shift` is an **absolute** addition to the diagonal, and `scale` is one **relative** to each
diagonal entry (Sorella's original prescription). A relative shift tracks the tensor's own
magnitude, which changes over the course of an optimization, so a value chosen at the start
does not become negligible or overwhelming later.

`diag_scale` defaults to zero — the absolute shift alone is the common case, and the relative
form is offered rather than imposed. Increasing either one damps the update.

!!! note "Regularization is not a substitute for a sensible step size"
    The geometric tensor is genuinely ill-conditioned — condition numbers of `1e15` are
    ordinary, since redundant parameter directions give exact zero modes — but that does not by
    itself make the update wrong. A run that stalls at a non-eigenstate while the gradient
    stays perfectly ordinary is usually **overshooting**: the preconditioned update can be two
    orders of magnitude larger than the gradient, so a learning rate tuned for plain descent is
    far too large for stochastic reconfiguration. Reach for a smaller step or a larger shift
    before concluding the conditioning is at fault.
"""
function _regularize(A::AbstractMatrix, shift::Real, scale::Real)
    iszero(scale) && return A
    return A + scale * Diagonal(diag(A))
end

"""
    _weighted_design(O, E, weights) -> (X, ε)

The real design matrix `X` and residual `ε` that both forms of the update are built from.

`X` is the centered log-derivative matrix scaled by `sqrt(p)`, with its real and imaginary parts
stacked into a `2N × P` **real** matrix; `ε` is the correspondingly stacked centered local
energy. Stacking rather than working in complex arithmetic is what makes `Re[X†X] = XᵀX`
literally true, so the kernel-trick identity applies without any special-casing.

Centering is not optional: it removes the direction corresponding to a global change of
normalization, which is unphysical and would otherwise be an exact zero mode of the geometric
tensor.
"""
function _weighted_design(O::AbstractMatrix, E::AbstractVector, weights)
    n = length(E)
    p = weights === nothing ? fill(1 / n, n) : weights ./ sum(weights)
    sqrt_p = sqrt.(p)

    Ō = O .- sum(p .* O; dims=1)                  # centered log-derivatives
    Ē = sum(p .* E)

    X = sqrt_p .* Ō
    ε = sqrt_p .* (E .- Ē)

    return vcat(real.(X), imag.(X)), vcat(real.(ε), imag.(ε))
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
    est = local_estimators(vs, operator; holomorphic=sr.holomorphic)
    X, ε = _weighted_design(est.O, est.E, est.weights)

    n_rows, n_params = size(X)
    mode = sr.mode === :auto ? (n_rows < n_params ? :minsr : :sr) : sr.mode

    gradient = 2 .* (transpose(X) * ε)

    δ = if mode === :sr
        A = _regularize(transpose(X) * X, sr.diag_shift, sr.diag_scale)
        solve(sr.solver, A, gradient, sr.diag_shift)
    else
        # (XᵀX + λI)⁻¹ Xᵀ b == Xᵀ (XXᵀ + λI)⁻¹ b, exactly.
        A = _regularize(X * transpose(X), sr.diag_shift, sr.diag_scale)
        transpose(X) * solve(sr.solver, A, 2 .* ε, sr.diag_shift)
    end

    stats = est.weights === nothing ?
            statistics(est.E) : weighted_statistics(est.E, est.weights)
    return stats, match_parameter_shape(δ, parameters(vs))
end

"""
    Identity() <: AbstractPreconditioner

No preconditioning: the update is the plain energy gradient.

Present so that a driver can treat plain gradient descent and stochastic reconfiguration
uniformly, and so that the two can be compared under otherwise identical conditions.
"""
struct Identity <: AbstractPreconditioner end

function NQSCore.precondition(::Identity, vs::AbstractVariationalState, operator)
    return expect_and_grad(vs, operator)
end

"""
    optimize!(state, operator, preconditioner; iterations, learning_rate, callback=nothing)
        -> Vector{Stats}

Run a variational optimization, returning the energy at every step.

This is the smallest useful driver: no logging, no checkpointing, no early stopping. Those
belong with the umbrella package; what is here is enough to test that a preconditioner actually
descends.

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
