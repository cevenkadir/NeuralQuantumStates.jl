"""
    flatten_parameters(θ) -> (flat, restore)

Flatten a parameter container to a plain vector, and return the function that rebuilds it.

Lux hands back parameters as nested `NamedTuple`s, while every piece of linear algebra
downstream — the quantum geometric tensor, the gradient, the optimizer step — wants a flat
vector. `ComponentArrays` provides the bijection; `restore` is what carries an updated flat
vector back into the nested shape the ansatz expects.
"""
function flatten_parameters(θ::AbstractVector)
    return θ, identity
end

function flatten_parameters(θ::NamedTuple)
    ca = ComponentArray(θ)
    ax = getaxes(ca)
    # Rebuild as a NamedTuple, not a ComponentArray. The promise downstream is that a gradient
    # has the same structure as the parameters it belongs to, so that `fmap` can walk the two
    # together; handing back a ComponentArray where a NamedTuple went in quietly breaks that.
    return getdata(ca), p -> NamedTuple(ComponentArray(p, ax))
end

function flatten_parameters(θ)
    ca = ComponentArray(θ)
    ax = getaxes(ca)
    return getdata(ca), p -> ComponentArray(p, ax)
end

"""
    log_derivatives(ansatz, parameters, x; backend, holomorphic=false) -> Matrix

The log-derivative matrix `O[s, k] = ∂ log ψ(x_s) / ∂θ_k`, of size `(batch, n_parameters)`.

`O` is the central object of variational Monte Carlo beyond plain gradient descent: the energy
gradient is a covariance between `O` and the local energies, and the quantum geometric tensor
that stochastic reconfiguration inverts is `S = ⟨O†O⟩ - ⟨O⟩†⟨O⟩`.

Rows are samples and columns are parameters, matching NetKet, so that `S = O'O` has the shape
one expects without transposing.

# Complex arithmetic

Two independent things can be complex here, and conflating them is the usual source of wrong
gradients.

**A complex log-amplitude with real parameters** — the common case for a neural network that
outputs a modulus and a phase. `O` is then complex even though `θ` is real, and both parts come
from a single automatic-differentiation pass over `[real(log ψ); imag(log ψ)]`. Doing it in one
pass rather than two halves the cost, since the two share all of their intermediate work.

**Complex parameters** are ambiguous until a convention is fixed:

- `holomorphic=true` assumes `ψ` is holomorphic in `θ`, so `∂/∂θ = ∂/∂θ_re`, and `O` keeps one
  column per parameter.
- `holomorphic=false` (the default) treats real and imaginary parts as `2n` independent real
  parameters, returning `[∂/∂θ_re  ∂/∂θ_im]` with `2n` columns. This is always valid; the
  holomorphic form is a shortcut that is silently wrong when the ansatz does not satisfy it, so
  it must be asked for explicitly.

# Backend

Differentiation goes through DifferentiationInterface.jl, so `backend` is any of its objects —
`AutoForwardDiff()`, `AutoZygote()`, `AutoEnzyme()`. None of those packages is a dependency
here; load the one you want and pass it.
"""
function log_derivatives(
    ansatz::AbstractAnsatz, θ, x::AbstractMatrix;
    backend, holomorphic::Bool=false
)
    flat, restore = flatten_parameters(θ)
    return _log_derivatives(ansatz, flat, restore, x, backend, holomorphic)
end

"""Real parameters: one pass over the stacked real and imaginary parts of `log ψ`."""
function _log_derivatives(
    ansatz, flat::AbstractVector{<:Real}, restore, x, backend, ::Bool
)
    batch = size(x, 2)
    function stacked(p)
        ψ = log_amplitude(ansatz, restore(p), x)
        return vcat(real.(ψ), imag.(ψ))
    end
    J = DifferentiationInterface.jacobian(stacked, backend, flat)
    return @views J[1:batch, :] .+ im .* J[(batch+1):(2batch), :]
end

"""Complex parameters: differentiate with respect to the real and imaginary parts separately."""
function _log_derivatives(
    ansatz, flat::AbstractVector{<:Complex}, restore, x, backend, holomorphic::Bool
)
    batch = size(x, 2)
    n = length(flat)
    split = vcat(real.(flat), imag.(flat))

    function stacked(v)
        p = @views v[1:n] .+ im .* v[(n+1):(2n)]
        ψ = log_amplitude(ansatz, restore(p), x)
        return vcat(real.(ψ), imag.(ψ))
    end
    J = DifferentiationInterface.jacobian(stacked, backend, split)

    # Blocks of ∂(re ψ, im ψ) / ∂(re θ, im θ).
    ∂reψ_∂reθ = @view J[1:batch, 1:n]
    ∂imψ_∂reθ = @view J[(batch+1):(2batch), 1:n]
    O_re = ∂reψ_∂reθ .+ im .* ∂imψ_∂reθ

    if holomorphic
        # For a holomorphic ψ, differentiating along the real axis is the full derivative.
        return O_re
    end

    ∂reψ_∂imθ = @view J[1:batch, (n+1):(2n)]
    ∂imψ_∂imθ = @view J[(batch+1):(2batch), (n+1):(2n)]
    O_im = ∂reψ_∂imθ .+ im .* ∂imψ_∂imθ
    return hcat(O_re, O_im)
end

"""
    match_parameter_shape(∇, θ) -> gradient

Put a raw gradient into the same shape as the parameters, so that `θ .- η .* ∇` is meaningful.

For complex parameters differentiated non-holomorphically, [`log_derivatives`](@ref) returns
`2n` columns — `∂/∂θ_re` followed by `∂/∂θ_im` — and a gradient built from it inherits that
length. Descending in the real parameterization means `θ_re -= η g_re` and `θ_im -= η g_im`
simultaneously, which is exactly `θ -= η (g_re + i g_im)`. Recombining here rather than at the
call site keeps the promise that a gradient always matches the parameters it belongs to,
whatever the parameterization: an optimizer should never have to ask how the ansatz was
parameterized.
"""
function match_parameter_shape(∇::AbstractVector, θ)
    flat, restore = flatten_parameters(θ)
    return _match_parameter_shape(∇, flat, restore)
end

function _match_parameter_shape(∇::AbstractVector, flat::AbstractVector, restore)
    n = length(flat)
    if eltype(flat) <: Complex && length(∇) == 2n
        return restore(@views ∇[1:n] .+ im .* ∇[(n+1):(2n)])
    end
    return restore(∇)
end

"""
    energy_gradient(ansatz, θ, x, E_loc, weights; backend) -> gradient

The variational energy gradient `2 Re[⟨O_k^* ΔE⟩]`, computed **without ever forming `O`**.

The gradient is a single contraction of the log-derivative matrix against the centered local
energies, and a contraction of a Jacobian is a job for one differentiation pass rather than for
`n_samples` of them. Writing

```math
L(θ) = 2 \\sum_s \\mathrm{Re}[\\overline{c_s} \\, \\log ψ(x_s)],
\\qquad c_s = p_s (E_s - \\bar{E})
```

makes that explicit: `c` does not depend on `θ`, so `∇L` *is* the energy gradient, and `L` is a
scalar. Differentiating a scalar is the one thing every backend does optimally — a reverse-mode
backend needs exactly one pass regardless of how many samples there are, where building the
full Jacobian would have cost it one pass per sample.

[`log_derivatives`](@ref) still builds `O` explicitly, because stochastic reconfiguration needs
the matrix itself and not just this one contraction of it.
"""
function energy_gradient(
    ansatz::AbstractAnsatz, θ, x::AbstractMatrix, E::AbstractVector,
    weights::Union{Nothing,AbstractVector}; backend
)
    flat, restore = flatten_parameters(θ)
    c = _gradient_cotangent(E, weights)
    ∇ = _energy_gradient(ansatz, flat, restore, x, c, backend)
    return _match_parameter_shape(∇, flat, restore)
end

"""
The cotangent `c_s = p_s (E_s - Ē)`.

Subtracting the mean is not cosmetic: without it the estimator has a non-vanishing variance
even at an exact eigenstate, where every local energy is the same number.
"""
function _gradient_cotangent(E::AbstractVector, weights::Union{Nothing,AbstractVector})
    p = weights === nothing ? fill(1 / length(E), length(E)) : weights ./ sum(weights)
    Ē = sum(p .* E)
    return p .* (E .- Ē)
end

"""Real parameters: differentiate the scalar loss directly."""
function _energy_gradient(ansatz, flat::AbstractVector{<:Real}, restore, x, c, backend)
    loss(p) = _gradient_loss(ansatz, restore(p), x, c)
    return DifferentiationInterface.gradient(loss, backend, flat)
end

"""
Complex parameters: differentiate with respect to real and imaginary parts as `2n` independent
real parameters, matching the non-holomorphic convention of [`log_derivatives`](@ref).
"""
function _energy_gradient(ansatz, flat::AbstractVector{<:Complex}, restore, x, c, backend)
    n = length(flat)
    split = vcat(real.(flat), imag.(flat))
    function loss(v)
        p = @views v[1:n] .+ im .* v[(n+1):(2n)]
        return _gradient_loss(ansatz, restore(p), x, c)
    end
    return DifferentiationInterface.gradient(loss, backend, split)
end

"""`L = 2 Σ_s Re[conj(c_s) log ψ(x_s)]`, the scalar whose gradient is the energy gradient."""
function _gradient_loss(ansatz, θ, x, c)
    ψ = log_amplitude(ansatz, θ, x)
    return 2 * sum(real.(conj.(c) .* ψ))
end

"""
    centered(O, weights=nothing) -> Matrix

Subtract the (optionally weighted) mean of each column of `O`.

Every use of `O` in variational Monte Carlo wants the centered version: the energy gradient is
a covariance, and the quantum geometric tensor is a covariance matrix. Centering also removes
the direction corresponding to a global change of normalization, which carries no physical
information and would otherwise make the geometric tensor singular.
"""
function centered(O::AbstractMatrix, weights::Union{Nothing,AbstractVector}=nothing)
    if weights === nothing
        return O .- sum(O; dims=1) ./ size(O, 1)
    end
    p = weights ./ sum(weights)
    return O .- sum(p .* O; dims=1)
end
