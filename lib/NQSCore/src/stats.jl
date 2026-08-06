"""
    Stats

The result of an [`expect`](@ref) call: an estimate together with everything needed to judge how
much to trust it.

Every state type returns this, so downstream code never has to branch on whether a number came
from exact summation or from a Markov chain.

# Fields
- `mean`: The estimate itself.
- `error_of_mean`: Standard error, corrected for autocorrelation. Zero for an exact result.
- `variance`: Variance of the sampled quantity. For an energy this is a physical quantity in its
  own right — it vanishes at an exact eigenstate, which makes it the sharpest available
  convergence diagnostic.
- `tau_corr`: Integrated autocorrelation time, in units of sampling steps. `1` means
  independent samples.
- `r_hat`: Split-R̂ convergence diagnostic. `1` means converged; anything above about `1.01`
  means the chains disagree and the error bar is not trustworthy.

# Why the error bar is not just `std/sqrt(n)`

Markov chain samples are correlated, so the naive formula understates the error by roughly
`sqrt(tau_corr)`. Reporting that number unadjusted is the single most common way to
"converge" to a wrong answer with confident-looking error bars, so the correction is applied
here rather than left to the caller.
"""
struct Stats{T}
    mean::T
    error_of_mean::Float64
    variance::Float64
    tau_corr::Float64
    r_hat::Float64
end

function Base.show(io::IO, s::Stats)
    print(io, round(real(s.mean), sigdigits=8))
    s.mean isa Complex && print(io, " + ", round(imag(s.mean), sigdigits=6), "im")
    s.error_of_mean > 0 && print(io, " ± ", round(s.error_of_mean, sigdigits=3))
    print(io, " (variance=", round(s.variance, sigdigits=4))
    isfinite(s.tau_corr) && print(io, ", τ=", round(s.tau_corr, sigdigits=3))
    isfinite(s.r_hat) && print(io, ", R̂=", round(s.r_hat, sigdigits=5))
    print(io, ")")
    return nothing
end

Base.isapprox(a::Stats, b::Stats; kwargs...) = isapprox(a.mean, b.mean; kwargs...)
Base.isapprox(a::Stats, b::Number; kwargs...) = isapprox(a.mean, b; kwargs...)
Base.isapprox(a::Number, b::Stats; kwargs...) = isapprox(a, b.mean; kwargs...)

"""
    exact_stats(value, variance=0.0) -> Stats

A `Stats` for a quantity computed exactly: zero error, no autocorrelation, converged by
construction.
"""
exact_stats(value, variance::Real=0.0) = Stats(value, 0.0, Float64(variance), 1.0, 1.0)

"""
    integrated_autocorrelation(chain) -> Float64

Integrated autocorrelation time of a single chain, by the initial-positive-sequence estimator.

The autocorrelation function is summed until it first goes non-positive, which is the standard
way of truncating a sum whose tail is pure noise. Summing the whole thing instead adds variance
without adding information, and can even produce a negative estimate.
"""
function integrated_autocorrelation(chain::AbstractVector{<:Real})
    n = length(chain)
    n < 4 && return 1.0

    x = chain .- mean(chain)
    denom = sum(abs2, x)
    iszero(denom) && return 1.0

    τ = 1.0
    for lag in 1:(n-1)
        ρ = sum(@views(x[1:(n-lag)]) .* @views(x[(1+lag):n])) / denom
        ρ <= 0 && break
        τ += 2ρ
    end
    return τ
end

"""
    split_rhat(chains) -> Float64

Split-R̂ over a matrix of chains, one chain per column.

Each chain is split in half before the comparison, so that a single chain that has drifted —
mean moving steadily through the run — is caught as disagreement between its own halves. R̂
computed without splitting is blind to exactly that failure.
"""
function split_rhat(chains::AbstractMatrix{<:Real})
    n_steps, n_chains = size(chains)
    n_steps < 4 && return 1.0

    half = n_steps ÷ 2
    pieces = [chains[1:half, c] for c in 1:n_chains]
    append!(pieces, [chains[(half+1):(2half), c] for c in 1:n_chains])

    m = length(pieces)
    means = mean.(pieces)
    variances = var.(pieces)

    W = mean(variances)                       # within-piece variance
    iszero(W) && return 1.0
    B = half * var(means)                     # between-piece variance
    V = (half - 1) / half * W + B / half      # marginal posterior variance estimate
    return sqrt(V / W)
end

"""
    statistics(values; n_chains=1) -> Stats

Reduce raw samples to a [`Stats`](@ref).

`values` may be a vector (one chain) or a `(steps, chains)` matrix. With more than one chain,
the error of the mean is estimated from the spread *between* chain means, which needs no
autocorrelation model at all; with a single chain it falls back to the autocorrelation-corrected
standard error.
"""
function statistics(values::AbstractMatrix)
    total = vec(values)
    μ = mean(total)
    σ² = var(total)

    n_steps, n_chains = size(values)
    real_values = real.(values)

    τ = mean(integrated_autocorrelation(real_values[:, c]) for c in 1:n_chains)
    r̂ = n_chains > 1 ? split_rhat(real_values) : split_rhat(reshape(real_values, :, 1))

    if n_chains > 1
        # Between-chain spread: assumes only that the chains are mutually independent.
        chain_means = [mean(values[:, c]) for c in 1:n_chains]
        err = sqrt(var(chain_means) / n_chains)
    else
        err = sqrt(max(σ², 0) * τ / length(total))
    end

    return Stats(μ, Float64(err), Float64(σ²), Float64(τ), Float64(r̂))
end

statistics(values::AbstractVector) = statistics(reshape(values, :, 1))

"""
    weighted_statistics(values, weights) -> Stats

Exact expectation of `values` under the normalized probability `weights`.

Used by [`FullSumState`](@ref), where the "samples" are the entire basis and the weights are
`|ψ(s)|²`. The reported variance is the physical variance of the quantity under that
distribution, and the error bar is zero because nothing was sampled.
"""
function weighted_statistics(values::AbstractVector, weights::AbstractVector{<:Real})
    total = sum(weights)
    p = weights ./ total
    μ = sum(p .* values)
    σ² = real(sum(p .* abs2.(values .- μ)))
    return exact_stats(μ, σ²)
end
