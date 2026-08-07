```@meta
CurrentModule = NQSCore
```

# Statistics

Every [`expect`](@ref) returns a [`Stats`](@ref), whichever state type produced it, so
downstream code never has to branch on whether a number came from exact summation or from a
Markov chain.

| Field | Meaning |
|---|---|
| `mean` | the estimate itself |
| `error_of_mean` | standard error, corrected for autocorrelation; zero for an exact result |
| `variance` | variance of the sampled quantity |
| `tau_corr` | integrated autocorrelation time, in sampling steps; `1` means independent samples |
| `r_hat` | split-R̂; `1` means converged, above about `1.01` means the chains disagree |

## Why the error bar is not `std/sqrt(n)`

Markov chain samples are correlated, so the naive formula understates the error by roughly
`sqrt(tau_corr)`. Reporting that unadjusted is the single most common way to converge
confidently to a wrong answer, so [`statistics`](@ref) applies the correction rather than
leaving it to the caller.

With more than one chain it does something better still: the error comes from the variance of
the per-chain means, which needs no autocorrelation model at all and assumes only that the
chains are mutually independent.

## The variance is a physical quantity

For an energy, `variance` is not merely a spread — it vanishes at an exact eigenstate. That
makes it the sharpest convergence diagnostic available, and a run whose energy has plateaued
while its variance has not is not converged, whatever the error bar says.

## Diagnosing a chain

[`integrated_autocorrelation`](@ref) uses the initial-positive-sequence estimator: the
autocorrelation function is summed until it first goes non-positive, which truncates a tail that
is pure noise. Summing the whole thing adds variance without adding information and can even
return a negative estimate.

[`split_rhat`](@ref) splits each chain in half before comparing. A single chain whose mean drifts
steadily through the run agrees with itself under an unsplit R̂ and disagrees under a split one,
which is exactly the failure that needs catching.

```@example stats
using NQSCore, Random

rng = Xoshiro(0)
φ, n = 0.8, 100_000
x = zeros(n)
for i in 2:n
    x[i] = φ * x[i-1] + randn(rng)
end

# An AR(1) chain with correlation φ has τ = (1+φ)/(1-φ) = 9.
integrated_autocorrelation(x)
```

```@example stats
# Four independent chains: R̂ near 1, and an error bar from their spread.
statistics(randn(Xoshiro(1), 2000, 4))
```

## Exact results

[`weighted_statistics`](@ref) is what a [`FullSumState`](@ref) uses: an exact expectation under
normalized weights, reporting the physical variance and a zero error bar because nothing was
sampled. [`exact_stats`](@ref) builds one directly.
