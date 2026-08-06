"""
    ExactSampler(basis, n_samples) <: AbstractSampler

Draws independent samples from the **exact** distribution `|ψ(s)|²`, by enumerating the basis
and computing every probability.

The counterpart of NetKet's `sampler.ExactSampler`, and it lives here alongside
[`FullSumState`](@ref) for the same reason: it is the reference implementation, not a
production tool. Enumerating the basis costs the dimension of the Hilbert space, so this is only
for small systems — but its samples are genuinely independent, with no burn-in, no
autocorrelation, and no proposal to tune.

That makes it the right instrument for a specific job: it isolates *Monte Carlo noise* from
*Markov chain pathology*. An `MCState` driven by this sampler must converge to the matching
`FullSumState`, and if it does not, the problem is in the estimator rather than in mixing.
`NQSSamplers` provides the Metropolis and autoregressive samplers meant for real use.

# Fields
- `basis`: the basis to enumerate.
- `n_samples::Int`: how many independent samples to draw.
"""
struct ExactSampler{B} <: AbstractSampler
    basis::B
    n_samples::Int
end

function sample(s::ExactSampler, a::AbstractAnsatz, θ, rng::AbstractRNG)
    states = s.basis.states
    x = ConnectedConfigs.configurations(a.dof, states, a.nsites)
    logψ = log_amplitude(a, θ, x)

    logp = 2 .* real.(logψ)
    logp .-= maximum(logp)              # shift before exponentiating, so nothing overflows
    p = exp.(logp)
    p ./= sum(p)

    # Inverse-CDF sampling over the enumerated distribution.
    cumulative = cumsum(p)
    out = similar(states, s.n_samples)
    for i in 1:s.n_samples
        r = rand(rng)
        out[i] = states[searchsortedfirst(cumulative, r)]
    end
    return out
end
