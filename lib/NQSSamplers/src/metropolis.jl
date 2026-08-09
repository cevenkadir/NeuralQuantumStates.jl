"""
    MetropolisSampler(rule, initial; n_chains, n_samples, burn_in, thinning, basis=nothing)

Metropolis–Hastings sampling of `|ψ(s)|²` over discrete configurations.

# Fields
- `rule`: an [`AbstractRule`](@ref) supplying proposals.
- `initial`: a configuration to start every chain from, or a vector of one per chain.
- `n_chains::Int`: chains advanced in parallel.
- `n_samples::Int`: samples **per chain** kept after burn-in and thinning.
- `burn_in::Int`: steps discarded at the start of each chain.
- `thinning::Int`: steps taken between kept samples.
- `basis`: if given, moves leaving this basis are rejected — how a sampler is confined to a
  symmetry sector.

# Why multiple chains

Not for speed. Several chains started from different configurations give a convergence
diagnostic no single chain can: split-R̂ compares them, and a chain stuck in one region is
invisible from the inside. The error bar also comes from the spread between chain means, which
assumes only that the chains are independent rather than relying on an autocorrelation model.

All chains advance together, so a step evaluates the ansatz once on a batch of `n_chains`
configurations rather than once per chain.

!!! tip "On a GPU, use hundreds of chains"
    A sweep costs one ansatz evaluation whatever its width, and on a device that is a handful of
    kernel launches whose latency does not depend on how much data they carry. Eight chains
    hands a GPU about a hundred numbers per launch and spends everything on overhead, while the
    per-sample cost on a host is flat in the width. The default of eight is right for a CPU and
    is left alone because more chains also means more memory.
"""
struct MetropolisSampler{R<:AbstractRule,S,B} <: AbstractSampler
    rule::R
    initial::Vector{S}
    n_chains::Int
    n_samples::Int
    burn_in::Int
    thinning::Int
    basis::B
end

function MetropolisSampler(
    rule::AbstractRule, initial;
    n_chains::Integer=8, n_samples::Integer=1000,
    burn_in::Integer=100, thinning::Integer=1, basis=nothing
)
    starts = initial isa AbstractVector ? collect(initial) : fill(initial, n_chains)
    length(starts) == n_chains || throw(ArgumentError(
        "got $(length(starts)) initial configurations for $n_chains chains"
    ))
    n_samples > 0 || throw(ArgumentError("n_samples must be positive"))
    thinning > 0 || throw(ArgumentError("thinning must be positive"))
    return MetropolisSampler(rule, starts, Int(n_chains), Int(n_samples),
        Int(burn_in), Int(thinning), basis)
end

"""Log of `|ψ|²` for a batch of packed configurations."""
function _log_prob(a::AbstractAnsatz, θ, states::AbstractVector)
    x = ConnectedBasisConfigurations.configurations(NQSCore.dof(a), states, NQSCore.n_sites(a))
    # To the host, because the acceptance test is scalar and sequential: against a device array
    # it would be a round trip per element. One transfer of `n_chains` numbers per step is the
    # cheap version of the same thing.
    return NQSCore.to_host(2 .* real.(log_amplitude(a, θ, x)))
end

"""Whether a proposed configuration is admissible: inside the basis, when one is given."""
_admissible(::Nothing, s) = true
_admissible(basis, s) = s in basis.states

# The sampler state is the configuration each chain finished on. Handing it back resumes those
# chains, which is also why the burn-in is skipped: they are already where it would take them.
function NQSCore.sample(
    sampler::MetropolisSampler, a::AbstractAnsatz, θ, rng::AbstractRNG, state=nothing
)
    warm = state !== nothing
    chains = warm ? copy(state) : copy(sampler.initial)
    length(chains) == sampler.n_chains || throw(ArgumentError(
        "got $(length(chains)) chain configurations for $(sampler.n_chains) chains"
    ))
    for s in chains
        _admissible(sampler.basis, s) || throw(ArgumentError(
            "an initial configuration is not in the sampler's basis"
        ))
    end

    logp = _log_prob(a, θ, chains)

    burn_in = warm ? 0 : sampler.burn_in
    total_steps = burn_in + sampler.n_samples * sampler.thinning
    out = Matrix{eltype(chains)}(undef, sampler.n_samples, sampler.n_chains)
    kept = 0
    accepted = 0
    proposed = 0

    # Hoisted out of the step loop; a Markov step should not allocate.
    proposals = similar(chains)
    corrections = zeros(Float64, sampler.n_chains)
    movable = falses(sampler.n_chains)

    spec, nsites = NQSCore.dof(a), NQSCore.n_sites(a)
    for step in 1:total_steps
        fill!(corrections, 0.0)
        fill!(movable, false)

        for c in 1:sampler.n_chains
            s′, correction = propose(sampler.rule, chains[c], spec, nsites, rng)
            if s′ == chains[c] || !_admissible(sampler.basis, s′)
                proposals[c] = chains[c]        # nothing to evaluate; the chain stays put
            else
                proposals[c] = s′
                corrections[c] = correction
                movable[c] = true
            end
        end

        # One batched ansatz evaluation per step, covering every chain that has a live proposal.
        if any(movable)
            logp′ = _log_prob(a, θ, proposals)
            for c in 1:sampler.n_chains
                movable[c] || continue
                proposed += 1
                if log(rand(rng)) < (logp′[c] - logp[c] + corrections[c])
                    chains[c] = proposals[c]
                    logp[c] = logp′[c]
                    accepted += 1
                end
            end
        end

        if step > burn_in && (step - burn_in) % sampler.thinning == 0
            kept += 1
            kept <= sampler.n_samples && (out[kept, :] .= chains)
        end
    end

    ACCEPTANCE[] = proposed == 0 ? 0.0 : accepted / proposed
    return out, chains
end

"""
    ACCEPTANCE

Acceptance rate of the most recent [`MetropolisSampler`](@ref) run.

A diagnostic, not part of the interface. Near zero means the chain is barely moving and the
samples are one configuration repeated; near one usually means the proposals are too timid to
explore. Both give error bars that look fine and mean nothing.
"""
const ACCEPTANCE = Ref(0.0)

"""
    random_configuration(dof, nsites, rng) -> BaseInt

A uniformly random configuration, for starting a chain.
"""
function random_configuration(dof, nsites::Integer, rng::AbstractRNG)
    values = local_values(dof)
    return ConnectedBasisConfigurations.packed(dof, [rand(rng, values) for _ in 1:nsites])
end

"""
    random_configurations(dof, nsites, n, rng) -> Vector

`n` independent random configurations — dispersed chain starts, which is what makes split-R̂
informative rather than vacuous.
"""
function random_configurations(dof, nsites::Integer, n::Integer, rng::AbstractRNG)
    return [random_configuration(dof, nsites, rng) for _ in 1:n]
end
