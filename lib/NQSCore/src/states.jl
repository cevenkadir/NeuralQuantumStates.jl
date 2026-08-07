"""
    born_probabilities(logψ) -> Vector{Float64}

Normalized Born probabilities `|ψ(s)|² / Σ|ψ|²` from a vector of log-amplitudes.

Computed by shifting the log-modulus by its maximum before exponentiating, so that a
wavefunction with a large log-modulus does not overflow on the way to a perfectly ordinary
probability.
"""
function born_probabilities(logψ::AbstractVector)
    logp = 2 .* real.(logψ)
    logp .-= maximum(logp)
    p = exp.(logp)
    return p ./ sum(p)
end

"""
    configurations_of(state, packed_states) -> Matrix

Unpack packed states into the `(nsites, batch)` numeric array an ansatz consumes.

This is the boundary between the two representations: packed integers are what the operator
kernel and the samplers work in, numeric arrays are what an ansatz consumes. Entry points that
need both the configurations and their log-amplitudes build this once and pass it along, rather
than rebuilding it for each consumer.
"""
function configurations_of(vs::AbstractVariationalState, states::AbstractArray)
    a = ansatz(vs)
    return ConnectedBasisConfigurations.configurations(dof(a), vec(states), n_sites(a))
end

"""
    log_amplitudes(state, packed_states) -> Vector

Log-amplitudes of `state`'s ansatz on a vector of **packed** configurations.
"""
log_amplitudes(vs::AbstractVariationalState, states::AbstractArray) =
    log_amplitude(ansatz(vs), parameters(vs), configurations_of(vs, states))

"""
    local_energy(state, operator, packed_states) -> AbstractArray

Local energies `E_loc(s) = Σ_{s'} ⟨s|Ô|s'⟩ ψ(s')/ψ(s)`.

The ratio is evaluated as `exp(log ψ(s') - log ψ(s))`, never as a quotient of amplitudes. That
is not a micro-optimization: amplitudes underflow to zero for any system worth studying, while
the *difference* of their logarithms stays perfectly well behaved.

Multi-chain samples arrive as a `(steps, chains)` array, and the result keeps that shape, so
the chain structure survives into [`statistics`](@ref) — which needs it for split-R̂ and for a
between-chain error bar. Flattening happens here rather than in the sampler, leaving the chain
layout the sampler's business.

`operator` is anything `ConnectedBasisConfigurations.connected_padded` accepts. Passing a
`compile`d operator skips recompiling it on every call, which is worth doing inside an
optimization loop and irrelevant for a single large batch, where the connected-configuration
kernel dominates.
"""
function local_energy(vs::AbstractVariationalState, operator, states::AbstractArray)
    E = local_energy(vs, operator, vec(states))
    return reshape(E, size(states))
end

local_energy(vs::AbstractVariationalState, operator, states::AbstractVector) =
    _local_energy(vs, operator, states, log_amplitudes(vs, states))

local_energy(vs::AbstractVariationalState, operator) =
    local_energy(vs, operator, samples(vs))

"""
The local-energy kernel, given the sample log-amplitudes the caller has already computed.

Every public entry point needs `log ψ` on the samples for something else as well — the Born
weights, the gradient — so it is computed once at the top and threaded down here rather than
recomputed behind each caller's back.
"""
function _local_energy(
    vs::AbstractVariationalState, operator, states::AbstractVector, logψ_s::AbstractVector
)
    res = ConnectedBasisConfigurations.connected_padded(operator, states)

    # One batched evaluation over every connected configuration of every sample, rather than
    # one call per sample: the ansatz is the expensive part, so it is called once.
    logψ_sp = reshape(log_amplitudes(vs, vec(res.configs)), size(res.configs))

    # The accumulator has to admit both factors: a real-valued ansatz is perfectly legitimate,
    # and its local energies are still complex whenever the operator's matrix elements are.
    T = promote_type(eltype(res.mels), eltype(logψ_sp), eltype(logψ_s))
    E = similar(logψ_s, T)
    @inbounds for b in eachindex(states)
        acc = zero(T)
        logψ_b = logψ_s[b]
        for j in 1:res.counts[b]
            acc += res.mels[j, b] * exp(logψ_sp[j, b] - logψ_b)
        end
        E[b] = acc
    end
    return E
end

"""
    local_estimators(state, operator; holomorphic=false) -> (; E, O, weights)

Everything an estimator built from this state needs, computed once from one set of samples.

- `E`: local energies, flattened.
- `O`: the log-derivative matrix, `(n_samples, n_parameters)`.
- `weights`: exact Born probabilities for a [`FullSumState`](@ref), `nothing` for an
  [`MCState`](@ref) whose samples are already distributed according to `|ψ|²`.

The plain energy gradient and the quantum geometric tensor are built from exactly these three
things, so exposing them keeps `NQSOptimisers` from having to recompute — and, more
importantly, guarantees that a preconditioned update and the energy it is derived from come
from the *same* samples. Recomputing would silently mix two sample sets, which shows up as an
optimizer that mysteriously fails to descend.

Note that `expect_and_grad` does **not** go through here: a plain gradient needs only one
contraction of `O`, which [`energy_gradient`](@ref) obtains without building the matrix at all.
`O` is materialized only for the preconditioners that genuinely need it.
"""
function local_estimators(vs::AbstractVariationalState, operator; holomorphic::Bool=false)
    states = samples(vs)
    a, θ = ansatz(vs), parameters(vs)

    x = configurations_of(vs, states)
    logψ = log_amplitude(a, θ, x)
    E = _local_energy(vs, operator, vec(states), logψ)
    O = log_derivatives(a, θ, x; backend=vs.backend, holomorphic=holomorphic)

    return (; E=E, O=O, weights=_sample_weights(vs, logψ))
end

"""
    sample_weights(state) -> Union{Nothing,Vector}

Probability weights attached to a state's samples, or `nothing` when they are already drawn
from `|ψ|²` and so carry equal weight.
"""
function sample_weights end

"""Weights from log-amplitudes the caller already has, avoiding a second pass over the basis."""
function _sample_weights end

# ---------------------------------------------------------------------------- FullSumState

"""
    FullSumState(ansatz, parameters; backend, basis=nothing) <: AbstractVariationalState

A variational state that sums **exactly** over the whole basis instead of sampling.

The counterpart of NetKet's `FullSumState`, and the reason it matters is testing: it satisfies
the same interface as [`MCState`](@ref) but has no sampling noise at all, so a disagreement
between the two is a bug rather than a fluctuation. Any model can be developed and debugged
against exact summation and then scaled up by swapping the state type.

Cost is the dimension of the Hilbert space, so this is for small systems only.

The basis is what the state sums over, so it belongs to the state rather than to the ansatz —
an ansatz is a functional form and has no opinion about which configurations exist. It defaults
to [`default_basis`](@ref) of the ansatz; pass `basis` explicitly to sum over a symmetry sector
instead.

# Fields
- `ansatz`, `parameters`: the wavefunction.
- `basis`: the configurations summed over.
- `backend`: the DifferentiationInterface backend used for [`log_derivatives`](@ref).
"""
mutable struct FullSumState{A<:AbstractAnsatz,P,S,B} <: AbstractVariationalState
    ansatz::A
    parameters::P
    basis::S
    backend::B
end

function FullSumState(ansatz::AbstractAnsatz, parameters; backend, basis=nothing)
    b = basis === nothing ? default_basis(ansatz) : basis
    return FullSumState(ansatz, parameters, b, backend)
end

"""
    default_basis(ansatz) -> basis

The basis a [`FullSumState`](@ref) should sum over when none is given.

There is deliberately no generic fallback. Constructing the full space from a degree-of-freedom
specification would mean depending on the basis library for one convenience method, and it
would be the wrong answer for any ansatz defined on a symmetry sector. Ansatze that do know
their basis define this; everything else passes `basis` explicitly.
"""
default_basis(a::AbstractAnsatz) = throw(ArgumentError(
    "no default basis is known for $(typeof(a)); pass `basis=...` when constructing the " *
    "state, or define `NQSCore.default_basis` for this ansatz"
))

"""A `LogStateVector` is indexed by a specific basis, so that is the one to sum over."""
default_basis(a::LogStateVector) = a.basis

ansatz(vs::FullSumState) = vs.ansatz
parameters(vs::FullSumState) = vs.parameters
setparameters!(vs::FullSumState, θ) = (vs.parameters = θ; vs)
samples(vs::FullSumState) = vs.basis.states
sample_weights(vs::FullSumState) = probabilities(vs)
_sample_weights(::FullSumState, logψ::AbstractVector) = born_probabilities(logψ)

"""
    probabilities(state) -> Vector{Float64}

The exact Born probabilities `|ψ(s)|² / Σ|ψ|²` over the whole basis.
"""
probabilities(vs::FullSumState) = born_probabilities(log_amplitudes(vs, samples(vs)))

function expect(vs::FullSumState, operator)
    states = samples(vs)
    logψ = log_amplitudes(vs, states)
    E = _local_energy(vs, operator, states, logψ)
    return weighted_statistics(E, born_probabilities(logψ))
end

function expect_and_grad(vs::FullSumState, operator)
    states = samples(vs)
    a, θ = ansatz(vs), parameters(vs)

    x = configurations_of(vs, states)
    logψ = log_amplitude(a, θ, x)
    E = _local_energy(vs, operator, states, logψ)
    p = born_probabilities(logψ)

    ∇ = energy_gradient(a, θ, x, E, p; backend=vs.backend)
    return weighted_statistics(E, p), ∇
end

# -------------------------------------------------------------------------------- MCState

"""
    MCState(ansatz, parameters, sampler; backend, rng) <: AbstractVariationalState

A variational state that estimates expectation values by Monte Carlo.

Samples are drawn when the state is built and reused until something invalidates them;
[`resample!`](@ref) draws a fresh set. Reuse is what makes `expect` and `expect_and_grad`
consistent with one another — both must be computed from the *same* samples, or the gradient
does not correspond to the reported energy.

Drawing eagerly, rather than on first use, is what lets the sample and sampler-state fields
have concrete types: their types are whatever the sampler returns, which is not knowable from
the sampler's type alone.

The sampler's own state survives a parameter change even though the samples do not, so a Markov
chain resumes from where it was instead of restarting. See [`sample`](@ref).

# Fields
- `ansatz`, `parameters`: the wavefunction.
- `sampler`: an [`AbstractSampler`](@ref).
- `backend`: the DifferentiationInterface backend for [`log_derivatives`](@ref).
- `rng`: the random source.
"""
mutable struct MCState{A<:AbstractAnsatz,P,S<:AbstractSampler,B,R<:AbstractRNG,C,T} <:
               AbstractVariationalState
    ansatz::A
    parameters::P
    sampler::S
    backend::B
    rng::R
    samples::C
    sampler_state::T
    stale::Bool
end

function MCState(
    ansatz::AbstractAnsatz, parameters, sampler::AbstractSampler;
    backend, rng::AbstractRNG=Random.default_rng()
)
    drawn, state = sample(sampler, ansatz, parameters, rng, nothing)
    return MCState(ansatz, parameters, sampler, backend, rng, drawn, state, false)
end

ansatz(vs::MCState) = vs.ansatz
parameters(vs::MCState) = vs.parameters

sample_weights(::MCState) = nothing
_sample_weights(::MCState, ::AbstractVector) = nothing

"""
    sampler_state(state)

The sampler's own state, carried between draws. `nothing` for a sampler that draws
independently.
"""
sampler_state(vs::MCState) = vs.sampler_state

"""
Changing the parameters invalidates the cached samples: they came from the old `|ψ|²`.

The *sampler* state is deliberately kept. A parameter update moves the distribution only
slightly, so the chain it describes is still very nearly equilibrated — throwing it away would
mean re-paying the burn-in on every optimization step.
"""
setparameters!(vs::MCState, θ) = (vs.parameters = θ; vs.stale = true; vs)

"""
    resample!(state) -> state

Draw a fresh set of Monte Carlo samples, resuming from the current sampler state.
"""
function resample!(vs::MCState)
    drawn, state = sample(vs.sampler, vs.ansatz, vs.parameters, vs.rng, vs.sampler_state)
    vs.samples = drawn
    vs.sampler_state = state
    vs.stale = false
    return vs
end

function samples(vs::MCState)
    vs.stale && resample!(vs)
    return vs.samples
end

function expect(vs::MCState, operator)
    return statistics(local_energy(vs, operator, samples(vs)))
end

function expect_and_grad(vs::MCState, operator)
    states = samples(vs)
    a, θ = ansatz(vs), parameters(vs)

    x = configurations_of(vs, states)
    logψ = log_amplitude(a, θ, x)
    E = _local_energy(vs, operator, vec(states), logψ)

    ∇ = energy_gradient(a, θ, x, E, nothing; backend=vs.backend)
    return statistics(reshape(E, size(states))), ∇
end
