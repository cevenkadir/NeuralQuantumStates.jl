"""
    log_amplitudes(state, packed_states) -> Vector

Log-amplitudes of `state`'s ansatz on a vector of **packed** configurations, unpacking them into
the numeric array the ansatz interface expects.

This is the boundary between the two representations: packed integers are what the operator
kernel and the samplers work in, numeric arrays are what an ansatz consumes.
"""
function log_amplitudes(vs::AbstractVariationalState, states::AbstractArray)
    a = ansatz(vs)
    x = ConnectedConfigs.configurations(a.dof, vec(states), a.nsites)
    return log_amplitude(a, parameters(vs), x)
end

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
"""
function local_energy(vs::AbstractVariationalState, operator, states::AbstractArray)
    E = local_energy(vs, operator, vec(states))
    return reshape(E, size(states))
end

function local_energy(vs::AbstractVariationalState, operator, states::AbstractVector)
    res = ConnectedConfigs.connected_padded(operator, states)
    a = ansatz(vs)

    logψ_s = log_amplitudes(vs, states)

    # One batched evaluation over every connected configuration of every sample, rather than
    # one call per sample: the ansatz is the expensive part, so it is called once.
    flat = vec(res.configs)
    logψ_sp = reshape(log_amplitudes(vs, flat), size(res.configs))

    E = similar(logψ_s)
    for b in eachindex(states)
        acc = zero(eltype(E))
        for j in 1:res.counts[b]
            acc += res.mels[j, b] * exp(logψ_sp[j, b] - logψ_s[b])
        end
        E[b] = acc
    end
    return E
end

local_energy(vs::AbstractVariationalState, operator) =
    local_energy(vs, operator, samples(vs))

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
"""
function local_estimators(vs::AbstractVariationalState, operator; holomorphic::Bool=false)
    states = samples(vs)
    E = local_energy(vs, operator, states)

    a = ansatz(vs)
    x = ConnectedConfigs.configurations(a.dof, vec(states), a.nsites)
    O = log_derivatives(a, parameters(vs), x; backend=vs.backend, holomorphic=holomorphic)

    return (; E=vec(E), O=O, weights=sample_weights(vs))
end

"""
    sample_weights(state) -> Union{Nothing,Vector}

Probability weights attached to a state's samples, or `nothing` when they are already drawn
from `|ψ|²` and so carry equal weight.
"""
function sample_weights end

"""
    _gradient(O, E_loc, weights) -> Vector

The variational energy gradient `2 Re[⟨O* E_loc⟩ - ⟨O*⟩⟨E_loc⟩]`.

`weights` is `nothing` for an unweighted Monte Carlo average, or the exact probabilities for a
full summation.
"""
function _gradient(O::AbstractMatrix, E::AbstractVector, weights::Union{Nothing,AbstractVector})
    p = weights === nothing ? fill(1 / length(E), length(E)) : weights ./ sum(weights)
    Ē = sum(p .* E)
    Ō = vec(sum(p .* O; dims=1))
    ŌE = vec(sum(p .* conj.(O) .* E; dims=1))
    return 2 .* real.(ŌE .- conj.(Ō) .* Ē)
end

"""
    _match_parameter_shape(∇, θ) -> gradient

Put a raw gradient into the same shape as the parameters, so that `θ .- η .* ∇` is meaningful.

For complex parameters differentiated non-holomorphically, `log_derivatives` returns `2n`
columns — `∂/∂θ_re` followed by `∂/∂θ_im` — and the gradient inherits that length. Descending
in the real parameterization means `θ_re -= η g_re` and `θ_im -= η g_im` simultaneously, which
is exactly `θ -= η (g_re + i g_im)`. Recombining here rather than at the call site keeps the
promise that the gradient always matches the parameters, whatever the parameterization: an
optimizer should never have to ask how the ansatz was parameterized.
"""
function _match_parameter_shape(∇::AbstractVector, θ)
    flat, restore = flatten_parameters(θ)
    n = length(flat)
    if eltype(flat) <: Complex && length(∇) == 2n
        return restore(∇[1:n] .+ im .* ∇[(n+1):(2n)])
    end
    return restore(∇)
end

# ---------------------------------------------------------------------------- FullSumState

"""
    FullSumState(ansatz, parameters; backend) <: AbstractVariationalState

A variational state that sums **exactly** over the whole basis instead of sampling.

The counterpart of NetKet's `FullSumState`, and the reason it matters is testing: it satisfies
the same interface as [`MCState`](@ref) but has no sampling noise at all, so a disagreement
between the two is a bug rather than a fluctuation. Any model can be developed and debugged
against exact summation and then scaled up by swapping the state type.

Cost is the dimension of the Hilbert space, so this is for small systems only.

# Fields
- `ansatz`, `parameters`: the wavefunction.
- `backend`: the DifferentiationInterface backend used for [`log_derivatives`](@ref).
"""
mutable struct FullSumState{A<:AbstractAnsatz,P,S,B} <: AbstractVariationalState
    ansatz::A
    parameters::P
    basis::S
    backend::B
end

"""
    FullSumState(ansatz, parameters, backend; basis=nothing)

The basis is what the state sums over, so it belongs to the state rather than to the ansatz —
an ansatz is a functional form and has no opinion about which configurations exist. It defaults
to the full space implied by the ansatz's `dof` and `nsites`; pass `basis` explicitly to sum
over a symmetry sector instead.
"""
function FullSumState(ansatz::AbstractAnsatz, parameters, backend; basis=nothing)
    b = basis === nothing ? default_basis(ansatz) : basis
    return FullSumState(ansatz, parameters, b, backend)
end

"""
    default_basis(ansatz) -> SymBasis.Basis

The full computational basis implied by an ansatz's degrees of freedom.
"""
default_basis(a::AbstractAnsatz) = SymBasis.Bases.basis(SymBasis.dof_object(a.dof), a.nsites)

"""A `LogStateVector` is indexed by a specific basis, so that is the one to sum over."""
default_basis(a::LogStateVector) = a.basis

ansatz(vs::FullSumState) = vs.ansatz
parameters(vs::FullSumState) = vs.parameters
setparameters!(vs::FullSumState, θ) = (vs.parameters = θ; vs)
samples(vs::FullSumState) = vs.basis.states
sample_weights(vs::FullSumState) = probabilities(vs)

"""
    probabilities(state) -> Vector{Float64}

The exact Born probabilities `|ψ(s)|² / Σ|ψ|²` over the whole basis.

Computed by shifting the log-amplitudes by their maximum before exponentiating, so that a
wavefunction with a large log-modulus does not overflow on the way to a perfectly ordinary
probability.
"""
function probabilities(vs::FullSumState)
    logψ = log_amplitudes(vs, samples(vs))
    logp = 2 .* real.(logψ)
    logp .-= maximum(logp)
    p = exp.(logp)
    return p ./ sum(p)
end

function expect(vs::FullSumState, operator)
    E = local_energy(vs, operator, samples(vs))
    return weighted_statistics(E, probabilities(vs))
end

function expect_and_grad(vs::FullSumState, operator)
    states = samples(vs)
    E = local_energy(vs, operator, states)
    p = probabilities(vs)

    a = ansatz(vs)
    x = ConnectedConfigs.configurations(a.dof, states, a.nsites)
    O = log_derivatives(a, parameters(vs), x; backend=vs.backend)

    return weighted_statistics(E, p), _match_parameter_shape(_gradient(O, E, p), parameters(vs))
end

# -------------------------------------------------------------------------------- MCState

"""
    MCState(ansatz, parameters, sampler; backend, rng) <: AbstractVariationalState

A variational state that estimates expectation values by Monte Carlo.

Samples are drawn once and cached; [`resample!`](@ref) draws a fresh set. Caching is what makes
`expect` and `expect_and_grad` consistent with one another — both must be computed from the
*same* samples, or the gradient does not correspond to the reported energy.

# Fields
- `ansatz`, `parameters`: the wavefunction.
- `sampler`: an [`AbstractSampler`](@ref).
- `backend`: the DifferentiationInterface backend for [`log_derivatives`](@ref).
- `rng`: the random source.
"""
mutable struct MCState{A<:AbstractAnsatz,P,S<:AbstractSampler,B,R<:AbstractRNG} <: AbstractVariationalState
    ansatz::A
    parameters::P
    sampler::S
    backend::B
    rng::R
    cache::Union{Nothing,AbstractArray}
end

function MCState(
    ansatz::AbstractAnsatz, parameters, sampler::AbstractSampler;
    backend, rng::AbstractRNG=Random.default_rng()
)
    return MCState(ansatz, parameters, sampler, backend, rng, nothing)
end

ansatz(vs::MCState) = vs.ansatz
parameters(vs::MCState) = vs.parameters

sample_weights(::MCState) = nothing

"""Changing the parameters invalidates the cached samples: they came from the old `|ψ|²`."""
setparameters!(vs::MCState, θ) = (vs.parameters = θ; vs.cache = nothing; vs)

"""
    resample!(state) -> state

Draw a fresh set of Monte Carlo samples.
"""
function resample!(vs::MCState)
    vs.cache = sample(vs.sampler, vs.ansatz, vs.parameters, vs.rng)
    return vs
end

function samples(vs::MCState)
    vs.cache === nothing && resample!(vs)
    return vs.cache
end

function expect(vs::MCState, operator)
    E = local_energy(vs, operator, samples(vs))
    return statistics(E)
end

function expect_and_grad(vs::MCState, operator)
    states = samples(vs)
    E = local_energy(vs, operator, states)

    a = ansatz(vs)
    x = ConnectedConfigs.configurations(a.dof, vec(states), a.nsites)
    O = log_derivatives(a, parameters(vs), x; backend=vs.backend)

    return statistics(E), _match_parameter_shape(_gradient(O, vec(E), nothing), parameters(vs))
end
