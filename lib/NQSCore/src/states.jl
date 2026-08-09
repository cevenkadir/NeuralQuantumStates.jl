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

Where it is built follows the *parameters*, not the samples: the array this returns is an input
to the ansatz, and an ansatz runs where its parameters are. On a host that is the loop over
digits below. On a device — with KernelAbstractions loaded — it is the same kernel the connected
configurations go through, which matters more than its size suggests, because this array is
consumed twice: once for `log ψ` and once inside the differentiated loss, where a host-to-device
conversion would sit in the middle of an automatic-differentiation pass.
"""
function configurations_of(vs::AbstractVariationalState, states::AbstractArray)
    flat = vec(states)
    reference = _reference_array(parameters(vs))
    reference isa AbstractArray || return _host_configurations(vs, flat)
    return _configurations(vs, flat, reference, _device_backend(reference))
end

_host_configurations(vs::AbstractVariationalState, states::AbstractVector) =
    ConnectedBasisConfigurations.configurations(
        dof(ansatz(vs)), states, n_sites(ansatz(vs))
    )

_configurations(vs::AbstractVariationalState, states::AbstractVector, ::AbstractArray, ::Nothing) =
    _host_configurations(vs, states)

"""
    NQSCore._reference_array(parameters) -> AbstractArray or nothing

Any array among the parameters, whose type says where the ansatz expects to be run.

Parameters are the only thing in a variational state that a caller deliberately places
somewhere. The samples are packed integers with no opinion, and the ansatz is a functional form
with none either, so asking the parameters is how anything here learns which memory it is
working in — without a device field to keep in sync, and without naming a GPU package.

`nothing` when there is no array to ask, which is a perfectly ordinary answer: it means the host.
"""
_reference_array(θ) = nothing
_reference_array(θ::NamedTuple) = isempty(θ) ? nothing : _reference_array(first(values(θ)))
_reference_array(θ::AbstractArray) = θ

"""
    log_amplitudes(state, packed_states) -> Vector

Log-amplitudes of `state`'s ansatz on a vector of **packed** configurations.
"""
log_amplitudes(vs::AbstractVariationalState, states::AbstractArray) =
    log_amplitude(ansatz(vs), parameters(vs), configurations_of(vs, states))

"""
Put `mels` on the same device as `like`, moving nothing when it is already there.

Only one of the four combinations needs work: host matrix elements against device
log-amplitudes, which is what the host connected-configuration kernel produces when the ansatz
runs on a device. `copyto!` into a `similar` of the log-amplitudes does that using nothing but
Base, which is why the local-energy kernel needs no GPU dependency to be GPU-ready.

The other three are deliberately left alone. On the host, broadcasting a real matrix-element
array against complex log-amplitudes promotes elementwise for free, and materializing a
converted copy would be a full `(max_conn, batch)` allocation bought for nothing. When *both*
are already on a device — which is what the device kernel gives — a copy would be the same
allocation, paid on the more expensive memory.
"""
_colocate(::Array, mels::Array) = mels
_colocate(like::AbstractArray, mels::Array) =
    copyto!(similar(like, eltype(mels), size(mels)), mels)
_colocate(::AbstractArray, mels::AbstractArray) = mels

"""
    NQSCore.to_host(x) -> AbstractArray

Bring `x` into host memory, leaving it alone when it is already there.

The mirror of `_colocate`, and it exists for the samplers. Accepting or rejecting a Metropolis
proposal, and searching an inverse cumulative distribution, are scalar sequential decisions —
run against a device array they are either an outright error or one round trip per element. The
quantities they branch on are small, one number per chain or per basis state, so fetching them
in a single transfer is the cheap half of that trade.

It does not make a sampler *fast* on a device: a Metropolis sweep still pays one transfer per
step, and the answer to that is a sampler that keeps its chains on the device rather than a
better transfer. It makes one work, and makes the cost measurable.
"""
to_host(x::Array) = x
to_host(x::AbstractArray) = Array(x)

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

When the ansatz's parameters live on a device and KernelAbstractions is loaded, the connected
configurations are computed there too — see [`connections`](@ref). That path re-uploads the
operator on every call unless it is already resident, so a loop can hand over
`to_backend(flatten(H), backend)` once, exactly as it would `compile` on the host. Measured, it
is worth 1.5%: seven small transfers are not much beside a millisecond of network.
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
    x, mels = connections(vs, operator, states, logψ_s)

    # One batched evaluation over every connected configuration of every sample, rather than
    # one call per sample: the ansatz is the expensive part, so it is called once.
    logψ_sp = reshape(log_amplitude(ansatz(vs), parameters(vs), x), size(mels))

    # A whole-column reduction rather than a loop bounded by each sample's connection count.
    # It needs no mask because the padding is already inert: a padded slot repeats the sample
    # itself with a zero matrix element, so it contributes `0 * exp(0) == 0` exactly — never
    # `0 * Inf`. Dropping the data-dependent trip count is also what makes this line run
    # unchanged on a GPU array.
    m = _colocate(logψ_sp, mels)
    return vec(sum(m .* exp.(logψ_sp .- transpose(logψ_s)); dims=1))
end

"""
    connections(state, operator, packed_states, like) -> (x, mels)

The connected configurations of `operator`, as the ansatz wants them, and their matrix elements.

`x` is the `(nsites, max_conn * batch)` numeric array [`log_amplitude`](@ref) consumes; `mels`
is `(max_conn, batch)`, so `size(mels)` is the shape to reshape the log-amplitudes back into.

This is the seam the device path attaches to, and it is a seam rather than a second copy of
[`local_energy`](@ref) because only *these two arrays* differ between host and device. The
reduction that follows is already written in terms that run unchanged on either, and
duplicating it would mean two versions of the one line where the physics is.

`like` says where the answer is wanted: it is the sample log-amplitudes, so it carries both the
memory space the ansatz put itself in and the float type it works in. Passing an array rather
than a device or an element type keeps `NQSCore` free of any notion of either — the host
implementation below ignores it entirely, and what a device *is* is the business of the
extension that KernelAbstractions activates.
"""
connections(vs::AbstractVariationalState, operator, states::AbstractVector, like) =
    _connections(vs, operator, states, like, _device_backend(like))

"""
The KernelAbstractions backend `like` lives on, or `nothing` for host memory.

Always `nothing` here, because without KernelAbstractions loaded there is no backend to name.
The extension replaces this with the real question, and still answers `nothing` for an `Array`:
loading KernelAbstractions must not silently divert host runs onto its CPU backend, which would
swap a tested serial kernel for a launch-per-batch one on the strength of an unrelated `using`.
"""
_device_backend(::Any) = nothing

function _connections(
    vs::AbstractVariationalState, operator, states::AbstractVector, ::AbstractArray, ::Nothing
)
    a = ansatz(vs)
    res = ConnectedBasisConfigurations.connected_padded(operator, states)
    x = ConnectedBasisConfigurations.configurations(dof(a), vec(res.configs), n_sites(a))
    return x, res.mels
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
function local_estimators(
    vs::AbstractVariationalState, operator; holomorphic::Bool=false, chunk_size=nothing
)
    states = samples(vs)
    a, θ = ansatz(vs), parameters(vs)

    x = configurations_of(vs, states)
    logψ = log_amplitude(a, θ, x)
    E = _local_energy(vs, operator, vec(states), logψ)
    O = log_derivatives(
        a, θ, x; backend=vs.backend, holomorphic=holomorphic, chunk_size=chunk_size
    )

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

function expect_and_grad(vs::FullSumState, operator; chunk_size=nothing)
    states = samples(vs)
    a, θ = ansatz(vs), parameters(vs)

    x = configurations_of(vs, states)
    logψ = log_amplitude(a, θ, x)
    E = _local_energy(vs, operator, states, logψ)
    p = born_probabilities(logψ)

    ∇ = energy_gradient(a, θ, x, E, p; backend=vs.backend, chunk_size=chunk_size)
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

function expect_and_grad(vs::MCState, operator; chunk_size=nothing)
    states = samples(vs)
    a, θ = ansatz(vs), parameters(vs)

    x = configurations_of(vs, states)
    logψ = log_amplitude(a, θ, x)
    E = _local_energy(vs, operator, vec(states), logψ)

    ∇ = energy_gradient(a, θ, x, E, nothing; backend=vs.backend, chunk_size=chunk_size)
    return statistics(reshape(E, size(states))), ∇
end
