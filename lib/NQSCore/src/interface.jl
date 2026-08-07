"""
    AbstractAnsatz

A parametrized wavefunction: something that maps a batch of configurations to log-amplitudes.

The single method every ansatz must provide is [`log_amplitude`](@ref). Keeping the interface
this narrow is deliberate — it is what lets `NQSAnsatze` build Lux models, `NQSSamplers` draw
from `|ψ|²`, and `NQSOptimisers` precondition gradients without any of them depending on each
other.

Parameters are *not* stored in the ansatz. An ansatz describes the functional form; the
parameters live in the [`AbstractVariationalState`](@ref) that owns it. That separation is what
makes it possible to evaluate the same ansatz at perturbed parameters — which is exactly what
computing a Jacobian does.
"""
abstract type AbstractAnsatz end

"""
    NQSCore.dof(ansatz)

The SymBasis degree-of-freedom specification the ansatz is defined over.

Defaults to the `dof` field, so an ansatz that stores one needs no method. Everything that
builds configurations goes through this rather than reaching for the field directly, which is
what lets an ansatz compute its degrees of freedom instead of storing them.

Deliberately **not exported**, and neither is [`NQSCore.n_sites`](@ref): `n_sites` is a lattice's
word as much as an ansatz's — `LatticeSpaceGroups` exports its own — and a package whose job is
to be depended on should not make that choice for everyone downstream. Extend and call them
qualified.
"""
dof(a::AbstractAnsatz) = a.dof

"""
    NQSCore.n_sites(ansatz) -> Int

Number of sites the ansatz is defined on. Defaults to the `nsites` field; see
[`NQSCore.dof`](@ref), which explains why neither is exported.
"""
n_sites(a::AbstractAnsatz) = a.nsites

"""
    log_amplitude(ansatz, parameters, x) -> AbstractVector

Log-amplitudes `log ψ(x)` for a batch of configurations.

`x` is a `(nsites, batch)` array of **physical local values** — magnetic quantum numbers for a
spin, occupation numbers for a boson — as produced by
`ConnectedBasisConfigurations.configurations`. The result is a length-`batch` vector, generally
complex: the real part carries the log-modulus and the imaginary part the phase.

Returning the *logarithm* rather than the amplitude is what makes this numerically usable at
all: amplitudes underflow catastrophically for any interesting system size, whereas the
differences of logarithms that local energies actually need stay bounded.
"""
function log_amplitude end

"""
    AbstractVariationalState

A variational wavefunction: an [`AbstractAnsatz`](@ref), a set of parameters, and a way of
estimating expectation values.

Two implementations ship here, and they differ only in how they average:

- [`FullSumState`](@ref) sums exactly over the whole basis. No sampling, no error bars, no
  Markov chain — feasible only for small systems, and indispensable for testing, since any
  disagreement with it is a bug rather than noise.
- [`MCState`](@ref) estimates the same quantities by Monte Carlo.

Both satisfy the same interface, so a model can be developed against exact summation and then
scaled up by swapping the state type and changing nothing else.

# Required methods
[`parameters`](@ref), [`setparameters!`](@ref), [`ansatz`](@ref), [`expect`](@ref),
[`expect_and_grad`](@ref), and `samples`.
"""
abstract type AbstractVariationalState end

"""
    parameters(state)

The variational parameters of `state`.
"""
function parameters end

"""
    setparameters!(state, θ) -> state

Replace the variational parameters of `state`, invalidating anything cached from the old ones.
"""
function setparameters! end

"""
    ansatz(state) -> AbstractAnsatz

The ansatz `state` is built on.
"""
function ansatz end

"""
    samples(state) -> AbstractVector

The configurations `state` currently estimates expectation values from, as packed states.

For an [`MCState`](@ref) these are the drawn Monte Carlo samples; for a [`FullSumState`](@ref)
they are every state in the basis.
"""
function samples end

"""
    expect(state, operator) -> Stats

Expectation value `⟨ψ|Ô|ψ⟩ / ⟨ψ|ψ⟩`, with an uncertainty estimate.

Always returns a [`Stats`](@ref), whichever state type produced it. A `FullSumState` reports a
zero error bar because its answer is exact; an `MCState` reports the standard error of the mean
along with the autocorrelation diagnostics needed to judge whether that error bar can be
believed.
"""
function expect end

"""
    expect_and_grad(state, operator) -> (Stats, gradient)

[`expect`](@ref) together with the gradient of that expectation value with respect to the
variational parameters.

The gradient is the usual variational-Monte-Carlo estimator

```math
\\partial_k \\langle E \\rangle = 2 \\, \\mathrm{Re}
    \\left[ \\langle O_k^* E_{\\mathrm{loc}} \\rangle
          - \\langle O_k^* \\rangle \\langle E_{\\mathrm{loc}} \\rangle \\right]
```

with `O_k` the log-derivatives from [`log_derivatives`](@ref). The subtraction of
`⟨O_k^*⟩⟨E_loc⟩` is not cosmetic: without it the estimator has a non-vanishing variance even at
an exact eigenstate.
"""
function expect_and_grad end

"""
    AbstractSampler

A way of drawing configurations distributed according to `|ψ|²`.

The interface is [`sample`](@ref). `NQSSamplers` provides Metropolis and autoregressive
samplers; the only one here is [`ExactSampler`](@ref), which enumerates the basis and is
meant for testing rather than for production use.
"""
abstract type AbstractSampler end

"""
    sample(sampler, ansatz, parameters, rng, state=nothing) -> (samples, sampler_state)

Draw configurations distributed according to `|ψ|²`, returned as packed states together with
the sampler's own state.

`state` is the `sampler_state` returned by a previous call, or `nothing` to start from scratch.
Handing it back is what lets a Markov chain **resume where it left off** instead of restarting
and re-paying its burn-in on every optimization step. Over a run that is the difference between
paying the equilibration cost once and paying it thousands of times — and consecutive steps
differ by one small parameter update, so the previous chain is already very nearly equilibrated
for the new parameters.

A sampler that draws independently — [`ExactSampler`](@ref) is the one here — has nothing to
carry between calls and returns `nothing`, which callers must accept.

Samples may be returned as a `(steps, chains)` matrix; that shape survives into
[`statistics`](@ref), which needs it for split-R̂ and for a between-chain error bar.
"""
function sample end

"""
    AbstractPreconditioner

A transformation applied to the raw gradient before it is handed to an optimizer — stochastic
reconfiguration, its kernel-trick variant, or plain identity.

`NQSOptimisers` provides the implementations; the type lives here so that a driver can accept
one without depending on that package.
"""
abstract type AbstractPreconditioner end

"""
    precondition(preconditioner, state, operator) -> (Stats, update)

Transform the raw energy gradient into the parameter update an optimizer should apply.

Declared here, next to [`AbstractPreconditioner`](@ref), so that the type and its verb travel
together: a driver written against this signature works with any implementation without
depending on the package that provides it. `NQSOptimisers` supplies the methods.
"""
function precondition end

"""
    local_energy(state, operator[, states]) -> AbstractVector

Local energies `E_loc(s) = Σ_{s'} ⟨s|Ô|s'⟩ ψ(s')/ψ(s)` for the given configurations, defaulting
to the state's current [`samples`](@ref).

This is where the variational method earns its keep: the sum runs only over configurations the
operator actually connects to `s`, which for a local Hamiltonian is a handful, not the whole
Hilbert space.
"""
function local_energy end
