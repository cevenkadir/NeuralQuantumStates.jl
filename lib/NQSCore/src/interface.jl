"""
    AbstractAnsatz

A parametrized wavefunction: something that maps a batch of configurations to log-amplitudes.
The one required method is [`log_amplitude`](@ref).

Parameters are *not* stored in the ansatz. An ansatz describes the functional form; the
parameters live in the [`AbstractVariationalState`](@ref) that owns it, which is what makes it
possible to evaluate the same ansatz at perturbed parameters.
"""
abstract type AbstractAnsatz end

"""
    NQSCore.dof(ansatz)

The SymBasis degree-of-freedom specification the ansatz is defined over. Defaults to the `dof`
field.

Not exported, and neither is [`NQSCore.n_sites`](@ref): `n_sites` is a lattice's word as much as
an ansatz's, and a package meant to be depended on should not claim it for everyone downstream.
Extend and call them qualified.
"""
dof(a::AbstractAnsatz) = a.dof

"""
    NQSCore.n_sites(ansatz) -> Int

Number of sites the ansatz is defined on. Defaults to the `nsites` field; see
[`NQSCore.dof`](@ref), which explains why neither is exported.
"""
n_sites(a::AbstractAnsatz) = a.nsites

"""
    NQSCore.input_type(ansatz, parameters) -> Union{Type,Nothing}

The element type this ansatz does its arithmetic in, or `nothing` for no preference.

`nothing` is the default and means "hand me the configurations as they naturally come" — exact
rationals for a spin. An ansatz that answers with a type gets a batch already in it.

A network with complex parameters should answer with them. Feeding it a real batch does not
avoid the promotion, it moves the promotion inside every matrix product, where BLAS has no
kernel for a mixed complex-real pair and the operation falls back to a generic one. That is
cheap in the forward direction, whose product is wide and shallow, and very expensive in the
pullback, whose product has few outputs over a long reduction.

Only [`configurations_of`](@ref) honours the answer, because that is the batch
`energy_gradient` and `log_derivatives` are handed. The connected configurations from
[`connections`](@ref) and the batches a sampler evaluates are never differentiated, so widening
them would be pure cost.
"""
input_type(::AbstractAnsatz, θ) = nothing

"""
    log_amplitude(ansatz, parameters, x) -> AbstractVector

Log-amplitudes `log ψ(x)` for a batch of configurations.

`x` is a `(nsites, batch)` array of physical local values, as produced by
`ConnectedBasisConfigurations.configurations`. The result is a length-`batch` vector, generally
complex: real part the log-modulus, imaginary part the phase.

Returning the logarithm is what makes this numerically usable — amplitudes underflow for any
interesting system size, while the differences of logarithms that local energies need stay
bounded.
"""
function log_amplitude end

"""
    AbstractVariationalState

An [`AbstractAnsatz`](@ref), a set of parameters, and a way of estimating expectation values.

Two implementations ship here and differ only in how they average: [`FullSumState`](@ref) sums
exactly over the whole basis, and [`MCState`](@ref) estimates by Monte Carlo. Both satisfy the
same interface, so a model can be developed against exact summation and scaled up by swapping
the state type.

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

The configurations `state` currently estimates expectation values from, as packed states: the
drawn samples for an [`MCState`](@ref), the whole basis for a [`FullSumState`](@ref).
"""
function samples end

"""
    expect(state, operator) -> Stats

Expectation value `⟨ψ|Ô|ψ⟩ / ⟨ψ|ψ⟩`, with an uncertainty estimate.

Always a [`Stats`](@ref), whichever state produced it. A `FullSumState` reports a zero error bar
because its answer is exact; an `MCState` reports the standard error of the mean along with the
diagnostics needed to judge whether that error bar can be believed.
"""
function expect end

"""
    expect_and_grad(state, operator) -> (Stats, gradient)

[`expect`](@ref) together with the gradient of that expectation value with respect to the
variational parameters:

```math
\\partial_k \\langle E \\rangle = 2 \\, \\mathrm{Re}
    \\left[ \\langle O_k^* E_{\\mathrm{loc}} \\rangle
          - \\langle O_k^* \\rangle \\langle E_{\\mathrm{loc}} \\rangle \\right]
```

with `O_k` the log-derivatives from [`log_derivatives`](@ref). Subtracting
`⟨O_k^*⟩⟨E_loc⟩` is not cosmetic: without it the estimator has a non-vanishing variance even at
an exact eigenstate.

`O` is never built — the gradient is one contraction of it, obtained as the gradient of a
scalar. Pass `chunk_size` to bound the memory of that differentiation pass; the answer does not
change.
"""
function expect_and_grad end

"""
    AbstractSampler

A way of drawing configurations distributed according to `|ψ|²`. The interface is
[`sample`](@ref). `NQSSamplers` provides the production samplers; the only one here is
[`ExactSampler`](@ref), which enumerates the basis and is meant for testing.
"""
abstract type AbstractSampler end

"""
    sample(sampler, ansatz, parameters, rng, state=nothing) -> (samples, sampler_state)

Draw configurations distributed according to `|ψ|²`, as packed states, together with the
sampler's own state.

`state` is the `sampler_state` from a previous call, or `nothing` to start from scratch.
Handing it back lets a Markov chain resume instead of re-paying its burn-in every optimization
step; consecutive steps differ by one small parameter update, so the previous chain is already
nearly equilibrated for the new parameters. A sampler that draws independently returns
`nothing`, which callers must accept.

Samples may be returned as a `(steps, chains)` matrix; that shape survives into
[`statistics`](@ref), which needs it for split-R̂ and for a between-chain error bar.
"""
function sample end

"""
    AbstractPreconditioner

A transformation applied to the raw gradient before it reaches an optimizer — stochastic
reconfiguration, its kernel-trick variant, or plain identity. `NQSOptimisers` provides the
implementations; the type lives here so a driver can accept one without depending on that
package.
"""
abstract type AbstractPreconditioner end

"""
    precondition(preconditioner, state, operator) -> (Stats, update)

Transform the raw energy gradient into the parameter update an optimizer should apply.
Declared next to [`AbstractPreconditioner`](@ref) so that the type and its verb travel together;
`NQSOptimisers` supplies the methods.
"""
function precondition end

function local_energy end
