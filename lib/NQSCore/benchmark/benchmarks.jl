"""
Benchmarks for the variational-state machinery: local energies, log-derivatives, energy
gradients, and the statistics reductions.

Run with

```
julia --project=lib/NQSCore/benchmark lib/NQSCore/benchmark/benchmarks.jl
```

These are not wired into CI. They exist so that a change to `expect`, `expect_and_grad` or the
sampling loop can be measured rather than argued about, and so that the cost of handing an
*uncompiled* operator to the local-energy kernel stays visible.

The script deliberately avoids the parts of the API that changed in the 0.1 refactor, going
through the two adapters below instead, so that the same file can be run against an older
checkout to produce a comparable baseline.
"""

using BenchmarkTools
using ConnectedBasisConfigurations
using DifferentiationInterface
using ForwardDiff
using NQSCore
using OperatorAlgebra
using Printf
using Random
using SymBasis

const BACKEND = AutoForwardDiff()

# ------------------------------------------------------------------------------- adapters

"""
`FullSumState` took its backend positionally before the 0.1 refactor and by keyword after it.
Detected once, by the shape of the two-positional-argument method.
"""
const BACKEND_IS_KEYWORD = hasmethod(FullSumState, Tuple{NQSCore.AbstractAnsatz,Any})

function fullsum(a, θ; basis=nothing)
    if BACKEND_IS_KEYWORD
        return FullSumState(a, θ; backend=BACKEND, basis=basis)
    else
        return FullSumState(a, θ, BACKEND; basis=basis)
    end
end

"""`sample` returned a bare vector before the refactor and `(samples, sampler_state)` after."""
drawn(x) = x isa Tuple ? first(x) : x

# -------------------------------------------------------------------------------- models

"""Transverse-field Ising on a periodic chain, as an `OpSum`."""
function tfi(nsites; J=1.0, h_x=1.0, h_z=0.0)
    ops = local_operators(Spin(1 // 2))
    σz, σx = 2 .* ops.sz, 2 .* ops.sx
    terms = AbstractOp[]
    for i in 1:nsites
        push!(terms, J * (Op(σz, i) * Op(σz, mod1(i + 1, nsites))))
        iszero(h_z) || push!(terms, h_z * Op(σz, i))
        iszero(h_x) || push!(terms, h_x * Op(σx, i))
    end
    return OpSum(terms)
end

"""
A restricted-Boltzmann-style ansatz whose parameters are a `NamedTuple` of real arrays.

Its point in this suite is the parameter container, not the physics: it is the shape Lux hands
back, so it exercises the `ComponentArrays` flattening path and the real-parameter branch of
[`log_derivatives`](@ref) — the two things `LogStateVector` never touches.
"""
struct ToyRBM{D} <: NQSCore.AbstractAnsatz
    dof::D
    nsites::Int
    nhidden::Int
end

function NQSCore.log_amplitude(a::ToyRBM, θ, x::AbstractMatrix)
    h = θ.W * x .+ θ.b
    logmod = vec(sum(log.(2 .* cosh.(h)); dims=1))
    phase = vec(transpose(x) * θ.v)
    return complex.(logmod, phase)
end

function toy_parameters(a::ToyRBM, rng)
    return (
        W=0.05 .* randn(rng, a.nhidden, a.nsites),
        b=0.05 .* randn(rng, a.nhidden),
        v=0.05 .* randn(rng, a.nsites),
    )
end

report(label, trial) = @printf("%-52s %10s\n", label, BenchmarkTools.prettytime(minimum(trial).time))

# The suite is dominated by a handful of slow automatic-differentiation entries; a shorter
# budget than the default keeps a full run to a couple of minutes without changing the minima,
# which is what is being reported.
BenchmarkTools.DEFAULT_PARAMETERS.seconds = 2

const SUITE = BenchmarkGroup()

# ------------------------------------------------ exact ansatz, full summation over the basis

# Small enough that the O(n_params × n_samples) Jacobian is affordable: the log-derivative
# entries below scale as the square of the basis dimension.
let rng = Xoshiro(0), nsites = 8
    dof = Spin(1 // 2)
    b = basis(dof_object(dof), nsites)
    a = LogStateVector(dof, nsites, b)
    θ = 0.3 .* randn(rng, ComplexF64, length(b.states))

    H = tfi(nsites; h_x=0.7, h_z=0.2)
    compiled = compile(H)
    vs = fullsum(a, θ)
    x = configurations(dof, b.states, nsites)

    SUITE["fullsum8"]["expect, compiling per call"] = @benchmarkable expect($vs, $H)
    SUITE["fullsum8"]["expect, precompiled"] = @benchmarkable expect($vs, $compiled)
    SUITE["fullsum8"]["expect_and_grad, compiling per call"] =
        @benchmarkable expect_and_grad($vs, $H)
    SUITE["fullsum8"]["expect_and_grad, precompiled"] =
        @benchmarkable expect_and_grad($vs, $compiled)
    SUITE["fullsum8"]["log_derivatives, complex theta"] =
        @benchmarkable log_derivatives($a, $θ, $x; backend=BACKEND, holomorphic=true)
    SUITE["fullsum8"]["log_derivatives, chunk_size=32"] = @benchmarkable log_derivatives(
        $a, $θ, $x; backend=BACKEND, holomorphic=true, chunk_size=32
    )
    SUITE["fullsum8"]["expect_and_grad, chunk_size=32"] =
        @benchmarkable expect_and_grad($vs, $compiled; chunk_size=32)
    SUITE["fullsum8"]["probabilities"] = @benchmarkable probabilities($vs)
end

# Bigger basis, no automatic differentiation: this is where compiling the operator once rather
# than per call is supposed to show up.
let rng = Xoshiro(5), nsites = 12
    dof = Spin(1 // 2)
    b = basis(dof_object(dof), nsites)
    a = LogStateVector(dof, nsites, b)
    θ = 0.3 .* randn(rng, ComplexF64, length(b.states))

    H = tfi(nsites; h_x=0.7, h_z=0.2)
    compiled = compile(H)
    vs = fullsum(a, θ)
    states = b.states

    SUITE["fullsum12"]["expect, compiling per call"] = @benchmarkable expect($vs, $H)
    SUITE["fullsum12"]["expect, precompiled"] = @benchmarkable expect($vs, $compiled)
    SUITE["fullsum12"]["local_energy, compiling per call"] =
        @benchmarkable local_energy($vs, $H, $states)
    SUITE["fullsum12"]["local_energy, precompiled"] =
        @benchmarkable local_energy($vs, $compiled, $states)
    SUITE["fullsum12"]["probabilities"] = @benchmarkable probabilities($vs)
end

# ---------------------------------------------------- NamedTuple parameters, real-theta path

let rng = Xoshiro(1), nsites = 8, nhidden = 8
    dof = Spin(1 // 2)
    b = basis(dof_object(dof), nsites)
    a = ToyRBM(dof, nsites, nhidden)
    θ = toy_parameters(a, rng)

    H = tfi(nsites; h_x=0.7)
    compiled = compile(H)
    vs = fullsum(a, θ; basis=b)
    x = configurations(dof, b.states, nsites)

    SUITE["rbm8"]["expect, precompiled"] = @benchmarkable expect($vs, $compiled)
    SUITE["rbm8"]["expect_and_grad, precompiled"] = @benchmarkable expect_and_grad($vs, $compiled)
    SUITE["rbm8"]["log_derivatives, NamedTuple theta"] =
        @benchmarkable log_derivatives($a, $θ, $x; backend=BACKEND)
    SUITE["rbm8"]["flatten_parameters"] = @benchmarkable flatten_parameters($θ)
end

# ----------------------------------------------------------------- Monte Carlo sampling loop

let nsites = 8
    dof = Spin(1 // 2)
    b = basis(dof_object(dof), nsites)
    a = LogStateVector(dof, nsites, b)
    θ = 0.3 .* randn(Xoshiro(2), ComplexF64, length(b.states))

    H = tfi(nsites; h_x=0.7)
    compiled = compile(H)
    sampler = ExactSampler(b, 4096)

    SUITE["mc8"]["sample"] =
        @benchmarkable drawn(NQSCore.sample($sampler, $a, $θ, $(Xoshiro(3))))
    SUITE["mc8"]["expect, precompiled"] = @benchmarkable expect(
        MCState($a, $θ, $sampler; backend=BACKEND, rng=Xoshiro(3)), $compiled
    )
    SUITE["mc8"]["expect_and_grad, precompiled"] = @benchmarkable expect_and_grad(
        MCState($a, $θ, $sampler; backend=BACKEND, rng=Xoshiro(3)), $compiled
    )
    SUITE["mc8"]["local_estimators, precompiled"] = @benchmarkable local_estimators(
        MCState($a, $θ, $sampler; backend=BACKEND, rng=Xoshiro(3)), $compiled; holomorphic=true
    )
end

# ------------------------------------------------------------------------ statistics kernels

let rng = Xoshiro(4)
    chains = randn(rng, 4096, 8)
    single = randn(rng, 100_000)
    weights = abs2.(randn(rng, 4096)); weights ./= sum(weights)
    values = randn(rng, ComplexF64, 4096)

    SUITE["stats"]["statistics, 4096x8"] = @benchmarkable statistics($chains)
    SUITE["stats"]["statistics, 4096 single chain"] = @benchmarkable statistics($(chains[:, 1]))
    SUITE["stats"]["integrated_autocorrelation, 1e5"] =
        @benchmarkable integrated_autocorrelation($single)
    SUITE["stats"]["split_rhat, 4096x8"] = @benchmarkable split_rhat($chains)
    SUITE["stats"]["weighted_statistics, 4096"] = @benchmarkable weighted_statistics($values, $weights)
end

if abspath(PROGRAM_FILE) == @__FILE__
    results = run(SUITE; verbose=false)
    for (group, trials) in sort(collect(results); by=first)
        println("\n", group)
        for (name, trial) in sort(collect(trials); by=first)
            report("  " * name, trial)
        end
    end
end
