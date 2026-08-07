using NQSSamplers

using ConnectedBasisConfigurations
using DifferentiationInterface
using ForwardDiff
using LinearAlgebra
using NQSCore
using OperatorAlgebra
using Random
using SymBasis
using Test

const BACKEND = AutoForwardDiff()

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

"""A fixed, non-uniform target wavefunction to sample from."""
function fixture(nsites=4; seed=0, symmetry=nothing)
    dof = Spin(1 // 2)
    dofo = dof_object(dof)
    b = symmetry === nothing ? basis(dofo, nsites) : basis(dofo, nsites, symmetry)
    a = LogStateVector(dof, nsites, b)
    θ = init_parameters(a, Xoshiro(seed); scale=0.6)
    return dof, b, a, θ
end

"""Empirical distribution of packed states over a basis."""
function empirical(states, b)
    index = Dict(s => i for (i, s) in pairs(b.states))
    counts = zeros(Int, length(b.states))
    for s in vec(states)
        counts[index[s]] += 1
    end
    return counts ./ sum(counts)
end

@testset "NQSSamplers.jl" begin
    @testset "dependency weight" begin
        loaded = Set(m.name for m in keys(Base.loaded_modules))
        for heavy in ("Lux", "Enzyme", "Reactant", "CUDA", "Metal")
            @test heavy ∉ loaded
        end
    end

    @testset "rules" begin
        dof, b, a, θ = fixture(4)
        rng = Xoshiro(1)

        @testset "LocalRule changes exactly one site" begin
            for s in b.states[1:8]
                s′, correction = propose(LocalRule(), s, dof, 4, rng)
                before = configurations(dof, s, 4)
                after = configurations(dof, s′, 4)
                @test count(before .!= after) == 1
                @test correction == 0.0              # symmetric proposal
            end
        end

        @testset "ExchangeRule conserves the total" begin
            for s in b.states
                s′, correction = propose(ExchangeRule(), s, dof, 4, rng)
                @test sum(configurations(dof, s′, 4)) == sum(configurations(dof, s, 4))
                @test correction == 0.0
            end
        end

        @testset "HamiltonianRule proposes only connected configurations" begin
            H = tfi(4; h_x=1.0)
            rule = HamiltonianRule(H)
            for s in b.states[1:8]
                s′, _ = propose(rule, s, dof, 4, rng)
                s′ == s && continue
                @test haskey(connected(H, s), s′)
            end
        end
    end

    @testset "Metropolis reproduces the exact distribution" begin
        # The decisive test: Metropolis must converge to |ψ|², which is known exactly here.
        # Anything wrong in the acceptance test, the proposal correction, or the burn-in shows
        # up as a mismatch between the sampled and exact distributions.
        dof, b, a, θ = fixture(4; seed=2)
        exact = NQSCore.probabilities(FullSumState(a, θ; backend=BACKEND))

        for rule in (LocalRule(), HamiltonianRule(tfi(4; h_x=1.0)))
            starts = random_configurations(dof, 4, 6, Xoshiro(3))
            sampler = MetropolisSampler(rule, starts;
                n_chains=6, n_samples=40_000, burn_in=2_000, thinning=2)
            drawn, _ = NQSCore.sample(sampler, a, θ, Xoshiro(4))

            tv = 0.5 * sum(abs.(empirical(drawn, b) .- exact))   # total variation distance
            @test tv < 0.02
        end
    end

    @testset "the sampler returns a (steps, chains) array" begin
        dof, b, a, θ = fixture(4)
        starts = random_configurations(dof, 4, 5, Xoshiro(0))
        sampler = MetropolisSampler(LocalRule(), starts;
            n_chains=5, n_samples=200, burn_in=50)
        @test size(first(NQSCore.sample(sampler, a, θ, Xoshiro(1)))) == (200, 5)
    end

    @testset "chains resume from the returned sampler state" begin
        dof, b, a, θ = fixture(4; seed=11)
        starts = random_configurations(dof, 4, 4, Xoshiro(0))
        sampler = MetropolisSampler(LocalRule(), starts;
            n_chains=4, n_samples=100, burn_in=500)

        drawn, state = NQSCore.sample(sampler, a, θ, Xoshiro(1))
        @test length(state) == 4
        @test state == vec(drawn[end, :])       # where each chain finished

        @testset "a resumed run picks up from there" begin
            # Handing the state back skips burn-in, so the first kept sample is one Metropolis
            # step from where the previous run ended rather than 500 steps from `starts`.
            again, state2 = NQSCore.sample(sampler, a, θ, Xoshiro(2), state)
            @test size(again) == (100, 4)
            @test length(state2) == 4
            for c in 1:4
                # One step can move a chain or leave it; either way it cannot have travelled
                # further than a single local flip from where it resumed.
                s0 = configurations(dof, state[c], 4)
                s1 = configurations(dof, again[1, c], 4)
                @test count(s0 .!= s1) <= 1
            end
        end

        @testset "a state of the wrong length is rejected" begin
            @test_throws ArgumentError NQSCore.sample(sampler, a, θ, Xoshiro(3), state[1:2])
        end

        @testset "an MCState keeps the chain warm across a parameter change" begin
            mc = MCState(a, θ, sampler; backend=BACKEND, rng=Xoshiro(4))
            first_state = NQSCore.sampler_state(mc)
            @test first_state !== nothing
            setparameters!(mc, θ .+ 0.01)
            @test NQSCore.sampler_state(mc) === first_state   # kept, not discarded
            samples(mc)                                       # redraw at the new parameters
            @test NQSCore.sampler_state(mc) !== first_state   # ...and it advanced
        end
    end

    @testset "MCState with Metropolis agrees with FullSumState" begin
        nsites = 4
        dof, b, a, θ = fixture(nsites; seed=5)
        H = tfi(nsites; J=1.0, h_x=0.8, h_z=0.2)

        exact = expect(FullSumState(a, θ; backend=BACKEND), H)

        starts = random_configurations(dof, nsites, 8, Xoshiro(6))
        sampler = MetropolisSampler(LocalRule(), starts;
            n_chains=8, n_samples=20_000, burn_in=2_000)
        sampled = expect(MCState(a, θ, sampler; backend=BACKEND, rng=Xoshiro(7)), H)

        @test sampled.error_of_mean > 0
        @test abs(real(sampled.mean - exact.mean)) < 5 * sampled.error_of_mean

        @testset "the chains agree with each other" begin
            @test sampled.r_hat < 1.05
        end

        @testset "acceptance is in a usable range" begin
            # Near 0 means a stuck chain, near 1 means proposals too timid to explore. Either
            # gives confident-looking error bars from samples carrying no information.
            @test 0.05 < ACCEPTANCE[] < 0.98
        end
    end

    @testset "sector-confined sampling" begin
        nsites = 6
        dof = Spin(1 // 2)
        dofo = dof_object(dof)
        b = basis(dofo, nsites, sym(TotalMagnetization(0 // 1, nsites), dofo))
        a = LogStateVector(dof, nsites, b)
        θ = init_parameters(a, Xoshiro(8); scale=0.5)
        starts = fill(first(b.states), 4)

        @testset "ExchangeRule stays in the sector and mixes" begin
            sampler = MetropolisSampler(ExchangeRule(), starts;
                n_chains=4, n_samples=20_000, burn_in=1_000, basis=b)
            drawn, _ = NQSCore.sample(sampler, a, θ, Xoshiro(9))

            @test all(s in b.states for s in vec(drawn))
            @test all(sum(configurations(dof, s, nsites)) == 0 for s in vec(drawn))
            # It must actually explore the sector, not merely stay legal within it.
            @test length(unique(vec(drawn))) == length(b.states)

            exact = NQSCore.probabilities(FullSumState(a, θ; backend=BACKEND))
            @test 0.5 * sum(abs.(empirical(drawn, b) .- exact)) < 0.02
        end

        @testset "LocalRule cannot move inside a conserved sector" begin
            # Documents the failure mode rather than leaving it to be rediscovered: every
            # single-site change alters the magnetization and is rejected by the basis check,
            # so the chain returns one configuration forever.
            sampler = MetropolisSampler(LocalRule(), starts;
                n_chains=4, n_samples=500, burn_in=100, basis=b)
            drawn, _ = NQSCore.sample(sampler, a, θ, Xoshiro(10))
            @test length(unique(vec(drawn))) == 1
        end
    end

    @testset "input validation" begin
        dof, b, a, θ = fixture(4)
        starts = random_configurations(dof, 4, 3, Xoshiro(0))
        @test_throws ArgumentError MetropolisSampler(LocalRule(), starts; n_chains=5)
        @test_throws ArgumentError MetropolisSampler(LocalRule(), starts; n_chains=3, n_samples=0)
        @test_throws ArgumentError MetropolisSampler(LocalRule(), starts; n_chains=3, thinning=0)
    end
end
