using NQSCore

using ConnectedBasisConfigurations
using DifferentiationInterface
using ForwardDiff
using LinearAlgebra
using OperatorAlgebra
using Random
using SymBasis
using Test

const BACKEND = AutoForwardDiff()

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

"""Dense Hamiltonian over the full basis, for exact reference answers."""
function dense(H, dof, nsites)
    b = basis(dof_object(dof), nsites)
    index = Dict(s => i for (i, s) in pairs(b.states))
    res = connected_padded(H, b.states)
    M = zeros(ComplexF64, length(b.states), length(b.states))
    for n in eachindex(b.states), j in 1:res.counts[n]
        M[index[res.configs[j, n]], n] += res.mels[j, n]
    end
    return M
end

"""A `FullSumState` on the exact ansatz, with the given parameters."""
function full_sum(dof, nsites, θ=nothing; rng=Xoshiro(0))
    b = basis(dof_object(dof), nsites)
    a = LogStateVector(dof, nsites, b)
    θ === nothing && (θ = init_parameters(a, rng; scale=0.3))
    return FullSumState(a, θ, BACKEND), a, b
end

@testset "NQSCore.jl" begin
    @testset "dependency weight" begin
        # AD backends and Lux are weak dependencies; loading NQSCore must not pull them in.
        loaded = Set(m.name for m in keys(Base.loaded_modules))
        for heavy in ("Lux", "Enzyme", "Reactant", "CUDA", "Metal")
            @test heavy ∉ loaded
        end
    end

    @testset "Stats" begin
        @testset "an exact result has no error bar" begin
            s = exact_stats(1.5, 0.25)
            @test s.mean == 1.5
            @test s.error_of_mean == 0.0
            @test s.r_hat == 1.0
        end

        @testset "independent samples have tau_corr near 1" begin
            x = randn(Xoshiro(42), 20_000)
            @test isapprox(integrated_autocorrelation(x), 1.0; atol=0.3)
        end

        @testset "a correlated chain has tau_corr above 1" begin
            # AR(1) with correlation φ has τ = (1+φ)/(1-φ); for φ=0.8 that is 9.
            rng = Xoshiro(7)
            φ, n = 0.8, 200_000
            x = zeros(n)
            for i in 2:n
                x[i] = φ * x[i-1] + randn(rng)
            end
            @test integrated_autocorrelation(x) > 4
            @test isapprox(integrated_autocorrelation(x), (1 + φ) / (1 - φ); rtol=0.3)
        end

        @testset "split-Rhat flags a drifting chain" begin
            converged = randn(Xoshiro(1), 2000, 4)
            @test isapprox(split_rhat(converged), 1.0; atol=0.05)

            drifting = converged .+ range(0, 8; length=2000)   # a steady trend in every chain
            @test split_rhat(drifting) > 1.1
        end

        @testset "the error of the mean shrinks like 1/sqrt(n)" begin
            rng = Xoshiro(3)
            e1 = statistics(randn(rng, 1000)).error_of_mean
            e2 = statistics(randn(rng, 4000)).error_of_mean
            @test isapprox(e2 / e1, 0.5; rtol=0.5)
        end

        @testset "weighted statistics reproduce a known average" begin
            s = weighted_statistics([1.0, 2.0, 3.0], [0.5, 0.25, 0.25])
            @test s.mean ≈ 1.75
            @test s.variance ≈ 0.5 * 0.75^2 + 0.25 * 0.25^2 + 0.25 * 1.25^2
            @test s.error_of_mean == 0.0
        end
    end

    @testset "log_derivatives" begin
        dof, nsites = Spin(1 // 2), 4
        vs, a, b = full_sum(dof, nsites)
        x = configurations(dof, b.states, nsites)
        n = n_parameters(a)

        @testset "matches the analytic answer for LogStateVector" begin
            # log ψ(s) = θ_{i(s)}, so O is exactly the identity indicator matrix.
            O = log_derivatives(a, parameters(vs), x; backend=BACKEND, holomorphic=true)
            @test size(O) == (length(b.states), n)
            @test O ≈ Matrix{ComplexF64}(I, length(b.states), n)
        end

        @testset "non-holomorphic mode splits real and imaginary parts" begin
            O = log_derivatives(a, parameters(vs), x; backend=BACKEND, holomorphic=false)
            @test size(O) == (length(b.states), 2n)
            @test O[:, 1:n] ≈ Matrix{ComplexF64}(I, length(b.states), n)
            @test O[:, (n+1):(2n)] ≈ im .* Matrix{ComplexF64}(I, length(b.states), n)
        end

        @testset "centering removes the column mean" begin
            O = log_derivatives(a, parameters(vs), x; backend=BACKEND, holomorphic=true)
            @test all(abs.(sum(centered(O); dims=1)) .< 1e-10)

            w = probabilities(vs)
            @test all(abs.(sum(w .* centered(O, w); dims=1)) .< 1e-10)
        end
    end

    @testset "FullSumState is exact" begin
        dof, nsites = Spin(1 // 2), 4
        H = tfi(nsites; J=1.0, h_x=0.7, h_z=0.2)
        M = dense(H, dof, nsites)

        @testset "expectation matches the Rayleigh quotient" begin
            vs, a, b = full_sum(dof, nsites)
            ψ = exp.(log_amplitude(a, parameters(vs), configurations(dof, b.states, nsites)))
            expected = real(dot(ψ, M * ψ) / dot(ψ, ψ))

            s = expect(vs, H)
            @test real(s.mean) ≈ expected
            @test s.error_of_mean == 0.0          # nothing was sampled
        end

        @testset "an eigenstate has zero variance" begin
            # The sharpest available check: E_loc is constant exactly on an eigenstate.
            vals, vecs = eigen(Hermitian(M))
            b = basis(dof_object(dof), nsites)
            a = LogStateVector(dof, nsites, b)
            vs = FullSumState(a, log.(complex.(vecs[:, 1])), BACKEND)

            s = expect(vs, H)
            @test real(s.mean) ≈ vals[1]
            @test s.variance < 1e-16
        end

        @testset "the gradient vanishes at the ground state" begin
            vals, vecs = eigen(Hermitian(M))
            b = basis(dof_object(dof), nsites)
            a = LogStateVector(dof, nsites, b)
            vs = FullSumState(a, log.(complex.(vecs[:, 1])), BACKEND)
            _, ∇ = expect_and_grad(vs, H)
            @test maximum(abs, ∇) < 1e-6
        end

        @testset "the gradient matches a finite difference" begin
            # The parameters are complex, so the gradient has a real and an imaginary part and
            # each must be checked against its own perturbation direction: `real(∇)` against a
            # perturbation of `θ_re`, `imag(∇)` against one of `θ_im`. Comparing the whole
            # complex gradient to a single real-direction difference conflates the two.
            vs, a, b = full_sum(dof, nsites)
            _, ∇ = expect_and_grad(vs, H)
            θ = copy(parameters(vs))
            E(t) = real(expect(FullSumState(a, t, BACKEND), H).mean)
            ε = 1e-6

            for k in (1, 5, 11)
                θp = copy(θ); θp[k] += ε
                θm = copy(θ); θm[k] -= ε
                @test isapprox(real(∇[k]), (E(θp) - E(θm)) / (2ε); atol=1e-5)

                θp = copy(θ); θp[k] += im * ε
                θm = copy(θ); θm[k] -= im * ε
                @test isapprox(imag(∇[k]), (E(θp) - E(θm)) / (2ε); atol=1e-5)
            end
        end
    end

    @testset "gradient descent reaches the exact ground state" begin
        # The end-to-end statement for this package: with an ansatz that can represent any
        # state, optimizing it must reproduce exact diagonalization. Anything less would mean
        # the energy, the gradient, or the log-derivatives are wrong.
        #
        # The transverse field is deliberately in the delocalized regime (h_x >= 2J). Plain
        # gradient descent on log-amplitudes stalls in the ordered regime, and not because of
        # a bug: the weighted gradient carries a factor of the Born probability p(s), which
        # goes to zero for exactly the configurations whose amplitude needs to grow. Raising
        # the learning rate helps, which is the signature of a vanishing gradient rather than
        # an overshooting one. Undoing that ill-conditioning is what stochastic reconfiguration
        # is for, and testing the ordered regime belongs with it in NQSOptimisers.
        dof, nsites = Spin(1 // 2), 4
        for (J, h_x) in ((1.0, 2.0), (0.5, 2.0), (1.0, 3.0))
            H = tfi(nsites; J=J, h_x=h_x)
            E_exact = minimum(real(eigvals(Hermitian(dense(H, dof, nsites)))))

            vs, _, _ = full_sum(dof, nsites; rng=Xoshiro(11))
            for _ in 1:5000
                _, ∇ = expect_and_grad(vs, H)
                setparameters!(vs, parameters(vs) .- 0.05 .* ∇)
            end

            s = expect(vs, H)
            # Converges to machine precision, not merely to a loose tolerance.
            @test isapprox(real(s.mean), E_exact; atol=1e-9)
            # ...and the energy variance vanishes, which is the definition of an eigenstate.
            @test s.variance < 1e-14
        end
    end

    @testset "MCState agrees with FullSumState" begin
        dof, nsites = Spin(1 // 2), 4
        H = tfi(nsites; J=1.0, h_x=0.7, h_z=0.2)
        b = basis(dof_object(dof), nsites)
        a = LogStateVector(dof, nsites, b)
        θ = init_parameters(a, Xoshiro(5); scale=0.3)

        exact = expect(FullSumState(a, θ, BACKEND), H)
        mc = MCState(a, θ, ExactSampler(b, 200_000); backend=BACKEND, rng=Xoshiro(5))
        sampled = expect(mc, H)

        @testset "the estimate agrees within its own error bar" begin
            @test sampled.error_of_mean > 0
            @test abs(real(sampled.mean - exact.mean)) < 5 * sampled.error_of_mean
        end

        @testset "independent samples show no autocorrelation" begin
            # ExactSampler draws independently, so much above 1 would be an estimator bug
            # rather than a property of the chain.
            @test sampled.tau_corr < 1.5
        end

        @testset "changing parameters invalidates the sample cache" begin
            n_before = length(samples(mc))
            setparameters!(mc, θ .+ 0.1)
            @test mc.cache === nothing
            @test length(samples(mc)) == n_before
        end
    end

    @testset "interface conformance" begin
        dof, nsites = Spin(1 // 2), 3
        vs, a, b = full_sum(dof, nsites)
        H = tfi(nsites)
        for state in (vs, MCState(a, parameters(vs), ExactSampler(b, 64); backend=BACKEND))
            @test ansatz(state) === a
            @test parameters(state) == parameters(vs)
            @test length(samples(state)) > 0
            @test expect(state, H) isa Stats
        end
    end
end
