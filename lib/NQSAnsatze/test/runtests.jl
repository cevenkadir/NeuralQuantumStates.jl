using NQSAnsatze

using ConnectedConfigs
using DifferentiationInterface
using LatticeSpaceGroups
using LinearAlgebra
using Lux
using NQSCore
using OperatorAlgebra
using Random
using SymBasis
using Test
using Zygote

const BACKEND = AutoZygote()

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

@testset "NQSAnsatze.jl" begin
    @testset "logtwocosh does not overflow" begin
        # A literal log(2cosh(z)) returns Inf past |Re z| ~ 710, which an RBM reaches once its
        # weights grow. One Inf poisons an entire batch.
        for x in (0.0, 1.0, 100.0, 800.0, 5000.0, -5000.0)
            @test isfinite(logtwocosh(x))
            @test isfinite(logtwocosh(x + 0.7im))
        end
        @test log(2cosh(800.0)) == Inf          # the naive form really does blow up
        @test logtwocosh(800.0) ≈ 800.0 atol = 1e-8
        for x in (0.0, 0.5, 2.0, 10.0)
            @test logtwocosh(x) ≈ log(2cosh(x))
        end
    end

    @testset "layers evaluate on a batch" begin
        rng = Xoshiro(0)
        nsites, batch = 6, 5
        x = randn(rng, Float64, nsites, batch)

        for layer in (RBM(nsites, 2), Jastrow(nsites),
                      SymmetricRBM([circshift(1:nsites, k) for k in 0:(nsites-1)], 2))
            ps = Lux.initialparameters(rng, layer)
            st = Lux.initialstates(rng, layer)
            y, _ = layer(x, ps, st)
            @test length(y) == batch
            @test eltype(y) <: Complex
            @test all(isfinite, y)
            @test Lux.parameterlength(layer) == sum(length, values(ps))
        end
    end

    @testset "SymmetricRBM is invariant under its group" begin
        # The defining property: permuting the input by a group element must leave log ψ
        # unchanged. If it does not, the weight sharing is wrong.
        rng = Xoshiro(1)
        nsites = 6
        perms = [collect(circshift(1:nsites, k)) for k in 0:(nsites-1)]
        layer = SymmetricRBM(perms, 2)
        ps = Lux.initialparameters(rng, layer)
        st = Lux.initialstates(rng, layer)

        x = randn(rng, Float64, nsites, 4)
        y, _ = layer(x, ps, st)
        for p in perms
            yp, _ = layer(x[p, :], ps, st)
            @test yp ≈ y
        end
    end

    @testset "an ordinary RBM is not invariant" begin
        # A control for the test above: without weight sharing the symmetry must be absent,
        # otherwise the invariance check proves nothing.
        rng = Xoshiro(2)
        nsites = 6
        layer = RBM(nsites, 2)
        ps = Lux.initialparameters(rng, layer)
        st = Lux.initialstates(rng, layer)
        x = randn(rng, Float64, nsites, 4)
        y, _ = layer(x, ps, st)
        yp, _ = layer(x[circshift(1:nsites, 1), :], ps, st)
        @test !isapprox(yp, y)
    end

    @testset "LuxAnsatz satisfies the NQSCore interface" begin
        rng = Xoshiro(3)
        dof, nsites = Spin(1 // 2), 4
        b = basis(dof_object(dof), nsites)
        a = LuxAnsatz(RBM(nsites, 2), dof, nsites; rng=rng)
        θ = init_parameters(a, rng)

        x = configurations(dof, b.states, nsites)
        logψ = log_amplitude(a, θ, x)
        @test length(logψ) == length(b.states)
        @test eltype(logψ) <: Complex
        @test all(isfinite, logψ)

        @testset "Rational inputs are converted for the network" begin
            @test eltype(x) <: Rational        # spins really do arrive as rationals
            @test all(isfinite, log_amplitude(a, θ, x))
        end

        @testset "a two-row model output is read as re/im" begin
            model = Chain(Dense(nsites => 4, tanh), Dense(4 => 2))
            ar = LuxAnsatz(model, dof, nsites; rng=rng)
            θr = init_parameters(ar, rng)
            y = log_amplitude(ar, θr, x)
            @test length(y) == length(b.states)
            @test eltype(y) <: Complex
        end
    end

    @testset "log_derivatives match finite differences" begin
        rng = Xoshiro(4)
        dof, nsites = Spin(1 // 2), 4
        b = basis(dof_object(dof), nsites)
        a = LuxAnsatz(RBM(nsites, 1), dof, nsites; rng=rng)
        θ = init_parameters(a, rng)
        x = configurations(dof, b.states[1:6], nsites)

        O = log_derivatives(a, θ, x; backend=BACKEND, holomorphic=true)
        flat, restore = flatten_parameters(θ)
        @test size(O) == (6, length(flat))

        ε = 1e-6
        for k in (1, 4, 9)
            fp = copy(flat); fp[k] += ε
            fm = copy(flat); fm[k] -= ε
            fd = (log_amplitude(a, restore(fp), x) .- log_amplitude(a, restore(fm), x)) ./ (2ε)
            @test isapprox(O[:, k], fd; atol=1e-5)
        end
    end

    @testset "an RBM optimizes to the exact ground state" begin
        # At this size an RBM with alpha=4 has ample capacity to represent the ground state
        # exactly, so the assertion is machine precision rather than a loose tolerance: a
        # weaker bound would pass even with a subtly wrong gradient.
        dof, nsites = Spin(1 // 2), 4
        H = tfi(nsites; J=1.0, h_x=2.0)
        E_exact = minimum(real(eigvals(Hermitian(dense(H, dof, nsites)))))

        rng = Xoshiro(5)
        a = LuxAnsatz(RBM(nsites, 4), dof, nsites; rng=rng)
        vs = FullSumState(a, init_parameters(a, rng), BACKEND)

        E_initial = real(expect(vs, H).mean)
        for _ in 1:3000
            _, ∇ = expect_and_grad(vs, H)
            setparameters!(vs, fmap((p, g) -> p .- 0.05 .* g, parameters(vs), ∇))
        end
        s = expect(vs, H)

        @test real(s.mean) < E_initial                  # it actually learned something
        @test real(s.mean) >= E_exact - 1e-9            # the variational bound holds
        @test isapprox(real(s.mean), E_exact; atol=1e-8)
        @test s.variance < 1e-8                         # ...and it is an eigenstate to that order
    end

    @testset "SymmetricRBM from a lattice" begin
        # The extension: permutations derived from geometry rather than written by hand.
        lat = build(Hypercube([6]; periodic=true))
        layer = SymmetricRBM(lat, 2)
        @test length(layer.permutations) == 6           # one per translation of the ring

        rng = Xoshiro(6)
        ps = Lux.initialparameters(rng, layer)
        st = Lux.initialstates(rng, layer)
        x = randn(rng, Float64, 6, 3)
        y, _ = layer(x, ps, st)
        @test length(y) == 3

        @testset "invariant under the lattice's translations" begin
            for p in layer.permutations
                yp, _ = layer(x[p, :], ps, st)
                @test yp ≈ y
            end
        end

        @test_throws ArgumentError SymmetricRBM(lat, 2; group=:nonsense)
    end

    @testset "input validation" begin
        @test_throws ArgumentError RBM(4, 0)
        @test_throws ArgumentError SymmetricRBM(Vector{Int}[], 1)
        @test_throws ArgumentError SymmetricRBM([[1, 2, 3], [1, 2]], 1)
    end
end
