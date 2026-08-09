using NQSAnsatze

using ConnectedBasisConfigurations
using DifferentiationInterface
using Functors: fmap
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

        @testset "a batch bound for a derivative is asked about separately" begin
            # Two different types, and the difference is the whole point. A forward pass wants
            # the narrowest faithful one, because a wider batch buys it nothing; a batch about
            # to be differentiated wants the network's own arithmetic type, because otherwise
            # every matrix product in the reverse pass is a mixed complex-real pair that no
            # BLAS has a kernel for. Promoting both was measured: 4.4x on the layer's reverse
            # pass, 13% worse on `expect` and 78% worse on a Metropolis sweep.
            @test eltype(θ.weight) <: Complex
            @test NQSAnsatze._input_type(θ) === real(eltype(θ.weight))
            @test NQSCore.input_type(a, θ) === eltype(θ.weight)

            # A real network asks for a real batch: this is about matching the parameters, not
            # about complex arithmetic being preferable.
            ar = LuxAnsatz(Chain(Dense(nsites => 2, tanh), Dense(2 => 2)), dof, nsites; rng=rng)
            θr = init_parameters(ar, rng)
            @test NQSCore.input_type(ar, θr) <: Real

            # Whatever the type, the answer is the same. A wider batch must pass through
            # untouched rather than being demoted back — which would throw, or silently drop
            # the imaginary part.
            @test log_amplitude(a, θ, ComplexF64.(x)) ≈ logψ
            @test log_amplitude(a, θ, Float64.(x)) ≈ logψ
            @test NQSAnsatze._as_input(Float64, ComplexF64.(x)) == ComplexF64.(x)
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

    @testset "colocation leaves host arrays alone" begin
        # The batch follows the parameters onto whatever device they are on. Everything here is
        # in host memory, so every one of these must be the identity -- returning a copy instead
        # would be a mutation, and reverse-mode AD refuses to differentiate through one. The
        # ComponentArray case is the one that matters: rebuilding parameters from a flat vector
        # hands back views, not `Array`s, and that is exactly what differentiating does.
        x = randn(Xoshiro(0), 3, 4)
        @test NQSAnsatze.colocate(nothing, x) === x
        @test NQSAnsatze.colocate(randn(2, 2), x) === x
        @test NQSAnsatze.colocate(view(randn(8), 2:5), x) === x
        @test NQSAnsatze.colocate(reshape(randn(8), 2, 4), x) === x

        θ = (W=randn(Xoshiro(1), 2, 3), b=randn(Xoshiro(2), 2))
        flat, restore = flatten_parameters(θ)
        rebuilt = restore(flat)
        @test !(rebuilt.W isa Array)                 # a view into a ComponentArray
        @test NQSAnsatze.colocate(rebuilt.W, x) === x
    end

    @testset "every layer differentiates correctly" begin
        # Jastrow assembles its coupling matrix by gathering from the parameter vector, and
        # SymmetricRBM gathers permuted copies of the input. Both are indexing operations that
        # reverse-mode AD has to carry a gradient back through, so each layer is checked against
        # a finite difference rather than only the one that happens to be simplest.
        rng = Xoshiro(12)
        dof, nsites = Spin(1 // 2), 4
        b = basis(dof_object(dof), nsites)
        perms = [circshift(1:nsites, s) for s in 0:(nsites-1)]

        @testset "$(nameof(typeof(model)))" for model in (
            RBM(nsites, 2), Jastrow(nsites), SymmetricRBM(perms, 2)
        )
            a = LuxAnsatz(model, dof, nsites; rng=rng)
            θ = init_parameters(a, Xoshiro(5))
            x = configurations(dof, b.states[1:6], nsites)

            O = log_derivatives(a, θ, x; backend=BACKEND, holomorphic=true)
            flat, restore = flatten_parameters(θ)
            @test size(O) == (6, length(flat))
            @test all(isfinite, O)

            ε = 1e-6
            for k in (1, min(3, length(flat)), length(flat))
                fp = copy(flat); fp[k] += ε
                fm = copy(flat); fm[k] -= ε
                fd = (log_amplitude(a, restore(fp), x) .- log_amplitude(a, restore(fm), x)) ./ (2ε)
                @test isapprox(O[:, k], fd; atol=1e-5)
            end
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
        vs = FullSumState(a, init_parameters(a, rng); backend=BACKEND)

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
        @test size(layer.permutations, 2) == 6       # one per translation of the ring

        rng = Xoshiro(6)
        ps = Lux.initialparameters(rng, layer)
        st = Lux.initialstates(rng, layer)
        x = randn(rng, Float64, 6, 3)
        y, _ = layer(x, ps, st)
        @test length(y) == 3

        @testset "invariant under the lattice's translations" begin
            for p in eachcol(layer.permutations)
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
