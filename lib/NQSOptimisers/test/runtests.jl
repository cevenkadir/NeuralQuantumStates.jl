using NQSOptimisers

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

function exact_ground_energy(H, dof, nsites)
    b = basis(dof_object(dof), nsites)
    index = Dict(s => i for (i, s) in pairs(b.states))
    res = connected_padded(H, b.states)
    M = zeros(ComplexF64, length(b.states), length(b.states))
    for n in eachindex(b.states), j in 1:res.counts[n]
        M[index[res.configs[j, n]], n] += res.mels[j, n]
    end
    return minimum(real(eigvals(Hermitian(M))))
end

function fresh_state(dof, nsites; seed=11, scale=0.3)
    b = basis(dof_object(dof), nsites)
    a = LogStateVector(dof, nsites, b)
    return FullSumState(a, init_parameters(a, Xoshiro(seed); scale=scale); backend=BACKEND)
end

@testset "NQSOptimisers.jl" begin
    @testset "dependency weight" begin
        loaded = Set(m.name for m in keys(Base.loaded_modules))
        for heavy in ("Lux", "Enzyme", "Reactant", "CUDA", "Metal")
            @test heavy ∉ loaded
        end
    end

    @testset "solvers" begin
        rng = Xoshiro(0)
        B = randn(rng, 12, 12)
        A = B' * B                       # symmetric positive semi-definite, as S always is
        b = randn(rng, 12)
        shift = 1e-3
        reference = (A + shift * I) \ b

        for solver in (CholeskySolver(), PseudoInverseSolver(; rtol=1e-14),
                       ConjugateGradientSolver(; tol=1e-12))
            @test solve(solver, A, b, shift) ≈ reference rtol = 1e-6
        end

        @testset "a singular matrix is handled, not silently mangled" begin
            # The geometric tensor genuinely has zero modes; a solver that returns garbage
            # there produces a wild update rather than no update.
            C = randn(rng, 12, 4)
            S = C * C'                   # rank 4 in 12 dimensions
            g = randn(rng, 12)
            for solver in (CholeskySolver(), PseudoInverseSolver(), ConjugateGradientSolver())
                x = solve(solver, S, g, 1e-6)
                @test all(isfinite, x)
            end
        end
    end

    @testset "SR and MinSR give the same update" begin
        # They are algebraically identical -- (XᵀX + λI)⁻¹Xᵀ = Xᵀ(XXᵀ + λI)⁻¹ -- so any
        # discrepancy beyond round-off is an implementation error, not an approximation.
        dof, nsites = Spin(1 // 2), 4
        H = tfi(nsites; J=1.0, h_x=0.9, h_z=0.1)

        for shift in (1e-2, 1e-4)
            vs = fresh_state(dof, nsites)
            _, δ_sr = precondition(
                StochasticReconfiguration(; diag_shift=shift, mode=:sr,
                    solver=PseudoInverseSolver(; rtol=1e-14)), vs, H)
            _, δ_min = precondition(
                StochasticReconfiguration(; diag_shift=shift, mode=:minsr,
                    solver=PseudoInverseSolver(; rtol=1e-14)), vs, H)
            @test δ_sr ≈ δ_min rtol = 1e-6
        end
    end

    @testset "the geometric tensor need not be built" begin
        X = randn(Xoshiro(21), 24, 9)
        v = randn(Xoshiro(22), 9)
        dense = transpose(X) * X

        @testset "it multiplies like the matrix it stands for" begin
            S = QuantumGeometricTensor(X)
            @test size(S) == (9, 9)
            @test size(S, 1) == 9 && size(S, 2) == 9
            @test eltype(S) === Float64
            @test to_dense(S) ≈ dense
            @test S * v ≈ dense * v
        end

        @testset "the relative shift lands on the diagonal" begin
            S = QuantumGeometricTensor(X, 0.25)
            expected = dense + 0.25 * Diagonal(diag(dense))
            @test to_dense(S) ≈ expected
            @test S * v ≈ expected * v
        end
    end

    @testset "matrix-free SR agrees with the built matrix" begin
        # :matrixfree never allocates the P x P tensor, so it must be checked against the mode
        # that does. Anything beyond solver tolerance is an implementation error.
        dof, nsites = Spin(1 // 2), 4
        H = tfi(nsites; J=1.0, h_x=0.9, h_z=0.1)

        for shift in (1e-2, 1e-4), scale in (0.0, 0.1)
            _, δ_sr = precondition(
                StochasticReconfiguration(; diag_shift=shift, diag_scale=scale, mode=:sr,
                    solver=PseudoInverseSolver(; rtol=1e-14)), fresh_state(dof, nsites), H)
            _, δ_mf = precondition(
                StochasticReconfiguration(; diag_shift=shift, diag_scale=scale,
                    mode=:matrixfree,
                    solver=ConjugateGradientSolver(; tol=1e-14, maxiter=5000)),
                fresh_state(dof, nsites), H)
            @test δ_sr ≈ δ_mf rtol = 1e-5
        end

        @testset "it reaches the ground state too" begin
            E_exact = exact_ground_energy(H, dof, nsites)
            vs = fresh_state(dof, nsites)
            optimize!(vs, H, StochasticReconfiguration(;
                    diag_shift=1e-3, mode=:matrixfree,
                    solver=ConjugateGradientSolver(; tol=1e-12));
                iterations=800, learning_rate=0.1)
            @test isapprox(real(expect(vs, H).mean), E_exact; atol=1e-5)
        end

        @testset "a direct solver has nothing to factorize" begin
            for solver in (CholeskySolver(), PseudoInverseSolver())
                @test_throws ArgumentError precondition(
                    StochasticReconfiguration(; mode=:matrixfree, solver=solver),
                    fresh_state(dof, nsites), H)
            end
        end
    end

    @testset "chunking does not change the update" begin
        dof, nsites = Spin(1 // 2), 4
        H = tfi(nsites; J=1.0, h_x=0.9, h_z=0.1)
        _, plain = precondition(
            StochasticReconfiguration(; diag_shift=1e-3), fresh_state(dof, nsites), H)
        for cs in (1, 5, 7)
            _, chunked = precondition(
                StochasticReconfiguration(; diag_shift=1e-3, chunk_size=cs),
                fresh_state(dof, nsites), H)
            @test chunked ≈ plain
        end
    end

    @testset "the update descends" begin
        # A preconditioned direction must still be a descent direction: its overlap with the
        # plain gradient has to be positive, or the "improvement" is an accident of step size.
        dof, nsites = Spin(1 // 2), 4
        H = tfi(nsites; J=1.0, h_x=0.9, h_z=0.1)
        vs = fresh_state(dof, nsites)

        _, ∇ = expect_and_grad(vs, H)
        _, δ = precondition(StochasticReconfiguration(; diag_shift=1e-3), vs, H)
        @test real(dot(vec(∇), vec(δ))) > 0
    end

    @testset "a large shift reduces to plain gradient descent" begin
        # With λ → ∞ the geometric tensor is negligible and S⁻¹∇ → ∇/λ. A limit the
        # implementation must respect, and a check that the two paths agree in scale.
        dof, nsites = Spin(1 // 2), 4
        H = tfi(nsites; J=1.0, h_x=0.9, h_z=0.1)
        vs = fresh_state(dof, nsites)

        _, ∇ = expect_and_grad(vs, H)
        _, δ = precondition(StochasticReconfiguration(; diag_shift=1e6), vs, H)
        @test vec(δ) .* 1e6 ≈ vec(∇) rtol = 1e-3
    end

    @testset "SR converges where plain gradient descent stalls" begin
        # The reason this package exists. NQSCore documents that plain descent stalls in the
        # ordered regime of the transverse-field Ising chain, settling near the classical
        # energy because the gradient carries a vanishing factor of p(s). SR divides that
        # factor out and must reach the true ground state.
        dof, nsites = Spin(1 // 2), 4
        H = tfi(nsites; J=1.0, h_x=0.9, h_z=0.1)
        E_exact = exact_ground_energy(H, dof, nsites)

        plain = fresh_state(dof, nsites)
        optimize!(plain, H, Identity(); iterations=5000, learning_rate=0.05)
        E_plain = real(expect(plain, H).mean)

        sr_state = fresh_state(dof, nsites)
        optimize!(sr_state, H, StochasticReconfiguration(; diag_shift=1e-3);
            iterations=800, learning_rate=0.1)
        E_sr = real(expect(sr_state, H).mean)

        # Plain descent really is stuck, well above the ground state...
        @test E_plain - E_exact > 0.5
        # ...and stochastic reconfiguration is not, in far fewer iterations.
        @test isapprox(E_sr, E_exact; atol=1e-6)
        @test E_sr < E_plain
        @test expect(sr_state, H).variance < 1e-8
    end

    @testset "both forms converge, and so does every solver" begin
        dof, nsites = Spin(1 // 2), 4
        H = tfi(nsites; J=1.0, h_x=0.9, h_z=0.1)
        E_exact = exact_ground_energy(H, dof, nsites)

        for mode in (:sr, :minsr), solver in (CholeskySolver(), PseudoInverseSolver())
            vs = fresh_state(dof, nsites)
            optimize!(vs, H, StochasticReconfiguration(;
                    diag_shift=1e-3, mode=mode, solver=solver);
                iterations=800, learning_rate=0.1)
            @test isapprox(real(expect(vs, H).mean), E_exact; atol=1e-5)
        end
    end

    @testset "optimize! reports a decreasing history" begin
        dof, nsites = Spin(1 // 2), 4
        H = tfi(nsites; J=1.0, h_x=0.9, h_z=0.1)
        vs = fresh_state(dof, nsites)
        history = optimize!(vs, H, StochasticReconfiguration(; diag_shift=1e-3);
            iterations=200, learning_rate=0.1)

        @test length(history) == 200
        @test all(h -> h isa Stats, history)
        @test real(history[end].mean) < real(history[1].mean)
    end

    @testset "input validation" begin
        @test_throws ArgumentError StochasticReconfiguration(; mode=:nonsense)
        @test StochasticReconfiguration(; mode=:matrixfree).mode === :matrixfree
        @test_throws ArgumentError StochasticReconfiguration(; diag_shift=-1.0)
    end
end
