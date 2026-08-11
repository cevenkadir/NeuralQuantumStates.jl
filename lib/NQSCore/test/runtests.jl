using NQSCore

# Recorded before anything else is loaded, so that the dependency-weight testset below sees
# what `using NQSCore` alone pulled in rather than what the rest of this file needs.
const LOADED_AFTER_NQSCORE = Set(m.name for m in keys(Base.loaded_modules))

using Aqua
using ConnectedBasisConfigurations
using DifferentiationInterface
using ForwardDiff
using KernelAbstractions
using LinearAlgebra
using OperatorAlgebra
using Random
using Statistics: mean
using SymBasis
using Test
using TOML

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
function dense(H, spec, nsites)
    b = basis(dof_object(spec), nsites)
    index = Dict(s => i for (i, s) in pairs(b.states))
    res = connected_padded(H, b.states)
    M = zeros(ComplexF64, length(b.states), length(b.states))
    for n in eachindex(b.states), j in 1:res.counts[n]
        M[index[res.configs[j, n]], n] += res.mels[j, n]
    end
    return M
end

"""A `FullSumState` on the exact ansatz, with the given parameters."""
function full_sum(spec, nsites, θ=nothing; rng=Xoshiro(0))
    b = basis(dof_object(spec), nsites)
    a = LogStateVector(spec, nsites, b)
    θ === nothing && (θ = init_parameters(a, rng; scale=0.3))
    return FullSumState(a, θ; backend=BACKEND), a, b
end

"""
An ansatz with **real** parameters held in a `NamedTuple`, and a complex log-amplitude.

That combination is what a Lux model produces and what `LogStateVector` never exercises: the
`ComponentArrays` flattening path, and the real-parameter branch of the gradient and of
[`log_derivatives`](@ref). Deliberately tiny, and deliberately not a good wavefunction.
"""
struct ToyAnsatz{D} <: NQSCore.AbstractAnsatz
    dof::D
    nsites::Int
    nhidden::Int
end

function NQSCore.log_amplitude(a::ToyAnsatz, θ, x::AbstractMatrix)
    h = θ.W * x .+ θ.b
    logmod = vec(sum(log.(2 .* cosh.(h)); dims=1))
    phase = vec(transpose(x) * θ.v)
    return complex.(logmod, phase)
end

toy_parameters(a::ToyAnsatz, rng) = (
    W=0.2 .* randn(rng, a.nhidden, a.nsites),
    b=0.2 .* randn(rng, a.nhidden),
    v=0.2 .* randn(rng, a.nsites),
)

"""The energy gradient built the long way, from an explicit log-derivative matrix."""
function gradient_from_jacobian(a, θ, x, E, weights)
    O = log_derivatives(a, θ, x; backend=BACKEND)
    p = weights === nothing ? fill(1 / length(E), length(E)) : weights ./ sum(weights)
    Ē = sum(p .* E)
    Ō = vec(sum(p .* O; dims=1))
    ŌE = vec(sum(p .* conj.(O) .* E; dims=1))
    return match_parameter_shape(2 .* real.(ŌE .- conj.(Ō) .* Ē), θ)
end

@testset "NQSCore.jl" begin
    @testset "dependency weight" begin
        # The direct dependency surface is part of the package's contract: this is an interface
        # package, and everything that depends on it inherits whatever it pulls in. Asserting
        # the set rather than counting it means an accidental `Pkg.add` fails here.
        project = TOML.parsefile(joinpath(pkgdir(NQSCore), "Project.toml"))
        @test Set(keys(project["deps"])) == Set([
            "ComponentArrays", "ConnectedBasisConfigurations", "DifferentiationInterface",
            "Random", "Statistics",
        ])

        # Weak dependencies are named too, but the assertion here is weaker on purpose: a
        # weakdep costs a caller nothing until they load it, so the list may grow. What must
        # not happen is one migrating into `[deps]`, and the exact set above is what stops that
        # — `Pkg.develop` has silently made that move twice in this repository's history.
        @test Set(keys(get(project, "weakdeps", Dict()))) ==
              Set(["Enzyme", "KernelAbstractions", "Reactant"])
        @test Set(keys(get(project, "extensions", Dict()))) ==
              Set(["NQSCoreKernelAbstractionsExt", "NQSCoreReactantExt"])

        # Every dependency is pinned, which registration requires and which nothing else checks.
        compat = project["compat"]
        for name in keys(project["deps"])
            @test haskey(compat, name)
        end
        for name in keys(get(project, "weakdeps", Dict()))
            @test haskey(compat, name)
        end
        @test haskey(compat, "julia")

        # Neural-network and GPU stacks are the caller's business, not this package's.
        for heavy in ("Lux", "Enzyme", "Reactant", "CUDA", "Metal", "Zygote", "ForwardDiff")
            @test heavy ∉ LOADED_AFTER_NQSCORE
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

        @testset "a short or constant chain is reported as uncorrelated" begin
            @test integrated_autocorrelation([1.0, 2.0, 3.0]) == 1.0
            @test integrated_autocorrelation(fill(2.5, 100)) == 1.0
        end

        @testset "split-Rhat flags a drifting chain" begin
            converged = randn(Xoshiro(1), 2000, 4)
            @test isapprox(split_rhat(converged), 1.0; atol=0.05)

            drifting = converged .+ range(0, 8; length=2000)   # a steady trend in every chain
            @test split_rhat(drifting) > 1.1

            @test split_rhat(reshape([1.0, 2.0, 3.0], :, 1)) == 1.0   # too short to judge
            @test split_rhat(fill(1.0, 100, 2)) == 1.0                # no within-chain variance
        end

        @testset "the error of the mean shrinks like 1/sqrt(n)" begin
            rng = Xoshiro(3)
            e1 = statistics(randn(rng, 1000)).error_of_mean
            e2 = statistics(randn(rng, 4000)).error_of_mean
            @test isapprox(e2 / e1, 0.5; rtol=0.5)
        end

        @testset "several chains give a between-chain error bar" begin
            # With more than one chain the error comes from the spread of the chain means, a
            # different code path from the single-chain autocorrelation estimate. Both must land
            # on the same answer for independent samples, since both estimate the same thing.
            rng = Xoshiro(9)
            many = randn(rng, 2000, 8)
            s = statistics(many)
            @test s.mean ≈ mean(vec(many))
            @test isapprox(s.error_of_mean, sqrt(1 / length(many)); rtol=0.4)
            @test isapprox(s.r_hat, 1.0; atol=0.05)

            # Chains that disagree with one another inflate the error bar and R̂ together.
            offset = copy(many)
            offset[:, 1] .+= 5
            @test statistics(offset).error_of_mean > 5 * s.error_of_mean
            @test statistics(offset).r_hat > 1.1
        end

        @testset "complex samples keep a complex mean and a real variance" begin
            v = ComplexF64[1+2im, 3-1im, 2+0im, 0+1im]
            s = statistics(v)
            @test s.mean ≈ mean(v)
            @test s.variance isa Float64
        end

        @testset "weighted statistics reproduce a known average" begin
            s = weighted_statistics([1.0, 2.0, 3.0], [0.5, 0.25, 0.25])
            @test s.mean ≈ 1.75
            @test s.variance ≈ 0.5 * 0.75^2 + 0.25 * 0.25^2 + 0.25 * 1.25^2
            @test s.error_of_mean == 0.0
        end

        @testset "data is fetched, not probed" begin
            # `statistics` walks lags sequentially, so it brings its input to the host once.
            # On the host that must cost nothing and change nothing — a copy here would be a
            # full sample array per call, and `to_host` returning something other than the
            # original would break the samplers that write back into what they were given.
            v = randn(Xoshiro(0), 64)
            m = reshape(v, :, 1)
            @test NQSCore.to_host(v) === v
            @test NQSCore.to_host(m) === m
            @test statistics(m).mean ≈ statistics(v).mean

            # A non-`Array` input is copied rather than indexed in place.
            sub = view(randn(Xoshiro(1), 128), 1:64)
            @test NQSCore.to_host(sub) isa Array
            @test NQSCore.to_host(sub) == collect(sub)
            @test statistics(reshape(collect(sub), :, 1)).mean ≈ mean(sub)
        end

        @testset "show and isapprox" begin
            s = Stats(1.25, 0.01, 0.5, 1.2, 1.001)
            str = sprint(show, s)
            @test occursin("1.25", str) && occursin("±", str) && occursin("R̂", str)
            @test occursin("im", sprint(show, Stats(1.0 + 2.0im, 0.0, 0.0, 1.0, 1.0)))

            @test s ≈ Stats(1.25 + 1e-12, 0.9, 0.9, 9.0, 9.0)   # compares means only
            @test s ≈ 1.25
            @test 1.25 ≈ s
        end
    end

    @testset "parameter flattening" begin
        @testset "a plain vector is its own flat form" begin
            θ = ComplexF64[1, 2, 3]
            flat, restore = flatten_parameters(θ)
            @test flat === θ
            @test restore(flat) === θ
        end

        @testset "a NamedTuple round-trips as a NamedTuple" begin
            # Not as a ComponentArray: the promise downstream is that a gradient has the same
            # structure as the parameters it belongs to, and `fmap` walks the two together.
            θ = (W=[1.0 2.0; 3.0 4.0], b=[5.0, 6.0])
            flat, restore = flatten_parameters(θ)
            @test flat isa AbstractVector{Float64}
            @test length(flat) == 6
            back = restore(flat)
            @test back isa NamedTuple
            @test back.W == θ.W && back.b == θ.b

            nested = (enc=(W=[1.0, 2.0],), dec=(W=[3.0], b=[4.0]))
            flat2, restore2 = flatten_parameters(nested)
            @test length(flat2) == 4
            @test restore2(flat2).enc.W == nested.enc.W
            @test restore2(flat2).dec.b == nested.dec.b
        end

        @testset "match_parameter_shape recombines a split complex gradient" begin
            θ = ComplexF64[1, 2, 3]
            ∇ = [1.0, 2.0, 3.0, 10.0, 20.0, 30.0]      # [∂/∂θ_re; ∂/∂θ_im]
            @test match_parameter_shape(∇, θ) == ComplexF64[1+10im, 2+20im, 3+30im]

            # A gradient that already matches the parameters is passed through untouched.
            @test match_parameter_shape([1.0, 2.0, 3.0], [1.0, 2.0, 3.0]) == [1.0, 2.0, 3.0]

            g = match_parameter_shape([1.0, 2.0, 3.0], (a=[1.0, 2.0], b=[3.0]))
            @test g isa NamedTuple && g.a == [1.0, 2.0] && g.b == [3.0]
        end
    end

    @testset "log_derivatives" begin
        spec, nsites = Spin(1 // 2),4
        vs, a, b = full_sum(spec, nsites)
        x = configurations(spec,b.states, nsites)
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

        @testset "real parameters give one column each" begin
            rng = Xoshiro(2)
            toy = ToyAnsatz(spec, nsites, 3)
            θ = toy_parameters(toy, rng)
            O = log_derivatives(toy, θ, x; backend=BACKEND)
            @test size(O) == (length(b.states), length(first(flatten_parameters(θ))))
            @test eltype(O) <: Complex          # complex ψ even though θ is real

            # Column k against a finite difference of log ψ in parameter k.
            flat, restore = flatten_parameters(θ)
            ε = 1e-6
            for k in (1, 4, length(flat))
                up = copy(flat); up[k] += ε
                dn = copy(flat); dn[k] -= ε
                fd = (log_amplitude(toy, restore(up), x) .- log_amplitude(toy, restore(dn), x)) ./ (2ε)
                @test isapprox(O[:, k], fd; atol=1e-6)
            end
        end

        @testset "centering removes the column mean" begin
            O = log_derivatives(a, parameters(vs), x; backend=BACKEND, holomorphic=true)
            @test all(abs.(sum(centered(O); dims=1)) .< 1e-10)

            w = probabilities(vs)
            @test all(abs.(sum(w .* centered(O, w); dims=1)) .< 1e-10)
        end
    end

    @testset "uniform weights live where the estimators do" begin
        # `_gradient_cotangent` built them with `fill`, which is a host `Vector` whatever `E` is.
        # On a device that is not a slow path, it is a compilation failure — "passing
        # non-bitstype argument", from the `Extruded` wrapper a host array arrives in — and it
        # is the `weights === nothing` branch, which is every `MCState`. No GPU is needed to see
        # it: the array type is the whole of the bug.
        E = ComplexF64[1.0, 2.0, 3.0, 4.0]
        c = NQSCore._gradient_cotangent(E, nothing)
        @test c isa Vector{ComplexF64}
        @test sum(c) ≈ 0 atol = 1e-12          # centred, which is the point of subtracting Ē

        # A weights vector is followed rather than replaced.
        w = [1.0, 1.0, 2.0, 4.0]
        @test NQSCore._gradient_cotangent(E, w) ≈ (w ./ sum(w)) .* (E .- sum((w ./ sum(w)) .* E))

        # `similar(E, ...)` is what carries the array type across, and that is the whole of
        # the fix: on a device `E` is a device array and the weights become one too.
        @test NQSCore._gradient_cotangent(view(E, :), nothing) isa AbstractVector
    end

    @testset "a compiler is a separate choice from a backend" begin
        # `backend` says which engine differentiates; `compiler` says whether anything compiles.
        # They were one field once, and that made `expect` and stochastic reconfiguration
        # disagree about what it meant: `log_derivatives` has no compiled form, so it would have
        # handed a compiler to DifferentiationInterface.
        spec, nsites = Spin(1 // 2), 4
        H = tfi(nsites; J=1.0, h_x=0.7, h_z=0.2)
        b = basis(dof_object(spec), nsites)
        a = LogStateVector(spec, nsites, b)
        θ = init_parameters(a, Xoshiro(0); scale=0.3)

        plain = FullSumState(a, θ; backend=BACKEND)
        asked = FullSumState(a, θ; backend=BACKEND, compiler=NQSCore.AutoReactant())
        @test plain.compiler === nothing
        @test asked.backend === BACKEND

        # Nothing provides a compiled step here, so asking for one is an error rather than a
        # silent hundredfold slowdown.
        @test_throws ArgumentError expect(asked, H)
        @test_throws ArgumentError expect_and_grad(asked, H)

        # Everything with no compiled form goes on using `backend`, untouched.
        @test last(expect_and_grad(asked, H; chunk_size=4)) ==
              last(expect_and_grad(plain, H; chunk_size=4))
        @test local_estimators(asked, H; holomorphic=true).O ==
              local_estimators(plain, H; holomorphic=true).O
    end

    @testset "precision follows the parameters" begin
        # Choosing single precision is choosing it for the arithmetic, and the whole device seam
        # is built on the rule that what a batch is and where it lives follow the parameters. The
        # matrix elements were the one array that did not: `local_operators` builds them
        # `Float64`, and a `Float64` array meeting a `ComplexF32` one promoted the reduction, so
        # `local_energy`, the cotangent and `expect` all came back double for a single-precision
        # ansatz — silently overriding the choice.
        spec, nsites = Spin(1 // 2), 4
        H = tfi(nsites; J=1.0, h_x=0.7, h_z=0.2)
        b = basis(dof_object(spec), nsites)

        @testset "$T" for T in (ComplexF64, ComplexF32)
            a = LogStateVector(spec, nsites, b)
            θ = T.(init_parameters(a, Xoshiro(0); scale=0.3))
            vs = FullSumState(a, θ; backend=BACKEND)

            logψ = log_amplitude(a, θ, configurations(spec, b.states, nsites))
            E = local_energy(vs, H, b.states)
            p = NQSCore.born_probabilities(logψ)

            @test eltype(E) === T
            @test eltype(NQSCore._gradient_cotangent(E, p)) === T
            @test typeof(expect(vs, H).mean) === T
        end

        @testset "the two precisions agree to single precision" begin
            # The point of the narrowing is speed, not a different answer.
            a = LogStateVector(spec, nsites, b)
            θ = init_parameters(a, Xoshiro(0); scale=0.3)
            e64 = expect(FullSumState(a, θ; backend=BACKEND), H).mean
            e32 = expect(FullSumState(a, ComplexF32.(θ); backend=BACKEND), H).mean
            @test real(e32) ≈ real(e64) rtol = 1e-5
        end
    end

    @testset "FullSumState is exact" begin
        spec, nsites = Spin(1 // 2),4
        H = tfi(nsites; J=1.0, h_x=0.7, h_z=0.2)
        M = dense(H, spec, nsites)

        @testset "expectation matches the Rayleigh quotient" begin
            vs, a, b = full_sum(spec, nsites)
            ψ = exp.(log_amplitude(a, parameters(vs), configurations(spec,b.states, nsites)))
            expected = real(dot(ψ, M * ψ) / dot(ψ, ψ))

            s = expect(vs, H)
            @test real(s.mean) ≈ expected
            @test s.error_of_mean == 0.0          # nothing was sampled
        end

        @testset "an eigenstate has zero variance" begin
            # The sharpest available check: E_loc is constant exactly on an eigenstate.
            vals, vecs = eigen(Hermitian(M))
            b = basis(dof_object(spec), nsites)
            a = LogStateVector(spec, nsites, b)
            vs = FullSumState(a, log.(complex.(vecs[:, 1])); backend=BACKEND)

            s = expect(vs, H)
            @test real(s.mean) ≈ vals[1]
            @test s.variance < 1e-16
        end

        @testset "the gradient vanishes at the ground state" begin
            vals, vecs = eigen(Hermitian(M))
            b = basis(dof_object(spec), nsites)
            a = LogStateVector(spec, nsites, b)
            vs = FullSumState(a, log.(complex.(vecs[:, 1])); backend=BACKEND)
            _, ∇ = expect_and_grad(vs, H)
            @test maximum(abs, ∇) < 1e-6
        end

        @testset "the gradient matches a finite difference" begin
            # The parameters are complex, so the gradient has a real and an imaginary part and
            # each must be checked against its own perturbation direction: `real(∇)` against a
            # perturbation of `θ_re`, `imag(∇)` against one of `θ_im`. Comparing the whole
            # complex gradient to a single real-direction difference conflates the two.
            vs, a, b = full_sum(spec, nsites)
            _, ∇ = expect_and_grad(vs, H)
            θ = copy(parameters(vs))
            E(t) = real(expect(FullSumState(a, t; backend=BACKEND), H).mean)
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

        @testset "the gradient agrees with the explicit Jacobian formula" begin
            # `expect_and_grad` contracts the Jacobian without building it. That shortcut has to
            # reproduce, exactly, what building the matrix and contracting it would have given —
            # for complex parameters and for a NamedTuple of real ones.
            vs, a, b = full_sum(spec, nsites)
            x = configurations(spec,b.states, nsites)
            _, ∇ = expect_and_grad(vs, H)
            E = local_energy(vs, H, b.states)
            @test ∇ ≈ gradient_from_jacobian(a, parameters(vs), x, E, probabilities(vs))

            toy = ToyAnsatz(spec, nsites, 3)
            θ = toy_parameters(toy, Xoshiro(4))
            tvs = FullSumState(toy, θ; backend=BACKEND, basis=b)
            _, ∇t = expect_and_grad(tvs, H)
            Et = local_energy(tvs, H, b.states)
            @test ∇t isa NamedTuple
            gold = gradient_from_jacobian(toy, θ, x, Et, probabilities(tvs))
            @test ∇t.W ≈ gold.W
            @test ∇t.b ≈ gold.b
            @test ∇t.v ≈ gold.v
        end

        @testset "a compiled operator gives the same answer" begin
            vs, _, _ = full_sum(spec, nsites)
            @test expect(vs, compile(H)).mean ≈ expect(vs, H).mean
            E1, ∇1 = expect_and_grad(vs, compile(H))
            E2, ∇2 = expect_and_grad(vs, H)
            @test E1.mean ≈ E2.mean
            @test ∇1 ≈ ∇2
        end

        @testset "an explicit basis restricts the sum to a sector" begin
            dofo = dof_object(spec)
            sector = basis(dofo, nsites, sym(TotalMagnetization(0 // 1, nsites), dofo))
            a = LogStateVector(spec, nsites, sector)
            vs = FullSumState(a, init_parameters(a, Xoshiro(6)); basis=sector, backend=BACKEND)
            @test length(samples(vs)) == length(sector.states)
            @test sum(probabilities(vs)) ≈ 1
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
        spec, nsites = Spin(1 // 2),4
        for (J, h_x) in ((1.0, 2.0), (0.5, 2.0), (1.0, 3.0))
            H = compile(tfi(nsites; J=J, h_x=h_x))
            E_exact = minimum(real(eigvals(Hermitian(dense(H, spec, nsites)))))

            vs, _, _ = full_sum(spec, nsites; rng=Xoshiro(11))
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
        spec, nsites = Spin(1 // 2),4
        H = tfi(nsites; J=1.0, h_x=0.7, h_z=0.2)
        b = basis(dof_object(spec), nsites)
        a = LogStateVector(spec, nsites, b)
        θ = init_parameters(a, Xoshiro(5); scale=0.3)

        exact = expect(FullSumState(a, θ; backend=BACKEND), H)
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

        @testset "the sampled gradient agrees with the exact one" begin
            _, ∇exact = expect_and_grad(FullSumState(a, θ; backend=BACKEND), H)
            _, ∇mc = expect_and_grad(mc, H)
            @test maximum(abs, ∇mc .- ∇exact) < 0.05
        end

        @testset "changing parameters invalidates the samples" begin
            n_before = length(samples(mc))
            setparameters!(mc, θ .+ 0.1)
            @test mc.stale
            @test length(samples(mc)) == n_before
            @test !mc.stale
        end

        @testset "resample! draws again" begin
            small = MCState(a, θ, ExactSampler(b, 64); backend=BACKEND, rng=Xoshiro(12))
            before = copy(samples(small))
            resample!(small)
            @test samples(small) != before      # a fresh draw, not the cached one
            @test length(samples(small)) == 64
        end
    end

    @testset "sampler contract" begin
        spec, nsites = Spin(1 // 2),3
        b = basis(dof_object(spec), nsites)
        a = LogStateVector(spec, nsites, b)
        θ = init_parameters(a, Xoshiro(8); scale=0.3)
        s = ExactSampler(b, 128)

        @testset "sample returns its own state alongside the draw" begin
            drawn, state = sample(s, a, θ, Xoshiro(0))
            @test length(drawn) == 128
            @test all(∈(b.states), drawn)
            @test state === nothing              # independent draws carry nothing forward
        end

        @testset "a handed-back state is accepted" begin
            drawn, state = sample(s, a, θ, Xoshiro(0))
            again, _ = sample(s, a, θ, Xoshiro(0), state)
            @test again == drawn                 # ExactSampler ignores it, deterministically
        end

        @testset "the state survives a parameter change" begin
            mc = MCState(a, θ, s; backend=BACKEND, rng=Xoshiro(0))
            @test sampler_state(mc) === nothing
            setparameters!(mc, θ .+ 0.1)
            @test sampler_state(mc) === nothing   # kept, not cleared
        end

        @testset "sampling follows the Born distribution" begin
            drawn, _ = sample(ExactSampler(b, 200_000), a, θ, Xoshiro(1))
            counts = Dict{eltype(b.states),Int}()
            for st in drawn
                counts[st] = get(counts, st, 0) + 1
            end
            vs = FullSumState(a, θ; backend=BACKEND)
            p = probabilities(vs)
            empirical = [get(counts, st, 0) / length(drawn) for st in b.states]
            @test maximum(abs, empirical .- p) < 0.01
        end
    end

    @testset "local estimators" begin
        spec, nsites = Spin(1 // 2),4
        H = tfi(nsites; J=1.0, h_x=0.7, h_z=0.2)
        vs, a, b = full_sum(spec, nsites)
        x = configurations(spec,b.states, nsites)

        @testset "the pieces agree with computing them separately" begin
            res = local_estimators(vs, H; holomorphic=true)
            @test res.E ≈ local_energy(vs, H, b.states)
            @test res.O ≈ log_derivatives(a, parameters(vs), x; backend=BACKEND, holomorphic=true)
            @test res.weights ≈ probabilities(vs)
            @test sum(res.weights) ≈ 1
        end

        @testset "non-holomorphic mode widens O" begin
            res = local_estimators(vs, H)
            @test size(res.O, 2) == 2 * n_parameters(a)
        end

        @testset "Monte Carlo samples carry no weights" begin
            mc = MCState(a, parameters(vs), ExactSampler(b, 256); backend=BACKEND, rng=Xoshiro(2))
            res = local_estimators(mc, H)
            @test res.weights === nothing
            @test sample_weights(mc) === nothing
            @test length(res.E) == 256
            @test size(res.O, 1) == 256
        end

        @testset "the gradient can be rebuilt from them" begin
            res = local_estimators(vs, H)
            _, ∇ = expect_and_grad(vs, H)
            p = res.weights ./ sum(res.weights)
            Ē = sum(p .* res.E)
            Ō = vec(sum(p .* res.O; dims=1))
            ŌE = vec(sum(p .* conj.(res.O) .* res.E; dims=1))
            rebuilt = match_parameter_shape(2 .* real.(ŌE .- conj.(Ō) .* Ē), parameters(vs))
            @test rebuilt ≈ ∇
        end
    end

    @testset "local energies" begin
        spec, nsites = Spin(1 // 2),4
        H = tfi(nsites; J=1.0, h_x=0.7, h_z=0.2)
        vs, a, b = full_sum(spec, nsites)

        @testset "the chain layout of the samples survives" begin
            # A (steps, chains) matrix must come back as one, or `statistics` loses the chain
            # structure it needs for R̂ and for a between-chain error bar.
            states = reshape(repeat(b.states, 2), :, 2)
            E = local_energy(vs, H, states)
            @test size(E) == size(states)
            @test vec(E[:, 1]) ≈ local_energy(vs, H, b.states)
        end

        @testset "defaults to the state's own samples" begin
            @test local_energy(vs, H) ≈ local_energy(vs, H, samples(vs))
        end

        @testset "padded slots contribute nothing" begin
            # The reduction runs over the full column rather than each sample's connection
            # count, which is only correct because a padded slot repeats the sample with a zero
            # matrix element. If that ever changed, this is what would catch it.
            res = ConnectedBasisConfigurations.connected_padded(H, b.states)
            for k in eachindex(b.states)
                for j in (res.counts[k]+1):size(res.mels, 1)
                    @test iszero(res.mels[j, k])
                    @test res.configs[j, k] == b.states[k]
                end
            end
        end

        @testset "a real-valued ansatz still gets complex local energies" begin
            # The accumulator has to admit the operator's matrix elements, not just the ansatz's
            # output type: a real ansatz with a complex operator is perfectly legitimate.
            Hc = OpSum([im * Op(2 .* local_operators(Spin(1 // 2)).sy, 1)])
            toy = ToyAnsatz(spec, nsites, 2)
            θ = toy_parameters(toy, Xoshiro(3))
            tvs = FullSumState(toy, θ; backend=BACKEND, basis=b)
            E = local_energy(tvs, Hc, b.states)
            @test eltype(E) <: Complex
            @test any(!iszero, E)
        end
    end

    @testset "the device local-energy path" begin
        # KernelAbstractions runs the same kernels on a CPU backend as on a GPU one, so the
        # whole device path can be checked here — against the host path, on a real model,
        # without any GPU present. What a device adds beyond this is the compiler and the
        # memory, not the algorithm.
        spec, nsites = Spin(1 // 2), 4
        H = tfi(nsites; J=1.0, h_x=0.7, h_z=0.2)
        vs, a, b = full_sum(spec, nsites)
        logψ = NQSCore.log_amplitudes(vs, b.states)

        @test Base.get_extension(NQSCore, :NQSCoreKernelAbstractionsExt) !== nothing

        @testset "loading KernelAbstractions does not divert a host run" begin
            # `get_backend` would happily call an `Array` a `CPU()` backend. Acting on that
            # would mean an unrelated `using` swapped the tested serial kernel for a
            # launch-per-batch one, on a machine with nothing to launch onto.
            @test NQSCore.device_backend(zeros(3)) === nothing
            @test NQSCore.device_backend(zeros(ComplexF64, 3)) === nothing
            @test NQSCore.device_backend(logψ) === nothing
        end

        host = NQSCore.connections(vs, H, b.states, logψ)
        device = NQSCore.device_connections(vs, H, b.states, logψ, CPU())

        @testset "the configurations are the same, as floats" begin
            # The host unpacking returns `Rational{Int64}` for a spin, which cannot live on a
            # GPU at all; the device kernel writes the float the network wants directly, which
            # is the second full-size array this path exists to avoid.
            @test eltype(host[1]) <: Rational
            @test eltype(device[1]) === Float64

            n = length(b.states)
            h_host, h_dev = size(host[2], 1), size(device[2], 1)
            @test h_dev >= h_host
            @test size(device[1]) == (nsites, h_dev * n)

            X_host = reshape(Float64.(host[1]), nsites, h_host, n)
            X_dev = reshape(device[1], nsites, h_dev, n)
            @test X_dev[:, 1:h_host, :] == X_host
            @test device[2][1:h_host, :] == host[2]

            # The device path keeps the operator's static bound rather than trimming to the
            # largest count it saw, so it can have extra rows. They must be inert padding —
            # the sample itself with a zero matrix element — or the reduction is wrong.
            samples_x = Float64.(configurations(spec, b.states, nsites))
            for j in (h_host+1):h_dev
                @test all(iszero, device[2][j, :])
                @test X_dev[:, j, :] == samples_x
            end
        end

        @testset "the local energies are identical" begin
            # A `view` is not an `Array`, so it takes the device branch, and KernelAbstractions
            # resolves its backend through the parent — which is how the whole path, dispatch
            # included, is reachable on a machine with no GPU.
            @test NQSCore.device_backend(view(logψ, :)) == CPU()
            E_device = NQSCore._local_energy(vs, H, b.states, view(logψ, :))
            @test E_device ≈ local_energy(vs, H, b.states)
        end

        @testset "the samples are unpacked where the parameters are" begin
            # The batch is an input to the ansatz, so it belongs wherever the ansatz runs, and
            # the parameters are the only thing in a variational state that a caller
            # deliberately placed somewhere.
            @test NQSCore._reference_array(parameters(vs)) === parameters(vs)
            @test NQSCore._reference_array((W=zeros(2, 2), b=zeros(2))) == zeros(2, 2)
            @test NQSCore._reference_array(NamedTuple()) === nothing
            @test NQSCore._reference_array(1.0) === nothing

            # Host parameters keep the host unpacking, exact rationals and all.
            @test NQSCore.configurations_of(vs, b.states) ==
                  configurations(spec, b.states, nsites)

            device_x = NQSCore.device_configurations(vs, b.states, nothing, parameters(vs), CPU())
            @test eltype(device_x) === Float64
            @test device_x == Float64.(configurations(spec, b.states, nsites))

            # An ansatz that names a type gets it written straight out, rather than a real
            # batch it then has to widen — which for a batch bound for a derivative is the
            # difference between a BLAS product and a generic one in every reverse pass.
            wide = NQSCore.device_configurations(vs, b.states, ComplexF64, parameters(vs), CPU())
            @test eltype(wide) === ComplexF64
            @test wide == ComplexF64.(configurations(spec, b.states, nsites))

            # `LogStateVector` has no opinion, and must not be given one: its parameters are
            # amplitudes rather than network weights, and it looks its configurations up in a
            # basis rather than doing arithmetic on them.
            @test NQSCore.input_type(a, parameters(vs)) === nothing
        end

        @testset "an operator already on the backend is not moved again" begin
            # The loop-friendly form: upload once, reuse. It has to give the same answer as
            # handing over the unflattened operator, or the optimization is paying for a
            # transfer per step to avoid one.
            resident = to_backend(flatten(H), CPU())
            x, mels = NQSCore.device_connections(vs, resident, b.states, logψ, CPU())
            @test x == device[1]
            @test mels == device[2]
        end
    end

    @testset "chunking" begin
        # Chunking bounds the memory of a differentiation pass. It must not change the answer,
        # for any block size — including ones that do not divide the batch, which is where an
        # off-by-one in the remainder would show up.
        spec, nsites = Spin(1 // 2), 4
        H = tfi(nsites; J=1.0, h_x=0.7, h_z=0.2)
        vs, a, b = full_sum(spec, nsites)
        x = configurations(spec, b.states, nsites)
        n = length(b.states)

        O = log_derivatives(a, parameters(vs), x; backend=BACKEND, holomorphic=true)
        _, ∇ = expect_and_grad(vs, H)

        @testset "block size $cs" for cs in (1, 3, 5, n - 1, n, n + 7)
            @test log_derivatives(
                a, parameters(vs), x; backend=BACKEND, holomorphic=true, chunk_size=cs
            ) == O
            _, ∇c = expect_and_grad(vs, H; chunk_size=cs)
            @test ∇c ≈ ∇
        end

        @testset "non-holomorphic chunking keeps the doubled width" begin
            wide = log_derivatives(a, parameters(vs), x; backend=BACKEND, chunk_size=3)
            @test size(wide) == (n, 2 * n_parameters(a))
            @test wide == log_derivatives(a, parameters(vs), x; backend=BACKEND)
        end

        @testset "a NamedTuple gradient chunks too" begin
            toy = ToyAnsatz(spec, nsites, 3)
            θ = toy_parameters(toy, Xoshiro(7))
            tvs = FullSumState(toy, θ; backend=BACKEND, basis=b)
            _, g = expect_and_grad(tvs, H)
            _, gc = expect_and_grad(tvs, H; chunk_size=5)
            @test gc.W ≈ g.W && gc.b ≈ g.b && gc.v ≈ g.v
        end

        @testset "local_estimators passes it through" begin
            @test local_estimators(vs, H; holomorphic=true, chunk_size=3).O ≈ O
        end

        @testset "a non-positive block size is rejected" begin
            @test_throws ArgumentError log_derivatives(
                a, parameters(vs), x; backend=BACKEND, chunk_size=0
            )
        end
    end

    @testset "ansatz interface" begin
        spec, nsites = Spin(1 // 2),3
        b = basis(dof_object(spec), nsites)
        a = LogStateVector(spec, nsites, b)

        @testset "accessors report the ansatz's own geometry" begin
            @test NQSCore.dof(a) === spec
            @test NQSCore.n_sites(a) == nsites
            @test NQSCore.dof(ToyAnsatz(spec, nsites, 2)) === spec
            @test NQSCore.n_sites(ToyAnsatz(spec, nsites, 2)) == nsites
        end

        @testset "the index is concretely typed" begin
            # An `Any` key type would box every packed state and dispatch hash/isequal
            # dynamically, once per configuration per batch, on the hottest path there is.
            @test isconcretetype(keytype(a.index))
            @test keytype(a.index) === eltype(b.states)
        end

        @testset "a configuration outside the basis is rejected" begin
            sector = basis(dof_object(spec), nsites, sym(TotalMagnetization(1 // 2, nsites), dof_object(spec)))
            restricted = LogStateVector(spec, nsites, sector)
            θ = init_parameters(restricted, Xoshiro(0))
            outside = configurations(spec,b.states, nsites)
            @test_throws ArgumentError log_amplitude(restricted, θ, outside)
        end

        @testset "init_parameters gives one small complex number per basis state" begin
            θ = init_parameters(a, Xoshiro(0); scale=0.01)
            @test length(θ) == n_parameters(a) == length(b.states)
            @test eltype(θ) === ComplexF64
            @test maximum(abs, θ) < 0.1
        end

        @testset "default_basis" begin
            @test default_basis(a) === b
            # Nothing else can guess: an ansatz on a sector would get the wrong answer, so the
            # generic case must ask rather than assume.
            @test_throws ArgumentError default_basis(ToyAnsatz(spec, nsites, 2))
            @test_throws ArgumentError FullSumState(ToyAnsatz(spec, nsites, 2), nothing; backend=BACKEND)
        end
    end

    @testset "interface conformance" begin
        spec, nsites = Spin(1 // 2),3
        vs, a, b = full_sum(spec, nsites)
        H = tfi(nsites)
        for state in (vs, MCState(a, parameters(vs), ExactSampler(b, 64); backend=BACKEND))
            @test ansatz(state) === a
            @test parameters(state) == parameters(vs)
            @test length(samples(state)) > 0
            @test expect(state, H) isa Stats
            @test expect_and_grad(state, H) isa Tuple{Stats,Any}
            @test setparameters!(state, parameters(state)) === state
        end
    end

    @testset "Aqua quality assurance" begin
        Aqua.test_all(NQSCore)
    end
end
