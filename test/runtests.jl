using NeuralQuantumStates

using LinearAlgebra
using Test

# Reference outputs from the pre-split `Operators.connected_basis_configs`, which has since
# been deleted (recover it with `git show 0774c25:archive/pre-split/`). Reproducing them via the new
# stack — lattice bonds, model definition, and local-energy kernel together — is what licensed
# that deletion, and this is the end-to-end statement of it.
# The data lives inside ConnectedBasisConfigurations so that package stays self-contained when
# published; the umbrella is only ever built from this repository, so reaching into `lib/` is
# safe here.
const GOLDEN_DIR =
    joinpath(@__DIR__, "..", "lib", "ConnectedBasisConfigurations", "test", "golden")
include(joinpath(GOLDEN_DIR, "tfi_chain8.jl"))
include(joinpath(GOLDEN_DIR, "bhm_chain16.jl"))

"""
    reference_dict(configs, mels) -> Dict

Reduce pre-split `(configs, mels)` to configuration => total matrix element, skipping the
`missing` padding, summing duplicates, and dropping exact zeros. See
`lib/ConnectedBasisConfigurations/test/golden.jl` for why the comparison is made this way
rather than element-wise.
"""
function reference_dict(configs::AbstractMatrix, mels::AbstractVector)
    out = Dict{Vector{eltype(configs)},Float64}()
    for c in axes(configs, 2)
        col = configs[:, c]
        (any(ismissing, col) || ismissing(mels[c])) && continue
        key = collect(skipmissing(col))
        out[key] = get(out, key, 0.0) + mels[c]
    end
    filter!(p -> !iszero(p.second), out)
    return out
end

"""Connected configurations of `model` for one configuration, in the same dictionary form."""
function model_dict(model, config)
    state = packed(model.dof, config)
    res = connected_padded(model.hamiltonian, [state])
    out = Dict{Vector{eltype(local_values(model.dof))},Float64}()
    for j in 1:res.counts[1]
        key = configurations(model.dof, res.configs[j, 1], model.nsites)
        out[key] = get(out, key, 0.0) + real(res.mels[j, 1])
    end
    filter!(p -> !iszero(p.second), out)
    return out
end

"""Dense Hamiltonian of `model`, built through the local-energy kernel."""
function dense(model)
    states = basis(model).states
    index_of = Dict(s => i for (i, s) in pairs(states))
    res = connected_padded(model.hamiltonian, states)
    M = zeros(ComplexF64, length(states), length(states))
    for n in eachindex(states), j in 1:res.counts[n]
        M[index_of[res.configs[j, n]], n] += res.mels[j, n]
    end
    return M
end

@testset "NeuralQuantumStates.jl" begin
    @testset "the stack is re-exported" begin
        # A user should need only `using NeuralQuantumStates` to reach the whole ecosystem.
        @test isdefined(@__MODULE__, :Hypercube)          # LatticeSpaceGroups
        @test isdefined(@__MODULE__, :connected_padded)   # ConnectedBasisConfigurations
        @test isdefined(@__MODULE__, :OpSum)              # OperatorAlgebra
        @test isdefined(@__MODULE__, :dof_object)         # SymBasis
    end

    @testset "dependency weight" begin
        # The umbrella deliberately pulls in the whole stack, Lux included -- that is what an
        # umbrella is for. The dependency-weight property belongs to the Tier 1 packages, and
        # is asserted in *their* test suites: `using LatticeSpaceGroups` or
        # `using ConnectedBasisConfigurations` must load none of this. Repeating that assertion
        # here would be testing the wrong package.
        loaded = Set(m.name for m in keys(Base.loaded_modules))
        @test "Lux" in loaded                       # via NQSAnsatze

        # A GPU backend is still not a dependency of anything in the stack.
        for gpu in ("CUDA", "Metal", "AMDGPU")
            @test gpu ∉ loaded
        end
    end

    @testset "lattices" begin
        lat = build(Hypercube([8]; periodic=true))
        @test n_sites(lat) == 8
        @test length(bonds(lat)) == 8                     # a ring has as many bonds as sites
        @test length(bonds(build(Hypercube([8]; periodic=false)))) == 7
        @test n_sites(build(Kagome([2, 2], 1.0))) == 12
        @test n_sites(build(Honeycomb([2, 2], 1.0))) == 8
    end

    @testset "models reproduce the pre-split implementation" begin
        @testset "transverse-field Ising" begin
            lat = build(Hypercube([8]; periodic=true))
            model = build(TransverseFieldIsing(lat; J=1.0, h_x=1.0, h_z=1.0))

            @test model.nsites == 8
            @test length(basis(model).states) == 2^8

            reference = reference_dict(tfi_single_configs, tfi_single_mels)
            @test model_dict(model, tfi_sample_vec) == reference
            @test length(reference) == 1 + model.nsites   # diagonal plus one flip per site
        end

        @testset "extended Bose-Hubbard" begin
            lat = build(Hypercube([16]; periodic=true))
            model = build(ExtendedBoseHubbard(lat, 5; J=1.0, U=1.0, V=1.0, μ=0.0))

            reference = reference_dict(bhm_single_configs, bhm_single_mels)
            @test model_dict(model, bhm_sample_vec) == reference
            @test length(reference) == 11                 # 1 diagonal + 10 hops
        end
    end

    @testset "models are Hermitian with sensible spectra" begin
        @testset "Ising" begin
            lat = build(Hypercube([6]; periodic=true))
            H = dense(build(TransverseFieldIsing(lat; J=1.0, h_x=0.5, h_z=0.2)))
            @test H ≈ H'
            # Every diagonal Ising energy lies within the bounds set by the couplings.
            @test all(abs.(eigvals(Hermitian(H))) .<= 6 * (1.0 + 0.5 + 0.2) + 1e-8)
        end

        @testset "Bose-Hubbard conserves particle number" begin
            lat = build(Hypercube([4]; periodic=true))
            model = build(ExtendedBoseHubbard(lat, 2; J=1.0, U=1.0, V=0.5, μ=0.0))
            states = basis(model).states
            res = connected_padded(model.hamiltonian, states)
            for n in eachindex(states)
                total = sum(configurations(model.dof, states[n], model.nsites))
                for j in 1:res.counts[n]
                    cfg = configurations(model.dof, res.configs[j, n], model.nsites)
                    @test sum(cfg) == total
                end
            end
        end
    end

    @testset "symmetry sectors partition the spectrum" begin
        # The end-to-end statement: lattice geometry produces the translation permutation,
        # SymBasis builds the sector, and the reduced kernel gives a block whose eigenvalues
        # are part of the full spectrum.
        nsites = 6
        lat = build(Hypercube([nsites]; periodic=true))
        model = build(TransverseFieldIsing(lat; J=1.0, h_x=0.7, h_z=0.0))
        dofo = dof_object(model.dof)

        full = sort(real(eigvals(Hermitian(dense(model)))))

        sector_eigs = Float64[]
        total_dim = 0
        for k in 0:(nsites-1)
            b = basis(model, sym(Translational(k, lat), dofo))
            isempty(b.states) && continue
            total_dim += length(b.states)

            index_of = Dict(s => i for (i, s) in pairs(b.states))
            res = connected_padded(model.hamiltonian, b.states, b)
            M = zeros(ComplexF64, length(b.states), length(b.states))
            for n in eachindex(b.states), j in 1:res.counts[n]
                M[index_of[res.configs[j, n]], n] += res.mels[j, n]
            end
            @test M ≈ M'
            append!(sector_eigs, real(eigvals(Hermitian(M))))
        end

        @test total_dim == 2^nsites
        @test sort(sector_eigs) ≈ full
    end
end

# ---------------------------------------------------------------------------------------
# End-to-end: the acceptance criterion for the whole package split.
#
# A full variational Monte Carlo run through the umbrella driver -- lattice geometry, model
# definition, operator kernel, ansatz, sampler, preconditioner, optimizer -- checked against
# exact diagonalization. Every package in the stack has to be right simultaneously for this
# to pass; it is the first end-to-end VMC result the project produces.
# ---------------------------------------------------------------------------------------

using DifferentiationInterface
using ForwardDiff
using Random

"""Exact ground-state energy by dense diagonalization."""
function exact_ground_energy(model)
    states = basis(model).states
    index = Dict(s => i for (i, s) in pairs(states))
    res = connected_padded(model.hamiltonian, states)
    M = zeros(ComplexF64, length(states), length(states))
    for n in eachindex(states), j in 1:res.counts[n]
        M[index[res.configs[j, n]], n] += res.mels[j, n]
    end
    return minimum(real(eigvals(Hermitian(M))))
end

@testset "end-to-end variational Monte Carlo" begin
    # Cost note: the exact ansatz has one parameter per basis state, so its Jacobian and
    # geometric tensor grow as the Hilbert space does. It is therefore used at N=6, where it
    # can be checked to machine precision, while N=10 uses a genuinely compressed RBM -- which
    # is the realistic case anyway.

    @testset "an exact ansatz reproduces exact diagonalization" begin
        nsites = 6
        lat = build(Hypercube([nsites]; periodic=true))
        model = build(TransverseFieldIsing(lat; J=1.0, h_x=1.0))
        E_exact = exact_ground_energy(model)

        b = basis(model)
        @test length(b.states) == 2^nsites

        a = LogStateVector(model.dof, nsites, b)
        vs = FullSumState(a, init_parameters(a, Xoshiro(0); scale=0.1), AutoForwardDiff())

        log = run!(VMC(vs, model.hamiltonian;
                preconditioner=StochasticReconfiguration(; diag_shift=1e-2),
                optimizer=Descent(0.05));
            iterations=2000, callbacks=(InvalidLossStopping(),))

        @test real(final_energy(log).mean) < real(log[1].mean)      # it descended
        @test isapprox(real(final_energy(log).mean), E_exact; atol=1e-6)
        @test final_energy(log).variance < 1e-8                     # ...to an eigenstate
    end

    @testset "a compressed RBM at N=10" begin
        nsites = 10
        lat = build(Hypercube([nsites]; periodic=true))
        model = build(TransverseFieldIsing(lat; J=1.0, h_x=2.0))
        E_exact = exact_ground_energy(model)
        b = basis(model)

        a = LuxAnsatz(RBM(nsites, 1), model.dof, nsites; rng=Xoshiro(3))
        # Genuinely compressed: far fewer parameters than the Hilbert space has dimensions.
        @test NQSCore.n_parameters(a) < length(b.states) ÷ 4

        vs = FullSumState(a, init_parameters(a, Xoshiro(3)), AutoForwardDiff())
        log = run!(VMC(vs, model.hamiltonian;
                preconditioner=StochasticReconfiguration(; diag_shift=1e-3),
                optimizer=Descent(0.05));
            iterations=300, callbacks=(InvalidLossStopping(),))

        E = real(final_energy(log).mean)
        @test E >= E_exact - 1e-8                                   # variational bound holds
        @test isapprox(E, E_exact; rtol=1e-3)

        @testset "sampling agrees with exact summation" begin
            # The same optimized parameters, measured by Metropolis instead of summed exactly.
            exact = expect(vs, model.hamiltonian)

            starts = random_configurations(model.dof, nsites, 8, Xoshiro(1))
            sampler = MetropolisSampler(LocalRule(), starts;
                n_chains=8, n_samples=8_000, burn_in=1_000)
            mc = MCState(a, parameters(vs), sampler;
                backend=AutoForwardDiff(), rng=Xoshiro(2))
            sampled = expect(mc, model.hamiltonian)

            @test sampled.error_of_mean > 0
            @test abs(real(sampled.mean - exact.mean)) < 5 * sampled.error_of_mean
            @test sampled.r_hat < 1.05
        end
    end
end

@testset "driver callbacks" begin
    nsites = 6
    lat = build(Hypercube([nsites]; periodic=true))
    model = build(TransverseFieldIsing(lat; J=1.0, h_x=2.0))
    b = basis(model)
    a = LogStateVector(model.dof, nsites, b)

    @testset "EarlyStopping halts a converged run" begin
        vs = FullSumState(a, init_parameters(a, Xoshiro(4); scale=0.1), AutoForwardDiff())
        log = run!(VMC(vs, model.hamiltonian;
                preconditioner=StochasticReconfiguration(; diag_shift=1e-3),
                optimizer=Descent(0.1));
            iterations=2000, callbacks=(EarlyStopping(; patience=20, min_delta=1e-10),))
        @test length(log) < 2000                             # it stopped early
    end

    @testset "InvalidLossStopping halts a diverged run" begin
        # A wildly oversized step makes the energy blow up; the run must stop rather than
        # propagate NaN through every parameter and report a finished run.
        vs = FullSumState(a, init_parameters(a, Xoshiro(5); scale=0.1), AutoForwardDiff())
        log = run!(VMC(vs, model.hamiltonian; optimizer=Descent(1e9));
            iterations=500, callbacks=(InvalidLossStopping(),))
        @test length(log) < 500
    end

    @testset "a user callback can stop the run" begin
        vs = FullSumState(a, init_parameters(a, Xoshiro(6); scale=0.1), AutoForwardDiff())
        seen = Int[]
        log = run!(VMC(vs, model.hamiltonian); iterations=100,
            callbacks=((it, stats, state) -> (push!(seen, it); it < 7),))
        @test length(log) == 7
        @test seen == 1:7
    end
end
