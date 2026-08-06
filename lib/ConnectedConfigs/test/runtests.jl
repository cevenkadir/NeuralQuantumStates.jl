using ConnectedConfigs
using OperatorAlgebra
using SymBasis
using Test

include("models.jl")

loaded_package_names() = Set(m.name for m in keys(Base.loaded_modules))

@testset "ConnectedConfigs.jl" begin
    @testset "dependency weight" begin
        # Tier 1: the local-energy kernel must stay usable by any VMC code, so it carries no
        # autodiff, Lux, or GPU dependency.
        loaded = loaded_package_names()
        for heavy in ("Lux", "Enzyme", "Reactant", "Zygote", "CUDA", "Metal", "ComponentArrays")
            @test heavy ∉ loaded
        end
    end

    @testset "local operators" begin
        @testset "match SymBasis digit ordering" begin
            # The whole point: digit d means local_values[d+1], so Sᶻ must be diag(-s..s) and
            # NOT OperatorAlgebra's PAULI_Z, which puts +1 on the first basis state.
            for s in (1 // 2, 1 // 1, 3 // 2)
                spec = Spin(s)
                ops = local_operators(spec)
                @test local_values(spec) == collect(dof_object(spec).ldof)
                @test [ops.sz[i, i] for i in 1:local_dimension(spec)] == local_values(spec)
            end
            @test local_values(Boson(3)) == 0:3
            @test [local_operators(Boson(3)).n[i, i] for i in 1:4] == 0:3
        end

        @testset "spin algebra" begin
            for s in (1 // 2, 1 // 1, 3 // 2)
                o = local_operators(Spin(s))
                d = local_dimension(Spin(s))
                # [Sˣ, Sʸ] = i Sᶻ, and S² = s(s+1) on every state.
                @test o.sx * o.sy - o.sy * o.sx ≈ im * o.sz
                @test o.sx^2 + o.sy^2 + o.sz^2 ≈ s * (s + 1) * o.id
                @test o.sp ≈ o.sx + im * o.sy
                @test o.sm ≈ o.sx - im * o.sy
            end
        end

        @testset "boson algebra" begin
            for n_max in (1, 3, 5)
                o = local_operators(Boson(n_max))
                # [a, a†] = 1 except on the truncated top state, where the ladder is cut.
                comm = o.a * o.adag - o.adag * o.a
                @test comm[1:n_max, 1:n_max] ≈ o.id[1:n_max, 1:n_max]
                @test o.adag * o.a ≈ o.n
            end
        end
    end

    @testset "packing round-trips" begin
        for (spec, nsites) in ((Spin(1 // 2), 6), (Spin(1 // 1), 4), (Boson(5), 5))
            states = basis(dof_object(spec), nsites).states
            for s in states[1:min(end, 40)]
                @test packed(spec, configurations(spec, s, nsites)) == s
            end
        end

        @testset "a value outside the local set is rejected" begin
            # Catches the classic mistake of writing spins as 0/1 instead of -1//2, 1//2.
            @test_throws ArgumentError packed(Spin(1 // 2), [0, 1, 0, 1])
            @test_throws ArgumentError packed(Boson(2), [0, 1, 7])
        end
    end

    @testset "kernel" begin
        spec, nsites = Spin(1 // 2), 4
        ops = local_operators(spec)

        @testset "a diagonal operator connects only to itself" begin
            H = OpSum([Op(ops.sz, i) for i in 1:nsites])
            for s in basis(dof_object(spec), nsites).states
                d = connected(H, s)
                total = sum(configurations(spec, s, nsites))
                if iszero(total)
                    @test isempty(d)                 # exact zeros are dropped
                else
                    @test collect(keys(d)) == [s]
                    @test only(values(d)) ≈ total
                end
            end
        end

        @testset "a single flip connects to exactly one configuration" begin
            H = Op(2 .* ops.sx, 2)
            s = first(basis(dof_object(spec), nsites).states)
            d = connected(H, s)
            @test length(d) == 1
            flipped = configurations(spec, only(keys(d)), nsites)
            original = configurations(spec, s, nsites)
            @test count(flipped .!= original) == 1
            @test flipped[2] != original[2]
            @test only(values(d)) ≈ 1.0
        end

        @testset "operators on the same site compose, not accumulate" begin
            # S⁺S⁻ is diagonal; if the intermediate index were not summed over correctly this
            # would produce off-diagonal entries.
            H = Op(ops.sp, 1) * Op(ops.sm, 1)
            for s in basis(dof_object(spec), nsites).states
                d = connected(H, s)
                @test all(k == s for k in keys(d))
            end
        end

        @testset "padding is inert" begin
            H = transverse_field_ising(nsites)
            states = basis(dof_object(spec), nsites).states
            res = connected_padded(H, states)
            for b in eachindex(states)
                # Every slot past `counts[b]` must carry a zero matrix element.
                @test all(iszero, res.mels[(res.counts[b]+1):end, b])
                @test size(res.configs) == size(res.mels)
            end
        end
    end
end

include("kernel.jl")
include("golden.jl")
include("exact.jl")
include("symmetry.jl")
