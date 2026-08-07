"""
The symmetry-reduced path.

`connected_padded(H, states, basis)` folds each connected configuration back onto the
representative of its symmetry orbit, rescaling by the symmetry character and the orbit-norm
ratio. Getting that rescaling wrong produces a matrix that still *looks* plausible — right
size, right sparsity, Hermitian — but has the wrong spectrum.

The decisive test is therefore spectral: block-diagonalizing a Hamiltonian by symmetry must
partition its spectrum exactly. Concatenating the eigenvalues of every sector has to reproduce
the eigenvalues of the full unreduced matrix, with multiplicities. Nothing weaker would catch a
wrong norm factor.
"""

using ConnectedBasisConfigurations
using LinearAlgebra
using SymBasis
using Test

"""Matrix of `H` restricted to a symmetry sector, via the reduced kernel."""
function sector_matrix(H, b::SymBasis.Bases.Basis)
    states = b.states
    index_of = Dict(s => i for (i, s) in pairs(states))
    res = connected_padded(H, states, b)

    M = zeros(ComplexF64, length(states), length(states))
    for n in eachindex(states), j in 1:res.counts[n]
        M[index_of[res.configs[j, n]], n] += res.mels[j, n]
    end
    return M
end

@testset "symmetry-reduced sectors" begin
    spec = Spin(1 // 2)
    dofo = dof_object(spec)

    @testset "momentum sectors partition the spectrum" begin
        for nsites in (4, 6)
            # A periodic chain's translation generator; LatticeSpaceGroups derives this from
            # geometry, but a bare chain needs no help and keeps this package dependency-free.
            perm = mod1.((1:nsites) .+ 1, nsites)
            H = transverse_field_ising(nsites; J=1.0, h_x=0.7, h_z=0.0)

            full = sort(real(eigvals(Hermitian(dense_from_connected(H, spec, nsites)))))

            sector_eigs = Float64[]
            total_dim = 0
            for k in 0:(nsites-1)
                b = basis(dofo, nsites, sym(Translational(k, perm), dofo))
                isempty(b.states) && continue
                total_dim += length(b.states)
                M = sector_matrix(H, b)
                @test M ≈ M'                       # each block must stay Hermitian
                append!(sector_eigs, real(eigvals(Hermitian(M))))
            end

            @test total_dim == 2^nsites            # the sectors tile the Hilbert space
            @test sort(sector_eigs) ≈ full         # ...and so does the spectrum
        end
    end

    @testset "magnetization sectors partition the spectrum" begin
        # A conserved quantity with no phase factors, isolating the orbit-norm bookkeeping
        # from the character bookkeeping.
        nsites = 6
        H = transverse_field_ising(nsites; J=1.0, h_x=0.0, h_z=0.4)
        full = sort(real(eigvals(Hermitian(dense_from_connected(H, spec, nsites)))))

        sector_eigs = Float64[]
        total_dim = 0
        for two_sz in (-nsites):2:nsites
            b = basis(dofo, nsites, sym(TotalMagnetization(two_sz // 2, nsites), dofo))
            isempty(b.states) && continue
            total_dim += length(b.states)
            append!(sector_eigs, real(eigvals(Hermitian(sector_matrix(H, b)))))
        end

        @test total_dim == 2^nsites
        @test sort(sector_eigs) ≈ full
    end

    @testset "combined magnetization and momentum" begin
        nsites = 6
        perm = mod1.((1:nsites) .+ 1, nsites)
        # Magnetization-conserving: a flip-flop term rather than a transverse field.
        ops = local_operators(spec)
        H = OpSum(vcat(
            [Op(ops.sp, i) * Op(ops.sm, mod1(i + 1, nsites)) for i in 1:nsites],
            [Op(ops.sm, i) * Op(ops.sp, mod1(i + 1, nsites)) for i in 1:nsites],
            [Op(ops.sz, i) * Op(ops.sz, mod1(i + 1, nsites)) for i in 1:nsites],
        ))

        full = sort(real(eigvals(Hermitian(dense_from_connected(H, spec, nsites)))))

        sector_eigs = Float64[]
        total_dim = 0
        for two_sz in (-nsites):2:nsites, k in 0:(nsites-1)
            sg = sym(TotalMagnetization(two_sz // 2, nsites), dofo) ∘
                 sym(Translational(k, perm), dofo)
            b = basis(dofo, nsites, sg)
            isempty(b.states) && continue
            total_dim += length(b.states)
            append!(sector_eigs, real(eigvals(Hermitian(sector_matrix(H, b)))))
        end

        @test total_dim == 2^nsites
        @test sort(sector_eigs) ≈ full
    end

    @testset "samples must be representatives of the basis" begin
        nsites = 4
        perm = mod1.((1:nsites) .+ 1, nsites)
        b = basis(dofo, nsites, sym(Translational(0, perm), dofo))
        H = transverse_field_ising(nsites)
        # A state that is not the representative of its orbit is a caller error, not something
        # to silently fold.
        outsider = packed(spec, [1 // 2, -1 // 2, -1 // 2, -1 // 2])
        if outsider ∉ b.states
            @test_throws ArgumentError connected_padded(H, [outsider], b)
        end
    end
end
