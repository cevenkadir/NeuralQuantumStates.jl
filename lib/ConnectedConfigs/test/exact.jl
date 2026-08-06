"""
Cross-check against OperatorAlgebra's own matrix construction.

The golden data proves this package reproduces the *pre-split* implementation. That is
necessary but not sufficient: both could be wrong in the same way. This file compares against
an entirely independent code path — OperatorAlgebra builds the full matrix by Kronecker
products from the operator's algebraic structure, never touching a packed state or a connected
configuration — so agreement between the two is real evidence, not a tautology.

The two use different basis orderings (OperatorAlgebra makes the first site the most
significant digit, SymBasis the least). Rather than translate indices, the comparison is made
on **eigenvalues**, which are invariant under any permutation of basis states.
"""

using ConnectedConfigs
using LinearAlgebra
using OperatorAlgebra
using SparseArrays
using SymBasis
using Test

"""
    dense_from_connected(H, spec, nsites) -> Matrix

Build the full matrix of `H` column by column, using only [`connected_padded`](@ref) — i.e.
through exactly the code path variational Monte Carlo uses.
"""
function dense_from_connected(H, spec, nsites::Integer)
    b = basis(dof_object(spec), nsites)
    states = b.states
    index_of = Dict(s => i for (i, s) in pairs(states))

    res = connected_padded(H, states)
    M = zeros(ComplexF64, length(states), length(states))
    for n in eachindex(states), j in 1:res.counts[n]
        M[index_of[res.configs[j, n]], n] += res.mels[j, n]
    end
    return M
end

"""Matrix of `H` via OperatorAlgebra's independent Kronecker-product construction."""
function dense_from_operator_algebra(H, spec, nsites::Integer)
    d = local_dimension(spec)
    bi = [i => d for i in 1:nsites]
    return Matrix(sparse(H, bi))
end

@testset "cross-check against OperatorAlgebra's matrix construction" begin
    cases = [
        ("TFI 6-site chain", transverse_field_ising(6; J=1.0, h_x=0.7, h_z=0.3), Spin(1 // 2), 6),
        ("TFI 4-site, no field", transverse_field_ising(4; J=1.0, h_x=0.0, h_z=0.0), Spin(1 // 2), 4),
        ("TFI 4-site, pure field", transverse_field_ising(4; J=0.0, h_x=1.0, h_z=0.0), Spin(1 // 2), 4),
        ("BHM 4-site, n_max=2", extended_bose_hubbard(4, 2; J=1.0, U=2.0, V=0.5, μ=0.3), Boson(2), 4),
        ("BHM 3-site, n_max=3", extended_bose_hubbard(3, 3; J=0.8, U=1.5, V=0.0, μ=1.0), Boson(3), 3),
    ]

    for (name, H, spec, nsites) in cases
        @testset "$name" begin
            mine = dense_from_connected(H, spec, nsites)
            theirs = dense_from_operator_algebra(H, spec, nsites)

            @test size(mine) == size(theirs)
            # Basis-ordering-independent invariants.
            @test eigvals(Hermitian(mine)) ≈ eigvals(Hermitian(theirs))
            @test tr(mine) ≈ tr(theirs)
            @test norm(mine) ≈ norm(theirs)
            # A Hermitian Hamiltonian must come out Hermitian.
            @test mine ≈ mine'
        end
    end
end

@testset "analytically known spectra" begin
    @testset "non-interacting spins in a longitudinal field" begin
        # H = h Σ σᶻ has eigenvalues h * (sum of ±1), i.e. h*(2k - N) with multiplicity C(N,k).
        nsites, h = 5, 0.75
        ops = local_operators(Spin(1 // 2))
        H = OpSum([h * Op(2 .* ops.sz, i) for i in 1:nsites])

        found = sort(real(eigvals(Hermitian(dense_from_connected(H, Spin(1 // 2), nsites)))))
        expected = sort(vcat([fill(h * (2k - nsites), binomial(nsites, k)) for k in 0:nsites]...))
        @test found ≈ expected
    end

    @testset "a single boson mode" begin
        # H = ω n on one site has eigenvalues 0, ω, 2ω, ...
        ω, n_max = 1.3, 4
        ops = local_operators(Boson(n_max))
        H = OpSum([ω * Op(ops.n, 1)])
        found = sort(real(eigvals(Hermitian(dense_from_connected(H, Boson(n_max), 1)))))
        @test found ≈ ω .* collect(0:n_max)
    end

    @testset "two-site transverse-field Ising by hand" begin
        # H = J σᶻ₁σᶻ₂ + h_x (σˣ₁ + σˣ₂) on two sites (the "ring" of 2 has a doubled bond,
        # so use an explicit single-bond Hamiltonian instead).
        J, hx = 1.0, 0.6
        ops = local_operators(Spin(1 // 2))
        σz, σx = 2 .* ops.sz, 2 .* ops.sx
        H = OpSum([J * (Op(σz, 1) * Op(σz, 2)), hx * Op(σx, 1), hx * Op(σx, 2)])

        found = sort(real(eigvals(Hermitian(dense_from_connected(H, Spin(1 // 2), 2)))))
        # Basis |↓↓⟩,|↑↓⟩,|↓↑⟩,|↑↑⟩: diagonal (J,-J,-J,J), σˣ couples within each parity block.
        expected = sort([-J, J, -sqrt(J^2 + 4hx^2), sqrt(J^2 + 4hx^2)])
        @test found ≈ expected
    end
end
