"""
Reference Hamiltonians as plain `OpSum`s, replacing the hand-written
`Operators.build(:TransverseFieldIsing, ...)` / `build(:ExtendedBoseHubbard, ...)` of the
pre-split package.

Both are written to match the pre-split conventions **exactly**, because the golden data was
generated from that code. Where those conventions are surprising, the surprise is documented
rather than quietly corrected.
"""

using ConnectedBasisConfigurations
using OperatorAlgebra
using SymBasis

"""Nearest-neighbour bonds of a periodic chain, each listed once."""
chain_bonds(n::Integer) = [(i, mod1(i + 1, n)) for i in 1:n]

"""
    transverse_field_ising(nsites; J, h_x, h_z)

``H = J \\sum_{\\langle ij \\rangle} \\sigma^z_i \\sigma^z_j + h_z \\sum_i \\sigma^z_i
     + h_x \\sum_i \\sigma^x_i`` on a periodic chain.

Note the factors of two. The pre-split code stored spins as `±1//2` and computed its diagonal
as `4J Σ zᵢzⱼ + 2h_z Σ zᵢ`, which is the same as writing the Hamiltonian in terms of the Pauli
matrices `σ = 2S`. Its off-diagonal element was `h_x` per single-site flip, again matching
`σˣ` rather than `Sˣ`.

`σᶻ` here is `diag(-1, +1)`: digit 0 is spin **down**, following SymBasis. This is the negative
of OperatorAlgebra's `PAULI_Z`, and using that constant instead would silently flip the sign of
the longitudinal field term.
"""
function transverse_field_ising(nsites::Integer; J=1.0, h_x=1.0, h_z=1.0)
    ops = local_operators(Spin(1 // 2))
    σz, σx = 2 .* ops.sz, 2 .* ops.sx

    terms = AbstractOp[]
    for (i, j) in chain_bonds(nsites)
        push!(terms, J * (Op(σz, i) * Op(σz, j)))
    end
    for i in 1:nsites
        push!(terms, h_z * Op(σz, i))
        push!(terms, h_x * Op(σx, i))
    end
    return OpSum(terms)
end

"""
    extended_bose_hubbard(nsites, n_max; J, U, V, μ)

``H = -J \\sum_{\\langle ij \\rangle} (b^\\dagger_i b_j + \\mathrm{h.c.})
     + \\frac{U}{2} \\sum_i n_i (n_i - 1)
     + V \\sum_{\\langle ij \\rangle} n_i n_j - \\mu \\sum_i n_i`` on a periodic chain.

The chemical-potential term enters with a **minus** sign, matching what the pre-split code
computed (`vals₀ .-= μ * sum(x)`) rather than the plus sign its own docstring advertised.
"""
function extended_bose_hubbard(nsites::Integer, n_max::Integer; J=1.0, U=1.0, V=1.0, μ=0.0)
    ops = local_operators(Boson(n_max))
    a, adag, n = ops.a, ops.adag, ops.n
    nn = n * n - n                      # n(n-1), for the on-site interaction

    terms = AbstractOp[]
    for (i, j) in chain_bonds(nsites)
        push!(terms, (-J) * (Op(adag, i) * Op(a, j)))
        push!(terms, (-J) * (Op(adag, j) * Op(a, i)))
        push!(terms, V * (Op(n, i) * Op(n, j)))
    end
    for i in 1:nsites
        push!(terms, (U / 2) * Op(nn, i))
        push!(terms, (-μ) * Op(n, i))
    end
    return OpSum(terms)
end
