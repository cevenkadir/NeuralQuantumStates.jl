"""
    RBM(nsites, alpha=1; T=ComplexF64) <: Lux.AbstractLuxLayer

A restricted Boltzmann machine, the canonical neural quantum state.

```math
\\log \\psi(x) = \\sum_i a_i x_i + \\sum_j \\log \\left( 2 \\cosh \\left( b_j
                 + \\sum_i W_{ji} x_i \\right) \\right)
```

The hidden units are summed analytically rather than sampled, which is what makes an RBM usable
as a wavefunction and why this is a bespoke layer rather than a composed `Dense`. `alpha` is the
hidden-unit density: there are `alpha * nsites` hidden units.

`T` defaults to `ComplexF64` because a wavefunction has a phase and a real RBM can only
represent a positive one; complex weights give modulus and phase from a single network.
"""
struct RBM <: Lux.AbstractLuxLayer
    nsites::Int
    nhidden::Int
    T::Type
end

function RBM(nsites::Integer, alpha::Real=1; T::Type=ComplexF64)
    nhidden = round(Int, alpha * nsites)
    nhidden > 0 || throw(ArgumentError("alpha * nsites must give at least one hidden unit"))
    return RBM(Int(nsites), nhidden, T)
end

"""
    logtwocosh(z)

`log(2 cosh z)`, computed without overflowing.

Factoring out the dominant exponential — `2cosh(z) = e^{a}(e^{z-a} + e^{-z-a})` with
`a = |Re z|` — keeps both exponentials bounded by 1 in modulus. A literal `log(2cosh(z))`
overflows to `Inf` once `|Re z|` passes about 710, which an RBM reaches routinely as its weights
grow, and an `Inf` log-amplitude poisons the whole batch rather than one sample.
"""
function logtwocosh(z::Number)
    a = abs(real(z))
    return a + log(exp(z - a) + exp(-z - a))
end

function Lux.initialparameters(rng::AbstractRNG, layer::RBM)
    T = layer.T
    scale = real(T)(0.01)
    return (
        visible=scale .* randn(rng, T, layer.nsites),
        hidden=scale .* randn(rng, T, layer.nhidden),
        weight=scale .* randn(rng, T, layer.nhidden, layer.nsites),
    )
end

Lux.initialstates(::AbstractRNG, ::RBM) = NamedTuple()
Lux.parameterlength(l::RBM) = l.nsites + l.nhidden + l.nhidden * l.nsites
Lux.statelength(::RBM) = 0

function (layer::RBM)(x::AbstractMatrix, ps, st)
    # x is (nsites, batch); θ = W x .+ b is (nhidden, batch).
    θ = ps.weight * x .+ ps.hidden
    # `reshape`, not `transpose(ps.visible)`: transposing a vector makes the reverse pass build
    # an `Adjoint{Transpose{Vector}}`, which cuBLAS has no method for, so LinearAlgebra falls
    # back to a generic matmul that indexes scalars — an error on a device array.
    logψ = reshape(ps.visible, 1, :) * x .+ sum(logtwocosh, θ; dims=1)
    return vec(logψ), st
end

"""
    Jastrow(nsites; T=ComplexF64) <: Lux.AbstractLuxLayer

A Jastrow factor: `log ψ(x) = Σ_{i<j} J_{ij} x_i x_j`.

Captures pair correlations exactly and nothing beyond them — a standard component to multiply
into a richer ansatz rather than to use alone.

Only the strict upper triangle is parameterized: `J_{ij}` and `J_{ji}` multiply the same product
`x_i x_j`, so carrying both would make the geometric tensor singular.
"""
struct Jastrow{I<:AbstractMatrix{Int}} <: Lux.AbstractLuxLayer
    nsites::Int
    T::Type
    # `index[i, j] == k + 1` means the `k`-th coupling; `1` means the padded zero. Precomputed
    # so assembling the matrix is one gather rather than a scatter loop.
    index::I
end

function Jastrow(nsites::Integer; T::Type=ComplexF64)
    n = Int(nsites)
    index = ones(Int, n, n)          # 1 selects the padded zero
    k = 0
    for i in 1:(n-1), j in (i+1):n
        k += 1
        index[i, j] = k + 1
    end
    return Jastrow(n, T, index)
end

function Lux.initialparameters(rng::AbstractRNG, layer::Jastrow)
    n = layer.nsites
    return (coupling=real(layer.T)(0.01) .* randn(rng, layer.T, n * (n - 1) ÷ 2),)
end
Lux.initialstates(::AbstractRNG, ::Jastrow) = NamedTuple()
Lux.parameterlength(l::Jastrow) = l.nsites * (l.nsites - 1) ÷ 2
Lux.statelength(::Jastrow) = 0

function (layer::Jastrow)(x::AbstractMatrix, ps, st)
    # `Σ_{i<j} J_ij x_i x_j` as the quadratic form `xᵀ J x` with `J` strictly upper triangular.
    # The loop over pairs it replaces would read parameters scalar-wise, which is an error on a
    # device array.
    padded = vcat(zero(eltype(ps.coupling)), ps.coupling)
    J = padded[colocate(padded, layer.index)]
    return vec(sum(x .* (J * x); dims=1)), st
end

"""
    SymmetricRBM(permutations, alpha; T=ComplexF64) <: Lux.AbstractLuxLayer

An RBM whose weights are shared across a symmetry group, so that `|ψ|` is invariant under it.

`permutations` is a vector of site permutations — exactly what
`LatticeSpaceGroups.site_permutation` produces — and each hidden filter is applied to every
permuted copy of the input. The resulting wavefunction respects the symmetry **by
construction**, rather than having to learn it.

The parameter count drops by roughly the order of the group, and the ansatz cannot waste
capacity on states the ground state is known not to occupy.
"""
struct SymmetricRBM{P<:AbstractMatrix{Int}} <: Lux.AbstractLuxLayer
    nsites::Int
    nfilters::Int
    # One permutation per column rather than a vector of vectors, so the whole table moves to a
    # device in one piece.
    permutations::P
    T::Type
end

function SymmetricRBM(permutations::AbstractVector, alpha::Real=1; T::Type=ComplexF64)
    isempty(permutations) && throw(ArgumentError("need at least one permutation"))
    nsites = length(first(permutations))
    all(length(p) == nsites for p in permutations) ||
        throw(ArgumentError("all permutations must have the same length"))
    nfilters = max(1, round(Int, alpha))
    table = Matrix{Int}(undef, nsites, length(permutations))
    for (k, p) in pairs(permutations)
        table[:, k] = p
    end
    return SymmetricRBM(nsites, nfilters, table, T)
end

function Lux.initialparameters(rng::AbstractRNG, layer::SymmetricRBM)
    scale = real(layer.T)(0.01)
    return (
        visible=scale .* randn(rng, layer.T, 1),
        hidden=scale .* randn(rng, layer.T, layer.nfilters),
        weight=scale .* randn(rng, layer.T, layer.nfilters, layer.nsites),
    )
end
Lux.initialstates(::AbstractRNG, ::SymmetricRBM) = NamedTuple()
Lux.parameterlength(l::SymmetricRBM) = 1 + l.nfilters + l.nfilters * l.nsites
Lux.statelength(::SymmetricRBM) = 0

function (layer::SymmetricRBM)(x::AbstractMatrix, ps, st)
    # A symmetric visible bias couples to the total, which is itself invariant. `ps.visible` is
    # a length-one array; broadcasting against it avoids `ps.visible[1]`, which would be a host
    # round-trip on a device array.
    out = vec(ps.visible .* sum(x; dims=1))
    perms = colocate(ps.weight, layer.permutations)      # once, not once per group element
    for k in axes(perms, 2)
        θ = ps.weight * x[view(perms, :, k), :] .+ ps.hidden
        out = out .+ vec(sum(logtwocosh, θ; dims=1))
    end
    return out, st
end
