"""
    RBM(nsites, alpha; T=ComplexF64, init=nothing) <: Lux.AbstractLuxLayer

A restricted Boltzmann machine, the canonical neural quantum state.

```math
\\log \\psi(x) = \\sum_i a_i x_i + \\sum_j \\log \\left( 2 \\cosh \\left( b_j
                 + \\sum_i W_{ji} x_i \\right) \\right)
```

The hidden units are summed analytically rather than sampled — that closed form is what makes an
RBM usable as a wavefunction at all, and it is why this is a bespoke layer rather than a `Dense`
composed with something.

`alpha` is the hidden-unit density: there are `alpha * nsites` hidden units.

# Complex parameters by default

`T` defaults to `ComplexF64` because a wavefunction has a phase, and a real RBM can only
represent a positive one. Complex weights give the modulus and the phase together from a single
network, rather than needing two.

`log(2cosh(z))` is evaluated as `logtwocosh` rather than literally, since `cosh` overflows for
|z| beyond about 710 and the layer would otherwise return `Inf` on perfectly ordinary inputs.
"""
struct RBM{F} <: Lux.AbstractLuxLayer
    nsites::Int
    nhidden::Int
    T::Type
    init::F
end

function RBM(nsites::Integer, alpha::Real=1; T::Type=ComplexF64, init=nothing)
    nhidden = round(Int, alpha * nsites)
    nhidden > 0 || throw(ArgumentError("alpha * nsites must give at least one hidden unit"))
    return RBM{typeof(init)}(Int(nsites), nhidden, T, init)
end

"""
    logtwocosh(z)

`log(2 cosh z)`, computed without overflowing.

Factoring out the dominant exponential,
`2cosh(z) = e^{z} + e^{-z} = e^{a}(e^{z-a} + e^{-z-a})` with `a = |Re z|`, keeps both
exponentials bounded by 1 in modulus. A literal `log(2cosh(z))` overflows to `Inf` once
`|Re z|` passes about 710, which an RBM reaches routinely as its weights grow during
optimization — and an `Inf` log-amplitude poisons the entire batch rather than one sample.
"""
function logtwocosh(z::Number)
    a = abs(real(z))
    return a + log(exp(z - a) + exp(-z - a))
end

function Lux.initialparameters(rng::AbstractRNG, layer::RBM)
    T = layer.T
    scale = real(T)(0.01)
    gen() = layer.init === nothing ? scale .* randn(rng, T, 1) : layer.init(rng, T, 1)
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
    # `reshape` rather than `transpose(ps.visible)`: transposing a *vector* makes the reverse
    # pass build an `Adjoint{Transpose{Vector}}`, a wrapper cuBLAS has no method for, so
    # LinearAlgebra silently falls back to its generic matmul — which indexes scalars and is
    # therefore an outright error on a device array. A row-shaped reshape is the same
    # arithmetic and stays in BLAS on both sides of the derivative.
    logψ = reshape(ps.visible, 1, :) * x .+ sum(logtwocosh, θ; dims=1)
    return vec(logψ), st
end

"""
    Jastrow(nsites; T=ComplexF64) <: Lux.AbstractLuxLayer

A Jastrow factor: `log ψ(x) = Σ_{i<j} J_{ij} x_i x_j`.

Captures pair correlations exactly and nothing beyond them. Cheap, interpretable, and a
standard component to multiply into a richer ansatz rather than to use alone.

Only the strict upper triangle is parameterized: `J_{ij}` and `J_{ji}` would multiply the same
product `x_i x_j`, so carrying both makes the parameterization redundant and the geometric
tensor singular.
"""
struct Jastrow{I<:AbstractMatrix{Int}} <: Lux.AbstractLuxLayer
    nsites::Int
    T::Type
    # Where each entry of the coupling matrix comes from: `index[i, j] == k + 1` means the
    # `k`-th coupling, and `1` means the padded zero. Precomputed so that assembling the matrix
    # is one gather rather than a scatter loop.
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
    # `Σ_{i<j} J_ij x_i x_j` written as the quadratic form `xᵀ J x` with `J` strictly upper
    # triangular, evaluated column-wise. That turns an `O(n²)` loop over pairs — each iteration
    # reading one parameter scalar-wise and allocating an accumulator — into one matrix
    # multiplication, which is both a large constant factor faster and the only version that
    # can run on a device array, where a scalar read of a parameter is an error.
    padded = vcat(zero(eltype(ps.coupling)), ps.coupling)
    # The index table has to sit beside the parameters it indexes into; on the CPU this is the
    # table itself and costs nothing.
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

That matters for two reasons. The parameter count drops by roughly the order of the group, and
more importantly the ansatz cannot waste capacity representing states the ground state is known
not to occupy. This is the payoff of deriving symmetry permutations from lattice geometry: the
same code gives a translation-invariant ansatz on a kagome torus as on a chain.
"""
struct SymmetricRBM{P<:AbstractMatrix{Int}} <: Lux.AbstractLuxLayer
    nsites::Int
    nfilters::Int
    # One permutation per **column**, rather than a vector of vectors: a single rectangular
    # array is what can be moved to a device in one piece, and `eachcol` recovers the old view.
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
    # a length-one array, and broadcasting against it gets its value without a scalar read —
    # `ps.visible[1]` would be a host round-trip on a device array. Starting the accumulator
    # from this term also removes the `zeros(...)`, which built a host array regardless of where
    # the parameters lived.
    out = vec(ps.visible .* sum(x; dims=1))
    # Moved once for the whole loop rather than once per group element; on the CPU this is the
    # table itself and the views below are free.
    perms = colocate(ps.weight, layer.permutations)
    for k in axes(perms, 2)
        θ = ps.weight * x[view(perms, :, k), :] .+ ps.hidden
        out = out .+ vec(sum(logtwocosh, θ; dims=1))
    end
    return out, st
end
