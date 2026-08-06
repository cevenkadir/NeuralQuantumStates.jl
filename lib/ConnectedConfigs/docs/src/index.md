```@meta
CurrentModule = ConnectedConfigs
```

# ConnectedConfigs.jl

*The batched local-energy kernel.*

Given an operator ``\hat{H}`` and a batch of basis configurations ``\vert s \rangle``, produce
every connected configuration ``\vert s' \rangle`` together with its matrix element
``\langle s' \vert \hat{H} \vert s \rangle``.

This is the hot loop of variational Monte Carlo — it runs once per sample per optimization step
— and the Julia counterpart of NetKet's `get_conn_padded`.

It depends only on [SymBasis.jl](https://github.com/cevenkadir/SymBasis.jl) and
[OperatorAlgebra.jl](https://github.com/h-mnzlr/OperatorAlgebra.jl): no autodiff, no Lux, no GPU
stack. That is deliberate, and asserted in the test suite — the kernel is useful to any
variational Monte Carlo code, neural or otherwise.

## Why this and not OperatorAlgebra alone

OperatorAlgebra builds operators, normal-orders them, handles fermionic signs, and converts to
sparse or dense matrices. What it does not provide is this: its `apply` consumes a single vector
or a `Dict` of states, not a *batch* of samples. Batching is the whole point here, so this is the
one piece that had to be written rather than reused, and nothing else from that package is
duplicated.

## Quick example

```@example index
using ConnectedConfigs, OperatorAlgebra, SymBasis

spec, nsites = Spin(1 // 2), 4
ops = local_operators(spec)

# Ising ring: σᶻσᶻ on each bond, plus a transverse field.
H = OpSum(vcat(
    [(2 .* ops.sz |> m -> Op(m, i) * Op(m, mod1(i + 1, nsites))) for i in 1:nsites],
    [Op(2 .* ops.sx, i) for i in 1:nsites],
))

states = basis(dof_object(spec), nsites).states
res = connected_padded(H, states[1:3])
res.mels
```

`res.configs` holds the connected configurations as packed integers,
[`configurations`](@ref) unpacks them into physical local values:

```@example index
configurations(spec, res.configs[:, 1], nsites)
```

## Padding

Columns are padded to a common height with a **zero** matrix element rather than `missing`. That
keeps `mels` a concretely typed numeric array, and it makes the local energy

```julia
E_loc(s) = sum(res.mels[:, b] .* exp.(logψ.(res.configs[:, b]) .- logψ(s)))
```

correct with no masking at all — a zero matrix element contributes nothing. `res.counts` records
how many entries of each column are real rather than padding.

## Two conventions worth knowing

Both are easy to get silently wrong.

**Digit position `i` is lattice site `i`**, counting from the least significant digit, which is
what `SymBasis.read` and `eachdigit` use. An operator's site identifier is that position
directly, so sites must be integers `1:nsites`.

**Digit `d` means local value `local_values(spec)[d+1]`** — `(-s, …, +s)` for a spin,
`0:max_occupancy` for a boson. This is *not* the convention behind OperatorAlgebra's `PAULI_Z`,
which puts `+1` on the first basis state where SymBasis puts spin-down. Building a Hamiltonian
from those constants and evaluating it on SymBasis states flips the sign of every ``S^z``, which
is invisible in an unbiased model and glaring the moment a longitudinal field is switched on.
[`local_operators`](@ref) exists to remove that trap: its matrices are derived from the
specification, so they agree with the digit ordering by construction.

## Symmetry sectors

`connected_padded` also takes a symmetry-reduced basis, folding each connected configuration
back onto its orbit representative with the appropriate character and orbit-norm factors:

```julia
b = basis(dofo, nsites, sym(Translational(0, perm), dofo))
connected_padded(H, b.states, b)
```

Getting those factors wrong yields a matrix that is the right size, sparse, and Hermitian, but
has the wrong spectrum — so the test suite checks that the sector spectra partition the full
one exactly.
