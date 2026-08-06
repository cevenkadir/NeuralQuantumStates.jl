```@meta
CurrentModule = ConnectedConfigs
```

# Backends

SymBasis and OperatorAlgebra are this package's defaults, not its requirement. Everything the
compiler and the kernel need from the outside world goes through a handful of generic
functions, and nothing else. Another library plugs in by adding methods to them — ordinary
multiple dispatch, with no package extension to load and no capability flag to set.

There are two independent seams. Implementing one does not oblige you to touch the other.

## The operator seam

Two functions, both consulted **once**, by [`compile`](@ref), and never on the hot path.

| Function | Contract |
|---|---|
| [`expand_terms`](@ref) | Flatten the operator into a plain sum of products of single-site factors |
| [`amplitude_type`](@ref) | The element type matrix elements should accumulate in |

`expand_terms` returns one entry per term; each term is a vector of `position => matrix` pairs
listing that term's factors in matrix-product order, with any scalar coefficient already folded
into a matrix. Because it runs only at compile time, an implementation is free to do arbitrary
work — the default one resolves Jordan–Wigner strings for fermionic sites, which is the
expensive part of the whole pipeline and exactly what compiling exists to hoist out.

The one obligation that is easy to miss: sums nested inside products must be **distributed**
here. `(A + B) * C` is two terms, and the kernel only ever sees terms.

A complete operator backend is therefore about ten lines:

```@example backends
using ConnectedConfigs, SymBasis

struct ToySum
    terms::Vector{Vector{Pair{Int,Matrix{Float64}}}}
end

ConnectedConfigs.expand_terms(op::ToySum) = op.terms
ConnectedConfigs.amplitude_type(::ToySum) = Float64
nothing # hide
```

That is enough to run the whole kernel:

```@example backends
spec, nsites = Spin(1 // 2), 4
ops = local_operators(spec)
σx, σz = Matrix(2 .* ops.sx), Matrix(2 .* ops.sz)

toy = ToySum([
    [[i => σz, mod1(i + 1, nsites) => σz] for i in 1:nsites]...,
    [[i => σx] for i in 1:nsites]...,
])

states = basis(dof_object(spec), nsites).states
res = connected_padded(toy, states)
res.mels[:, 1]
```

An empty factor list is the identity, which is how a term that cancels to a scalar is
expressed.

## The state seam

| Function | Contract |
|---|---|
| [`read_digit`](@ref) | The zero-based local digit at a position |
| [`write_digit`](@ref) | A copy of the configuration with one digit replaced |

These two *are* on the hot path — they run in the innermost loop — so an implementation should
be allocation-free and inlineable. Packed states are treated as immutable values, so
`write_digit` returns a new one rather than mutating.

## The symmetry seam

Only needed if you want the symmetry-reduced path with a basis type of your own.

| Function | Contract |
|---|---|
| [`fold_state`](@ref) | Map a configuration onto its orbit representative, with the phase |
| [`state_lookup`](@ref) | Build, once, whatever finds a state's index |
| [`lookup_index`](@ref) | A state's index, or `0` when it is absent |
| [`orbit_norms`](@ref) | The orbit norms, indexed as `lookup_index` indexes them |

`lookup_index` has a default for any `AbstractDict`, so a backend whose `state_lookup` returns
a dictionary implements three functions, not four. Hoisting `state_lookup` out of the kernel is
the entire point of its existence: the pre-compilation version of this package rebuilt an
equivalent dictionary over the whole basis on every single call.

## Where the defaults live

`src/backends/symbasis.jl` and `src/backends/operatoralgebra.jl`, deliberately separate from
the kernel so the seam is visible in the file layout rather than only in the documentation. The
compiler and the kernel consume nothing but the functions above, which is what makes them
backend-agnostic by construction rather than by good intentions — and the test suite runs a toy
operator backend, defined without any OperatorAlgebra type, to keep it that way.
