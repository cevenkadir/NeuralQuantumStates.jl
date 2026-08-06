```@meta
CurrentModule = ConnectedConfigs
```

# Local operators and configurations

Two conventions govern how a packed configuration is read, and both are easy to get silently
wrong. [`local_operators`](@ref) exists to remove the second one as a source of error.

## Digit position `i` is lattice site `i`

Positions count from the **least significant** digit, which is what `SymBasis.read` and
`eachdigit` use. An operator's site identifier is used directly as that position, so sites must
be integers `1:nsites`.

!!! warning "OperatorAlgebra numbers digits the other way"
    OperatorAlgebra's `sparse` and `apply` build their Kronecker products with site 1 as the
    **most** significant digit. The two matrices are therefore related by a permutation of the
    basis, not equal. Chain models that are symmetric under reversing the site order hide the
    difference completely, so a comparison that happens to pass proves less than it looks like
    it does — compare eigenvalues, traces and norms, or translate the index explicitly.

## Digit `d` means local value `local_values(spec)[d + 1]`

For a spin these are the magnetic quantum numbers ``(-s, \dots, +s)``; for a boson the
occupation numbers `0:max_occupancy`.

```@example ops
using ConnectedConfigs, SymBasis

local_values(Spin(1 // 2)), local_values(Boson(3))
```

This is *not* the convention behind OperatorAlgebra's `PAULI_Z`, which puts `+1` on the first
basis state where SymBasis puts spin **down**. Building a Hamiltonian from those constants and
evaluating it on SymBasis states flips the sign of every ``S^z`` — invisible in an unbiased
model, and glaring the moment a longitudinal field is switched on.

[`local_operators`](@ref) derives its matrices from the specification itself, so they agree
with the digit ordering by construction:

```@example ops
ops = local_operators(Spin(1 // 2))
ops.sz
```

Note the sign: ``S^z`` is `diag(-1/2, +1/2)`, the negative of `PAULI_Z / 2`.

The constants are also 2×2 only, so they cannot describe a boson with `max_occupancy > 1` or a
spin above 1/2 at all:

```@example ops
local_operators(Boson(3)).n
```

## What each specification returns

| Specification | Fields |
|---|---|
| `Spin` | `sx`, `sy`, `sz`, `sp` (``S^+``), `sm` (``S^-``), `id` |
| `Boson` | `a`, `adag`, `n`, `id` |

Matrices follow the usual `mat[out, in]` convention, so `Op(ops.adag, 1) * Op(ops.a, 2)` hops a
particle from site 2 to site 1.

## Packing and unpacking

[`configurations`](@ref) turns packed integers into physical local values, and
[`packed`](@ref) turns them back. The degree-of-freedom axis comes first throughout, because
Julia is column-major and that puts each configuration in contiguous memory — the layout a
batched forward pass wants.

```@example ops
spec, nsites = Boson(3), 4
s = packed(spec, [0, 2, 1, 3])
configurations(spec, s, nsites)
```

`packed` rejects a value that is not one of the specification's local values, which catches the
classic mistake of writing spins as `0, 1` rather than `-1//2, 1//2`:

```@example ops
try
    packed(Spin(1 // 2), [0, 1, 0, 1])
catch err
    err
end
```
