```@meta
CurrentModule = ConnectedConfigs
```

# ConnectedConfigs.jl

*The batched connected-configurations kernel for quantum lattice models.*

Variational Monte Carlo spends nearly all its time answering one question, several thousand
times per optimization step: given a configuration ``\vert s \rangle`` and a Hamiltonian
``\hat{H}``, which configurations ``\vert s' \rangle`` does ``\hat{H}`` connect it to, and with
what matrix element ``\langle s' \vert \hat{H} \vert s \rangle``? This package answers it for a
whole **batch** at once, returning padded arrays rather than a dictionary per sample — the
Julia counterpart of NetKet's `get_conn_padded`.

Its only dependencies are [SymBasis.jl](https://github.com/cevenkadir/SymBasis.jl) and
[OperatorAlgebra.jl](https://github.com/cevenkadir/OperatorAlgebra.jl), and even those are
replaceable: the kernel reaches them through a small documented interface, described in
[Backends](@ref).

## Key features

- **Compiled operators** — [`compile`](@ref) hoists every configuration-independent decision
  out of the hot loop: Jordan–Wigner strings are resolved, sums are distributed over products,
  same-site factors are multiplied out, and every local matrix is column-compressed. What is
  left runs with no hash tables and no per-sample allocation.
- **Batches of any shape** — [`connected_padded`](@ref) takes an array of configurations with
  any leading dimensions and prepends the connection axis, so a `(steps, chains)` block of
  Monte Carlo samples needs no reshaping on the way in or out.
- **Padding you do not have to mask** — unused slots carry a zero matrix element and repeat the
  sample itself, so a local energy is a plain sum over the connection axis.
- **An in-place form** — [`connected_padded!`](@ref) writes into caller-owned buffers sized by
  [`max_conn_size`](@ref), so an optimization loop stops reallocating its output every step.
- **Symmetry sectors** — pass a symmetry-reduced basis and matrix elements come back folded
  onto orbit representatives, rescaled by the character and the orbit norms.
- **Swappable backends** — SymBasis and OperatorAlgebra are the defaults, not the requirement.
- **Almost no dependencies** — no autodiff, no Lux, no GPU stack, so the kernel stays usable by
  any Monte Carlo or exact-diagonalization code. The test suite asserts it.

## Installation

**Requirements**: Julia 1.11 or later.

To install the latest stable version, use the Julia package manager. Either use the Julia REPL
package mode (by pressing `]`):

```julia
pkg> add ConnectedConfigs
```

or open the Julia REPL and run the following command:

```julia
julia> import Pkg; Pkg.add("ConnectedConfigs")
```

## Quick example

Build a transverse-field Ising chain and compile it once:

```@example index
using ConnectedConfigs, OperatorAlgebra, SymBasis

spec, nsites = Spin(1 // 2), 4
ops = local_operators(spec)          # matrices in SymBasis's digit ordering

H = OpSum(vcat(
    [Op(2 .* ops.sz, i) * Op(2 .* ops.sz, mod1(i + 1, nsites)) for i in 1:nsites],
    [Op(2 .* ops.sx, i) for i in 1:nsites],
))

compiled = compile(H)
```

Ask what a batch of configurations connects to. The diagonal comes first, then one spin flip
per site:

```@example index
states = basis(dof_object(spec), nsites).states
res = connected_padded(compiled, states)

size(res.configs), res.mels[:, 1]
```

The batch axis follows the input, so a block of Monte Carlo samples goes straight in:

```@example index
block = reshape(states, 4, 4)
r = connected_padded(compiled, block)

size(r.configs), size(r.counts)
```

And packed configurations turn into the numeric array a neural network wants, with the
degree-of-freedom axis first:

```@example index
size(configurations(spec, res.configs, nsites))
```

## Why compiling matters

Every call that receives a bare operator has to flatten it first. Over a batch that cost is
amortized; over the two single-configuration queries a Hamiltonian-driven Metropolis rule makes
per step, it is the entire runtime. Compiling once, outside the loop, is the difference:

| Case | Bare operator | Compiled |
|---|---|---|
| 16-site Ising, batch of 1024 | 318 µs | 278 µs |
| 10⁴ single-configuration queries | 484 ms | 5.8 ms |

## Scope

Building operators — algebra, normal ordering, fermionic signs, sparse and dense conversion —
belongs upstream in OperatorAlgebra.jl, and this package deliberately duplicates none of it.
What it adds is the one thing an operator library has no reason to provide: evaluation over a
*batch of samples* rather than over a state vector.

Sampling lives in NQSSamplers.jl and wavefunctions in NQSCore.jl; neither is a dependency, and
neither needs to be for this package to be useful on its own.

## Part of a larger ecosystem

ConnectedConfigs.jl is developed in the
[NeuralQuantumStates.jl](https://github.com/cevenkadir/NeuralQuantumStates.jl) monorepo, and is
usable entirely on its own — it has no idea the neural-network stack exists, and its test suite
asserts as much.

## Important notice

This project is still under active development. While it includes an extensive test suite and
is developed with high scientific rigor, you should always benchmark your own code. Please
report any issues you encounter via the
[GitHub issue tracker](https://github.com/cevenkadir/NeuralQuantumStates.jl/issues/new).

## Module

What `?ConnectedConfigs` shows at the REPL:

```@docs
ConnectedConfigs
```
