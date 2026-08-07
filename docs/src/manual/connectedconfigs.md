```@meta
CurrentModule = NeuralQuantumStates
```

# Connected configurations

*Provided by `ConnectedConfigs`, which has [its own
site](https://cevenkadir.github.io/NeuralQuantumStates.jl/ConnectedConfigs/dev/).*

Variational Monte Carlo spends nearly all its time answering one question, several thousand
times per optimization step: given a configuration ``\vert s \rangle`` and a Hamiltonian
``\hat{H}``, which ``\vert s' \rangle`` does it connect to, and with what matrix element
``\langle s' \vert \hat{H} \vert s \rangle``? That is the local energy,

```math
E_{\mathrm{loc}}(s) = \sum_{s'} \langle s' \vert \hat{H} \vert s \rangle \,
                      \frac{\psi(s')}{\psi(s)}
```

and `ConnectedConfigs` answers it for a whole **batch** at once, returning padded arrays rather
than a dictionary per sample — the Julia counterpart of NetKet's `get_conn_padded`.

It is re-exported, so everything below is available from `using NeuralQuantumStates`.

```@example cc
using NeuralQuantumStates

lat = build(Hypercube([6]; periodic=true))
model = build(TransverseFieldIsing(lat; J=1.0, h_x=1.0))
states = basis(model).states
nothing # hide
```

## Compile once

`compile` flattens an operator into the form the kernel runs on: Jordan–Wigner strings
resolved, sums distributed over products, same-site factors multiplied out, every local matrix
column-compressed. Everything that does not depend on the configuration is decided here, once.

```@example cc
compiled = compile(model.hamiltonian)
```

Passing the Hamiltonian straight to `connected_padded` works too — it just compiles on every
call. Over a large batch that is amortized; over the two single-configuration queries a
Hamiltonian-driven Metropolis rule makes per step, it *is* the runtime. Compile whenever the
operator outlives the call, which in an optimization loop it always does.

## The batched call

```@example cc
res = connected_padded(compiled, states)
size(res.configs), size(res.mels), size(res.counts)
```

`configs` holds the connected configurations as packed integers, `mels` their matrix elements,
and `counts` how many entries per sample are real rather than padding. The connection axis
comes **first**: NetKet puts it last, but Julia is column-major, and this is what keeps one
sample's connections contiguous.

`states` may be an array of any shape — the connection axis is prepended and the batch shape
preserved — so a `(steps, chains)` block of Monte Carlo samples goes in without reshaping.

## Two conventions that matter

**Padding is inert.** Unused slots carry a matrix element of zero *and* repeat the sample
itself. The zero makes the local energy correct with no masking; the repeat means every padded
row is still a valid configuration, so the wavefunction can be evaluated on the whole array
without special-casing.

```@example cc
b = argmin(res.counts)
res.counts[b], res.mels[(res.counts[b]+1):end, b]
```

**Repeated configurations are not merged.** Two terms reaching the same configuration produce
two rows rather than one summed row, which is what keeps the kernel free of a per-sample hash
table. Every consumer sums over the connection axis, so results are unchanged; only the row
count differs. Use `connected` when you want the summed, deduplicated matrix elements as a
dictionary.

## Unpacking for a network

Packed integers suit the kernel — compact, hashable, cheap to permute. A network wants an
explicit numeric array, and `configurations` is the boundary. The degree-of-freedom axis comes
first, so each configuration is contiguous:

```@example cc
size(configurations(model.dof, res.configs, model.nsites))
```

## Full documentation

`ConnectedConfigs` is registered as a package in its own right and has a complete site:

- [Connected configurations](https://cevenkadir.github.io/NeuralQuantumStates.jl/ConnectedConfigs/dev/manual/connected_states/)
  — compiling, batch shapes, padding, and reusing buffers with `connected_padded!`
- [Local operators](https://cevenkadir.github.io/NeuralQuantumStates.jl/ConnectedConfigs/dev/manual/local_operators/)
  — the digit ordering, and why `PAULI_Z` is the wrong constant to build a Hamiltonian from
- [Symmetry sectors](https://cevenkadir.github.io/NeuralQuantumStates.jl/ConnectedConfigs/dev/manual/symmetry_sectors/)
  — folding onto orbit representatives, characters and norm factors
- [Backends](https://cevenkadir.github.io/NeuralQuantumStates.jl/ConnectedConfigs/dev/manual/backends/)
  — swapping SymBasis or OperatorAlgebra for another library
- [Exact diagonalization](https://cevenkadir.github.io/NeuralQuantumStates.jl/ConnectedConfigs/dev/examples/exact_diagonalization/)
- [API reference](https://cevenkadir.github.io/NeuralQuantumStates.jl/ConnectedConfigs/dev/api/connected/)
