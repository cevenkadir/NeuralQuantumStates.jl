```@meta
CurrentModule = ConnectedBasisConfigurations
```

# Connected configurations

Given an operator ``\hat{H}`` and a configuration ``\vert s \rangle``, the *connected*
configurations are those ``\vert s' \rangle`` for which
``\langle s' \vert \hat{H} \vert s \rangle`` is non-zero. Everything in this package is built
to produce them for a whole batch at once, because that is what a variational Monte Carlo local
energy needs:

```math
E_{\mathrm{loc}}(s) = \sum_{s'} \langle s' \vert \hat{H} \vert s \rangle \,
                      \frac{\psi(s')}{\psi(s)}
```

## Compile first

[`compile`](@ref) flattens an operator once into the form the kernel runs on. It resolves
Jordan–Wigner strings, distributes sums nested inside products, multiplies out factors acting
on the same site, promotes every matrix element to one concrete element type, and
column-compresses each local matrix so the kernel visits only the entries it needs.

```@example connected
using ConnectedBasisConfigurations, OperatorAlgebra, SymBasis

spec, nsites = Spin(1 // 2), 6
ops = local_operators(spec)

H = OpSum(vcat(
    [Op(2 .* ops.sz, i) * Op(2 .* ops.sz, mod1(i + 1, nsites)) for i in 1:nsites],
    [0.7 * Op(2 .* ops.sx, i) for i in 1:nsites],
))

compiled = compile(H)
```

The result reports how the terms split. Terms that cannot change a configuration at all —
``\sigma^z \sigma^z``, ``n_i n_j``, an uncancelled Jordan–Wigner tail — are separated out and
evaluated by a branch-free loop that never touches a state, which in a typical lattice
Hamiltonian is most of them.

Passing the operator itself to [`connected_padded`](@ref) works too; it just compiles on every
call. Over a large batch that is amortized, but over single-configuration queries it dominates
the runtime, so hoist it whenever the operator outlives the call.

## The batched call

```@example connected
states = basis(dof_object(spec), nsites).states
res = connected_padded(compiled, states)

size(res.configs), size(res.mels), size(res.counts)
```

- `configs` holds the connected configurations as packed integers.
- `mels` holds the matrix elements.
- `counts` says how many entries of each sample are real rather than padding.

The connection axis comes **first**. NetKet puts it last, but Julia is column-major, and this
is the ordering that keeps one sample's connections contiguous in memory — which is what the
inner loop of a local energy walks.

## Batches of any shape

`states` may be an array with any leading dimensions. The connection axis is prepended and the
batch shape is preserved, so a `(steps, chains)` block of Monte Carlo samples goes in directly:

```@example connected
block = reshape(states[1:24], 4, 6)
r = connected_padded(compiled, block)

size(r.configs), size(r.counts)
```

## Padding

Different configurations have different numbers of connections, so the connection axis is
padded to a common length. Padded slots carry a matrix element of **zero** and repeat the
sample itself.

Both halves of that convention earn their place. A zero matrix element makes the local energy
correct with no masking at all:

```julia
E_loc = sum(res.mels[:, b] .* exp.(logψ.(res.configs[:, b]) .- logψ(s)))
```

and repeating the sample means every padded row is still a *valid* configuration, so the
wavefunction can be evaluated on the whole array without special-casing. Marking padding
`missing` instead would turn `mels` into a `Union{T,Missing}` array and put a branch in the
innermost loop.

```@example connected
b = argmin(res.counts)
res.counts[b], res.mels[(res.counts[b]+1):end, b]
```

## Repeated configurations are not merged

Two terms reaching the same configuration produce **two** rows rather than one summed row. This
follows NetKet, and it is what keeps the kernel free of a per-sample hash table.

Every consumer sums over the connection axis, so the result is unchanged; only the row count
differs. When you want the summed, deduplicated matrix elements — to build a matrix column, or
to count how many distinct configurations an operator reaches — use [`connected`](@ref), which
returns a dictionary:

```@example connected
d = connected(compiled, first(states))
length(d), sum(abs, values(d))
```

Entries are ordered by term: the diagonal first when it is non-zero, then each off-diagonal
term's contribution. That ordering is deterministic, which hash iteration order is not.

## Reusing buffers

[`connected_padded!`](@ref) writes into arrays you own, so an optimization loop stops
allocating a fresh output every step. The buffers are sized by [`max_conn_size`](@ref), a bound
known from the operator alone without looking at a single sample:

```@example connected
height = max_conn_size(compiled)
configs = Matrix{eltype(states)}(undef, height, length(states))
mels = Matrix{Float64}(undef, height, length(states))
counts = Vector{Int}(undef, length(states))

connected_padded!(configs, mels, counts, compiled, states)
height, maximum(counts)
```

Unlike the allocating form, the connection axis is not trimmed to the batch's actual maximum:
the buffers keep their full height, with every slot above a sample's count padded inert. That
is what lets the same buffers serve calls whose connection counts differ.

## Unpacking for a neural network

Packed integers are the right representation for the kernel — compact, hashable, cheap to
permute. A network wants an explicit numeric array, and [`configurations`](@ref) is the
boundary between the two. The degree-of-freedom axis comes first, so each configuration is
contiguous:

```@example connected
size(configurations(spec, res.configs, nsites))
```
