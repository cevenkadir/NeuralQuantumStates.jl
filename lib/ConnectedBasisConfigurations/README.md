<div align="center">

# ConnectedBasisConfigurations.jl

*The batched connected-basis-configuration kernel for quantum lattice models, in Julia*

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://cevenkadir.github.io/NeuralQuantumStates.jl/ConnectedBasisConfigurations/stable/) [![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://cevenkadir.github.io/NeuralQuantumStates.jl/ConnectedBasisConfigurations/dev/) [![Build Status](https://github.com/cevenkadir/NeuralQuantumStates.jl/actions/workflows/CI-ConnectedBasisConfigurations.yml/badge.svg?branch=main)](https://github.com/cevenkadir/NeuralQuantumStates.jl/actions/workflows/CI-ConnectedBasisConfigurations.yml?query=branch%3Amain) [![Coverage](https://codecov.io/gh/cevenkadir/NeuralQuantumStates.jl/branch/main/graph/badge.svg?flag=ConnectedBasisConfigurations)](https://codecov.io/gh/cevenkadir/NeuralQuantumStates.jl/flags?flag=ConnectedBasisConfigurations) [![Downloads](https://img.shields.io/badge/dynamic/json?url=http%3A%2F%2Fjuliapkgstats.com%2Fapi%2Fv1%2Fmonthly_downloads%2FConnectedBasisConfigurations&query=total_requests&suffix=%2Fmonth&label=Downloads)](https://juliapkgstats.com/pkg/ConnectedBasisConfigurations) [![Total Downloads](https://img.shields.io/badge/dynamic/json?url=http%3A%2F%2Fjuliapkgstats.com%2Fapi%2Fv1%2Ftotal_downloads%2FConnectedBasisConfigurations&query=total_requests&label=Total%20Downloads)](https://juliapkgstats.com/pkg/ConnectedBasisConfigurations)
</div>

Variational Monte Carlo spends nearly all its time answering one question, several thousand times per optimization step: given a basis configuration `|s⟩` and a Hamiltonian `Ĥ`, which basis configurations `|s′⟩` does `Ĥ` connect it to, and with what matrix element `⟨s′|Ĥ|s⟩`? This package answers it for a whole **batch** at once, returning padded arrays rather than a dictionary per sample — the Julia counterpart of NetKet's `get_conn_padded`. The speed comes from a `compile` step that flattens the operator once, resolving Jordan–Wigner strings, distributing sums over products and column-compressing every local matrix, after which the per-sample loop is pure array indexing with no hash tables and no allocation. Its only dependencies are [SymBasis.jl](https://github.com/cevenkadir/SymBasis.jl) and [OperatorAlgebra.jl](https://github.com/cevenkadir/OperatorAlgebra.jl), and even those are replaceable — the kernel talks to them through a small documented interface.

## Key features
- **Compiled operators**: `compile` hoists every configuration-independent decision out of the hot loop. For a 16-site transverse-field Ising chain over a batch of 1024 configurations, that is the difference between 37 ms and 0.28 ms per call.
- **Batches of any shape**: `connected_padded` accepts an array of configurations with any leading dimensions and returns results with the connection axis prepended, so a `(chains, steps)` block of Monte Carlo samples needs no reshaping.
- **Padding you do not have to mask**: unused slots carry a zero matrix element and repeat the sample itself, so a local energy is a plain sum over the connection axis and the wavefunction can be evaluated on every row, real or padded.
- **An in-place form**: `connected_padded!` writes into caller-owned buffers sized by `max_conn_size`, so an optimization loop stops reallocating its output every step.
- **Symmetry sectors**: pass a symmetry-reduced basis and matrix elements come back folded onto orbit representatives, rescaled by the character and the orbit norms — the ingredients of a sector-resolved exact diagonalization.
- **Swappable backends**: SymBasis and OperatorAlgebra are the defaults, not the requirement. A different basis or operator library plugs in by implementing `expand_terms`, `amplitude_type`, `read_digit` and `write_digit` — ordinary methods on ordinary generic functions, with nothing to load.
- **Almost no dependencies**: no autodiff, no Lux, no GPU stack, so the kernel stays usable by any Monte Carlo or exact-diagonalization code, neural or otherwise. The test suite asserts it.

## Installation
**Requirements**: Julia 1.11 or later.

To install the latest stable version, use the Julia package manager. Either use the Julia REPL package mode (by pressing `]`):
```julia
pkg> add ConnectedBasisConfigurations
```
or open the Julia REPL and run the following command:
```julia
julia> import Pkg; Pkg.add("ConnectedBasisConfigurations")
```

## Documentation
For detailed information on using this package, check out the [stable documentation](https://cevenkadir.github.io/NeuralQuantumStates.jl/ConnectedBasisConfigurations/stable/).

## Manual
- [Connected basis configurations](https://cevenkadir.github.io/NeuralQuantumStates.jl/ConnectedBasisConfigurations/stable/manual/connected_states/)
- [Local operators](https://cevenkadir.github.io/NeuralQuantumStates.jl/ConnectedBasisConfigurations/stable/manual/local_operators/)
- [Symmetry sectors](https://cevenkadir.github.io/NeuralQuantumStates.jl/ConnectedBasisConfigurations/stable/manual/symmetry_sectors/)
- [Backends](https://cevenkadir.github.io/NeuralQuantumStates.jl/ConnectedBasisConfigurations/stable/manual/backends/)

## Examples
- [Exact diagonalization](https://cevenkadir.github.io/NeuralQuantumStates.jl/ConnectedBasisConfigurations/stable/examples/exact_diagonalization/)

## Quick example
Build a transverse-field Ising chain, compile it once, and ask what a batch of configurations connects to:
```julia
julia> using ConnectedBasisConfigurations, OperatorAlgebra, SymBasis

julia> spec, nsites = Spin(1 // 2), 4;

julia> ops = local_operators(spec);   # matrices in SymBasis's digit ordering

julia> H = OpSum(vcat(
           [Op(2 .* ops.sz, i) * Op(2 .* ops.sz, mod1(i + 1, nsites)) for i in 1:nsites],
           [Op(2 .* ops.sx, i) for i in 1:nsites],
       ));

julia> compiled = compile(H)
CompiledOperator{Float64}(4 diagonal + 4 off-diagonal terms, max_conn 5)

julia> states = basis(dof_object(spec), nsites).states;

julia> res = connected_padded(compiled, states);

julia> size(res.configs)         # (connections, batch)
(5, 16)

julia> res.mels[:, 1]            # the diagonal first, then one flip per site
5-element Vector{Float64}:
 4.0
 1.0
 1.0
 1.0
 1.0
```

The batch axis follows the input, so nothing needs reshaping on the way in or out:
```julia
julia> block = reshape(states, 4, 4);      # e.g. (steps, chains) of Monte Carlo samples

julia> r = connected_padded(compiled, block);

julia> size(r.configs), size(r.counts)
((5, 4, 4), (4, 4))
```

Turn packed configurations into the numeric array a neural network wants, degree-of-freedom axis first:
```julia
julia> size(configurations(spec, res.configs, nsites))
(4, 5, 16)
```

## Scope
Building operators — algebra, normal ordering, fermionic signs, sparse and dense conversion — belongs upstream in OperatorAlgebra.jl, and this package deliberately duplicates none of it. What it adds is the one thing an operator library has no reason to provide: evaluation over a *batch of samples* rather than over a state vector. Sampling lives in NQSSamplers.jl and wavefunctions in NQSCore.jl; neither is a dependency, and neither needs to be for this package to be useful on its own.

## Part of a larger ecosystem
ConnectedBasisConfigurations.jl is developed in the [NeuralQuantumStates.jl](https://github.com/cevenkadir/NeuralQuantumStates.jl) monorepo, and is usable entirely on its own — it has no idea the neural-network stack exists, and its test suite asserts as much.

## Important notice
This project is still under active development. While it includes an extensive test suite and is developed with high scientific rigor, you should always benchmark your own code. Please report any issues you encounter via the [GitHub issue tracker](https://github.com/cevenkadir/NeuralQuantumStates.jl/issues/new).

## Citation
If you use this package in your work, we would appreciate the following reference as in [CITATION.bib](https://github.com/cevenkadir/NeuralQuantumStates.jl/blob/main/CITATION.bib).
