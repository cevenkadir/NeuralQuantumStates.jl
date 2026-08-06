# Salvaged scratch files

These were untracked working files in the repo root and `test/` before the package split. They
are preserved here because they were never committed and contain the only implementations of
work that has not yet landed in a package. **Do not treat anything here as maintained code** —
it does not run against the post-split API and is kept purely as a source to port from.

| File | Contains | Ports into |
|---|---|---|
| `JACOBIAN.jl` | `chunked_jacobian` for `ComponentArray` parameters, including the real/imaginary splitting path for complex parameters, and Enzyme `chunk=Val(n)` usage | `NQSCore.log_derivatives` (the `O_k` matrix) |
| `JACOBIAN copy.jl`, `JACOBIAN_with_Reactant.jl` | Variants of the above; the Reactant one has the `@compile` wiring | `NQSCore` Reactant extension |
| `runtest_network.jl`, `runtest_network copy*.jl` | Lux + Reactant + Enzyme wiring, `DataLoader` batching, per-sample VJP accumulation via `fmap`, the `ComplexF32[1 im] * model(x)` log-ψ trick | `NQSCore` extensions, `NQSAnsatze` log-ψ layer |
| `runtest_ising.jl` | Concrete `connected_basis_configs` calls covering **all four input ranks** (single vector, vector-of-vectors, matrix, 3-D reshape) with fixed spin configurations | Golden-test inputs for `ConnectedConfigs` — see `test/golden/` |
| `batch_test.jl` | `Lux.batched_jacobian` with `AutoForwardDiff`/`AutoZygote` | `NQSCore` backend comparison |
| `Project copy.toml` | An alternative dependency set that additionally pinned Metal, Zygote and MLUtils | Reference only; the split moves these to weakdeps |
| `kkk.ipynb` | Empty (0 bytes) | — |

Once the corresponding package has landed the functionality and its tests pass, the matching
row here can be deleted.
