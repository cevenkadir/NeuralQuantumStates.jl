# Golden reference data

These files were generated from the **pre-split** `Operators.connected_basis_configs`
implementation (`src/operators/predefined_operators/{ising,bosehubbard}.jl`, deleted as part of
the package split) at commit `5ca3f22` plus the uncommitted working-tree changes. They exist so
the new `ConnectedConfigs` package can be proved to reproduce the old behaviour before the old
code is removed.

They are Julia source files containing literal arrays — `include` them and compare. No
serialization format, no version fragility, and they diff readably.

## Contents

| File | Model |
|---|---|
| `tfi_chain8.jl` | Transverse-field Ising, 8-site periodic chain, `J=1.0, h_x=1.0, h_z=1.0`, spin-1/2, unconstrained |
| `bhm_chain16.jl` | Extended Bose-Hubbard, 16-site periodic chain, `J=1.0, U=1.0, V=1.0, μ=0.0`, Fock `n_max=5`, total number fixed to 5 |

Each defines, for the single-sample and batched-array call paths, `*_configs` and `*_mels`,
plus the `*_sample_*` inputs that produced them, and `*_vov_*` for the vector-of-vectors path.

## Two conventions are baked in — do not "fix" them silently

The pre-split code was **inconsistent between the two models**, and the golden data records
that faithfully:

- **Ising array method** asserts `size(samples, ndims(samples)) == N` — the DoF axis is **last** —
  and returns `configs` shaped `(M_max, batch..., N)`, `mels` shaped `(M_max, batch...)`.
- **Bose-Hubbard array method** asserts `size(samples, 1) == N` — the DoF axis is **first** —
  and returns `configs` shaped `(N, M_max, batch...)`, `mels` shaped `(1, M_max, batch...)`
  (note the extra leading singleton dimension that the Ising path does not have).

`ConnectedConfigs` should settle on **one** convention. When it does, the comparison against
this data must transpose/reshape explicitly and say so, rather than the test being quietly
written to whatever the new code happens to produce. Picking a convention is a deliberate API
decision; silently inheriting one is not.

## Verified properties

Both datasets were checked against analytically known values at capture time, not merely
recorded:

- **TFI**: column 1 is the unflipped configuration; the remaining 8 differ in exactly one site
  each with matrix element `h_x = 1.0`; the diagonal element equals
  `4J Σ⟨ij⟩ zᵢzⱼ + 2h_z Σᵢ zᵢ = -4.0` for the stored sample.
- **BHM**: 11 connected configurations (1 diagonal + 10 hops) for 5 particles on a 16-site ring;
  total particle number is conserved in every connected configuration; off-diagonal elements are
  `-J√(nⱼ)√(nᵢ+1)`, giving the `-√2 ≈ -1.4142` entries where a particle hops onto an
  already-occupied site.

## What this data cannot check

`tfi_sample_vec` has **total magnetization zero**, so the `h_z Σ σᶻ` term contributes nothing
to it. Reproducing this data therefore does *not* pin down the sign convention of `σᶻ`:
substituting OperatorAlgebra's `PAULI_Z` (which puts `+1` on the first basis state) for
`2Sᶻ` (which follows SymBasis's `ldof = (-1//2, 1//2)`, putting `-1` there) passes every
comparison against these arrays. This was verified by mutation, not assumed.

`lib/ConnectedConfigs/test/golden.jl` therefore carries a separate testset using **magnetized**
configurations, where the two conventions differ by the entire field term. Do not treat
agreement with the arrays here as evidence that a model's spin convention is correct.

## Regenerating

The capture script is not kept in the repo because it depends on deleted modules. To regenerate,
check out a commit that still has `src/operators/`, and run the script recorded in the split
commit message. In normal use these files are **immutable inputs** — regenerating them to make a
failing test pass defeats their entire purpose.
