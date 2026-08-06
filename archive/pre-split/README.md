# Retired pre-split modules

`Hilberts`, `Operators`, and `Extras`, as they stood immediately before being removed from
`src/`. They are kept here for one reason: **this is the code the reference data in
`test/golden/` was generated from**, so it is the provenance of the equivalence proof that
licensed deleting them.

They also carried uncommitted working-tree changes at the time of removal, which existed
nowhere else — the committed versions are in git history, those edits were not.

| Module | Replaced by |
|---|---|
| `hilberts/` | [SymBasis.jl](https://github.com/cevenkadir/SymBasis.jl) — packed integer states, a far wider range of symmetries |
| `operators/` | [OperatorAlgebra.jl](https://github.com/h-mnzlr/OperatorAlgebra.jl) for the algebra, `lib/ConnectedConfigs` for the batched local-energy kernel |
| `extras/` | Nothing — `get_array_type_params` existed only to prop up the `T_Array` type-parameter gymnastics of the other two |

**This is not maintained code.** It does not load against the current package and is not
tested. Once the split has settled and the golden data has served its purpose, this directory
can be deleted outright; git history keeps it.
