"""
    ConnectedBasisConfigurations

The batched local-energy kernel: given an operator `Ĥ` and a *batch* of basis configurations
`|s⟩`, produce every connected basis configuration `|s′⟩` together with its matrix element
`⟨s′|Ĥ|s⟩`.

This is the Julia counterpart of NetKet's `operator.get_conn_padded`, and the one piece of the
operator layer OperatorAlgebra.jl does not provide: its `apply` consumes a vector or a `Dict`
of states, not a batch of samples. Building `Ĥ`, normal ordering, fermionic signs and sparse
conversion all belong upstream and are not duplicated here.

The kernel is the hot loop of variational Monte Carlo. It depends only on SymBasis and
OperatorAlgebra — no Lux, autodiff or GPU — and reaches even those through the generic functions
in `interface.jl`, so another basis or operator library plugs in by adding methods.

# Entry points
- [`compile`](@ref) — flatten an operator once, ahead of the loop
- [`connected_padded`](@ref) — the batched kernel, with and without a symmetry-reduced basis
- [`connected_padded!`](@ref) — the same, into buffers the caller owns
- [`connected`](@ref) — a single configuration, as a `Dict`
- [`configurations`](@ref) / [`packed`](@ref) — convert between packed `BaseInt` states and
  arrays of physical local values, and [`configurations!`](@ref) for the same on a device
- [`local_operators`](@ref) — single-site matrices in SymBasis's digit ordering

# Conventions

Three, all easy to get silently wrong:

1. **Digit position `i` is lattice site `i`**, counting from the least significant digit, which
   is what `SymBasis.read`/`eachdigit` use. An `Op`'s site identifier is used directly as that
   position, so sites must be integers `1:nsites`.
2. **Digit `d` means local value `local_values(spec)[d+1]`** — `(-s, …, +s)` for a spin,
   `0:max_occupancy` for a boson. This is *not* the convention behind OperatorAlgebra's `PAULI_Z`,
   which puts `+1` on the first basis state; see [`local_operators`](@ref).
3. **Two terms reaching the same configuration produce two entries**, not one summed entry, in
   the padded output. Consumers sum over the connection axis, so the result is unchanged; only
   the row count differs. [`connected`](@ref) gives the summed, deduplicated form.

!!! warning "OperatorAlgebra numbers digits the other way"
    Convention 1 is SymBasis's. OperatorAlgebra's `sparse` and `apply` build their Kronecker
    products with site 1 as the **most** significant digit, so their matrix and this package's
    are related by a permutation of the basis rather than being equal. Chain models symmetric
    under site reversal hide the difference completely.

# Example
```julia
using ConnectedBasisConfigurations, OperatorAlgebra, SymBasis

spec, nsites = Spin(1 // 2), 4
ops = local_operators(spec)
H = OpSum([Op(ops.sz, i) * Op(ops.sz, mod1(i + 1, nsites)) for i in 1:nsites])

compiled = compile(H)                    # once, outside the loop
states = basis(dof_object(spec), nsites).states
res = connected_padded(compiled, states)
res.mels        # (max_connections, batch)
```
"""
module ConnectedBasisConfigurations

using LinearAlgebra: I

using OperatorAlgebra
using OperatorAlgebra: AbstractOp, Op, OpChain, OpSum, basis_info

using SymBasis
using SymBasis: Boson, Spin, representative
using SymBasis.DigitBase: BaseInt, read, write

include("interface.jl")
include("backends/symbasis.jl")
include("backends/operatoralgebra.jl")
include("local_operators.jl")
include("compile.jl")
include("connected.jl")
include("unpack.jl")
include("flat.jl")

export compile, CompiledOperator, CompiledSector, max_conn_size
export flatten, FlatOperator, to_backend
export connected, connected_padded, connected_padded!
export configurations, configurations!, packed
export local_operators, local_values, local_dimension
export expand_terms, amplitude_type, read_digit, write_digit
export fold_state, state_lookup, lookup_index, orbit_norms

end # module ConnectedBasisConfigurations
