"""
    ConnectedConfigs

The batched local-energy kernel: given an operator `Ĥ` and a *batch* of basis configurations
`|s⟩`, produce every connected configuration `|s′⟩` together with its matrix element
`⟨s′|Ĥ|s⟩`.

This is the Julia counterpart of NetKet's `operator.get_conn_padded`, and it is the one piece
of the operator layer that OperatorAlgebra.jl does not already provide: its `apply` consumes a
vector or a `Dict` of states, not a batch of samples. Everything else — building `Ĥ`, normal
ordering, fermionic signs, sparse and dense conversion — belongs upstream in OperatorAlgebra,
and this package deliberately does not duplicate any of it.

The kernel is the hot loop of variational Monte Carlo: it runs once per sample per optimization
step. It is kept in its own package, depending only on SymBasis and OperatorAlgebra, so that it
carries **no** Lux, autodiff, or GPU dependency and stays usable by any VMC code — neural or
otherwise.

# Entry points
- [`connected_padded`](@ref) — the batched kernel, with and without a symmetry-reduced basis
- [`connected`](@ref) — a single configuration, as a `Dict`
- [`configurations`](@ref) / [`packed`](@ref) — convert between packed `BaseInt` states and
  arrays of physical local values
- [`local_operators`](@ref) — single-site matrices in SymBasis's digit ordering

# Conventions

Two, both easy to get silently wrong:

1. **Digit position `i` is lattice site `i`**, counting from the least significant digit, which
   is what `SymBasis.read`/`eachdigit` use. An `Op`'s site identifier is used directly as that
   position, so sites must be integers `1:nsites`.
2. **Digit `d` means local value `local_values(spec)[d+1]`** — `(-s, …, +s)` for a spin,
   `0:max_occupancy` for a boson. This is *not* the convention behind OperatorAlgebra's `PAULI_Z`,
   which puts `+1` on the first basis state; see [`local_operators`](@ref).

# Example
```julia
using ConnectedConfigs, OperatorAlgebra, SymBasis

spec, nsites = Spin(1 // 2), 4
ops = local_operators(spec)
H = OpSum([Op(ops.sz, i) * Op(ops.sz, mod1(i + 1, nsites)) for i in 1:nsites])

states = basis(dof_object(spec), nsites).states
res = connected_padded(H, states)
res.mels        # (max_connections, batch)
```
"""
module ConnectedConfigs

using LinearAlgebra: I

using OperatorAlgebra
using OperatorAlgebra: AbstractOp, Op, OpChain, OpSum, basis_info

using SymBasis
using SymBasis: Boson, Spin, representative
using SymBasis.DigitBase: BaseInt, read, write

include("local_operators.jl")
include("connected.jl")
include("unpack.jl")

export connected, connected_padded
export configurations, packed
export local_operators, local_values, local_dimension

end # module ConnectedConfigs
