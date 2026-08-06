"""
    local_operators(spec) -> NamedTuple

Single-site operator matrices for a SymBasis degree-of-freedom specification, in the **digit
ordering SymBasis actually uses**.

# Why this exists rather than OperatorAlgebra's constants

OperatorAlgebra ships `PAULI_X`, `PAULI_Z`, `RAISE`, and friends, but they are 2×2 only, so
they cannot describe a boson with `max_occupancy > 1` or a spin above 1/2.

More importantly, they use the **opposite spin convention**. `PAULI_Z = [1 0; 0 -1]` maps the
first basis state to `+1`, i.e. digit `0` is spin *up*. SymBasis orders local states as
`dofo.ldof = (-s, …, +s)`, so digit `0` is spin *down*. Building a Hamiltonian from
`PAULI_Z` and evaluating it on SymBasis states therefore flips the sign of every ``S^z``
— which is invisible in the ground-state energy of an unbiased model and glaring the moment a
longitudinal field is switched on.

The operators returned here are constructed directly from `spec`, so they agree with the digit
ordering by construction. Matrices follow the usual convention `mat[out, in]`.

# Returned fields
- `Spin`: `sx`, `sy`, `sz`, `sp` (``S^+``), `sm` (``S^-``), `id`
- `Boson`: `a`, `adag`, `n`, `id`

# Example
```julia
ops = local_operators(Boson(3))
ops.n                       # diag(0, 1, 2, 3)
Op(ops.adag, 1) * Op(ops.a, 2)   # hopping from site 2 to site 1
```
"""
function local_operators end

"""
    local_dimension(spec) -> Int

Number of local states of `spec`, i.e. the base of the packed integer representation.
"""
function local_dimension end

local_dimension(spec::Spin) = Int(2 * spec.s + 1)
local_dimension(spec::Boson) = Int(spec.max_occupancy) + 1

"""
    local_values(spec) -> Vector

The physical value of each local state, indexed by digit + 1 — the same ordering as
`dof_object(spec).ldof`. For a spin these are the magnetic quantum numbers `-s:s`; for a boson
the occupation numbers `0:max_occupancy`.
"""
function local_values end

local_values(spec::Spin) = collect(-spec.s:spec.s)
local_values(spec::Boson) = collect(0:Int(spec.max_occupancy))

function local_operators(spec::Spin)
    s = spec.s
    d = local_dimension(spec)
    m = local_values(spec)          # m[i] is the magnetic quantum number of digit i-1

    sz = zeros(Float64, d, d)
    for i in 1:d
        sz[i, i] = m[i]
    end

    # S⁺|m⟩ = sqrt(s(s+1) - m(m+1)) |m+1⟩, which raises the digit by one.
    sp = zeros(Float64, d, d)
    for i in 1:(d-1)
        sp[i+1, i] = sqrt(s * (s + 1) - m[i] * (m[i] + 1))
    end
    sm = collect(transpose(sp))

    sx = (sp + sm) / 2
    sy = (sp - sm) / (2im)

    return (; sx=sx, sy=sy, sz=sz, sp=sp, sm=sm, id=Matrix{Float64}(I, d, d))
end

function local_operators(spec::Boson)
    d = local_dimension(spec)

    # a|n⟩ = sqrt(n)|n-1⟩ lowers the digit by one; a†|n⟩ = sqrt(n+1)|n+1⟩ raises it.
    a = zeros(Float64, d, d)
    for n in 1:(d-1)
        a[n, n+1] = sqrt(n)
    end
    adag = collect(transpose(a))

    num = zeros(Float64, d, d)
    for n in 0:(d-1)
        num[n+1, n+1] = n
    end

    return (; a=a, adag=adag, n=num, id=Matrix{Float64}(I, d, d))
end
