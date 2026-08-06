"""
    AbstractModelSpec

A specification of a predefined many-body model.

Like a lattice specification, a model spec *describes* a Hamiltonian — its lattice, its local
degrees of freedom, and its couplings — without building anything. Pass one to
`build` to get a [`Model`](@ref):

```julia
lat = build(Hypercube([8]; periodic=true))
model = build(TransverseFieldIsing(lat; J=1.0, h_x=1.0))
```

Concrete specs are [`TransverseFieldIsing`](@ref) and [`ExtendedBoseHubbard`](@ref).

This replaces the pre-split `Operators.build(:TransverseFieldIsing, hilbert, lattice; ...)`.
The Hamiltonian is now an ordinary OperatorAlgebra `OpSum` rather than a bespoke type, so it
composes with everything that package offers — `sparse`, `LinearMap`, `commutator`, ITensor
conversion — none of which the old hand-written operators supported.
"""
abstract type AbstractModelSpec end

"""
    Model

A Hamiltonian together with the space it acts on.

# Fields
- `hamiltonian::OpSum`: The Hamiltonian, as an OperatorAlgebra operator.
- `dof`: The SymBasis degree-of-freedom specification (`Spin`, `Boson`, ...).
- `nsites::Int`: Number of sites.
- `lattice`: The lattice it was built on.

`dof` and `nsites` are carried alongside the operator because the operator alone does not
determine them: an `OpSum` knows which sites it touches, but not that a site is a spin-1/2
rather than a spin-1, nor that sites it happens not to act on still exist.
"""
struct Model{O,S,L}
    hamiltonian::O
    dof::S
    nsites::Int
    lattice::L
end

"""
    basis(model::Model, symmetry...) -> SymBasis.Basis

The computational basis of `model`, optionally restricted to a symmetry sector.

```julia
basis(model)                                        # full space
basis(model, sym(Translational(0, lat), dof_object(model.dof)))
```
"""
SymBasis.Bases.basis(model::Model) = basis(dof_object(model.dof), model.nsites)
SymBasis.Bases.basis(model::Model, symmetry) =
    basis(dof_object(model.dof), model.nsites, symmetry)

"""
    TransverseFieldIsing{L,T} <: AbstractModelSpec

``H = J \\sum_{\\langle ij \\rangle} \\sigma^z_i \\sigma^z_j
     + h_x \\sum_i \\sigma^x_i + h_z \\sum_i \\sigma^z_i``

# Constructor
    TransverseFieldIsing(lattice; J=1.0, h_x=1.0, h_z=0.0)

!!! note "Pauli matrices, not spin operators"
    The couplings multiply ``\\sigma`` rather than ``S = \\sigma/2``, matching the pre-split
    package. `σᶻ` is `diag(-1, +1)`: digit 0 is spin **down**, following SymBasis's
    `ldof = (-1//2, 1//2)`. This is the negative of OperatorAlgebra's `PAULI_Z` constant — see
    `ConnectedConfigs.local_operators` for why that distinction matters.
"""
struct TransverseFieldIsing{L,T<:Real} <: AbstractModelSpec
    lattice::L
    J::T
    h_x::T
    h_z::T

    function TransverseFieldIsing(lattice::L; J=1.0, h_x=1.0, h_z=0.0) where {L}
        T = promote_type(typeof(J), typeof(h_x), typeof(h_z), Float64)
        return new{L,T}(lattice, T(J), T(h_x), T(h_z))
    end
end

function build(spec::TransverseFieldIsing)
    dof = Spin(1 // 2)
    ops = local_operators(dof)
    σz, σx = 2 .* ops.sz, 2 .* ops.sx
    nsites = n_sites(spec.lattice)

    terms = AbstractOp[]
    for (i, j) in bonds(spec.lattice)
        iszero(spec.J) || push!(terms, spec.J * (Op(σz, i) * Op(σz, j)))
    end
    for i in 1:nsites
        iszero(spec.h_z) || push!(terms, spec.h_z * Op(σz, i))
        iszero(spec.h_x) || push!(terms, spec.h_x * Op(σx, i))
    end
    return Model(OpSum(terms), dof, nsites, spec.lattice)
end

"""
    ExtendedBoseHubbard{L,T} <: AbstractModelSpec

``H = -J \\sum_{\\langle ij \\rangle} (b^\\dagger_i b_j + \\mathrm{h.c.})
     + \\frac{U}{2} \\sum_i n_i (n_i - 1)
     + V \\sum_{\\langle ij \\rangle} n_i n_j - \\mu \\sum_i n_i``

# Constructor
    ExtendedBoseHubbard(lattice, n_max; J=1.0, U=1.0, V=0.0, μ=0.0)

`n_max` is the maximum occupation per site, which truncates the bosonic ladder.

!!! note "Sign of the chemical potential"
    ``\\mu`` enters with a **minus** sign, matching what the pre-split package computed. Its
    own docstring advertised a plus sign, which its code did not implement.
"""
struct ExtendedBoseHubbard{L,T<:Real} <: AbstractModelSpec
    lattice::L
    n_max::Int
    J::T
    U::T
    V::T
    μ::T

    function ExtendedBoseHubbard(lattice::L, n_max::Integer; J=1.0, U=1.0, V=0.0, μ=0.0) where {L}
        n_max > 0 || throw(ArgumentError("n_max must be positive"))
        T = promote_type(typeof(J), typeof(U), typeof(V), typeof(μ), Float64)
        return new{L,T}(lattice, Int(n_max), T(J), T(U), T(V), T(μ))
    end
end

function build(spec::ExtendedBoseHubbard)
    dof = Boson(spec.n_max)
    ops = local_operators(dof)
    a, adag, n = ops.a, ops.adag, ops.n
    nn = n * n - n                      # n(n-1), for the on-site interaction
    nsites = n_sites(spec.lattice)

    terms = AbstractOp[]
    for (i, j) in bonds(spec.lattice)
        if !iszero(spec.J)
            push!(terms, (-spec.J) * (Op(adag, i) * Op(a, j)))
            push!(terms, (-spec.J) * (Op(adag, j) * Op(a, i)))
        end
        iszero(spec.V) || push!(terms, spec.V * (Op(n, i) * Op(n, j)))
    end
    for i in 1:nsites
        iszero(spec.U) || push!(terms, (spec.U / 2) * Op(nn, i))
        iszero(spec.μ) || push!(terms, (-spec.μ) * Op(n, i))
    end
    return Model(OpSum(terms), dof, nsites, spec.lattice)
end
