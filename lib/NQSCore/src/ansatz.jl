"""
    LogStateVector(dof, nsites, basis) <: AbstractAnsatz

The exact ansatz: one complex parameter per basis state, holding that state's log-amplitude
directly.

This is NetKet's `models.LogStateVector`, and it exists for the same reason. It has as many
parameters as the Hilbert space has dimensions, so it is useless for anything large — but it can
represent *any* state exactly, which makes it the reference every approximate ansatz is measured
against. Optimizing it must reproduce exact diagonalization; if it does not, the bug is in the
optimization machinery rather than in the ansatz.

It is also the one ansatz whose log-derivatives are known in closed form: since
`log ψ(s) = θ_{i(s)}`, the matrix `O` is a one-hot indicator of which basis state each sample
is. That makes it the natural test of [`log_derivatives`](@ref) itself.

# Fields
- `dof`: The SymBasis degree-of-freedom specification.
- `nsites::Int`: Number of sites.
- `basis`: The basis whose states the parameters are indexed by.
"""
struct LogStateVector{D,B} <: AbstractAnsatz
    dof::D
    nsites::Int
    basis::B
    index::Dict{Any,Int}

    function LogStateVector(dof::D, nsites::Integer, basis::B) where {D,B}
        index = Dict{Any,Int}(s => i for (i, s) in pairs(basis.states))
        return new{D,B}(dof, Int(nsites), basis, index)
    end
end

"""Number of parameters of a `LogStateVector`: the dimension of its basis."""
n_parameters(a::LogStateVector) = length(a.basis.states)

function log_amplitude(a::LogStateVector, θ, x::AbstractMatrix)
    out = similar(θ, size(x, 2))
    for j in axes(x, 2)
        s = ConnectedBasisConfigurations.packed(a.dof, @view x[:, j])
        i = get(a.index, s, 0)
        i == 0 && throw(ArgumentError(
            "configuration $(collect(x[:, j])) is not in this ansatz's basis"
        ))
        out[j] = θ[i]
    end
    return out
end

"""
    init_parameters(ansatz, rng; scale=0.01) -> Vector

Random initial parameters for `ansatz`.

Small values keep the initial state close to uniform, which is a deliberately uninformative
starting point rather than an accidental one.
"""
function init_parameters(a::LogStateVector, rng::AbstractRNG=Random.default_rng(); scale::Real=0.01)
    n = n_parameters(a)
    return scale .* (randn(rng, ComplexF64, n))
end
