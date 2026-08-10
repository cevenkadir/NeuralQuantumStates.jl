"""
    LuxAnsatz(model, dof, nsites; rng, states=nothing) <: AbstractAnsatz

Wraps any Lux model as a variational wavefunction.

`NQSCore` asks for `log_amplitude(ansatz, θ, x)` and knows nothing about Lux; Lux models want
`Lux.apply(model, x, ps, st)`. Everything Lux offers is available through this one type, and
`NQSCore` keeps no dependency on it.

# Output convention

The model must return either

- a length-`batch` vector of log-amplitudes, which the layers here do; or
- a `(2, batch)` array, whose two rows are read as the real and imaginary parts of `log ψ`.

The second form is how a real-valued network represents a complex wavefunction.

Configurations arrive as physical local values — `Rational` for spins — which no network can
consume, so they are converted once at this boundary and no layer has to think about it.

# Fields
- `model`: the Lux layer or chain.
- `states`: Lux layer states (empty for the layers here, which are all stateless).
- `dof`, `nsites`: the space the wavefunction is defined on.
"""
struct LuxAnsatz{M,S,D} <: AbstractAnsatz
    model::M
    states::S
    dof::D
    nsites::Int
end

function LuxAnsatz(
    model, dof, nsites::Integer;
    rng::AbstractRNG=Random.default_rng(), states=nothing
)
    st = states === nothing ? Lux.initialstates(rng, model) : states
    return LuxAnsatz(model, st, dof, Int(nsites))
end

"""
    init_parameters(ansatz, rng) -> NamedTuple

Initial parameters for a [`LuxAnsatz`](@ref), from the model's own initialization.
"""
function NQSCore.init_parameters(a::LuxAnsatz, rng::AbstractRNG=Random.default_rng())
    return Lux.initialparameters(rng, a.model)
end

NQSCore.n_parameters(a::LuxAnsatz) = Lux.parameterlength(a.model)

"""
The full space implied by the ansatz's degrees of freedom.

A network is defined on every configuration of its input space, so that is what a `FullSumState`
built on one sums over unless told otherwise. It lives here rather than in `NQSCore`, which has
no reason to depend on a basis library in order to guess.
"""
NQSCore.default_basis(a::LuxAnsatz) =
    SymBasis.Bases.basis(SymBasis.dof_object(a.dof), a.nsites)

"""
The type this network computes in, usually complex for a wavefunction.

Answering lets `NQSCore` build a to-be-differentiated batch in this type directly, so that the
mixed complex-real matrix products in the reverse pass — which no BLAS has a kernel for — never
arise.
"""
NQSCore.input_type(::LuxAnsatz, θ) =
    (r = _reference_array(θ); r === nothing ? nothing : eltype(r))

"""
    colocate(reference, input)

Put `input` wherever `reference` lives, when the two are not already on the same side of the
host/device boundary.

The parameters decide the device — `Lux` moves them, and the batch has to follow or the first
matrix multiplication mixes host and device memory — so there is no device field to keep in
sync, and on the CPU this does nothing at all.

Membership is decided by unwrapping rather than by testing for `Array`. Parameters do not always
arrive as plain arrays even on the CPU: differentiating through `NQSCore.flatten_parameters`
hands back views into a `ComponentArray`, and copying one of those is a mutation reverse-mode AD
refuses to differentiate.
"""
colocate(reference, input::AbstractArray) =
    _is_host(reference) == _is_host(input) ? input :
    copyto!(similar(reference, eltype(input), size(input)), input)

# The transfer is a `copyto!`, which reverse-mode AD refuses to differentiate. It does not have
# to: configurations are data, not parameters, and nothing needs a gradient with respect to
# them. Saying so is what lets the gradient reach `θ` instead of stopping here.
ChainRulesCore.@non_differentiable colocate(::Any, ::Any)

_is_host(::Array) = true
function _is_host(a::AbstractArray)
    p = parent(a)
    return p === a ? false : _is_host(p)
end
_is_host(_) = true

function NQSCore.log_amplitude(a::LuxAnsatz, θ, x::AbstractMatrix)
    reference = _reference_array(θ)

    # Promote, never convert. A configuration is real, so a forward pass wants the narrowest
    # real type the parameters imply — but `NQSCore` builds the batch a derivative follows in
    # the network's own complex type, and narrowing that would drop the imaginary part.
    narrowest = reference === nothing ? Float64 : real(eltype(reference))
    T = promote_type(eltype(x), narrowest)
    input = colocate(reference, T === eltype(x) ? x : T.(x))

    y, _ = Lux.apply(a.model, input, θ, a.states)
    return _as_log_amplitude(y)
end

"""Interpret a model's output as a vector of complex log-amplitudes."""
_as_log_amplitude(y::AbstractVector) = y
function _as_log_amplitude(y::AbstractMatrix)
    if size(y, 1) == 1
        return vec(y)
    elseif size(y, 1) == 2
        # A real network encoding a complex wavefunction: modulus and phase in two rows.
        return @views complex.(y[1, :], y[2, :])
    end
    throw(ArgumentError(
        "an ansatz model must return 1 or 2 rows, got $(size(y, 1)); a single row is a " *
        "complex log-amplitude, two rows are its real and imaginary parts"
    ))
end
