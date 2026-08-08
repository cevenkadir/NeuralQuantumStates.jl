"""
    LuxAnsatz(model, dof, nsites; rng, states=nothing) <: AbstractAnsatz

Wraps any Lux model as a variational wavefunction.

This is the bridge between the two worlds: `NQSCore` asks for `log_amplitude(ansatz, θ, x)` and
knows nothing about Lux, while Lux models want `Lux.apply(model, x, ps, st)`. Everything Lux
offers — layer composition, GPU movement, the whole ecosystem — is available through this one
type, and `NQSCore` keeps no dependency on it.

# Output convention

The model must return either

- a length-`batch` vector of log-amplitudes, which the layers here do; or
- a `(2, batch)` array, whose two rows are read as the real and imaginary parts of `log ψ`.

The second form is how a real-valued network represents a complex wavefunction, and it replaces
the `ComplexF32[1 im] * model(x)` idiom that was scattered through the pre-split scratch code.

# Input conversion

Configurations arrive as physical local values — `Rational` for spins — which no neural network
can consume. They are converted to the real element type matching the parameters, once, at this
boundary. Doing it here rather than inside each layer means a layer never has to think about it.

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

A network is defined on every configuration of its input space, so that is what a
`FullSumState` built on one should sum over unless told otherwise. The default lives here rather
than in `NQSCore`, which has no way to know what an arbitrary ansatz spans — and no reason to
depend on a basis library in order to guess.
"""
NQSCore.default_basis(a::LuxAnsatz) =
    SymBasis.Bases.basis(SymBasis.dof_object(a.dof), a.nsites)

"""Real element type to feed the network, matching whatever the parameters are made of."""
_input_type(θ) = Float64
_input_type(θ::NamedTuple) = isempty(θ) ? Float64 : _input_type(first(values(θ)))
_input_type(θ::AbstractArray) = real(eltype(θ))

"""Any array among the parameters, whose type says where the network expects to be run."""
_reference_array(θ) = nothing
_reference_array(θ::NamedTuple) = isempty(θ) ? nothing : _reference_array(first(values(θ)))
_reference_array(θ::AbstractArray) = θ

"""
    colocate(reference, input)

Put `input` wherever `reference` lives.

The parameters decide the device: `Lux` moves them, and the batch has to follow or the first
matrix multiplication mixes host and device memory. Inferring it from the parameters rather than
storing a device in the ansatz means there is no new field, no new dependency and nothing to
keep in sync — and on the CPU it does nothing at all.

Whether the reference is in host memory is decided by unwrapping it rather than by testing for
`Array` directly. Parameters do not always arrive as plain arrays even on the CPU: rebuilding
them from a flat vector, which is exactly what differentiating through
`NQSCore.flatten_parameters` does, hands back views into a `ComponentArray`. Treating one of
those as foreign would copy it, and a copy is a mutation that reverse-mode AD refuses to
differentiate through.
"""
colocate(reference, input::AbstractArray) =
    _is_host(reference) ? input : copyto!(similar(reference, eltype(input), size(input)), input)

_is_host(::Nothing) = true
_is_host(::Array) = true
function _is_host(a::AbstractArray)
    p = parent(a)
    return p === a ? false : _is_host(p)
end
_is_host(_) = true

function NQSCore.log_amplitude(a::LuxAnsatz, θ, x::AbstractMatrix)
    # Configurations arrive as exact rationals, which no network and no accelerator wants, so
    # the conversion to a float type has to happen anyway; doing the transfer in the same step
    # keeps the host/device boundary to this one line.
    T = _input_type(θ)
    input = colocate(_reference_array(θ), T.(x))
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
        return @views y[1, :] .+ im .* y[2, :]
    end
    throw(ArgumentError(
        "an ansatz model must return 1 or 2 rows, got $(size(y, 1)); a single row is a " *
        "complex log-amplitude, two rows are its real and imaginary parts"
    ))
end
