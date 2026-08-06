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

"""Real element type to feed the network, matching whatever the parameters are made of."""
_input_type(θ) = Float64
_input_type(θ::NamedTuple) = isempty(θ) ? Float64 : _input_type(first(values(θ)))
_input_type(θ::AbstractArray) = real(eltype(θ))

function NQSCore.log_amplitude(a::LuxAnsatz, θ, x::AbstractMatrix)
    T = _input_type(θ)
    input = T.(x)
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
