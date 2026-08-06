using Lux
using Random;
using Enzyme
using ComponentArrays
using BenchmarkTools

add_dim(x::Array) = reshape(x, (1, size(x)...))

rng = Xoshiro(42);

input, output = 6, 2
model = Chain(Dense(input => input^2, tanh), Dense(input^2 => output));
ps, st = Lux.setup(rng, model);

θ = ComponentArray(ps)
ax = only(getaxes(θ))

n = 100
σₛ = rand(rng, Float32, input, n)


@benchmark model(σₛ, θ, st)

model(σₛ, θ, st)[1]

function rmodel(_σₛ, _raw_full_θ, _st)
    return model(_σₛ, ComponentArray(_raw_full_θ, ax), _st)
end

@benchmark rmodel(σₛ, θ |> getdata, st)
@benchmark model(σₛ, θ, st)

#! working 2?

# Reverse
# ReverseHolomorphic
# ReverseWithPrimal
# ReverseHolomorphicWithPrimal

# function chunked_jacobian(
#     mode::Union{typeof(Reverse),typeof(ReverseWithPrimal)},
#     sfun,
#     σₛ,
#     θ::ComponentArray{T_p};
#     chunk_size::Integer=16
# ) where {T_p<:AbstractFloat}
#     ax_θ = getaxes(θ)
#     raw_θ = getdata(θ)

#     function f(_raw_θ)
#         _θ = ComponentArray(_raw_θ, ax_θ)
#         return Base.Fix1(sfun, σₛ)(_θ)
#     end

#     raw_O = Enzyme.jacobian(mode, f, raw_θ; chunk=Val(chunk_size)) |> only

#     return ComponentArray.(
#         eachslice(
#             raw_O;
#             dims=tuple((collect(1:length(size(raw_O))-1))...)
#         ),
#         ax_θ
#     )
# end

function chunked_jacobian(
    mode::Union{typeof(Reverse),typeof(ReverseWithPrimal)},
    sfun,
    σₛ,
    θ::ComponentArray{T_p};
    chunk_size::Integer=16
) where {T_p<:AbstractFloat}
    ax_θ = getaxes(θ)
    raw_θ = getdata(θ)

    function f(_raw_θ)
        _θ = ComponentArray(_raw_θ, ax_θ)
        return Base.Fix1(sfun, σₛ)(_θ)
    end

    raw_O = Enzyme.jacobian(mode, f, raw_θ; chunk=Val(chunk_size)) |> only

    return ComponentArray.(
        eachslice(
            raw_O;
            dims=tuple((collect(1:length(size(raw_O))-1))...)
        ),
        ax_θ
    )
end

function chunked_jacobian(
    mode::Union{typeof(Reverse),typeof(ReverseWithPrimal)},
    sfun,
    σₛ,
    θ::ComponentArray{T_p};
    chunk_size::Integer=16
) where {T_p<:Complex{<:AbstractFloat}}
    #! might requires further optimization

    full_θ = ComponentArray(real=real(θ), imag=imag(θ))
    ax_full_θ = getaxes(full_θ)
    raw_full_θ = getdata(full_θ)

    function f(_raw_full_θ)
        _full_θ = ComponentArray(_raw_full_θ, ax_full_θ)
        _θ = _full_θ.real .+ im .* _full_θ.imag
        return Base.Fix1(sfun, σₛ)(_θ)
    end

    raw_O = Enzyme.jacobian(mode, f, raw_full_θ; chunk=Val(chunk_size)) |> only

    O = ComponentArray.(
        eachslice(
            raw_O;
            dims=tuple((collect(1:length(size(raw_O))-1))...)
        ),
        ax_full_θ
    )

    return map(_θ -> _θ.real .+ im .* _θ.imag, O)
end

sfun = StatefulLuxLayer{true}(model, θ, st)

sfun(σₛ, θ)[1, :]

logpsi = (_σₛ, _θ) -> sfun(_σₛ, _θ)[1, :] .+ im .* sfun(_σₛ, _θ)[2, :]

logpsi(σₛ, θ)

logpsi_real = (_σₛ, _θ) -> sfun(_σₛ, _θ)[1, :]

logpsi_real(σₛ, θ)

a = Enzyme.jacobian(Enzyme.Reverse, logpsi_real, Enzyme.Const(σₛ), ps)

a = Enzyme.jacobian(Enzyme.Reverse, logpsi_real, Enzyme.Const(σₛ), θ)


@benchmark Oₖ = chunked_jacobian(Enzyme.Reverse, sfun, σₛ, θ) #|> only


(first ∘ Lux.apply)(model, σₛ, θ, st)
(first ∘ rmodel)(σₛ, θ |> getdata, st)

# unsafe to copy error (probably related to ComponentArrays)
Enzyme.jacobian(Enzyme.Reverse, first ∘ Lux.apply, Const(model), Const(σₛ), θ, Const(st))


@benchmark _, _, a, _ = Enzyme.jacobian(Enzyme.Reverse, first ∘ Lux.apply, Const(model), Const(σₛ), ps, Const(st))


@benchmark chunked_jacobian(Enzyme.Reverse, sfun, σₛ, θ)

@benchmark begin
    _, a = Enzyme.jacobian(Enzyme.Reverse, sfun, Enzyme.Const(σₛ), ps)

    #ComponentArray.(a)
end


chunked_jacobian_compiled = @compile chunked_jacobian(Enzyme.Reverse, sfun, σₛ, θ)

@benchmark chunked_jacobian(Enzyme.Reverse, sfun, σₛ, θ) #|> only

cθ = @views Oₖ[1, :] .+ im .* Oₖ[2, :] |> sum

cθ |> imag

sfun(σₛ, cθ)[1, :]

@benchmark chunked_jacobian(Enzyme.Reverse, (_σₛ, _θ) -> real(logpsi(_σₛ, _θ)), σₛ, θ)
@benchmark chunked_jacobian(Enzyme.Reverse, (_σₛ, _θ) -> real(logpsi(_σₛ, _θ)), σₛ, cθ)

chunked_jacobian(Enzyme.Reverse, (_σₛ, _θ) -> real(logpsi(_σₛ, _θ)), σₛ, cθ) |> sum;


@benchmark ComponentArray(real=real(cθ), imag=imag(cθ))


hcat(θ, θ)

θ

hcat(θ, θ)[:layer_2, 2]
