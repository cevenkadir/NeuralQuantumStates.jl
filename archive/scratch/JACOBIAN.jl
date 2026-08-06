using Lux
using Random;
using Enzyme, Reactant
using ComponentArrays
using BenchmarkTools

const dev = reactant_device(; force=true)

add_dim(x::Array) = reshape(x, (1, size(x)...))

rng = Xoshiro(42);

input, output = 6, 2
model = Chain(Dense(input => input^2, tanh), Dense(input^2 => output));
ps, st = Lux.setup(rng, model) |> dev;
#st = st |> dev

θ = ComponentArray(ps) |> dev
ax = only(getaxes(θ))

n = 100
σₛ = rand(rng, Float32, input, n) |> dev


model_compiled = @compile model(σₛ, θ, Lux.testmode(st))

@benchmark model_compiled(σₛ, θ, st)
@benchmark model(σₛ, θ, st)

model_compiled(σₛ, θ, st)[1] .- model(σₛ, θ, st)[1]

#! working 2?

# Reverse
# ReverseHolomorphic
# ReverseWithPrimal
# ReverseHolomorphicWithPrimal

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

function chunked_jacobian_raw(
    mode::Union{typeof(Reverse),typeof(ReverseWithPrimal)},
    sfun,
    σₛ,
    raw_θ::AbstractVector;
    chunk_size::Integer=16, ax_θ
) where {T_p<:AbstractFloat}
    #ax_θ = getaxes(θ)
    #raw_θ = getdata(θ)

    function f(_raw_θ)
        _θ = ComponentArray(_raw_θ, ax_θ)
        return Base.Fix1(sfun, σₛ)(_θ)
    end

    return Enzyme.jacobian(mode, f, raw_θ; chunk=Val(chunk_size)) |> only

    # return ComponentArray.(
    #     eachslice(
    #         raw_O;
    #         dims=tuple((collect(1:length(size(raw_O))-1))...)
    #     ),
    #     ax_θ
    # )
end


sfun = StatefulLuxLayer{true}(model, θ, st)
ssfun = StatefulLuxLayer{true}(model, ps, st)

sfun(σₛ, θ)[1, :]

logpsi = (_σₛ, _θ) -> sfun(_σₛ, _θ)[1, :] .+ im .* sfun(_σₛ, _θ)[2, :]

chunked_jacobian_raw(Enzyme.Reverse, sfun, σₛ, θ |> getdata; ax_θ=getaxes(θ))


Enzyme.jacobian(Enzyme.Reverse, Base.Fix1(ssfun, σₛ), ps)

Base.Fix1(sfun, σₛ)(θ)

@compile sfun(σₛ, θ)


@compile first(Lux.apply(model, σₛ, θ, st))

@compile Enzyme.jacobian(Enzyme.Reverse, sfun, Const(σₛ), θ)

dθ = Enzyme.make_zero(θ)
Enzyme.jacobian(Enzyme.Reverse, sfun, Const(σₛ), Duplicated(θ, dθ))

@compile Enzyme.jacobian(Enzyme.Reverse, mmodel, σₛ |> Const, θ, st |> Const)

# crashes
chunked_jacobian(Enzyme.Reverse, logpsi_real, σₛ, θ)
Oₖ_real = chunked_jacobian(Enzyme.Reverse, (_σₛ, _θ) -> real(logpsi(_σₛ, _θ)), σₛ, θ)
Oₖ_imag = chunked_jacobian(Enzyme.Reverse, (_σₛ, _θ) -> imag(logpsi(_σₛ, _θ)), σₛ, θ)

Oₖ = chunked_jacobian(Enzyme.Reverse, sfun, σₛ, θ) #|> only


chunked_jacobian_compiled = @compile chunked_jacobian(Enzyme.Reverse, sfun, σₛ, θ)

@benchmark chunked_jacobian_compiled(Enzyme.Reverse, sfun, σₛ, θ)

@benchmark Oₖ = chunked_jacobian(Enzyme.Reverse, sfun, σₛ, θ) #|> only

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
