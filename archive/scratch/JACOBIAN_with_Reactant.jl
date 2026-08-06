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
ps = ps |> dev
st = st |> dev


n = 100
σₛ = rand(rng, Float32, input, n) |> dev

σₛ = rand(rng, Float32, input) |> dev


# it seems CompenentArrays creates an issue together with Reactant.jl
# so let's stick to NamedTuples (for now)

Enzyme.jacobian(Enzyme.Reverse, first ∘ Lux.apply,
    Enzyme.Const(model), Enzyme.Const(σₛ), ps, Enzyme.Const(st)
)

Enzyme.jacobian(Enzyme.Reverse, first ∘ Lux.apply,
    Enzyme.Const(model), σₛ, Enzyme.Const(ps), Enzyme.Const(st)
)

(first ∘ Lux.apply)(model, σₛ, ps, st)


function f(model, x, _ps, _st)
    return sum((first ∘ Lux.apply)(model, x, _ps, _st); dims=1)
end

f(model, σₛ, ps, st)

Enzyme.jacobian(Enzyme.Reverse, f,
    Enzyme.Const(model), Enzyme.Const(σₛ), ps, Enzyme.Const(st)
)

Enzyme.jacobian(Enzyme.Reverse, f,
    Enzyme.Const(model), σₛ, Enzyme.Const(ps), Enzyme.Const(st)
)


Enzyme.jacobian(Enzyme.Reverse, (x, y) -> x .* y,
    rand(3), Const(rand(3))
)


Enzyme.jacobian(Enzyme.Reverse,
    Base.Fix1(StatefulLuxLayer{true}(model, ps, st), σₛ),
    ps |> ComponentArray
)
