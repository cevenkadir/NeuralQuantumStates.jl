using NeuralQuantumStates: Lattices, Hilberts, Operators
using Lux, Zygote#Reactant
using Optimisers
using MLUtils: splitobs, DataLoader, shuffleobs, load_iris
using Random
using ComponentArrays
using LinearAlgebra
using Functors

#using Metal
using BenchmarkTools
#Metal.allowscalar(false)

const dev = reactant_device()

opt = Adam(0.03f0);
vjp_rule = AutoEnzyme();


lat = Lattices.build(:Hypercube, [32], 1.0; periodic=[true]);
hil = Hilberts.build(:Fock, 5, Lattices.nv(lat); ∑n=5.0);
ham = Operators.build(:ExtendedBoseHubbard, hil, lat; J=1.0, U=1.0, V=1.0, μ=0.0);


α = 2
model = Chain(
    Dense(Lattices.nv(lat) => α * Lattices.nv(lat), tanh),
    Dense(α * Lattices.nv(lat) => 1)
)

batchsize = 40;


test_set₀ = Int32[
    1 1 1 0 1 0 0 1 0 0 0 0 0 0 0 0 1 1 1 0 1 0 0 1 0 0 0 0 0 0 0 0;
    1 1 0 1 0 0 0 1 0 0 0 1 0 0 0 0 1 1 0 1 0 0 0 1 0 0 0 1 0 0 0 0;
    0 1 0 0 0 0 0 1 0 0 0 0 0 1 1 1 0 1 0 0 0 0 0 1 0 0 0 0 0 1 1 1;
    1 0 0 1 2 0 0 1 0 0 0 0 0 0 0 0 1 0 0 1 2 0 0 1 0 0 0 0 0 0 0 0;
    0 1 0 0 0 0 0 1 0 0 0 0 0 1 1 1 0 1 0 0 0 0 0 1 0 0 0 0 0 1 1 1;
    1 0 0 1 2 0 0 1 0 0 0 0 0 0 0 0 1 0 0 1 2 0 0 1 0 0 0 0 0 0 0 0
]' |> collect;
test_set = convert(Matrix{Float32}, repeat(test_set₀, inner=(1, 20)));
test_set_device = test_set |> dev;
test_set_device_2 = test_set₀ |> dev;


rng = Random.default_rng()
Random.seed!(rng, 0)

ps, st = Lux.setup(rng, model) |> dev
cps = ComponentArray(ps)

f = StatefulLuxLayer{true}(model, cps, st)
ff = Base.Fix1(f, test_set_device_2)

ff(cps)

rr = only(Zygote.jacobian(ff, cps))

ax = only(getaxes(cps))
@benchmark er_rr = rr |> eachrow
@benchmark Ø = ComponentArray.(er_rr, ax |> Ref)
