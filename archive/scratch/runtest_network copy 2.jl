using NeuralQuantumStates: Lattices, Hilberts, Operators
using Lux, Reactant, Enzyme
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
    Dense(α * Lattices.nv(lat) => 2)
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

@assert size(test_set, 2) % batchsize == 0
train_dataloader = DataLoader(test_set, batchsize=batchsize, shuffle=false) |> dev

rng = Random.default_rng()
Random.seed!(rng, 0)

ps, st = Lux.setup(rng, model) |> dev

train_state = Lux.Training.TrainState(model, ps, st, opt)


function log_ψ(m, x, ps, st)
    #st_model = StatefulLuxLayer{true}(model, ps, st)
    return ComplexF32[1.0 im] * first(Lux.apply(m, x, ps, st)) #sum(model(x, ps, st) |> first, dims=1)
    #return st_model(x)
end

function loss_function_batched(model, x, ps, st)
    # Make it a stateful layer
    smodel = StatefulLuxLayer{true}(model, ps, st)

    J = batched_jacobian(smodel, AutoZygote(), x)

    return J
end

loss_function_batched(model, test_set_device_2, ps, st)

function main(tstate::Training.TrainState, vjp, data, epochs)
    for epoch in 1:epochs
        for (x, y) in data
            _, loss, _, tstate = Training.single_train_step!(vjp, loss_function, (x, y), tstate)
            # if epoch % 50 == 1 || epoch == epochs
            #     @printf "Epoch: %3d \t Loss: %.5g\n" epoch loss
            # end
        end
    end
    return tstate
end




epochs = 10;




g_batched = Flux.gradient((f, x) -> sum(abs2, f(x)), model, test_set)

@btime for epoch = 1:epochs
    g_batched = Flux.gradient((f, x) -> sum(abs2, f(x)), model, test_set)
    #Flux.Optimise.update!(opt_cpu, Flux.params(model), g_batched)
end


g_batched = Flux.gradient((f, x) -> sum(abs2, f(x)), m_model, m_test_set)

@btime for epoch = 1:epochs
    g_batched = Flux.gradient((f, x) -> sum(abs2, f(x)), m_model, m_test_set)
    #Flux.Optimise.update!(opt_cpu, Flux.params(model), g_batched)
end
