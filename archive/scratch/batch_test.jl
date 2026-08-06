using Lux
using Random;
using Zygote
rng = Xoshiro(42);

input, output = 6, 2
model = Chain(Dense(input => input^2, tanh), Dense(input^2 => output));
ps, st = Lux.setup(rng, model);

n = 100
x = rand(rng, Float32, input, n);
f = StatefulLuxLayer{true}(model, ps, st)

@code_warntype f(x) # Type stable, return type `Matrix{Float32}`
#const backend = AutoForwardDiff()
@code_warntype batched_jacobian(f, AutoForwardDiff(; chunksize=6), x)

batched_jacobian(f, AutoForwardDiff(; chunksize=6), x)
batched_jacobian(f, AutoZygote(), x)
