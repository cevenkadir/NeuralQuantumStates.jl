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

train_state = Lux.Training.TrainState(model, ps, st, opt)

function log_ψ(m, x, ps, st)
    return Lux.apply(m, x, ps, st)
end

ps.layer_1.weight

fmap(vcat, ps, ps).layer_1.weight

a = log_ψ(model, test_set_device, ps, st)

log_ψ_c = @compile log_ψ(model, test_set_device, ps, Lux.testmode(st))



log_ψ_c(model, test_set_device, ps, st)[1] .- a[1]

log_ψ(model, test_set[:, 1], ps, st)[1]

log_ψ_c(model, test_set[:, 1], ps, st)[1]

# add MODE
function jacobian(
    f::Function,
    backend::Lux.AbstractADType,
    params,
    samples::AbstractArray,
    model::Lux.AbstractLuxLayer,
    layer_states;
    batch_size::Integer=100, shuffle::Bool=false,
    dev::Lux.MLDataDevices.AbstractDevice=Lux.cpu_device()
)
    @assert ndims(samples) >= 2

    batches = DataLoader(samples, batchsize=batch_size, shuffle=shuffle) |> dev

    #ADD compiled function via Reactant.jl

    #@compile f = f₀(model, samples[:, 1], params, Lux.testmode(layer_states))
    #f = f₀

    _f = (ps, x) -> f(model, x, ps, layer_states) |> first

    @views y₀, = _f(params, samples[:, 1] |> dev)

    function vjp(sample)
        return vector_jacobian_product(
            Base.Fix2(_f, sample),
            backend,
            params,
            ones(size(y₀))
        )
    end

    batch_results = [
        fmap(
            (x...) -> begin
                n_dims = x[1] |> ndims

                c_arr = cat(x...; dims=n_dims + 1)

                permutedims(c_arr, (n_dims + 1, (1:n_dims)...))
            end,
            vjp.(eachcol(batchᵢ))...
        )
        for batchᵢ in batches
    ]

    return fmap(vcat, batch_results...)
end

function jjacobian(
    f::Function,
    backend::Lux.AbstractADType,
    params,
    batchᵢ,
    model::Lux.AbstractLuxLayer,
    layer_states
)
    @assert ndims(batchᵢ) >= 2

    _f = (ps, x) -> f(model, x, ps, layer_states) |> first

    #@views y₀, = _f(params, batchᵢ[:, 1])

    function vjp(sample)
        return vector_jacobian_product(
            Base.Fix2(_f, sample),
            backend,
            params,
            ones((1,))
            #ones(size(y₀))
        )
    end

    # return fmap(
    #     (x...) -> begin
    #         n_dims = x[1] |> ndims

    #         c_arr = cat(x...; dims=n_dims + 1)

    #         permutedims(c_arr, (n_dims + 1, (1:n_dims)...))
    #     end,
    #     fmap(vjp, eachcol(batchᵢ))...
    # )

    #return 1.0#fmap(vjp, eachcol(batchᵢ))
    return vjp.(eachcol(batchᵢ))
end

function jjj(model, samples, params, layer_states)
    return jjacobian(log_ψ, AutoZygote(), params, samples, model, layer_states)
end

@benchmark qq = jacobian(log_ψ, AutoZygote(), ps, test_set, model, st;
    batch_size=16, dev=dev)

function merge_namedtuples(vec::Vector{<:NamedTuple})
    keys = fieldnames(typeof(vec[1]))
    merged = NamedTuple{keys}(map(k -> [getfield(nt, k) for nt in vec], keys))
    return merged
end

function stack_namedtuples(tuples::Vector{<:NamedTuple})
    keys = fieldnames(typeof(tuples[1]))
    stacked = NamedTuple{keys}(
        map(k -> stack([getfield(t, k) for t in tuples]), keys)
    )
    return stacked
end

# Base method for arrays or other non-struct fields
_stackfields(xs::Vector) = stack(xs)

# Recursive method for NamedTuples
function _stackfields(xs::Vector{<:NamedTuple})
    #println(xs[1] |> typeof)
    xkeys = keys(xs[1])
    println(NamedTuple{Tuple(xkeys)})
    return NamedTuple{Tuple(xkeys)}((
        k => _stackfields(getfield.(xs, k)) for k in xkeys
    ))
end

# Public API
stack_namedtuples_nested(vec::Vector) = _stackfields(vec)

function stack_layer()

end

function stack_namedtuples_nestedd(pss, st)
    layer_keys = st |> keys

    return NamedTuple{Tuple(layer_keys)}(
        (lk => stack_layer(getfield.(pss, lk))
         for lk in layer_keys)
    )
end

batchess = DataLoader(test_set, batchsize=2) |> dev
for batchᵢ in batchess
    (jjacobian(log_ψ, AutoZygote(), ps, batchᵢ, model, st) |> stack_namedtuples_nested).layer_1# |> println
    break
end

[jjj(model, batchᵢ, ps, st) for batchᵢ in batchess][1]

for batchᵢ in batchess
    ###jjacobian(log_ψ, AutoZygote(), ps, batchᵢ, model, st)[2].layer_1.weight |> size |> println
    jjj(model, batchᵢ, ps, st)#.layer_1.weight |> size |> println
end

for batchᵢ in batchess
    #jjacobian(log_ψ, AutoZygote(), ps, batchᵢ, model, st).layer_1.weight |> size |> println
    @compile jjj(model, batchᵢ, ps, st)#.layer_1.weight |> size |> println
end

@compile c_jjacobian = jjacobian(log_ψ, AutoZygote(), ps, test_set_device, model, Lux.testmode(st))

@compile log_ψ(model, test_set_device, ps, Lux.testmode(st))


fmap

qq.layer_1.weight |> size

hj = @compile jjacobian(
    model,
    test_set₀,
    ps,
    Lux.testmode(st)
)

@compile jjj(model, test_set_device, ps, Lux.testmode(st))





bbatches = DataLoader(test_set, batchsize=20) |> dev

for batchᵢ in bbatches
    println(size(batchᵢ |> eachcol))
end

stack(rand(3, 10) for _ in 1:2; dims=1)

fmap(x -> sum(x, dims=1), [[1 2; 3 8], [1 2; 3 8]])

#!
f = (tsd, params) -> first(log_ψ(model, tsd, params, st))
wq, = Zygote.gradient(f1, ps)

Base.Fix1()

f(ps |> ComponentArray, test_set_device_2)

#!






#!

y, st = log_ψ(model, test_set_device_2, c_ps, st)

f = x -> first(log_ψ(model, test_set_device_2, x, st))

f(ps)


a = vector_jacobian_product(f, AutoZygote(), c_ps, ones(size(y)))

a.layer_1.weight
ps.layer_1.weight

#!

batched_jacobian(x -> first(log_ψ(model, x, c_ps, st)), AutoZygote(), test_set_device_2)

Zygote.jacobian(x -> log_ψ(model, test_set_device_2[:, 1], x, st), c_ps)

Zygote.jacobian(Base.Fix1(smodel, x), ps)

ndims(y)

size(y, ndims(y))
size(ps, ndims(ps))

batched_jacobian(x -> log_ψ(model, test_set_device_2, x, st), AutoZygote(), ps)

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
