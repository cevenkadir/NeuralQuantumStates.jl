"""
Benchmarks for the connected-configurations kernel.

Run with

```
julia --project=lib/ConnectedBasisConfigurations/benchmark lib/ConnectedBasisConfigurations/benchmark/benchmarks.jl
```

These are not wired into CI. They exist so that a change to the kernel can be measured against
something rather than argued about, and so the cost of *not* compiling an operator up front
stays visible.
"""

using BenchmarkTools
using ConnectedBasisConfigurations
using OperatorAlgebra
using Printf
using Random
using SymBasis

chain_bonds(n) = [(i, mod1(i + 1, n)) for i in 1:n]

function transverse_field_ising(nsites; J=1.0, h_x=1.0, h_z=0.5)
    ops = local_operators(Spin(1 // 2))
    σz, σx = 2 .* ops.sz, 2 .* ops.sx
    terms = AbstractOp[]
    for (i, j) in chain_bonds(nsites)
        push!(terms, J * (Op(σz, i) * Op(σz, j)))
    end
    for i in 1:nsites
        push!(terms, h_z * Op(σz, i))
        push!(terms, h_x * Op(σx, i))
    end
    return OpSum(terms)
end

function extended_bose_hubbard(nsites, n_max; J=1.0, U=1.0, V=1.0, μ=0.5)
    ops = local_operators(Boson(n_max))
    a, adag, n = ops.a, ops.adag, ops.n
    terms = AbstractOp[]
    for (i, j) in chain_bonds(nsites)
        push!(terms, (-J) * (Op(adag, i) * Op(a, j)))
        push!(terms, (-J) * (Op(adag, j) * Op(a, i)))
        push!(terms, V * (Op(n, i) * Op(n, j)))
    end
    for i in 1:nsites
        push!(terms, (U / 2) * Op(n * n - n, i))
        push!(terms, (-μ) * Op(n, i))
    end
    return OpSum(terms)
end

random_states(spec, nsites, count, rng) =
    [packed(spec, [rand(rng, local_values(spec)) for _ in 1:nsites]) for _ in 1:count]

report(label, trial) = @printf("%-52s %10s\n", label, BenchmarkTools.prettytime(minimum(trial).time))

const SUITE = BenchmarkGroup()

let rng = Xoshiro(0), nsites = 16
    spec = Spin(1 // 2)
    H = transverse_field_ising(nsites)
    compiled = compile(H)
    states = random_states(spec, nsites, 1024, rng)

    height = max_conn_size(compiled)
    configs = Matrix{eltype(states)}(undef, height, length(states))
    mels = Matrix{Float64}(undef, height, length(states))
    counts = Vector{Int}(undef, length(states))

    SUITE["tfi16"]["compile"] = @benchmarkable compile($H)
    SUITE["tfi16"]["padded, compiling per call"] = @benchmarkable connected_padded($H, $states)
    SUITE["tfi16"]["padded, precompiled"] = @benchmarkable connected_padded($compiled, $states)
    SUITE["tfi16"]["padded, in place"] =
        @benchmarkable connected_padded!($configs, $mels, $counts, $compiled, $states)
    SUITE["tfi16"]["connected, single state"] =
        @benchmarkable connected($compiled, $(first(states)))
end

let rng = Xoshiro(1), nsites = 10, n_max = 3
    spec = Boson(n_max)
    H = extended_bose_hubbard(nsites, n_max)
    compiled = compile(H)
    states = random_states(spec, nsites, 512, rng)

    SUITE["bhm10"]["padded, compiling per call"] = @benchmarkable connected_padded($H, $states)
    SUITE["bhm10"]["padded, precompiled"] = @benchmarkable connected_padded($compiled, $states)
end

let nsites = 12
    spec = Spin(1 // 2)
    dofo = dof_object(spec)
    b = basis(dofo, nsites, sym(TotalMagnetization(0 // 1, nsites), dofo))
    H = transverse_field_ising(nsites; h_z=0.0)
    sector = compile(H, b)

    SUITE["sector12"]["padded, compiling per call"] = @benchmarkable connected_padded($H, $(b.states), $b)
    SUITE["sector12"]["padded, precompiled"] = @benchmarkable connected_padded($sector, $(b.states))
end

if abspath(PROGRAM_FILE) == @__FILE__
    results = run(SUITE; verbose=false)
    for (group, trials) in sort(collect(results); by=first)
        println("\n", group)
        for (name, trial) in sort(collect(trials); by=first)
            report("  " * name, trial)
        end
    end
end
