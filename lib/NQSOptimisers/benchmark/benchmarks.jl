"""
Benchmarks for the preconditioners and the linear solvers behind them.

Run with

```
julia --project=lib/NQSOptimisers/benchmark lib/NQSOptimisers/benchmark/benchmarks.jl
```

These are not wired into CI. They exist so that the cost of *forming* the quantum geometric
tensor stays visible next to the cost of merely multiplying by it — which is the whole argument
for `mode=:matrixfree`, and the thing that decides whether a run with a real network fits in
memory at all.

The `qgt` group deliberately uses a synthetic design matrix rather than a wavefunction. The
shapes are what matter here, and a `LogStateVector` cannot reach the parameter counts where the
`P × P` tensor becomes the problem.
"""

using BenchmarkTools
using ConnectedBasisConfigurations
using DifferentiationInterface
using ForwardDiff
using LinearAlgebra
using NQSCore
using NQSOptimisers
using OperatorAlgebra
using Printf
using Random
using SymBasis

const BACKEND = AutoForwardDiff()

BenchmarkTools.DEFAULT_PARAMETERS.seconds = 2

"""Transverse-field Ising on a periodic chain, as an `OpSum`."""
function tfi(nsites; J=1.0, h_x=1.0, h_z=0.0)
    ops = local_operators(Spin(1 // 2))
    σz, σx = 2 .* ops.sz, 2 .* ops.sx
    terms = AbstractOp[]
    for i in 1:nsites
        push!(terms, J * (Op(σz, i) * Op(σz, mod1(i + 1, nsites))))
        iszero(h_z) || push!(terms, h_z * Op(σz, i))
        iszero(h_x) || push!(terms, h_x * Op(σx, i))
    end
    return OpSum(terms)
end

function fresh_state(nsites; seed=11, scale=0.3)
    spec = Spin(1 // 2)
    b = basis(dof_object(spec), nsites)
    a = LogStateVector(spec, nsites, b)
    return FullSumState(a, init_parameters(a, Xoshiro(seed); scale=scale); backend=BACKEND)
end

report(label, trial) = @printf("%-52s %10s\n", label, BenchmarkTools.prettytime(minimum(trial).time))

const SUITE = BenchmarkGroup()

# --------------------------------------------------- the three modes on a real wavefunction

let nsites = 8
    H = compile(tfi(nsites; h_x=0.9, h_z=0.1))
    vs = fresh_state(nsites)
    cg = ConjugateGradientSolver(; tol=1e-10)

    SUITE["sr8"][":sr, cholesky"] = @benchmarkable precondition(
        StochasticReconfiguration(; diag_shift=1e-3, mode=:sr, solver=CholeskySolver()), $vs, $H
    )
    SUITE["sr8"][":minsr, cholesky"] = @benchmarkable precondition(
        StochasticReconfiguration(; diag_shift=1e-3, mode=:minsr, solver=CholeskySolver()),
        $vs, $H
    )
    SUITE["sr8"][":sr, conjugate gradient"] = @benchmarkable precondition(
        StochasticReconfiguration(; diag_shift=1e-3, mode=:sr, solver=$cg), $vs, $H
    )
    SUITE["sr8"][":matrixfree, conjugate gradient"] = @benchmarkable precondition(
        StochasticReconfiguration(; diag_shift=1e-3, mode=:matrixfree, solver=$cg), $vs, $H
    )
end

# ------------------------------------ forming the tensor versus multiplying by it, by shape

# (rows, params). The last two are the regime a network actually lives in: far more parameters
# than samples, where the P x P tensor is the thing that does not fit.
for (n, p) in ((512, 128), (512, 1024), (1024, 4096))
    X = randn(Xoshiro(0), n, p)
    v = randn(Xoshiro(1), p)
    S = QuantumGeometricTensor(X)
    group = "qgt $(n)x$(p)"

    SUITE[group]["form XᵀX"] = @benchmarkable transpose($X) * $X
    SUITE[group]["one product, matrix-free"] = @benchmarkable $S * $v
    SUITE[group]["50 products, matrix-free"] =
        @benchmarkable (for _ in 1:50; $S * $v; end)
end

# --------------------------------------------------------------- the solvers on their own

# A realistic geometric tensor rather than a random dense one: rank-deficient, because a
# geometric tensor always is — redundant parameter directions and unexplored directions both
# give exact zero modes, which is what the shift exists to cover.
let rng = Xoshiro(2), n = 256, p = 512
    X = randn(rng, n, p)
    A = transpose(X) * X
    b = randn(rng, p)

    SUITE["solvers"]["cholesky, 512 (rank 256)"] =
        @benchmarkable solve(CholeskySolver(), $A, $b, 1e-3)
    SUITE["solvers"]["pinv, 512 (rank 256)"] =
        @benchmarkable solve(PseudoInverseSolver(), $A, $b, 1e-3)
    SUITE["solvers"]["conjugate gradient, 512 (rank 256)"] =
        @benchmarkable solve(ConjugateGradientSolver(; tol=1e-8), $A, $b, 1e-3)
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
