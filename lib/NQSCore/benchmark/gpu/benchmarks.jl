"""
CPU versus CUDA, for the parts of a variational run that could plausibly move to a device.

Run with

```
julia --project=lib/NQSCore/benchmark/gpu lib/NQSCore/benchmark/gpu/benchmarks.jl
```

On a machine without a working CUDA device it prints why and exits cleanly, so it is harmless
to run anywhere. No package in the stack gains a GPU dependency from this directory: the
environment here is separate, and CUDA is a weak dependency of `NQSOptimisers` only.

# What to look at first

The **transfer** section, not the speedups. Connected configurations are computed on the host —
the compiled operator is a nested structure that cannot be uploaded as it stands — so every step
ships the samples out and the connected configurations and matrix elements back. The return leg
is `max_conn` times larger than the outbound one, and that asymmetry is the whole question:
if it dominates, porting the connected-configuration kernel to the device is worth the work,
and if it does not, the split architecture is fine as it is. NetKet shipped the same split for
years before writing device-side operators.

The `max_conn / mean(n_conn)` ratio printed alongside is the second half of the answer. At 1 the
padding is free and a device kernel would waste nothing; well above 1 it means most of what
would be transferred, and most of what a device kernel would evaluate, is padding.

Each measurement is guarded: a failure prints and the run continues, because most of the device
paths here have never executed on hardware and one broken step should not hide the rest.
"""

using BenchmarkTools
using CUDA
using ConnectedBasisConfigurations
using cuDNN                    # Lux needs it alongside CUDA, or gpu_device() silently returns a CPU
using DifferentiationInterface
using LinearAlgebra
using Functors: fmap
using Lux
using MLDataDevices: AbstractGPUDevice
using NQSAnsatze
using NQSCore
using NQSOptimisers
using OperatorAlgebra
using Printf
using Random
using Statistics: mean
using SymBasis
using Zygote

"""
Whether there is a device this run can actually use, and what to say if not.

`CUDA.functional()` is necessary but not sufficient: it can return true while the first attempt
to touch the device fails — because the toolkit is too new for the card, or because the card
itself is reporting a fault. Both have to be caught here rather than surfacing as a stacktrace
from the middle of a benchmark.
"""
function device_status()
    CUDA.functional() || return (false, "CUDA is not available on this machine")
    try
        dev = CUDA.device()
        return (true, "$(CUDA.name(dev)) (compute capability $(CUDA.capability(dev)))" *
                      "  |  CUDA $(CUDA.runtime_version())")
    catch err
        return (false, "a device is present but unusable — " *
                       first(split(sprint(showerror, err), '\n')))
    end
end

let (ok, message) = device_status()
    println()
    println(ok ? "device: $message" : "No usable CUDA device: $message")
    if !ok
        println()
        println("The CPU benchmarks need none of this and run anywhere:")
        println()
        println("    julia --project=lib/NQSCore/benchmark lib/NQSCore/benchmark/benchmarks.jl")
        println()
        if CUDA.functional()
            # A device was found and then refused. That is worth distinguishing from having no
            # GPU at all, because it is usually fixable and always specific.
            println("Two things commonly cause this, and they look alike from here:")
            println()
            println("  * The toolkit is newer than the card. CUDA 13 dropped compute")
            println("    capability below 7.5, so a V100 (7.0), a P100 (6.0) or a T4 in a")
            println("    CUDA 13 environment will refuse. Select an older runtime:")
            println("        CUDA.set_runtime_version!(v\"12.6\")")
            println("    or point at a CUDA 12 module with local_toolkit=true.")
            println()
            println("  * The card is reporting a fault — ERROR_ECC_UNCORRECTABLE and friends")
            println("    are hardware, not software. Check `nvidia-smi -q | grep -i -A3 ecc`;")
            println("    a retained uncorrectable error usually needs a GPU reset")
            println("    (`nvidia-smi -r`, root) or simply another node.")
            println()
            exit(1)
        end
        exit(0)
    end
end

BenchmarkTools.DEFAULT_PARAMETERS.seconds = 3

const BACKEND = AutoZygote()

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

report(label, t) = @printf("  %-46s %12s\n", label, BenchmarkTools.prettytime(t * 1e9))

"""Time `f`, reporting rather than raising if the device path is not there yet."""
function timed(label, f)
    try
        t = @belapsed $f()
        report(label, t)
        return t
    catch err
        @printf("  %-46s %12s\n", label, "FAILED")
        println("      ", first(split(sprint(showerror, err), '\n')))
        return nothing
    end
end

speedup(cpu, gpu) = (cpu === nothing || gpu === nothing) ? nothing : cpu / gpu
function compare(label, cpu, gpu)
    s = speedup(cpu, gpu)
    s === nothing ? println("  $label: not comparable") :
    @printf("  %-46s %11.1fx\n", label, s)
end

# ============================================================ the network, CPU versus device

const NSITES = 12
const ALPHA = 4

let spec = Spin(1 // 2), nsites = NSITES
    b = basis(dof_object(spec), nsites)
    H = compile(tfi(nsites; h_x=0.9, h_z=0.1))
    model = RBM(nsites, ALPHA)
    a = LuxAnsatz(model, spec, nsites; rng=Xoshiro(0))
    θ_cpu = init_parameters(a, Xoshiro(0))

    # `gpu_device()` falls back to a CPU device with only a warning when its trigger packages
    # are missing, and every "CUDA" measurement below would then quietly run on the host and
    # report a speedup of one. A silently wrong benchmark is worse than none, so stop here.
    dev = gpu_device()
    if !(dev isa AbstractGPUDevice)
        println()
        println("Lux resolved a $(typeof(dev)) rather than a GPU device, so every measurement")
        println("below would run on the host and compare it against itself.")
        println()
        println("CUDA alone is not enough for Lux — cuDNN has to be loaded too. It is declared")
        println("in this environment, so this usually means it failed to install:")
        println()
        println("    julia --project=lib/NQSCore/benchmark/gpu -e 'import Pkg; Pkg.instantiate()'")
        println()
        exit(1)
    end

    println("\nmodel: RBM($nsites, $ALPHA) over $(length(b.states)) configurations, ",
            "$(sum(length, values(θ_cpu))) parameters")

    println("\nforward pass")
    x = configurations(spec, b.states, nsites)
    t_cpu = timed("log_amplitude, CPU", () -> log_amplitude(a, θ_cpu, x))
    θ_gpu = nothing
    t_gpu = try
        # `fmap(CuArray, ...)`, not `gpu_device()`: the latter demotes ComplexF64 to
        # ComplexF32, which would make every comparison below a precision comparison as much
        # as a hardware one. Single precision is the faster way to run on a GPU and worth
        # measuring — but as a separate number, not silently folded into this one.
        θ_gpu = fmap(CuArray, θ_cpu)
        timed("log_amplitude, CUDA", () -> CUDA.@sync log_amplitude(a, θ_gpu, x))
    catch err
        println("      moving parameters to the device failed: ",
                first(split(sprint(showerror, err), '\n')))
        nothing
    end
    compare("speedup", t_cpu, t_gpu)

    θ_32 = fmap(y -> CuArray(ComplexF32.(y)), θ_cpu)
    t_32 = timed("log_amplitude, CUDA (ComplexF32)",
                 () -> CUDA.@sync log_amplitude(a, θ_32, x))
    compare("speedup, single precision", t_cpu, t_32)

    println("\nlocal energy and gradient")
    vs_cpu = FullSumState(a, θ_cpu; backend=BACKEND, basis=b)
    e_cpu = timed("expect, CPU", () -> expect(vs_cpu, H))
    g_cpu = timed("expect_and_grad, CPU", () -> expect_and_grad(vs_cpu, H))
    if θ_gpu !== nothing
        vs_gpu = FullSumState(a, θ_gpu; backend=BACKEND, basis=b)
        e_gpu = timed("expect, CUDA", () -> CUDA.@sync expect(vs_gpu, H))
        g_gpu = timed("expect_and_grad, CUDA", () -> CUDA.@sync expect_and_grad(vs_gpu, H))
        compare("expect speedup", e_cpu, e_gpu)
        compare("expect_and_grad speedup", g_cpu, g_gpu)

        # Agreement matters more than speed: a device result that disagrees is not a result.
        try
            # Tolerance follows the arithmetic: comparing a double-precision host result
            # against a single-precision device one at 1e-8 fails on precision alone and says
            # nothing about correctness.
            tol = real(eltype(θ_gpu.weight)) === Float32 ? 1e-5 : 1e-10
            ec = real(expect(vs_cpu, H).mean)
            eg = real(expect(vs_gpu, H).mean)
            @printf("  CPU %.12f  vs  CUDA %.12f   (relative %.2e, tolerance %.0e)\n",
                    ec, eg, abs(ec - eg) / abs(ec), tol)
            println("  CPU and CUDA energies agree: ", isapprox(ec, eg; rtol=tol))
        catch err
            println("  energy comparison failed: ", first(split(sprint(showerror, err), '\n')))
        end
    end

    # ===================================================== the host/device boundary itself

    println("\ntransfer, per optimization step")
    res = connected_padded(H, b.states)
    height, batch = size(res.configs)
    xp = configurations(spec, vec(res.configs), nsites)

    t_out = timed("samples out  ($(size(x, 1))x$(size(x, 2)))",
                  () -> CUDA.@sync CuArray(Float64.(x)))
    t_back = timed("connected configs back  ($(size(xp, 1))x$(size(xp, 2)))",
                   () -> CUDA.@sync CuArray(Float64.(xp)))
    t_mels = timed("matrix elements back  ($(height)x$(batch))",
                   () -> CUDA.@sync CuArray(res.mels))

    total = sum(t for t in (t_out, t_back, t_mels) if t !== nothing; init=0.0)
    @printf("  %-46s %12s\n", "total per step", BenchmarkTools.prettytime(total * 1e9))
    if t_gpu !== nothing
        @printf("  %-46s %11.1f%%\n", "as a fraction of one forward pass", 100 * total / t_gpu)
    end

    counts = res.counts
    @printf("\n  %-46s %12d\n", "max_conn (static bound)", max_conn_size(H))
    @printf("  %-46s %12.2f\n", "mean connections per sample", mean(counts))
    @printf("  %-46s %12.2f\n", "padding ratio max_conn / mean(n_conn)",
            max_conn_size(H) / mean(counts))
    println("""
      A ratio near 1 means a device-side connected-configuration kernel would waste nothing;
      well above 1 means most of the transfer above, and most of what such a kernel would
      evaluate, is padding.""")
end

# ================================================= the linear algebra behind the SR solve

println("\nstochastic reconfiguration, solve only")
let n = 2048, p = 4096
    X_cpu = randn(Xoshiro(0), n, p)
    g_cpu = randn(Xoshiro(1), p)
    cg = ConjugateGradientSolver(; tol=1e-8, maxiter=200)

    t_form_cpu = timed("form XᵀX, CPU  ($(p)x$(p))", () -> transpose(X_cpu) * X_cpu)
    t_mf_cpu = timed("matrix-free solve, CPU",
                     () -> solve(cg, QuantumGeometricTensor(X_cpu), g_cpu, 1e-3))

    t_form_gpu, t_mf_gpu = try
        X_gpu = CuArray(X_cpu)
        g_gpu = CuArray(g_cpu)
        (timed("form XᵀX, CUDA", () -> CUDA.@sync transpose(X_gpu) * X_gpu),
         timed("matrix-free solve, CUDA",
               () -> CUDA.@sync solve(cg, QuantumGeometricTensor(X_gpu), g_gpu, 1e-3)))
    catch err
        println("      device solve unavailable: ", first(split(sprint(showerror, err), '\n')))
        (nothing, nothing)
    end
    compare("forming the tensor: speedup", t_form_cpu, t_form_gpu)
    compare("matrix-free solve: speedup", t_mf_cpu, t_mf_gpu)

    @printf("\n  the tensor this avoids allocating: %.0f MB\n", p * p * 8 / 2^20)
end

println()
