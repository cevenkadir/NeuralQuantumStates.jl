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

`expect, CUDA` against the **transfer** section below it. Connected configurations used to be
computed on the host — a compiled operator is a nested structure that cannot be uploaded as it
stands — so every step shipped the samples out and the connected configurations and matrix
elements back, the return leg being `max_conn` times larger than the outbound one. That
asymmetry is what motivated `FlatOperator` and the device kernel, and the transfer section is
kept as the standing measurement of what they removed: those numbers are no longer part of an
`expect`, and the comparison says whether the work paid for itself.

The `max_conn / mean(n_conn)` ratio printed alongside is now the device path's own overhead,
since it keeps the operator's static bound rather than trimming to the largest connection count
it saw. At 1 the padding is free; well above 1 it means most of what the kernel and the network
evaluate is padding.

Each measurement is guarded: a failure prints and the run continues, because most of the device
paths here have never executed on hardware and one broken step should not hide the rest.
"""

using BenchmarkTools
using CUDA
using ConnectedBasisConfigurations
using cuDNN                    # Lux needs it alongside CUDA, or gpu_device() silently returns a CPU
using DifferentiationInterface
using KernelAbstractions
using LinearAlgebra
using Functors: fmap
using Lux
using MLDataDevices: AbstractGPUDevice
using NQSAnsatze
using NQSCore
using NQSOptimisers
using NQSSamplers
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
    e_gpu = nothing              # referenced by the transfer section, which runs either way
    if θ_gpu !== nothing
        vs_gpu = FullSumState(a, θ_gpu; backend=BACKEND, basis=b)
        # With parameters on the device and KernelAbstractions loaded, the connected
        # configurations are computed there too, so this covers the device kernel and not just
        # the network. The operator is uploaded on each call in this form.
        e_gpu = timed("expect, CUDA", () -> CUDA.@sync expect(vs_gpu, H))
        g_gpu = timed("expect_and_grad, CUDA", () -> CUDA.@sync expect_and_grad(vs_gpu, H))
        compare("expect speedup", e_cpu, e_gpu)
        compare("expect_and_grad speedup", g_cpu, g_gpu)

        # The loop-friendly form: upload the operator once instead of on every call. Seven
        # small transfers is not much against a millisecond, but an optimization run pays them
        # once per step for as many steps as it takes, and the fix is one line at the call site.
        H_dev = nothing
        e_res = try
            H_dev = to_backend(flatten(H), CUDABackend())
            timed("expect, CUDA (operator uploaded once)",
                  () -> CUDA.@sync expect(vs_gpu, H_dev))
        catch err
            println("      uploading the operator failed: ",
                    first(split(sprint(showerror, err), '\n')))
            nothing
        end
        compare("expect speedup, operator uploaded once", e_cpu, e_res)

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
            if H_dev !== nothing
                # The uploaded operator is a different code path through `_connections`, so it
                # gets its own comparison rather than being assumed equivalent.
                er = real(expect(vs_gpu, H_dev).mean)
                println("  the uploaded operator gives the same energy: ",
                        isapprox(ec, er; rtol=tol))
            end
        catch err
            println("  energy comparison failed: ", first(split(sprint(showerror, err), '\n')))
        end
    end

    # ============================================ where `expect_and_grad` actually spends

    # The gradient costs 16% on top of `expect` on a host and far more than that on a device,
    # and that asymmetry is the question this section exists to answer. It should not be there:
    # `energy_gradient` differentiates a scalar loss over the *samples*, the same batch the
    # forward pass above covers, so one reverse pass ought to be a small multiple of one
    # forward pass. Anything much larger is overhead — a host round trip inside the
    # differentiated function, or a parameter flattening that does not like device memory —
    # rather than arithmetic, and the breakdown says which.
    if θ_gpu !== nothing
        println("\ngradient, broken down")
        vs_gpu = FullSumState(a, θ_gpu; backend=BACKEND, basis=b)
        xs = NQSCore.configurations_of(vs_gpu, b.states)
        # Every diagnostic below that is meant to show what a *real* batch costs takes this
        # one, never `xs`. `configurations_of` decides `xs`'s element type, so a measurement
        # written against it silently changes what it measures the moment that decision
        # changes — which is how `complex matmul, reverse` once appeared to regress twelvefold
        # while nothing had regressed at all.
        xs_real = Float64.(real.(xs))
        @printf("  %-46s %12s\n", "sample batch is on the device",
                string(!(xs isa Array)))
        @printf("  %-46s %12s\n", "sample batch element type", string(eltype(xs)))

        lg = timed("log_amplitude on the samples (forward)",
                   () -> CUDA.@sync log_amplitude(a, θ_gpu, xs))
        logψ = log_amplitude(a, θ_gpu, xs)
        p = NQSCore.born_probabilities(logψ)
        E = NQSCore.local_energy(vs_gpu, H, b.states)

        eg = timed("energy_gradient alone (reverse)",
                   () -> CUDA.@sync NQSCore.energy_gradient(
                       a, θ_gpu, xs, E, p; backend=BACKEND))
        le = timed("local_energy alone", () -> CUDA.@sync NQSCore.local_energy(vs_gpu, H, b.states))
        if lg !== nothing && eg !== nothing
            @printf("  %-46s %11.1fx\n", "reverse / forward", eg / lg)
            println("""
      A reverse pass is normally two to three times a forward one. Well above that is
      overhead in the gradient path rather than the arithmetic of differentiating.""")
        end
        if le !== nothing && eg !== nothing && g_gpu !== nothing
            @printf("  %-46s %11.1f%%\n", "the two together, of expect_and_grad",
                    100 * (le + eg) / g_gpu)
        end

        # The layer itself is not the problem. Measured on this device, the reverse of the
        # `logtwocosh` reduction is 4.6x its forward and the complex matmul's is 2.9x — 529 us
        # between them — and `flatten_parameters` with its `restore` costs 85 us standing
        # still. That is 614 us against a 4.17 ms `energy_gradient`, so the cost is in neither,
        # and guessing a third time is worse than bisecting.
        #
        # These two stay as a standing check that the layer's reverse pass remains cheap: they
        # are what would move if a Zygote or CUDA upgrade dropped the complex broadcast onto a
        # slower path, and that would otherwise show up only as a slower run with no cause.
        #
        # The reduction's 380 us is not `logtwocosh` being expensive to differentiate. Giving
        # it an exact rule — `@scalar_rule logtwocosh(z) tanh(z)`, verified to reproduce
        # Zygote's own gradient to 1.4e-14 — was measured here and came back at **1.0x**,
        # 383.9 us against 382.1 us. The cost is Zygote's per-element pullback machinery for a
        # complex broadcast, which a scalar rule does not remove; it only replaces what runs
        # inside each closure. Do not try it again without changing that structure.
        println("\n  the two halves of the layer, and the flattening")
        nh = ALPHA * NSITES
        try
            u = CuArray(randn(Xoshiro(0), ComplexF64, nh, length(b.states)))
            V = CuArray(randn(Xoshiro(1), ComplexF64, nh, NSITES))
            f_red(z) = sum(real, sum(NQSAnsatze.logtwocosh, z; dims=1))
            f_mm(M) = sum(real, M * xs_real)

            for (label, f, arg) in (("logtwocosh reduction", f_red, u),
                                    ("complex matmul", f_mm, V))
                fw = timed("  $label, forward", () -> CUDA.@sync f(arg))
                rv = timed("  $label, reverse",
                           () -> CUDA.@sync Zygote.gradient(f, arg))
                fw === nothing || rv === nothing ||
                    @printf("  %-46s %11.1fx\n", "    reverse / forward", rv / fw)
            end
        catch err
            println("      the split measurement failed: ",
                    first(split(sprint(showerror, err), '\n')))
        end

        # The bisection. `energy_gradient` is a ladder of four wrappers around the model, and
        # each rung below removes exactly one of them, so a jump between two consecutive rungs
        # names the wrapper responsible rather than suggesting one.
        #
        #   split      the real/imag reparameterization, `v[1:n] .+ im .* v[n+1:2n]`
        #   restore    rebuilding a NamedTuple from a ComponentArray, differentiated
        #   cotangent  `2 sum(real(conj(c) * psi))`
        #   model      `log_amplitude` itself, whose forward is 109 us
        println("\n  bisecting energy_gradient")
        try
            flat, restore = NQSCore.flatten_parameters(θ_gpu)
            n = length(flat)
            c = NQSCore._gradient_cotangent(E, p)
            split = vcat(real.(flat), imag.(flat))

            function full(v)
                q = @views v[1:n] .+ im .* v[(n+1):(2n)]
                return NQSCore._gradient_loss(a, restore(q), xs, c)
            end
            no_split(q) = NQSCore._gradient_loss(a, restore(q), xs, c)
            no_restore(θ) = NQSCore._gradient_loss(a, θ, xs, c)
            model_only(θ) = sum(real, log_amplitude(a, θ, xs))

            t_fwd = timed("  the whole closure, forward only", () -> CUDA.@sync full(split))
            t_full = timed("  + split + restore + cotangent + model",
                           () -> CUDA.@sync Zygote.gradient(full, split))
            t_nosplit = timed("  - split", () -> CUDA.@sync Zygote.gradient(no_split, flat))
            t_norestore = timed("  - split - restore",
                                () -> CUDA.@sync Zygote.gradient(no_restore, θ_gpu))
            t_model = timed("  - split - restore - cotangent",
                            () -> CUDA.@sync Zygote.gradient(model_only, θ_gpu))

            rungs = (("the real/imag split", t_full, t_nosplit),
                     ("restore, differentiated", t_nosplit, t_norestore),
                     ("the cotangent", t_norestore, t_model))
            println("\n  what each rung costs")
            for (label, upper, lower) in rungs
                upper === nothing || lower === nothing ||
                    @printf("  %-46s %12s\n", label,
                            BenchmarkTools.prettytime((upper - lower) * 1e9))
            end
            t_model === nothing ||
                @printf("  %-46s %12s\n", "the model's own reverse pass",
                        BenchmarkTools.prettytime(t_model * 1e9))
            t_fwd === nothing || t_full === nothing ||
                @printf("  %-46s %11.1fx\n", "closure reverse / closure forward",
                        t_full / t_fwd)
        catch err
            println("      the bisection failed: ", first(split(sprint(showerror, err), '\n')))
        end

        # The ladder above put 3.516 ms of a 4.167 ms gradient inside the model, and the two
        # operations the model is made of measured 530 us between them in isolation. Composing
        # them therefore costs about seven times what they cost apart, and only on a device:
        # on a host the isolated pieces add up to the whole (3.41 ms and 1.83 ms against
        # 5.85 ms). So the second ladder walks into the model the same way the first walked
        # into `energy_gradient`.
        println("\n  bisecting the model")
        try
            layer, st = a.model, a.states
            lt = NQSAnsatze.logtwocosh
            # `xs` is whatever `configurations_of` decided on, which is the point of the top
            # rung. The rungs below take the deliberately real batch, so that they keep showing
            # the cost the type choice avoids rather than quietly agreeing with it once the
            # choice is right.
            r_ansatz(θ) = sum(real, log_amplitude(a, θ, xs))
            r_apply(θ) = sum(real, first(Lux.apply(layer, xs_real, θ, st)))
            r_layer(θ) = sum(real, first(layer(xs_real, θ, st)))
            r_body(θ) = sum(real, reshape(θ.visible, 1, :) * xs_real .+
                                  sum(lt, θ.weight * xs_real .+ θ.hidden; dims=1))
            r_hidden(θ) = sum(real, sum(lt, θ.weight * xs_real .+ θ.hidden; dims=1))
            r_mm(θ) = sum(real, θ.weight * xs_real .+ θ.hidden)

            for (label, f) in (("through LuxAnsatz, batch as built", r_ansatz),
                               ("through Lux.apply, real batch", r_apply),
                               ("the layer directly, real batch", r_layer),
                               ("the layer body inlined, real batch", r_body),
                               ("...without the visible term", r_hidden),
                               ("...matmul and bias only", r_mm))
                timed("  $label", () -> CUDA.@sync Zygote.gradient(f, θ_gpu))
            end

            # The one type question in the whole path. Parameters are `ComplexF64` and the
            # batch is `Float64`, because `real(eltype(θ))` is what a configuration naturally
            # is. cuBLAS has no mixed complex-real `gemm`, so such a product falls to a generic
            # kernel — the same class of mistake as the `Adjoint` wrapper already commented in
            # the layer, and invisible on a host where generic matmul is merely somewhat
            # slower rather than dramatically so.
            #
            # The *direction* is what makes it expensive, and it is worth being explicit about
            # because measuring the wrong one says the opposite. The forward product is
            # `(48x12)(12x4096)`: two and a half million outputs over a reduction of length 12,
            # so a generic kernel has plenty to be parallel over and is no slower than cuBLAS.
            # Its pullback is `(48x4096)(4096x12)`: 576 outputs over a reduction of length
            # 4096, which a generic kernel runs on 576 threads that each loop four thousand
            # times. Both are timed below, and only the second should differ.
            println("\n  does the matmul reach cuBLAS?")
            xc = ComplexF64.(xs_real)
            Δ = CuArray(randn(Xoshiro(2), ComplexF64, nh, length(b.states)))

            f_mixed = timed("  forward, ComplexF64 x Float64",
                            () -> CUDA.@sync θ_gpu.weight * xs_real)
            f_same = timed("  forward, ComplexF64 x ComplexF64", () -> CUDA.@sync θ_gpu.weight * xc)
            compare("  converting the batch would be", f_mixed, f_same)

            r_mixed = timed("  pullback, ComplexF64 x adjoint Float64",
                            () -> CUDA.@sync Δ * xs_real')
            r_same = timed("  pullback, ComplexF64 x adjoint ComplexF64",
                           () -> CUDA.@sync Δ * xc')
            compare("  converting the batch would be", r_mixed, r_same)

            # And the same question asked end to end, which is the one that decides it: the
            # whole layer differentiated against a complex batch, against the 3.5 ms it costs
            # against a real one. If this is small, the fix is one line in `_input_type` —
            # promote the batch once, on a small array, rather than forcing a promotion inside
            # every BLAS call on both sides of the derivative.
            r_layer_c(θ) = sum(real, first(layer(xc, θ, st)))
            timed("  the layer differentiated, ComplexF64 batch",
                  () -> CUDA.@sync Zygote.gradient(r_layer_c, θ_gpu))
        catch err
            println("      the model bisection failed: ",
                    first(split(sprint(showerror, err), '\n')))
        end
    end

    # ===================================================== the host/device boundary itself

    # This section measured the case for the device kernel, and now measures what it removed:
    # with the connections computed on the device, none of these transfers happens during an
    # `expect`. They are kept because they are the standing answer to "was that worth it" —
    # the numbers below, against the `expect, CUDA` above, are the whole argument.
    println("\ntransfer per step, as it was before the device kernel")
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
    # The samples still cross: `expect` unpacks them on the host to get `log ψ(s)`, and only the
    # connected configurations — the array `max_conn` times larger — now stay on the device.
    if e_gpu !== nothing
        @printf("  %-46s %11.1f%%\n", "against one device expect today", 100 * total / e_gpu)
    end

    counts = res.counts
    @printf("\n  %-46s %12d\n", "max_conn (static bound)", max_conn_size(H))
    @printf("  %-46s %12.2f\n", "mean connections per sample", mean(counts))
    @printf("  %-46s %12.2f\n", "padding ratio max_conn / mean(n_conn)",
            max_conn_size(H) / mean(counts))
    println("""
      A ratio near 1 means the device kernel wastes nothing; well above 1 means most of what
      it evaluates, and most of what the transfers above carried, is padding. The device path
      keeps the full static bound rather than trimming to the largest count it saw, so this
      ratio is exactly its overhead.""")
end

# ======================================== connected configurations, host versus device kernel

# The measurement this whole exercise was for. Computing connections on the host costs the
# kernel plus the transfer of its results; computing them on the device costs a kernel launch
# and nothing else. The comparison is only meaningful if the two agree, so that is checked
# first and the timings are skipped if they do not.
println("\nconnected configurations")

let spec = Spin(1 // 2), nsites = NSITES
    b = basis(dof_object(spec), nsites)
    states = b.states
    H = tfi(nsites; h_x=0.9, h_z=0.1)
    flat = flatten(H)
    h, n = max_conn_size(flat), length(states)

    host_c = Matrix{eltype(states)}(undef, h, n)
    host_m = Matrix{Float64}(undef, h, n)
    host_k = Vector{Int}(undef, n)
    connected_padded!(host_c, host_m, host_k, flat, states)

    # The kernel lives in an extension, so the first question is whether it is there at all.
    # A missing extension and a broken kernel look alike from the error message otherwise.
    let ext = Base.get_extension(
            ConnectedBasisConfigurations, :ConnectedBasisConfigurationsKernelAbstractionsExt
        )
        println("  kernel extension loaded: ", ext !== nothing)
        ext === nothing && println("      `using KernelAbstractions` should activate it; if it " *
                                   "does not, the extension failed to precompile")
    end

    ok = false
    dev_flat = nothing
    dev_states = dev_c = dev_m = dev_k = nothing
    try
        backend = CUDABackend()
        dev_flat = to_backend(flat, backend)
        dev_states = CuArray(states)
        dev_c = CUDA.similar(dev_states, h, n)
        dev_m = CuArray{Float64}(undef, h, n)
        dev_k = CuArray{Int}(undef, n)
        connected_padded!(dev_c, dev_m, dev_k, dev_flat, dev_states, backend)
        ok = Array(dev_k) == host_k && Array(dev_c) == host_c && Array(dev_m) == host_m
        println("  device kernel matches the host kernel: ", ok)
    catch err
        println("  device kernel FAILED: ", first(split(sprint(showerror, err), '\n')))
    end

    # The unpacking that turns packed states into the network's input. On the host it produces
    # `Rational`s, which no accelerator can hold, so the float conversion is a second full-size
    # array before anything is transferred; the device kernel writes the float directly.
    host_x = configurations(spec, vec(host_c), nsites)
    dev_x = dev_values = ok_x = nothing
    try
        dev_values = CuArray(collect(Float64, local_values(spec)))
        dev_x = CuArray{Float64}(undef, nsites, h * n)
        configurations!(dev_x, dev_values, dev_c, nsites, CUDABackend())
        ok_x = Array(dev_x) == Float64.(host_x)
        println("  device unpacking matches the host unpacking: ", ok_x)
    catch err
        println("  device unpacking FAILED: ", first(split(sprint(showerror, err), '\n')))
    end

    t_host = timed("connections on host (kernel only)",
                   () -> connected_padded!(host_c, host_m, host_k, flat, states))
    t_move = timed("...plus unpacking and moving the results",
                   () -> CUDA.@sync (CuArray(Float64.(configurations(spec, vec(host_c), nsites)));
                                     CuArray(host_m)))
    if ok
        t_dev = timed("connections on device",
                      () -> CUDA.@sync connected_padded!(
                          dev_c, dev_m, dev_k, dev_flat, dev_states, CUDABackend()))
        t_unpack = ok_x === true ?
                   timed("...plus unpacking on device",
                         () -> CUDA.@sync begin
                             connected_padded!(dev_c, dev_m, dev_k, dev_flat, dev_states,
                                               CUDABackend())
                             configurations!(dev_x, dev_values, dev_c, nsites, CUDABackend())
                         end) : nothing
        if t_host !== nothing && t_move !== nothing && t_dev !== nothing
            compare("against host kernel alone", t_host, t_dev)
            compare("against host kernel plus transfer", t_host + t_move, t_dev)
            t_unpack === nothing ||
                compare("whole path, host versus device", t_host + t_move, t_unpack)
        end
    end
end

# ====================================================== sampling, which is the real loop

# `FullSumState` is a testing instrument. A production run samples, and a Metropolis sampler
# evaluates the ansatz **once per sweep step** — so with parameters on a device and samples on
# the host, every step is a host/device round trip. Those are latency-bound, they do not
# amortize, and there are `burn_in + n_samples * thinning` of them. This section is here to
# find out whether that swamps the speedups measured above.
println("\nsampling")

let spec = Spin(1 // 2), nsites = NSITES
    b = basis(dof_object(spec), nsites)
    H = compile(tfi(nsites; h_x=0.9, h_z=0.1))
    a = LuxAnsatz(RBM(nsites, ALPHA), spec, nsites; rng=Xoshiro(0))
    θ_cpu = init_parameters(a, Xoshiro(0))
    θ_gpu = fmap(CuArray, θ_cpu)

    n_chains, n_samples, burn_in = 8, 200, 50
    steps = burn_in + n_samples
    starts = random_configurations(spec, nsites, n_chains, Xoshiro(1))
    sampler = MetropolisSampler(LocalRule(), starts;
        n_chains=n_chains, n_samples=n_samples, burn_in=burn_in)

    @printf("  %d chains x %d samples, burn-in %d — %d ansatz evaluations of %d configurations\n",
            n_chains, n_samples, burn_in, steps, n_chains)

    s_cpu = timed("Metropolis, parameters on host",
                  () -> NQSCore.sample(sampler, a, θ_cpu, Xoshiro(2)))
    s_gpu = timed("Metropolis, parameters on device",
                  () -> CUDA.@sync NQSCore.sample(sampler, a, θ_gpu, Xoshiro(2)))
    compare("speedup", s_cpu, s_gpu)

    if s_cpu !== nothing && s_gpu !== nothing
        @printf("  %-46s %11.1f µs\n", "per sweep step, host", 1e6 * s_cpu / steps)
        @printf("  %-46s %11.1f µs\n", "per sweep step, device", 1e6 * s_gpu / steps)
        if s_gpu > s_cpu
            println("""
      The device is slower, and that is the point: a sweep evaluates the ansatz on only
      $(n_chains) configurations, which is far too little work to cover a round trip. The
      fix is not a faster kernel but a sampler that keeps the chains on the device.""")
        end
    end

    # Is the sampler's problem its design, or its width? A sweep costs one ansatz evaluation
    # however many chains it advances, and on a device that evaluation is a handful of kernel
    # launches whose latency does not depend on how much data they carry. If that is what the
    # per-step cost is made of, then widening the sweep should leave the per-step time almost
    # unchanged while dividing the per-sample time by the chain count — and the fix is a
    # configuration rather than a rewrite. If instead per-step time grows with the width, the
    # cost is real work and only a device-resident chain would help.
    println("\n  per sample, by sweep width")
    @printf("  %-14s %14s %14s %10s\n", "chains", "host", "device", "speedup")
    for nc in (8, 64, 512, 2048)
        s = MetropolisSampler(LocalRule(), random_configurations(spec, nsites, nc, Xoshiro(1));
            n_chains=nc, n_samples=10, burn_in=5)
        drawn = 10 * nc
        h = try
            @belapsed NQSCore.sample($s, $a, $θ_cpu, Xoshiro(2))
        catch
            nothing
        end
        d = try
            @belapsed CUDA.@sync NQSCore.sample($s, $a, $θ_gpu, Xoshiro(2))
        catch
            nothing
        end
        if h === nothing || d === nothing
            @printf("  %-14d %14s %14s %10s\n", nc, "-", "-", "-")
        else
            @printf("  %-14d %11.3f ns %11.3f ns %9.1fx\n",
                    nc, 1e9 * h / drawn, 1e9 * d / drawn, h / d)
        end
    end

    # For contrast: one batched draw, where the ansatz sees the whole basis at once and there
    # is exactly one transfer rather than one per step.
    exact = ExactSampler(b, n_chains * n_samples)
    e_cpu = timed("ExactSampler, parameters on host",
                  () -> NQSCore.sample(exact, a, θ_cpu, Xoshiro(2)))
    e_gpu = timed("ExactSampler, parameters on device",
                  () -> CUDA.@sync NQSCore.sample(exact, a, θ_gpu, Xoshiro(2)))
    compare("speedup", e_cpu, e_gpu)

    # And the whole loop, which is what a run actually pays.
    m_cpu = timed("MCState expect, host",
                  () -> expect(MCState(a, θ_cpu, sampler; backend=BACKEND, rng=Xoshiro(3)), H))
    m_gpu = timed("MCState expect, device",
                  () -> CUDA.@sync expect(
                      MCState(a, θ_gpu, sampler; backend=BACKEND, rng=Xoshiro(3)), H))
    compare("speedup", m_cpu, m_gpu)
end

# ================================================= the linear algebra behind the SR solve

println("\nstochastic reconfiguration, solve only")
let n = 2048, p = 4096
    X_cpu = randn(Xoshiro(0), n, p)
    g_cpu = randn(Xoshiro(1), p)
    # A shift and tolerance that actually converge. With a tighter pair the solver hits its
    # iteration cap on both sides and the comparison measures non-convergence rather than a
    # solve — and a real run regularizes for exactly this reason.
    shift = 1e-2
    cg = ConjugateGradientSolver(; tol=1e-6, maxiter=100)

    t_form_cpu = timed("form XᵀX, CPU  ($(p)x$(p))", () -> transpose(X_cpu) * X_cpu)
    t_mf_cpu = timed("matrix-free solve, CPU",
                     () -> solve(cg, QuantumGeometricTensor(X_cpu), g_cpu, shift))

    t_form_gpu, t_mf_gpu = try
        X_gpu = CuArray(X_cpu)
        g_gpu = CuArray(g_cpu)
        (timed("form XᵀX, CUDA", () -> CUDA.@sync transpose(X_gpu) * X_gpu),
         timed("matrix-free solve, CUDA",
               () -> CUDA.@sync solve(cg, QuantumGeometricTensor(X_gpu), g_gpu, shift)))
    catch err
        println("      device solve unavailable: ", first(split(sprint(showerror, err), '\n')))
        (nothing, nothing)
    end
    compare("forming the tensor: speedup", t_form_cpu, t_form_gpu)
    compare("matrix-free solve: speedup", t_mf_cpu, t_mf_gpu)

    @printf("\n  the tensor this avoids allocating: %.0f MB\n", p * p * 8 / 2^20)
end

println()
