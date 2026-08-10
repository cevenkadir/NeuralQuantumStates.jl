"""
The whole variational step in XLA, against the same step on CUDA.jl.

Run with

```
XLA_REACTANT_GPU_PREALLOCATE=false NQS_REACTANT_BACKEND=gpu \\
  julia --project=lib/NQSCore/benchmark/kernels lib/NQSCore/benchmark/kernels/step_probe.jl
```

# The number this exists for

`expect_and_grad` under CUDA.jl is 5.245 ms on an RTX 5000 Ada — `local_energy` 2.609 ms plus
`energy_gradient` 2.471 ms, 96.9% between them. Reactant does the gradient alone in 355 µs, seven
times faster, but only if the connected configurations are already where it can reach them; left
on the host they cost 3.8 ms and swallow the win whole.

Everything needed to close that is now in place. The kernel runs under Reactant over the raw
integer a `BaseInt` wraps; `complex.` replaced the `Complex{Bool}` the reparameterization used to
build; and the parameters need no flattening, because differentiating the `NamedTuple` directly
gives the same gradient — Zygote to 0.0 and Reactant to 2.206e-15 across all three rungs.

So this measures the step end to end, in three compiled regions, and says whether the projection
of roughly 0.8–1.3 ms survives contact.

# Why three regions and not one

The kernel does not raise — `cannot raise op to stablehlo` on its `scf.for`, its loops being
data-dependent — so XLA cannot fuse across it. That costs fusion, not the ability to run, and
nothing forces the step into a single region: only the middle region's output crosses to the
gradient. The three are timed separately and together, so the cost of not fusing is visible
rather than assumed.

    prepare   the kernel, then both batches unpacked      (no derivative)
    energy    two network evaluations, the reduction,
              the Born weights and the cotangent          (no derivative)
    gradient  the one region a derivative runs through
"""

using BenchmarkTools
using CUDA
using ConnectedBasisConfigurations
using DifferentiationInterface
using Enzyme
using Functors: fmap
using KernelAbstractions
using LinearAlgebra
using Lux
using NQSAnsatze
using NQSCore
using OperatorAlgebra
using Printf
using Random
using Reactant
using SymBasis
using cuDNN

let requested = get(ENV, "NQS_REACTANT_BACKEND", "")
    isempty(requested) || Reactant.set_default_backend(requested)
end

BenchmarkTools.DEFAULT_PARAMETERS.seconds = 2

# Bounded rather than left to the clock: a `ConcretePJRTArray` looks like a few hundred bytes of
# host object to Julia's collector, which has no notion of the XLA memory behind it.
const SAMPLES = 500

println()
println("Reactant ", pkgversion(Reactant), "  |  Enzyme ", pkgversion(Enzyme),
        "  |  CUDA ", pkgversion(CUDA), "  |  Lux ", pkgversion(Lux))
println("CUDA.jl: ", CUDA.functional() ?
        "$(CUDA.name(CUDA.device())), CUDA $(CUDA.runtime_version())" : "not functional")
println("XLA devices: ", try
    length(Reactant.devices())
catch err
    "unavailable — " * first(split(sprint(showerror, err), '\n'))
end)
println("""
Both allocators are live here. XLA_REACTANT_GPU_PREALLOCATE=false is what leaves CUDA.jl room;
without it XLA takes the better part of the card before CUDA.jl allocates anything.""")

report(label, t) = @printf("  %-46s %12s\n", label, BenchmarkTools.prettytime(t * 1e9))
note(label, what) = @printf("  %-46s %12s\n", label, what)

const FAILURES = String[]

function failed(what, err)
    text = sprint(showerror, err)
    lines = split(text, '\n')
    push!(FAILURES, "$what: $(first(lines))")
    println("  ", what, " FAILED")
    for l in lines[1:min(8, length(lines))]
        println("      ", length(l) > 300 ? l[1:prevind(l, 300)] * " …" : l)
    end
    length(lines) > 8 && println("      … ", length(lines) - 8, " more lines")
    return nothing
end

function timed(label, f)
    try
        t = minimum((@benchmark $f() evals = 1 samples = SAMPLES)).time / 1e9
        GC.gc()
        report(label, t)
        return t
    catch err
        failed(strip(label), err)
        return nothing
    end
end

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

# ============================================================================ the regions

# Top-level functions of explicit arguments, which is the form Reactant caches on. Buffers are
# allocated inside with `similar`, the way the kernel tutorial writes one, rather than passed in.

"""Connected configurations, and both batches the network consumes, from packed states."""
function prepare_region(op, states, values_real, values_complex, height::Int, nsites::Int,
                        base::Val)
    n = length(states)
    backend = KernelAbstractions.get_backend(states)

    configs = similar(states, height, n)
    mels = similar(op.vals, height, n)
    counts = similar(states, Int, n)
    ConnectedBasisConfigurations.connected_padded!(
        configs, mels, counts, op, states, backend; base=base
    )

    # The connected configurations are unpacked real and the samples complex, which is what
    # `NQSCore` does and for the reason it does it: only the sample batch is followed by a
    # derivative, and only there does a mixed complex-real product cost anything.
    x_conn = similar(values_real, nsites, height * n)
    ConnectedBasisConfigurations.configurations!(
        x_conn, values_real, configs, nsites, backend; base=base
    )
    xs = similar(values_complex, nsites, n)
    ConnectedBasisConfigurations.configurations!(
        xs, values_complex, states, nsites, backend; base=base
    )
    return xs, x_conn, mels
end

"""Local energies, Born weights and the cotangent — everything the gradient needs but `θ`."""
function energy_region(a, θ, xs, x_conn, mels)
    logψ_s = log_amplitude(a, θ, xs)
    logψ_sp = reshape(log_amplitude(a, θ, x_conn), size(mels))
    # The column reduces without a mask because the padding is inert: a padded slot repeats the
    # sample with a zero matrix element.
    E = vec(sum(mels .* exp.(logψ_sp .- transpose(logψ_s)); dims=1))
    p = NQSCore.born_probabilities(logψ_s)
    return E, p, NQSCore._gradient_cotangent(E, p)
end

"""The scalar whose gradient is the energy gradient, with the parameters last."""
step_loss(a, x, c, θ) = NQSCore._gradient_loss(a, θ, x, c)

"""
The gradient with respect to `θ`, everything else held constant.

`Const` on all but the last, and the parameters as a `NamedTuple` rather than a flat vector:
Enzyme returns a structure gradient, and that is measured to be the same gradient the library's
real/imag split produces — to 2.206e-15 under Reactant — so the flattening the host path does is
not needed here at all.
"""
step_gradient(a, x, c, θ) =
    Enzyme.gradient(Enzyme.Reverse, Const(step_loss), Const(a), Const(x), Const(c), θ)[end]

"""Compile `f(args...)`, reporting what that cost. `sync=true` is the barrier XLA needs."""
function compiled(label, f::F, args) where {F}
    t = @elapsed thunk = Reactant.compile(f, args; sync=true)
    report(label, t)
    return thunk
end

"""A gradient in one comparable vector, whatever shape the engine returned it in."""
host_gradient(g::NamedTuple) = first(NQSCore.flatten_parameters(map(Array, g)))
host_gradient(g) = first(NQSCore.flatten_parameters(g))

const NSITES = 12
const ALPHA = 4

let spec = Spin(1 // 2), nsites = NSITES
    b = basis(dof_object(spec), nsites)
    states = b.states
    H = compile(tfi(nsites; h_x=0.9, h_z=0.1))
    op = flatten(tfi(nsites; h_x=0.9, h_z=0.1))
    a = LuxAnsatz(RBM(nsites, ALPHA), spec, nsites; rng=Xoshiro(0))
    θ = init_parameters(a, Xoshiro(0))

    height, n = max_conn_size(op), length(states)
    S = eltype(states)
    V, B = S.parameters[1], S.parameters[3]
    base = Val(B)

    println("\nmodel: RBM($nsites, $ALPHA) over $n configurations, ",
            "$(sum(length, values(θ))) parameters; max_conn $height")
    note("connected configurations per step", string(height * n))

    # ------------------------------------------------------ the baseline, on CUDA.jl

    println("\nthe step on CUDA.jl, with Zygote")
    t_cuda_step = t_cuda_le = t_cuda_grad = nothing
    E_cuda = ∇_cuda = nothing
    try
        θ_gpu = fmap(CuArray, θ)
        vs = FullSumState(a, θ_gpu; backend=AutoZygote(), basis=b)
        t_cuda_step = timed("expect_and_grad, CUDA.jl", () -> CUDA.@sync expect_and_grad(vs, H))
        t_cuda_le = timed("  local_energy alone",
                          () -> CUDA.@sync NQSCore.local_energy(vs, H, states))

        xs_gpu = NQSCore.configurations_of(vs, states)
        E_g = NQSCore.local_energy(vs, H, states)
        p_g = NQSCore.born_probabilities(log_amplitude(a, θ_gpu, xs_gpu))
        t_cuda_grad = timed("  energy_gradient alone",
                            () -> CUDA.@sync NQSCore.energy_gradient(
                                a, θ_gpu, xs_gpu, E_g, p_g; backend=AutoZygote()))
        E_cuda = Array(E_g)
        ∇_cuda = host_gradient(NQSCore.energy_gradient(
            a, θ_gpu, xs_gpu, E_g, p_g; backend=AutoZygote()))
    catch err
        failed("the CUDA.jl baseline", err)
    end

    # ------------------------------------------------------------- the step in XLA

    println("\nthe step in XLA, three compiled regions")
    t1 = t2 = t3 = nothing
    E_xla = ∇_xla = nothing
    try
        op_ra = Reactant.to_rarray(op)
        states_ra = Reactant.to_rarray(collect(reinterpret(V, states)))
        values_real = Reactant.to_rarray(collect(Float64, local_values(spec)))
        values_complex = Reactant.to_rarray(collect(ComplexF64, local_values(spec)))
        θ_ra = Reactant.to_rarray(θ)
        note("parameter element type", string(eltype(θ_ra.weight)))

        prep_args = (op_ra, states_ra, values_real, values_complex, height, nsites, base)
        prep = compiled("prepare, compiling (once)", prepare_region, prep_args)
        xs, x_conn, mels = prep(prep_args...)
        note("sample batch", string(size(xs)))
        note("connected batch", string(size(x_conn)))
        t1 = timed("prepare (kernel + both unpackings)", () -> prep(prep_args...))

        energy_args = (a, θ_ra, xs, x_conn, mels)
        energy = compiled("energy, compiling (once)", energy_region, energy_args)
        E, p, c = energy(energy_args...)
        t2 = timed("energy (two forwards + reduction)", () -> energy(energy_args...))
        E_xla = Array(E)

        grad_args = (a, xs, c, θ_ra)
        grad = compiled("gradient, compiling (once)", step_gradient, grad_args)
        t3 = timed("gradient", () -> grad(grad_args...))
        ∇_xla = host_gradient(grad(grad_args...))
    catch err
        failed("the XLA step", err)
    end

    # --------------------------------------------------------------------- the verdict

    println("\nagreement")
    if E_cuda === nothing || E_xla === nothing
        println("  one of the two paths produced no energies")
    else
        rel = norm(E_xla .- E_cuda) / norm(E_cuda)
        @printf("  %-46s %12.3e %s\n", "local energies, XLA vs CUDA.jl", rel,
                rel <= 1e-8 ? "agrees" : "DISAGREES")
        rel <= 1e-8 || push!(FAILURES, "local energies disagree between XLA and CUDA.jl")
    end
    if ∇_cuda === nothing || ∇_xla === nothing
        println("  one of the two paths produced no gradient")
    else
        rel = norm(∇_xla .- ∇_cuda) / norm(∇_cuda)
        @printf("  %-46s %12.3e %s\n", "the gradient, XLA vs CUDA.jl", rel,
                rel <= 1e-8 ? "agrees" : "DISAGREES")
        rel <= 1e-8 || push!(FAILURES, "gradients disagree between XLA and CUDA.jl")
    end

    println("\nthe step, end to end")
    total = any(isnothing, (t1, t2, t3)) ? nothing : t1 + t2 + t3
    total === nothing || report("XLA, three regions together", total)
    t_cuda_step === nothing || report("CUDA.jl with Zygote", t_cuda_step)
    if total !== nothing && t_cuda_step !== nothing
        @printf("  %-46s %11.2fx\n", "the whole step", t_cuda_step / total)
        println("""
      Against the projection this probe was written to test: roughly 0.8-1.3 ms, from the
      gradient's measured 7x and the assumption that the network forward over the connected
      configurations goes at a similar factor. The `energy` region is where that assumption
      lives, so it is the line to read if the total misses.""")
    end
end

println()
if isempty(FAILURES)
    println("no failures")
else
    println(length(FAILURES), " failure", length(FAILURES) == 1 ? "" : "s", ":")
    foreach(f -> println("  * ", f), FAILURES)
end
println()
