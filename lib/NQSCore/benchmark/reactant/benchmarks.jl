"""
Zygote against Enzyme, and against Reactant's XLA compilation of Enzyme, on this stack's
gradient.

Run with

```
julia --project=lib/NQSCore/benchmark/reactant lib/NQSCore/benchmark/reactant/benchmarks.jl
```

It needs no accelerator: Reactant runs on the CPU and picks up an XLA GPU when one is there.
`NQS_REACTANT_BACKEND=cpu` or `=gpu` forces the choice; unset, Reactant picks the best it finds.

# The question

Lux recommends Reactant + Enzyme over Zygote for both CPU and GPU. This stack's differentiated
region is not the one that recommendation was written for: the parameters are `ComplexF64`, the
closure that reaches the network is real-to-real with the complex reparameterization inside it,
and a `ComponentArray` is rebuilt within the derivative. Any of those can be the thing that
decides it, so the comparison is made rung by rung rather than end to end — `energy_gradient` is
four wrappers around the model, and a backend can win on one and lose on another.

The rungs are the same four the CUDA suite bisects, so the two reports read side by side:

    loss_split    the real/imag reparameterization, `restore`, the cotangent, and the model
    loss_flat     without the reparameterization
    loss_theta    without `restore` as well — parameters straight from the ansatz
    model_loss    the model alone

# Why Reactant is not measured end to end

`energy_gradient` calls `DifferentiationInterface.gradient` directly, and DifferentiationInterface
has no Reactant integration — Reactant is a compiler, not an ADTypes backend. Diverting an
XLA-compiled region out of `NQSCore` needs a seam that does not exist yet, and inventing one
before there is a number to justify it is the wrong order. So Reactant is measured on the same
closures, and `expect_and_grad` is reported for the two backends that reach it today.

# A failure here is the result

Reverse mode over `ComplexF64` is the risk. Every measurement is guarded and the run ends with a
tally, because a wall of `FAILED` lines is easy to skim past and is exactly what this suite might
legitimately produce.
"""

using BenchmarkTools
using ConnectedBasisConfigurations   # `compile`, which turns an OpSum into the operator form
using DifferentiationInterface
using Enzyme
using EnzymeCore
using KernelAbstractions   # loads NQSCore's extension, which is what makes the seam below needed
using LinearAlgebra
using Lux
using NQSAnsatze
using NQSCore
using OperatorAlgebra
using Printf
using Random
using Reactant
using SymBasis
using Zygote

# The backend has to be chosen before Reactant initializes its client, so this comes before
# anything that could touch one.
let requested = get(ENV, "NQS_REACTANT_BACKEND", "")
    isempty(requested) || Reactant.set_default_backend(requested)
end

println()
println("Reactant ", pkgversion(Reactant), "  |  Enzyme ", pkgversion(Enzyme),
        "  |  Zygote ", pkgversion(Zygote), "  |  Lux ", pkgversion(Lux))
# Which hardware this actually ran on, said here and not only in the README. The rung labels
# below are deliberately identical to the CUDA suite's, which is what makes the two reports
# readable side by side — and also exactly what invites someone to set a number here against
# `expect, CUDA` and conclude something false. Unset, Reactant picks the best device it finds,
# which on a GPU node is a GPU: the choice is worth printing rather than assuming.
let requested = get(ENV, "NQS_REACTANT_BACKEND", "")
    n = try
        length(Reactant.devices())
    catch err
        println("XLA devices unavailable — ", first(split(sprint(showerror, err), '\n')))
        0
    end
    println("XLA target: ", isempty(requested) ? "chosen by Reactant" : requested,
            ", $n device", n == 1 ? "" : "s",
            " — see the `service.cc` lines above for the platform and the cards.")
    println("Rung labels below match ../gpu deliberately. The hardware need not: do not diff a\n",
            "number here against one there without checking both headers first.")
    # XLA's BFC allocator takes the better part of each card the moment the client initializes.
    # On a shared node that is antisocial, and in the same process as CUDA.jl it is fatal — which
    # is why the CUDA suite is a separate environment. `NQS_REACTANT_BACKEND=cpu` avoids it.
    n > 0 && isempty(requested) &&
        println("Reactant has preallocated most of each device. Do not run ../gpu concurrently.")
end

BenchmarkTools.DEFAULT_PARAMETERS.seconds = 2

# ------------------------------------------------------------------------------- reporting

"""
Each failure's label, its first line, and its whole text.

The whole text is kept because Enzyme's failures are not one-liners: the actionable part — the
instruction it could not handle, or the type it could not infer — sits in the middle, and line
one is often just an exception name. Everything is dumped after the tally, so a run that fails
says why on the first attempt rather than costing another trip to the machine it ran on.
"""
const FAILURES = Tuple{String,String,String}[]

report(label, t) = @printf("  %-46s %12s\n", label, BenchmarkTools.prettytime(t * 1e9))

"""Time `f`, recording rather than raising when a path is not there."""
function timed(label, f)
    try
        trial = @benchmark $f()
        t = minimum(trial).time / 1e9
        report(label, t)
        # A handful of evaluations is a compile time wearing a benchmark's clothes. Saying so is
        # cheaper than someone reading a confident number off a single sample.
        length(trial.times) < 5 &&
            @printf("  %-46s %12d\n", "    samples (few — treat as indicative)",
                    length(trial.times))
        return t
    catch err
        text = sprint(showerror, err)
        push!(FAILURES, (strip(label), first(split(text, '\n')), text))
        @printf("  %-46s %12s\n", label, "FAILED")
        println("      ", first(split(text, '\n')))
        return nothing
    end
end

"""Record a failure from a guarded block that is not a timing, and print its first line."""
function failed(label, err)
    text = sprint(showerror, err)
    push!(FAILURES, (label, first(split(text, '\n')), text))
    println("      ", first(split(text, '\n')))
end

"""
The head and tail of a long error, with the middle elided.

Both engines bury the sentence that matters. Enzyme prints the whole LLVM function it could not
differentiate and then says why underneath — "Complex inputs not yet supported in reverse mode
for BLAS calls" arrived two hundred lines below its own headline. Reactant prints the entire
StableHLO module and then the XLA diagnostic. Neither middle has ever been the answer, and a
report nobody scrolls to the end of is not a report.
"""
function abridged(text, head::Int=4, tail::Int=25, width::Int=400)
    # Long lines are clipped before long files are, because the two failure modes are different
    # and only one of them is fixed by taking fewer lines. Reactant reports a rejected XLA option
    # by listing every option it *would* have accepted — one line, twenty thousand characters —
    # where the name that was rejected sits in the first eighty.
    clip(l) = length(l) <= width ? l : l[1:prevind(l, width)] * " … [$(length(l)) chars]"
    lines = map(clip, split(text, '\n'))
    length(lines) <= head + tail + 1 && return join(lines, '\n')
    elided = length(lines) - head - tail
    return join(vcat(lines[1:head], ["", "  … $elided lines elided …", ""],
                     lines[(end-tail+1):end]), '\n')
end

speedup(a, b) = (a === nothing || b === nothing) ? nothing : a / b
function compare(label, a, b)
    s = speedup(a, b)
    s === nothing ? println("  $label: not comparable") :
    @printf("  %-46s %11.2fx\n", label, s)
end

# ------------------------------------------------------------------- the two seam methods

# Both belong in package extensions if this evaluation says Reactant is worth adopting —
# `NQSCoreReactantExt` and `NQSAnsatzeEnzymeCoreExt`. They live here while it is still a
# question, so that no package in the stack gains a weak dependency on the strength of a
# benchmark that has not run yet.

# `ConcreteRArray` is not an `Array`, so the KernelAbstractions extension's
# `device_backend(x::AbstractArray) = KernelAbstractions.get_backend(x)` claims it and asks a
# GPU backend about an XLA buffer. `nothing` means "host path", which is always correct and
# merely slower: connected configurations are computed on the CPU, as they were before the
# device kernel existed.
NQSCore.device_backend(::Reactant.AnyConcreteRArray) = nothing

# `colocate` is protected by `ChainRulesCore.@non_differentiable`, which Enzyme does not consult.
# Without this the derivative stops at a `copyto!` of the batch instead of reaching the
# parameters — configurations are data, and nothing wants a gradient with respect to them.
EnzymeCore.EnzymeRules.inactive(::typeof(NQSAnsatze.colocate), args...) = nothing

# ---------------------------------------------------------------------------- the engines

"""
The gradient with respect to `f`'s **last** argument, everything before it held constant.

`Enzyme.gradient` returns a tuple aligned with the arguments and `nothing` wherever one was
`Const`, so putting the active argument last makes `[end]` the answer however many constants
precede it. That is what lets one function serve every rung, and — because Reactant compiles it
like any other call — serve the compiled column too.
"""
enzyme_gradient(f::F, args...) where {F} = Enzyme.gradient(
    Enzyme.Reverse, Const(f), map(Const, Base.front(args))..., last(args)
)[end]

"""The same shape for Zygote, which takes a closure and returns a one-tuple."""
zygote_gradient(f::F, args...) where {F} =
    Zygote.gradient(w -> f(Base.front(args)..., w), last(args))[1]

"""
Compile `f(args...)` and return a thunk, reporting what compilation cost.

Compilation is reported on its own line rather than folded into a measurement: a caller pays it
once, and a benchmark that includes it is measuring the compiler. The function form of `@compile`
is used because the rungs are walked in a loop; `sync=true` is what makes the returned thunk
block until XLA has finished, which is the barrier this suite needs and which no `CUDA.@sync`
would provide.
"""
function compile_thunk(label, f::F, args) where {F}
    t = @elapsed thunk = Reactant.compile(f, args; sync=true)
    report(label, t)
    return thunk
end

# ------------------------------------------------------------------------ the rungs, named

# Top-level functions rather than closures, and the active argument last. A closure over `x` and
# `c` would hand Reactant the batch as a traced constant, which recompiles on every call and
# reports a compile time as though it were a gradient — the one failure here that looks like a
# plausible number rather than an error.

"""`L(θ) = 2 Σ Re[conj(c) log ψ]` with the parameters in the shape the ansatz uses."""
loss_theta(a, x, c, θ) = NQSCore._gradient_loss(a, θ, x, c)

"""...with the parameters flattened, so `restore` is inside the derivative."""
loss_flat(a, restore, x, c, q) = NQSCore._gradient_loss(a, restore(q), x, c)

"""...and with the real/imag reparameterization on top: the closure `energy_gradient` builds."""
function loss_split(a, restore, x, c, n, v)
    # `complex.(re, im)` and not `re .+ im .* im_part`, which is what NQSCore actually writes.
    # `im` is `Complex{Bool}`, and Reactant rejects that outright: "in RNumber, in T, expected
    # T<:Union{…, ComplexF64, ComplexF32}, got Type{Complex{Bool}}". The two are identical in
    # value and in what Zygote and Enzyme do with them, so this rung still measures the library's
    # closure — but the substitution is itself a finding, and it is the one-line change
    # `_energy_gradient` and `_log_derivatives` would need before Reactant could reach them.
    q = @views complex.(v[1:n], v[(n+1):(2n)])
    return NQSCore._gradient_loss(a, restore(q), x, c)
end

"""The model alone, with no cotangent contraction to carry."""
model_loss(a, x, θ) = sum(real, log_amplitude(a, θ, x))

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

# The same system the CUDA suite measures, so that the rung timings there and here are the same
# quantity on different hardware rather than two different benchmarks that resemble each other.
const NSITES = 12
const ALPHA = 4

let spec = Spin(1 // 2), nsites = NSITES
    b = basis(dof_object(spec), nsites)
    H = compile(tfi(nsites; h_x=0.9, h_z=0.1))
    a = LuxAnsatz(RBM(nsites, ALPHA), spec, nsites; rng=Xoshiro(0))
    θ = init_parameters(a, Xoshiro(0))

    println("\nmodel: RBM($nsites, $ALPHA) over $(length(b.states)) configurations, ",
            "$(sum(length, values(θ))) parameters")

    # `to_rarray` rather than `reactant_device()`. The CUDA suite avoids `gpu_device()` for the
    # same reason: Lux's device machinery demotes ComplexF64 to ComplexF32, which would make
    # every comparison below a precision comparison as much as an engine one. Asserted rather
    # than assumed, because a silent demotion is exactly what it would look like if it happened.
    θ_ra = nothing
    try
        θ_ra = Reactant.to_rarray(θ)
        @printf("  %-46s %12s\n", "parameter element type, host", string(eltype(θ.weight)))
        @printf("  %-46s %12s\n", "parameter element type, Reactant",
                string(eltype(θ_ra.weight)))
        eltype(θ_ra.weight) === eltype(θ.weight) ||
            println("      the conversion changed precision; every Reactant number below is a\n" *
                    "      precision comparison as well as an engine one")
    catch err
        println("      moving the parameters to Reactant failed:")
        failed("to_rarray(θ)", err)
    end

    # ------------------------------------------------------------ what the derivative sees

    vs = FullSumState(a, θ; backend=AutoZygote(), basis=b)
    xs = NQSCore.configurations_of(vs, b.states)
    logψ = log_amplitude(a, θ, xs)
    p = NQSCore.born_probabilities(logψ)
    E = NQSCore.local_energy(vs, H, b.states)
    c = NQSCore._gradient_cotangent(E, p)
    flat, restore = NQSCore.flatten_parameters(θ)
    n = length(flat)
    v = vcat(real.(flat), imag.(flat))
    # The batch `configurations_of` builds is complex, because `input_type` answers with
    # `eltype(θ)` so that a device derivative reaches cuBLAS `zgemm` rather than a generic
    # mixed complex-real kernel. That choice is load-bearing for CUDA and is what the last rung
    # below exists to price against, since a complex batch is also what puts `zgemm` inside the
    # derivative — and Enzyme has no reverse rule for a complex one.
    xs_real = Float64.(real.(xs))

    @printf("\n  %-46s %12s\n", "sample batch element type", string(eltype(xs)))
    @printf("  %-46s %12d\n", "parameters (complex)", n)

    println("\nforward pass")
    t_fwd = timed("log_amplitude, host", () -> log_amplitude(a, θ, xs))

    # The cost of *not* compiling, measured once rather than benchmarked. Eager Reactant traces
    # and compiles each operation as it meets it, so this is not a slow path but a different
    # activity altogether, and `@belapsed` on it would take a coffee break. It is here because
    # it is the number that says why the compiled region has to exist.
    if θ_ra !== nothing
        try
            x_ra = Reactant.to_rarray(xs)
            t = @elapsed Reactant.synchronize(log_amplitude(a, θ_ra, x_ra))
            report("log_amplitude, Reactant uncompiled (once)", t)
            t_fwd === nothing || compare("  the price of not compiling", t, t_fwd)
        catch err
            println("      eager Reactant forward failed:")
            failed("log_amplitude, uncompiled", err)
        end
    end

    # ------------------------------------------------------- the rungs, engine by engine

    # `flatten_parameters` runs *outside* every differentiated region, on whatever the parameters
    # are. On Reactant arrays its ComponentArrays concatenation is a sequence of traced
    # operations with no compiled region around them, so it gets its own line: it would otherwise
    # hide inside every rung below and be attributed to the derivative.
    println("\nflattening, which no engine differentiates")
    timed("flatten_parameters, host", () -> NQSCore.flatten_parameters(θ))
    if θ_ra !== nothing
        # Once, not benchmarked, for the same reason as the eager forward above: outside a
        # compiled region every one of these operations is a separate trace-and-compile, and
        # `@belapsed` would sit through several hundred of them.
        try
            t = @elapsed NQSCore.flatten_parameters(θ_ra)
            report("flatten_parameters, Reactant (once)", t)
        catch err
            println("      flattening Reactant parameters failed:")
            failed("flatten_parameters on Reactant arrays", err)
        end
    end

    rungs = (
        ("split + restore + cotangent + model", loss_split, (a, restore, xs, c, n), v),
        ("  - split", loss_flat, (a, restore, xs, c), flat),
        ("  - split - restore", loss_theta, (a, xs, c), θ),
        ("  - split - restore - cotangent", model_loss, (a, xs), θ),
        # Not a rung of `energy_gradient` — the control that says whether an engine's failure is
        # about this model or about one operation in it. The batch is real here, so the layer's
        # product is complex x real and falls to a generic kernel instead of `zgemm`. An engine
        # that dies on the rung above and lives on this one is telling us that the complex batch,
        # not the network, is what it cannot differentiate.
        ("  model only, real batch (control)", model_loss, (a, xs_real), θ),
    )

    println("\nthe gradient, rung by rung")
    times = Dict{String,Float64}()
    gradients = Dict{String,Any}()
    for (label, f, constants, active) in rungs
        println("\n  ", label)
        for (engine, grad) in (("Zygote", zygote_gradient), ("Enzyme", enzyme_gradient))
            t = timed("  $engine", () -> grad(f, constants..., active))
            t === nothing && continue
            times["$label|$engine"] = t
            # One more call, for the agreement check below. It cannot throw: the timed run just
            # made several hundred of them.
            gradients["$label|$engine"] = grad(f, constants..., active)
        end

        if θ_ra !== nothing
            try
                # The active argument and every array constant become Reactant arrays. `n` and
                # `restore` are neither, and stay as they are — traced into the compiled region
                # as constants, which is what they are.
                ra = map(z -> z isa AbstractArray || z isa NamedTuple ?
                              Reactant.to_rarray(z) : z, constants)
                args = (f, ra..., Reactant.to_rarray(active))
                thunk = compile_thunk("  Reactant, compiling (once)", enzyme_gradient, args)
                t = timed("  Reactant + Enzyme", () -> thunk(args...))
                if t !== nothing
                    times["$label|Reactant"] = t
                    gradients["$label|Reactant"] = thunk(args...)
                end
            catch err
                @printf("  %-46s %12s\n", "  Reactant + Enzyme", "FAILED")
                failed("$(strip(label)) under Reactant", err)
            end
        end

        for engine in ("Enzyme", "Reactant")
            compare("    $engine against Zygote", get(times, "$label|Zygote", nothing),
                    get(times, "$label|$engine", nothing))
        end
    end

    # ------------------------------------------------------------------------- agreement

    # A faster wrong gradient is not a result, and this suite is read from a file by someone who
    # did not watch it run, so it has to say so itself. The comparison is on the rung the library
    # actually differentiates, whose gradient is a plain real vector under every engine.
    println("\nagreement on the rung the library uses")
    reference = get(gradients, "split + restore + cotangent + model|Zygote", nothing)
    if reference === nothing
        println("  Zygote produced no gradient to compare against")
    else
        for engine in ("Enzyme", "Reactant")
            g = get(gradients, "split + restore + cotangent + model|$engine", nothing)
            if g === nothing
                println("  $engine: no gradient")
                continue
            end
            # Reactant's tolerance is looser on purpose: XLA reassociates, and Lux's own
            # documentation reports differences around 1e-8 against the eager path.
            tol = engine == "Reactant" ? 1e-5 : 1e-10
            h = Array(g)
            rel = norm(h .- reference) / norm(reference)
            @printf("  %-30s relative %.3e   (tolerance %.0e)  %s\n",
                    engine, rel, tol, isapprox(h, reference; rtol=tol) ? "agrees" : "DISAGREES")
            isapprox(h, reference; rtol=tol) ||
                push!(FAILURES, ("$engine vs Zygote, split rung",
                                 "gradients disagree beyond the tolerance", ""))
        end
    end

    # --------------------------------------------------- the two halves of the layer

    # The most interesting probe in the suite. On the GV100 the `logtwocosh` reduction's Zygote
    # reverse is 4.5x its forward, and an exact `@scalar_rule` for it measured at 1.0x — the cost
    # is Zygote's per-element pullback machinery for a complex broadcast, not the transcendental.
    # Enzyme and XLA use entirely different mechanisms, so this number either improves a lot or
    # fails outright, and either answer decides more than the rungs above do.
    println("\nthe two halves of the layer, on ComplexF64")
    let nh = ALPHA * NSITES, xs_real = Float64.(real.(xs))
        u = randn(Xoshiro(0), ComplexF64, nh, length(b.states))
        V = randn(Xoshiro(1), ComplexF64, nh, NSITES)

        reduction(z) = sum(real, sum(NQSAnsatze.logtwocosh, z; dims=1))
        matmul(x, M) = sum(real, M * x)

        for (label, f, constants, active) in (("logtwocosh reduction", reduction, (), u),
                                              ("complex matmul", matmul, (xs_real,), V))
            println("\n  ", label)
            fw = timed("  forward", () -> f(constants..., active))
            for (engine, grad) in (("Zygote", zygote_gradient), ("Enzyme", enzyme_gradient))
                rv = timed("  $engine reverse", () -> grad(f, constants..., active))
                fw === nothing || rv === nothing ||
                    @printf("  %-46s %11.2fx\n", "    reverse / forward", rv / fw)
            end
            θ_ra === nothing && continue
            try
                args = (f, map(Reactant.to_rarray, constants)..., Reactant.to_rarray(active))
                thunk = compile_thunk("  Reactant, compiling (once)", enzyme_gradient, args)
                rv = timed("  Reactant + Enzyme reverse", () -> thunk(args...))
                fw === nothing || rv === nothing ||
                    @printf("  %-46s %11.2fx\n", "    reverse / forward", rv / fw)
            catch err
                @printf("  %-46s %12s\n", "  Reactant + Enzyme reverse", "FAILED")
                failed("$label under Reactant", err)
            end
        end
    end

    # ------------------------------------------------------ end to end, where it reaches

    # Only the two engines DifferentiationInterface can drive. `AutoEnzyme()` needs no change to
    # any package — `backend` is a plain ADTypes object on the state — so this measures the whole
    # step including `local_energy`, which the rungs above deliberately exclude.
    println("\nend to end, through DifferentiationInterface")
    for (engine, backend) in (("Zygote", AutoZygote()), ("Enzyme", AutoEnzyme()))
        state = FullSumState(a, θ; backend=backend, basis=b)
        timed("expect_and_grad, $engine", () -> expect_and_grad(state, H))
    end
    timed("expect (no derivative)", () -> expect(vs, H))
    timed("local_energy alone", () -> NQSCore.local_energy(vs, H, b.states))
end

# ------------------------------------------------------------------------------ the tally

println()
if isempty(FAILURES)
    println("no failures")
else
    println(length(FAILURES), " failure", length(FAILURES) == 1 ? "" : "s", ":")
    for (label, headline, _) in FAILURES
        println("  * ", label, ": ", headline)
    end
    # And then the whole of each one. Enzyme puts the actionable part — the instruction it could
    # not handle, or the type it could not infer — several lines in, so a summary that stops at
    # the first line reports that something failed without reporting what.
    for (label, _, text) in FAILURES
        isempty(text) && continue
        println("\n", "-"^78, "\n", label, "\n")
        println(abridged(text))
    end
end
println()
