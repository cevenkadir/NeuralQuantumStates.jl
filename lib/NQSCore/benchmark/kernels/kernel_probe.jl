"""
Can Reactant run the connected-configuration kernel, and how fast?

Run with

```
XLA_REACTANT_GPU_PREALLOCATE=false NQS_REACTANT_BACKEND=cpu \\
  julia --project=lib/NQSCore/benchmark/kernels lib/NQSCore/benchmark/kernels/kernel_probe.jl
```

`NQS_REACTANT_BACKEND` picks Reactant's target and `XLA_REACTANT_GPU_MEM_FRACTION` caps its share
of a card (`XLA_REACTANT_GPU_PREALLOCATE=false` stops it taking one up front). Those are the names
Reactant reads — *not* the `XLA_PYTHON_CLIENT_*` ones, which belong to JAX and are ignored here.

Reactant on the CPU with CUDA.jl on the device is a perfectly good configuration for this probe:
it answers whether the kernel runs and whether it raises, which is what decides feasibility, and
it leaves the CUDA.jl timing untouched. Only the Reactant-versus-CUDA.jl comparison needs both on
the same hardware.

`NQS_REACTANT_CUDA_DIR=/path/to/cuda-12` points Reactant's XLA at a toolkit other than the one it
bundles. That is the lever for a card its own toolkit refuses — CUDA 13 dropped compute capability
below 7.5, so a V100 or a GV100 gets `ptxas too old` and `BlasLt is unavailable` from a stock
install.

# Why this exists

On a GPU, Reactant is a net loss today. It cannot share arrays with CUDA.jl, so the connected
configurations fall back to the host — 3.834 ms against the 74.8 µs the device kernel already
achieves — while the gradient it speeds up is only worth 1.5 ms of a 3.336 ms step. Reactant
therefore only pays off on a device if the kernel comes along.

Two documented facts make that cheaper than rewriting the algorithm in StableHLO:
KernelAbstractions kernels *run* inside a `@compile` region, and raising — turning a kernel into
tensor operations — is needed for differentiation and fusion, not for running. Nothing
differentiates connected configurations; they are data. So a kernel that merely runs is enough.

The suspected obstacle is the element type. `configs` holds `BaseInt`, and XLA tensors need
primitive element types. `BaseInt{T,Ti,B}` is a single-field isbits wrapper around `value::T`, and
the kernel touches it only through the two unchecked accessors and one equality — so carrying the
raw integer instead should be the whole of the change.

This probe tests exactly that, in two steps, and changes no package. Probe 1 runs the kernel as it
stands and reports what Reactant says. Probe 2 defines a copy of it over raw integers, here in
this file, and asks three questions: does it run, does it raise, and how does it compare with
CUDA.jl on the same data in the same process.

# What decides the port

Probe 2's time against the CUDA.jl number printed beside it. Close to it, carrying the kernel
into Reactant is worth doing; far above, CUDA.jl keeps the GPU and Reactant stays a CPU story.
"""

using BenchmarkTools
using CUDA
using ConnectedBasisConfigurations
using KernelAbstractions
using OperatorAlgebra
using Printf
using Random
using Reactant
using SymBasis
using SymBasis.DigitBase: BaseInt

let requested = get(ENV, "NQS_REACTANT_BACKEND", "")
    isempty(requested) || Reactant.set_default_backend(requested)
end

# The toolkit XLA compiles kernels with. Reactant defaults to the one inside its own artifact,
# which on a stock install is CUDA 13 — and CUDA 13 dropped compute capability below 7.5, so a
# GV100 (7.0) gets `ptxas too old` and `BlasLt is unavailable` before any of this code runs.
# Pointing it at the CUDA 12 toolkit that CUDA.jl is already using on the same machine is the
# cheapest thing to try, and this is the only knob for it.
let dir = get(ENV, "NQS_REACTANT_CUDA_DIR", "")
    isempty(dir) || (Reactant.XLA.CUDA_DATA_DIR[] = dir)
end

BenchmarkTools.DEFAULT_PARAMETERS.seconds = 2

println()
println("Reactant ", pkgversion(Reactant), "  |  CUDA ", pkgversion(CUDA),
        "  |  KernelAbstractions ", pkgversion(KernelAbstractions))
println("CUDA.jl functional: ", CUDA.functional(),
        CUDA.functional() ? "  ($(CUDA.name(CUDA.device())), CUDA $(CUDA.runtime_version()))" : "")
println("XLA devices: ", try
    length(Reactant.devices())
catch err
    "unavailable — " * first(split(sprint(showerror, err), '\n'))
end)

# Printed because it is the thing most likely to need changing, and because a stock Reactant and
# a working CUDA.jl on the same machine can be using two different CUDA versions without saying so.
println("Reactant's CUDA toolkit: ", Reactant.XLA.CUDA_DATA_DIR[])
println("CUDA.jl's CUDA toolkit:  ", try
    @eval(Main, import CUDA_Runtime_jll)
    Main.CUDA_Runtime_jll.artifact_dir
catch
    "not resolvable here — `CUDA.versioninfo()` names the version"
end)
println("""
If Reactant's XLA refuses this card, NQS_REACTANT_CUDA_DIR=<the second path> is the one lever.
Both allocators are live in this process; XLA_REACTANT_GPU_MEM_FRACTION and
XLA_REACTANT_GPU_PREALLOCATE cap XLA's share if CUDA.jl reports out of memory below.""")

report(label, t) = @printf("  %-46s %12s\n", label, BenchmarkTools.prettytime(t * 1e9))
note(label, what) = @printf("  %-46s %12s\n", label, what)

"""Print the first line of a failure and keep going, since a failure here is a result."""
function failed(what, err)
    text = sprint(showerror, err)
    lines = split(text, '\n')
    println("  ", what, " FAILED")
    for l in lines[1:min(6, length(lines))]
        println("      ", l)
    end
    length(lines) > 6 && println("      … ", length(lines) - 6, " more lines")
    return nothing
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

# ------------------------------------------------------ the kernel again, over raw integers

# A copy of the extension's kernel with `BaseInt{V,Ti,B}` replaced by `V` and the base carried as
# a `Val`. The body is otherwise line for line the same, deliberately: the question is whether the
# element type is the only thing standing between this kernel and Reactant, and any other
# difference would blur the answer. If it is, this is the shape the extension would take.

@inline function raw_read(value::V, position::Integer, ::Val{B}) where {V,B}
    if ispow2(B)
        bits = trailing_zeros(B)
        return Int((value >> ((position - 1) * bits)) & V(B - 1))
    else
        return Int(rem(div(value, V(B)^(position - 1)), V(B)))
    end
end

@inline function raw_write(value::V, position::Integer, digit::Integer, ::Val{B}) where {V,B}
    if ispow2(B)
        bits = trailing_zeros(B)
        shift = (position - 1) * bits
        mask = V(B - 1) << shift
        return (value & ~mask) | (V(digit) << shift)
    else
        p = V(B)^(position - 1)
        old = rem(div(value, p), V(B))
        return value + (V(digit) - old) * p
    end
end

@kernel function raw_connected_kernel!(
    configs, mels, counts,
    @Const(term_start), @Const(factor_position), @Const(factor_col_start),
    @Const(colptr), @Const(outs), @Const(vals), @Const(states),
    scratch_states, scratch_vals,
    n_diagonal::Int, n_terms::Int, height::Int, base::Val,
)
    b = @index(Global)
    T = eltype(mels)

    @inbounds begin
        state = states[b]

        diagonal = zero(T)
        for t in 1:n_diagonal
            v = one(T)
            alive = true
            for fi in term_start[t]:(term_start[t+1]-1)
                d = raw_read(state, factor_position[fi], base)
                cs = factor_col_start[fi]
                lo = colptr[cs+d]
                if lo >= colptr[cs+d+1]
                    alive = false
                    break
                end
                v *= vals[lo]
            end
            alive && (diagonal += v)
        end

        k = 1
        for t in (n_diagonal+1):n_terms
            n = 1
            src, dst = 1, 2
            scratch_states[1, src, b] = state
            scratch_vals[1, src, b] = one(T)

            alive = true
            for fi in term_start[t]:(term_start[t+1]-1)
                pos = factor_position[fi]
                cs = factor_col_start[fi]
                m = 0
                for j in 1:n
                    s = scratch_states[j, src, b]
                    amp = scratch_vals[j, src, b]
                    d = raw_read(s, pos, base)
                    for p in colptr[cs+d]:(colptr[cs+d+1]-1)
                        m += 1
                        scratch_states[m, dst, b] = raw_write(s, pos, outs[p], base)
                        scratch_vals[m, dst, b] = amp * vals[p]
                    end
                end
                if m == 0
                    alive = false
                    break
                end
                n = m
                src, dst = dst, src
            end

            if alive
                for j in 1:n
                    v = scratch_vals[j, src, b]
                    iszero(v) && continue
                    s′ = scratch_states[j, src, b]
                    if s′ == state
                        diagonal += v
                    else
                        k += 1
                        configs[k, b] = s′
                        mels[k, b] = v
                    end
                end
            end
        end

        if iszero(diagonal)
            for j in 2:k
                configs[j-1, b] = configs[j, b]
                mels[j-1, b] = mels[j, b]
            end
            k -= 1
        else
            configs[1, b] = state
            mels[1, b] = diagonal
        end
        counts[b] = k

        for j in (k+1):height
            configs[j, b] = state
            mels[j, b] = zero(T)
        end
    end
end

"""
The function Reactant compiles: allocate, launch, return.

Written the way the Reactant kernel tutorial writes one — buffers allocated inside from `similar`,
the backend asked of the data with `get_backend` — rather than taking preallocated buffers, since
that is the form known to trace.
"""
function raw_connected(
    term_start, factor_position, factor_col_start, colptr, outs, vals, states,
    n_diagonal::Int, n_terms::Int, height::Int, width::Int, base::Val,
)
    n = length(states)
    configs = similar(states, height, n)
    mels = similar(vals, height, n)
    counts = similar(states, Int, n)
    scratch_states = similar(states, width, 2, n)
    scratch_vals = similar(vals, width, 2, n)

    backend = KernelAbstractions.get_backend(states)
    kernel = raw_connected_kernel!(backend, 64)
    kernel(
        configs, mels, counts,
        term_start, factor_position, factor_col_start, colptr, outs, vals, states,
        scratch_states, scratch_vals,
        n_diagonal, n_terms, height, base;
        ndrange=n,
    )
    return configs, mels, counts
end

# ================================================================================ the system

const NSITES = 12

let spec = Spin(1 // 2), nsites = NSITES
    b = basis(dof_object(spec), nsites)
    states = b.states
    op = flatten(tfi(nsites; h_x=0.9, h_z=0.1))
    height, n = max_conn_size(op), length(states)
    width = max(op.max_branch, 1)

    S = eltype(states)
    V, B = S.parameters[1], S.parameters[3]
    base = Val(B)

    println("\nsystem: TFI on $nsites sites, $n states, max_conn $height")
    note("packed state type", string(S))
    note("raw integer type", string(V))
    note("base", string(B))

    # The reference every variant is checked against.
    host_c = Matrix{S}(undef, height, n)
    host_m = Matrix{eltype(op)}(undef, height, n)
    host_k = Vector{Int}(undef, n)
    connected_padded!(host_c, host_m, host_k, op, states)
    raw_states = collect(reinterpret(V, states))
    raw_host_c = collect(reinterpret(V, host_c))

    # ------------------------------------------------------- probe 1, the kernel as it stands

    println("\nprobe 1 — the kernel as it stands, on Reactant arrays")
    try
        rstates = Reactant.to_rarray(states)
        note("to_rarray(states::Vector{BaseInt})", "ok")
        println("      unexpected: BaseInt was accepted, so the element type is not the obstacle")
        println("      type is ", typeof(rstates))
    catch err
        failed("to_rarray(states::Vector{BaseInt})", err)
        println("""
      This is the expected outcome and the reason probe 2 exists. What matters is whether the
      message names the element type and nothing else.""")
    end

    # -------------------------------------------------- probe 2, the same kernel over integers

    println("\nprobe 2 — the same kernel over raw integers")

    # First on the host, so a wrong answer is caught before any compiler is involved.
    cpu_c, cpu_m, cpu_k = raw_connected(
        op.term_start, op.factor_position, op.factor_col_start,
        op.colptr, op.outs, op.vals, raw_states,
        op.n_diagonal, op.n_terms, height, width, base,
    )
    agrees = cpu_k == host_k && cpu_c == raw_host_c && cpu_m == host_m
    note("matches the reference kernel, on host", string(agrees))
    agrees || println("""
      The raw-integer rewrite is wrong, so every number below is meaningless. Fix this first;
      the two accessors are where to look.""")

    # CUDA.jl, the bar to beat. Same kernel, same data, same process — so this is a comparison
    # rather than a number remembered from another run on another machine.
    t_cuda = nothing
    if CUDA.functional()
        try
            cu = (CuArray(op.term_start), CuArray(op.factor_position),
                  CuArray(op.factor_col_start), CuArray(op.colptr), CuArray(op.outs),
                  CuArray(op.vals), CuArray(raw_states))
            gpu_c, gpu_m, gpu_k = raw_connected(
                cu..., op.n_diagonal, op.n_terms, height, width, base
            )
            CUDA.@sync nothing
            ok = Array(gpu_k) == host_k && Array(gpu_c) == raw_host_c && Array(gpu_m) == host_m
            note("matches the reference kernel, on CUDA", string(ok))
            t_cuda = @belapsed CUDA.@sync raw_connected(
                $cu..., $(op.n_diagonal), $(op.n_terms), $height, $width, $base
            )
            report("connections on CUDA.jl", t_cuda)
        catch err
            failed("the CUDA.jl run", err)
        end
    else
        println("  no CUDA device, so there is nothing to compare against")
    end

    # And Reactant.
    try
        # Bound to names rather than splatted from a tuple: `@code_hlo` below parses the call
        # expression it is given, and a `...` inside it is not something it has to understand.
        ts = Reactant.to_rarray(op.term_start)
        fp = Reactant.to_rarray(op.factor_position)
        fc = Reactant.to_rarray(op.factor_col_start)
        cp = Reactant.to_rarray(op.colptr)
        ou = Reactant.to_rarray(op.outs)
        va = Reactant.to_rarray(op.vals)
        st = Reactant.to_rarray(raw_states)
        nd, nt = op.n_diagonal, op.n_terms
        args = (ts, fp, fc, cp, ou, va, st, nd, nt, height, width, base)

        # Did it raise? A raised kernel lowers to plain StableHLO; an unraised one keeps a kernel
        # call in the module. Both are results — raising buys fusion with the network that
        # consumes this, running is already enough — but which one it is decides how wide the
        # compiled region in the real thing can be.
        hlo = try
            sprint(show, Reactant.@code_hlo raise = true raw_connected(
                ts, fp, fc, cp, ou, va, st, nd, nt, height, width, base
            ))
        catch err
            failed("@code_hlo", err)
            ""
        end
        if !isempty(hlo)
            launched = occursin("kernel_call", hlo) || occursin("custom_call", hlo)
            note("raised to plain StableHLO", string(!launched))
            note("HLO module size (lines)", string(count(==('\n'), hlo) + 1))
        end

        t0 = @elapsed thunk = Reactant.compile(raw_connected, args; sync=true, raise=true)
        report("Reactant, compiling (once)", t0)

        r_c, r_m, r_k = thunk(args...)
        ok = Array(r_k) == host_k && Array(r_c) == raw_host_c && Array(r_m) == host_m
        note("matches the reference kernel, on Reactant", string(ok))

        t_reactant = @belapsed $thunk($args...)
        report("connections on Reactant", t_reactant)
        t_cuda === nothing ||
            @printf("  %-46s %11.2fx\n", "Reactant / CUDA.jl", t_reactant / t_cuda)
    catch err
        failed("the Reactant run", err)
    end
end

println("""

What this decides: the Reactant time against the CUDA.jl one. Within a small factor, carrying the
kernel into Reactant is worth doing and the whole step can live in one compiled region. Far above,
CUDA.jl keeps the GPU and Reactant stays a CPU story.
""")
