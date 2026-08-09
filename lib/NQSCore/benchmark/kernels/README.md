# Kernel probe

Can Reactant run the connected-configuration kernel, and how fast against CUDA.jl?

This is a decision probe, not a benchmark suite. It answers one question and is meant to be
deleted once it has.

## Why

On a GPU, Reactant is a net loss today. It cannot share arrays with CUDA.jl, so under Reactant the
connected configurations fall back to the host — **3.834 ms** against the **74.8 µs** the device
kernel already achieves — while the gradient it accelerates is worth only ~1.5 ms of a 3.336 ms
step. Reactant pays off on a device only if the kernel comes along.

Two documented facts make that cheaper than a StableHLO rewrite:

- **KernelAbstractions kernels run inside a `@compile` region.** They do not have to be rewritten.
- **Raising is for differentiation and fusion, not for running.** Nothing differentiates connected
  configurations — they are data — so a kernel that merely runs is already enough.

The suspected obstacle is the element type: `configs` holds `BaseInt`, and XLA tensors need
primitive element types. `BaseInt{T,Ti,B}` is a single-field isbits wrapper around `value::T`, and
the kernel touches it only through the two unchecked accessors and one equality.

## Setup

```
julia --project=lib/NQSCore/benchmark/kernels -e 'using Pkg; Pkg.instantiate()'
XLA_PYTHON_CLIENT_MEM_FRACTION=0.25 \
  julia --project=lib/NQSCore/benchmark/kernels lib/NQSCore/benchmark/kernels/kernel_probe.jl
```

`NQS_REACTANT_BACKEND=cpu` or `=gpu` pins the target; unset, Reactant picks.

**Cap XLA's memory.** CUDA.jl and XLA are both live in this process — Reactant needs CUDA.jl
loaded to lower a KernelAbstractions kernel at all, which is the whole subject here — and XLA
otherwise takes the better part of every visible card at startup, leaving CUDA.jl nothing. If the
CUDA.jl half reports out of memory, that fraction is the knob.

This is also why the probe has its own environment rather than living in `../reactant`: that one
holds the validated 15–18× AD measurement, and adding CUDA to it would re-resolve it.

## What it does

**Probe 1** hands the kernel's own `Vector{BaseInt}` to `Reactant.to_rarray` and reports what
happens. The expected outcome is a failure naming the element type — what matters is whether it
names *only* that.

**Probe 2** defines a copy of the kernel over raw integers, in the probe file, changing no package.
The body is line for line the extension's, with `BaseInt{V,Ti,B}` replaced by `V` and the base
carried as a `Val`; any other difference would blur the answer. It then reports:

- whether the rewrite is correct at all, checked on the host against `connected_padded!` slot for
  slot, the way `lib/ConnectedBasisConfigurations/test/kernelabstractions.jl` does;
- whether Reactant raises it to plain StableHLO or keeps a kernel call in the module;
- its time against **CUDA.jl running the same kernel on the same data in the same process** — a
  comparison rather than a number remembered from another machine.

## Reading the result

The Reactant time against the CUDA.jl one, on the last line.

Within a small factor, carrying the kernel into Reactant is worth doing, and the whole variational
step can live in one compiled region. Far above, CUDA.jl keeps the GPU and Reactant stays a CPU
story — where it is already a measured 15–18× on the gradient, and needs no kernel at all.

Whether it *raised* is the secondary question. Unraised but running is a success; raised also buys
fusion with the network evaluation that consumes the kernel's output, which is what would let the
compiled region span `local_energy` and the gradient together.
