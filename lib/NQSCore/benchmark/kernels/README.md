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
XLA_REACTANT_GPU_PREALLOCATE=false NQS_REACTANT_BACKEND=cpu \
  julia --project=lib/NQSCore/benchmark/kernels lib/NQSCore/benchmark/kernels/kernel_probe.jl
```

| variable | what it does |
|---|---|
| `NQS_REACTANT_BACKEND` | `cpu` or `gpu`; unset, Reactant picks |
| `XLA_REACTANT_GPU_PREALLOCATE` | `false` stops XLA taking the card up front |
| `XLA_REACTANT_GPU_MEM_FRACTION` | XLA's share of a card, default `0.75` |
| `NQS_REACTANT_CUDA_DIR` | a CUDA toolkit for XLA other than the one Reactant bundles |

Those are the names **Reactant** reads. The `XLA_PYTHON_CLIENT_*` variables are JAX's and are
ignored here.

**Reactant on CPU with CUDA.jl on the device is a valid configuration for this probe**, and the
right one when Reactant's XLA refuses the card. It answers whether the kernel runs and whether it
raises — which is what decides feasibility — and leaves the CUDA.jl timing untouched. Only the
head-to-head timing needs both on the same hardware.

**When Reactant's XLA refuses the card.** It compiles kernels with the toolkit inside its own
artifact, which on a stock install is CUDA 13 — and CUDA 13 dropped compute capability below 7.5,
so a V100 or GV100 (7.0) fails with `ptxas too old` and `BlasLt is unavailable` before any of this
code runs. `NQS_REACTANT_CUDA_DIR` pointing at the CUDA 12 toolkit CUDA.jl is already using on the
same machine is the one lever; the probe prints both paths so you can see whether they differ.

**Cap XLA's memory when both are on the device.** CUDA.jl and XLA are both live in this process —
Reactant needs CUDA.jl loaded to lower a KernelAbstractions kernel at all, which is the whole
subject here — and XLA otherwise takes 75% of every visible card at startup.

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

**Does it run, and is the answer right.** That is the whole question, and the bar is lower than it
looks.

Measured on an RTX 5000 Ada, `expect_and_grad` under CUDA.jl is 5.245 ms — `local_energy` 2.609 ms
plus `energy_gradient` 2.471 ms, 96.9% between them — and the connected-configuration kernel is
**38.8 µs** of that. So the kernel is 1.5% of `local_energy`, and porting it is not about kernel
speed at all. It is the enabler: the bulk of `local_energy` is the network forward over 53,248
connected configurations, ~2.57 ms, and that can only move into XLA if the data it consumes is
already there. Reactant does the gradient in 355 µs against Zygote-on-CUDA's 2.471 ms; if the
forward goes at a similar factor the step lands near 0.75 ms.

A kernel several times slower than CUDA.jl's is therefore still worth having. Only a catastrophic
ratio — or a wrong answer — argues against.

**Raising is the secondary question, and it already has an answer: no.** The loops are
data-dependent, and Reactant says `cannot raise op to stablehlo` on the `scf.for`. That costs
fusion with the network evaluation downstream, not the ability to run, so the probe compiles
unraised *first* and treats the raised attempt as a bonus. Asking only for the raised form takes
the whole compilation down with it, which is how an earlier version of this probe managed to leave
the real question untested.
