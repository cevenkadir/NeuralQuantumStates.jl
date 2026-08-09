# GPU benchmarks

CPU versus CUDA for the parts of a variational run that could plausibly move to a device.

This environment is deliberately separate from every package in the repository, so that nothing
in the stack acquires a GPU dependency. `CUDA` is a weak dependency of `NQSOptimisers` and a
direct dependency of nothing.

## One AD engine, and why

Zygote only. The engine is still an axis in the source (`ENGINES`), so adding a second one is a
single line, but both candidates are out for reasons worth recording.

**Enzyme was tried here and removed.** On CUDA arrays its failed reverse-mode compilations leave
the GPUArrays allocator double-releasing buffers. BenchmarkTools runs `gcscrub` before every
trial, so the wreckage surfaces as `ArgumentError("Attempt to release freed data")` inside
whichever *later, unrelated* measurement happens to trigger the collection — including Zygote
ones. That corrupts the baseline this file exists to hold, which is too high a price for a column
whose answer is already known.

[`../reactant`](../reactant) answered it on the CPU, where nothing else was at risk, using a
control rung to isolate the cause: Enzyme has no reverse rule for a complex `zgemm`, which the
differentiated region always reaches because `input_type` answers `eltype(θ)`. Given a real batch
Enzyme works and is 2.2× *slower* than Zygote. Nothing was lost.

**Reactant is absent for a different reason.** XLA preallocates the bulk of the card when it
initializes and CUDA.jl keeps its own pool, so the two cannot share a process. That column lives
in [`../reactant`](../reactant), which prints the same rung labels so the two reports read side
by side — but on whatever hardware each one ran, which its header states.

## Setup

Two commands, from the repository root. There is nothing to configure — the environment declares
its own path dependencies, so it does not matter whether any other environment in the repository
has been instantiated.

```
julia --project=lib/NQSCore/benchmark/gpu -e 'using Pkg; Pkg.instantiate()'
julia --project=lib/NQSCore/benchmark/gpu lib/NQSCore/benchmark/gpu/benchmarks.jl
```

On a machine with no usable CUDA device the second command prints why and exits 0, so both are
safe to run anywhere.

### What the machine needs

Only an NVIDIA driver. CUDA.jl ships the toolkit itself as an artifact and downloads it on first
use — expect a one-off download of a gigabyte or so the first time `using CUDA` runs, which is
part of `instantiate` here. You do not need a system CUDA installation, and you should not set
`JULIA_CUDA_USE_BINARYBUILDER` or point CUDA.jl at one unless you have a specific reason.

To check what CUDA.jl found:

```
julia --project=lib/NQSCore/benchmark/gpu -e 'using CUDA; CUDA.versioninfo()'
```

If that reports a device but the suite still skips, `CUDA.functional()` is the thing returning
false — usually a driver older than the toolkit CUDA.jl selected.

With more than one GPU, `CUDA.device!(i)` before the run picks one; the suite uses whichever
device is current.

## Reading the output

The speedups are the least interesting numbers here.

Connected configurations are computed on the **host**, because the compiled operator is a nested
structure that cannot be uploaded as it stands. Every optimization step therefore ships the
samples out and the connected configurations and matrix elements back, and the return leg is
`max_conn` times larger than the outbound one. The suite reports that transfer cost explicitly,
including as a fraction of one forward pass, alongside the padding ratio
`max_conn / mean(n_conn)`.

Those two numbers answer a question that was deliberately left open rather than guessed at:
whether porting the connected-configuration kernel to the device is worth the work, or whether
the split architecture is good enough. NetKet shipped the same split for years before writing
device-side operators, and says of its device-side ones that the benefit is composability rather
than speed. A padding ratio near 1 means such a kernel would waste nothing; well above 1 means
most of the transfer, and most of what the kernel would evaluate, is padding.

## Expect failures

**No part of this stack has run on a GPU.** The device paths are written and their dispatch is
tested, but the first real execution is this suite. Every measurement is individually guarded:
one unavailable path prints its error and the run continues, rather than hiding the rest.

A failure here is a result, not a problem — it says which of the untested paths needs work.
