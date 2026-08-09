# Reactant benchmarks

Zygote against Enzyme, and against Reactant's XLA compilation of Enzyme, on this stack's
gradient.

Lux recommends Reactant + Enzyme over Zygote for both CPU and GPU. This suite exists to find out
whether that holds here, because the differentiated region in this stack is not the one that
recommendation was written for: the parameters are `ComplexF64`, the closure that reaches the
network is real-to-real with the complex reparameterization inside it, and a `ComponentArray` is
rebuilt within the derivative.

Nothing in the repository depends on Reactant or Enzyme. This environment is separate, as
`../gpu` is, and the two are deliberately **not** merged — see below.

## Setup

Two commands, from the repository root.

```
julia --project=lib/NQSCore/benchmark/reactant -e 'using Pkg; Pkg.instantiate()'
julia --project=lib/NQSCore/benchmark/reactant lib/NQSCore/benchmark/reactant/benchmarks.jl
```

No accelerator is required, but do not assume it stays off one. **Unset, Reactant picks the best
device it finds, which on a GPU node is a GPU** — and XLA's BFC allocator then takes the better
part of every card it can see (on a three-GV100 node, 23.8 GiB each) the moment the client
initializes. That is antisocial on a shared machine and fatal in the same process as CUDA.jl.

`NQS_REACTANT_BACKEND=cpu` forces the host; `=gpu` asks for a device explicitly. The run prints
which it got, and warns when it has preallocated. **Never run this and `../gpu` at the same
time.**

Expect the first run to be slow. Reactant compiles six regions and reports each compile time on
its own line; Enzyme compiles too, on first call. The measurements exclude all of it.

## Why this is not in `../gpu`

XLA preallocates the bulk of a GPU's memory when Reactant initializes, and CUDA.jl keeps its own
pool. In one process the two fight over the card and neither set of numbers survives it. So
`../gpu` compares Zygote and Enzyme, both native CUDA.jl, and this directory adds the Reactant
column in its own process.

The two scripts print the same rung labels, so the reports read side by side. That is how the
Zygote-to-Reactant comparison is actually made.

## Reading the output

`energy_gradient` is four wrappers around the model, and an engine can win on one and lose on
another, so the comparison is made rung by rung rather than end to end:

| rung | what it adds |
|---|---|
| `split + restore + cotangent + model` | the closure `energy_gradient` actually differentiates |
| `- split` | without the real/imag reparameterization |
| `- split - restore` | without the `ComponentArray` rebuild as well |
| `- split - restore - cotangent` | the model alone |

A jump between two consecutive rungs names the wrapper responsible instead of suggesting one.
These are the same rungs `../gpu` bisects.

Two other sections carry more weight than their size suggests:

- **the two halves of the layer.** On the GV100 the `logtwocosh` reduction's Zygote reverse is
  4.5× its forward, and an exact `@scalar_rule` for it measured at 1.0× — the cost is Zygote's
  per-element pullback machinery for a complex broadcast, not the transcendental. Enzyme and XLA
  use entirely different mechanisms, so this number either improves a lot or fails outright.
- **agreement.** A faster wrong gradient is not a result. Each engine's gradient on the rung the
  library uses is compared against Zygote's — `1e-10` for Enzyme, `1e-5` for Reactant, since XLA
  reassociates and Lux's own documentation reports differences around `1e-8`.

## Reactant is not measured end to end

`energy_gradient` calls `DifferentiationInterface.gradient` directly, and
DifferentiationInterface has no Reactant integration — Reactant is a compiler, not an ADTypes
backend. Diverting an XLA-compiled region out of `NQSCore` needs a seam that does not exist, and
inventing one before there is a number to justify it is the wrong order. Reactant is therefore
measured on the same closures, and `expect_and_grad` is reported for the two engines that reach
it today.

`AutoEnzyme()` needs no change to any package: `backend` is a plain ADTypes object on the
variational state, and `NQSCore` reaches AD only through DifferentiationInterface.

## The two provisional methods

Near the top of `benchmarks.jl` are two method definitions that belong in package extensions if
this evaluation says Reactant is worth adopting — `NQSCoreReactantExt` and
`NQSAnsatzeEnzymeCoreExt`. They live in the script while it is still a question, so that no
package in the stack gains a weak dependency on the strength of a benchmark that has not run.

- `NQSCore.device_backend(::Reactant.AnyConcreteRArray) = nothing`. A `ConcreteRArray` is not an
  `Array`, so the KernelAbstractions extension's `device_backend(x::AbstractArray)` claims it and
  asks a GPU backend about an XLA buffer. `nothing` means "host path", which is always correct
  and merely slower.
- `EnzymeCore.EnzymeRules.inactive(::typeof(NQSAnsatze.colocate), args...)`. `colocate` is
  protected by `ChainRulesCore.@non_differentiable`, which Enzyme does not consult.

`KernelAbstractions` is loaded on purpose: without it the generic `device_backend(::Any) =
nothing` would handle Reactant arrays by itself and the first of those two methods would look
unnecessary. Loading it makes the run representative of the setting that matters.

## Expect failures

Reverse mode over `ComplexF64` is the risk this suite exists to price. Every measurement is
guarded and the run ends with an explicit tally, because a wall of `FAILED` lines is easy to skim
past and is exactly what this comparison might legitimately produce.

A failure here is a result. It says which engine cannot do what this stack needs, which is worth
knowing before any of it is built into the packages.
