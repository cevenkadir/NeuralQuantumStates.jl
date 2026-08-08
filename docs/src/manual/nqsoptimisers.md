```@meta
CurrentModule = NQSOptimisers
```

# NQSOptimisers.jl

*Stochastic reconfiguration and natural-gradient preconditioning.*

Plain gradient descent follows the steepest direction in *parameter* space. For a wavefunction
that is the wrong geometry: two parameters can be rescaled against each other without changing
the state at all, while completely changing the gradient.
[`StochasticReconfiguration`](@ref) follows the steepest direction in *state* space instead,
using the quantum geometric tensor

```math
S_{kk'} = \mathrm{Re} \left[ \langle O_k^* O_{k'} \rangle
          - \langle O_k^* \rangle \langle O_{k'} \rangle \right]
```

as the metric, and solving ``S \delta = \nabla E``.

## Why this is not optional

In the ordered regime of a transverse-field Ising chain, plain descent **stalls** at the
classical energy. The reason is structural rather than a matter of tuning: the gradient carries
a factor of the Born probability ``p(s)``, which vanishes for exactly the configurations whose
amplitude needs to grow. The geometric tensor carries the same factor, and dividing by it undoes
the suppression.

On a 4-site chain at ``J = 1``, ``h_x = 0.9``:

| | iterations | gap to the exact ground state |
|---|---|---|
| plain gradient descent | 5000 | 0.98 |
| stochastic reconfiguration | 200 | 0.0 (machine precision) |

## Two forms of the same update

`:sr` builds the ``P \times P`` matrix ``X^T X``; `:minsr` builds the ``N \times N`` matrix
``X X^T``. The identity

```math
(X^T X + \lambda I)^{-1} X^T = X^T (X X^T + \lambda I)^{-1}
```

means they give **exactly** the same update — not approximately — so the choice is purely one
of cost. For a network with far more parameters than samples, which is the usual case, the
second is dramatically cheaper. That is the kernel trick, known here as MinSR or SRt. `:auto`
picks whichever matrix is smaller.

### A third form, which builds neither

Both of those allocate a square matrix, and for a large enough network that matrix — not the
sampling — is what runs out of memory first. `mode=:matrixfree` allocates neither: it wraps the
design matrix in a [`QuantumGeometricTensor`](@ref), whose only operation is multiplication, and
hands that to [`ConjugateGradientSolver`](@ref), which needs nothing else.

The cost moves from one factorization to one matrix-vector product per solver iteration. Which
is cheaper depends on the shape and on how quickly the solve converges: for a `1024 × 4096`
design matrix, forming ``X^T X`` costs about as much as fifty matrix-free products, so the
iterative route wins outright if the solver converges in fewer than that — and it never
allocates the 134 MB the tensor itself would need. For a small model, forming the matrix is
cheaper, and `:sr` remains the right default.

A direct solver has nothing to factorize in this mode and is rejected with an error rather than
silently materializing the matrix behind your back.

### On a GPU

Loading `CUDA` alongside this package brings in an extension that changes exactly two things,
because everything else already works on a device array unchanged:

- `CholeskySolver`'s indefinite fallback uses conjugate gradients instead of Bunch–Kaufman,
  which cuSOLVER does not provide. The Cholesky path itself needs no help.
- `PseudoInverseSolver` refuses to run. `LinearAlgebra.pinv` has no CUDA method, and the generic
  fallback would process singular values by scalar indexing — a host round-trip per element, or
  an outright error. Better to say so than to appear to work.

`ConjugateGradientSolver` needs nothing: it only multiplies, and with `:matrixfree` there is no
square matrix to keep on the device in the first place.

!!! warning "Not yet run on hardware"
    The device paths are written and their dispatch is tested, but no part of this stack has
    executed on a GPU. Treat the first run as an experiment; `lib/NQSCore/benchmark/gpu` exists
    to make it a measured one.

## Regularization, and what it is not for

The geometric tensor is routinely singular: redundant parameter directions and directions no
sample explores both give exact zero modes, and condition numbers of `1e15` are ordinary. Some
regularization is therefore mandatory — `diag_shift` adds an absolute amount to the diagonal,
`diag_scale` an amount relative to each diagonal entry.

!!! warning "Use a smaller step than for plain descent"
    The preconditioned update is routinely one to two orders of magnitude larger than the raw
    gradient. A learning rate that is fine for plain gradient descent will overshoot here and
    can leave the run plateaued at a non-eigenstate — with a perfectly ordinary gradient, which
    makes it look like ill-conditioning when it is not. Start near `0.05` with
    `diag_shift = 0.01`, and reduce the step if the energy plateaus above the ground state
    while its variance stays large.

## Solvers

Which one matters, because the system being solved is genuinely singular.

| Solver | Behaviour on zero modes |
|---|---|
| [`CholeskySolver`](@ref) | Fast; falls back to an indefinite factorization if the shift is too small |
| [`PseudoInverseSolver`](@ref) | **Projects out** small singular directions rather than inflating them |
| [`ConjugateGradientSolver`](@ref) | Iterative; needs only matrix-vector products, so it scales |

The [`AbstractLinearSolver`](@ref) interface is also the seam where a GPU solver attaches.
[`ConjugateGradientSolver`](@ref) needs only matrix-vector products, so pairing it with
`mode=:matrixfree` — which wraps the design matrix in a [`QuantumGeometricTensor`](@ref) instead
of forming `S` — means neither the `P × P` nor the `2N × 2N` matrix is ever allocated. For a
network with many parameters that is the difference between a run that fits in memory and one
that does not.

## Not yet here

Time evolution (TDVP) uses the same geometric tensor with a real or imaginary time step. It is
left out rather than shipped half-done.

## Quick example

```@example nqsoptimisers
using NQSOptimisers, NQSCore, ConnectedBasisConfigurations, SymBasis, OperatorAlgebra
using DifferentiationInterface, ForwardDiff, Random, LinearAlgebra

spec, nsites = Spin(1 // 2), 4
b = basis(dof_object(spec), nsites)
ops = local_operators(spec)
H = OpSum(vcat(
    [Op(2 .* ops.sz, i) * Op(2 .* ops.sz, mod1(i + 1, nsites)) for i in 1:nsites],
    [0.9 * Op(2 .* ops.sx, i) for i in 1:nsites],
))

a = LogStateVector(spec, nsites, b)
vs = FullSumState(a, init_parameters(a, Xoshiro(11); scale=0.3); backend=AutoForwardDiff())

optimize!(vs, H, StochasticReconfiguration(; diag_shift=1e-2);
    iterations=600, learning_rate=0.05)
expect(vs, H)
```
