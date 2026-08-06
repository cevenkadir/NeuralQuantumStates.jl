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

The [`AbstractLinearSolver`](@ref) interface is also the seam where a batched GPU solver such as
BatchSolve.jl would attach. That is deliberately not wired up: whether a batched path beats a
plain Cholesky or CG on a realistic geometric tensor is an open question, and batching is only a
win if it pays for itself here.

## Not yet here

Time evolution (TDVP) uses the same geometric tensor with a real or imaginary time step. It is
left out rather than shipped half-done.

## Quick example

```@example index
using NQSOptimisers, NQSCore, ConnectedConfigs, SymBasis, OperatorAlgebra
using DifferentiationInterface, ForwardDiff, Random, LinearAlgebra

spec, nsites = Spin(1 // 2), 4
b = basis(dof_object(spec), nsites)
ops = local_operators(spec)
H = OpSum(vcat(
    [Op(2 .* ops.sz, i) * Op(2 .* ops.sz, mod1(i + 1, nsites)) for i in 1:nsites],
    [0.9 * Op(2 .* ops.sx, i) for i in 1:nsites],
))

a = LogStateVector(spec, nsites, b)
vs = FullSumState(a, init_parameters(a, Xoshiro(11); scale=0.3), AutoForwardDiff())

optimize!(vs, H, StochasticReconfiguration(; diag_shift=1e-2);
    iterations=600, learning_rate=0.05)
expect(vs, H)
```
