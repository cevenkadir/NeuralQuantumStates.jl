```@meta
CurrentModule = NQSCore
```

# Log-derivatives and gradients

[`log_derivatives`](@ref) computes

```math
O_{sk} = \frac{\partial \log \psi(x_s)}{\partial \theta_k}
```

the `(n_samples, n_parameters)` matrix behind both the energy gradient — a covariance between
`O` and the local energies — and the quantum geometric tensor that stochastic reconfiguration
inverts. Rows are samples and columns are parameters, matching NetKet, so `S = O'O` has the
shape one expects without transposing.

## Two independent kinds of complexity

Conflating them is the usual source of gradients that almost work.

**A complex log-amplitude with real parameters** is the common case for a network emitting a
modulus and a phase. `O` is complex even though `θ` is real, and both parts come from a single
differentiation pass over `[real(log ψ); imag(log ψ)]` — one pass rather than two, since the two
share all of their intermediate work.

**Complex parameters** are ambiguous until a convention is fixed:

- `holomorphic=true` assumes `ψ` is holomorphic in `θ`, so `∂/∂θ = ∂/∂θ_re`, and `O` keeps one
  column per parameter.
- `holomorphic=false`, the default, treats real and imaginary parts as `2n` independent real
  parameters and returns `[∂/∂θ_re  ∂/∂θ_im]`. This is always valid; the holomorphic form is a
  shortcut that is silently wrong when the ansatz does not satisfy it, so it must be asked for
  explicitly.

Whatever the parameterization, [`match_parameter_shape`](@ref) puts a raw gradient back into the
shape of the parameters it belongs to, so `θ .- η .* ∇` is always meaningful and an optimizer
never has to ask how the ansatz was parameterized.

## The gradient does not build `O`

The variational energy gradient is one contraction of the Jacobian:

```math
\partial_k \langle E \rangle = 2 \, \mathrm{Re}
    \left[ \langle O_k^* \, \Delta E_{\mathrm{loc}} \rangle \right],
\qquad \Delta E_s = E_s - \bar{E}
```

Building the whole matrix in order to contract it is wasteful, and on a reverse-mode backend it
is dramatically so: a Jacobian costs one pass *per sample*. Writing the same quantity as the
gradient of a scalar

```math
L(\theta) = 2 \sum_s \mathrm{Re}[\overline{c_s} \log \psi(x_s)],
\qquad c_s = p_s (E_s - \bar{E})
```

makes it a single differentiation of a scalar-valued function, which every backend does
optimally — one reverse pass regardless of the sample count. `c` does not depend on `θ`, so
`∇L` *is* the energy gradient. This is what [`expect_and_grad`](@ref) does. For a small RBM under
Zygote it is about 17× faster than forming the Jacobian, and the gap widens with the sample
count.

Subtracting the mean is not cosmetic: without it the estimator has a non-vanishing variance even
at an exact eigenstate, where every local energy is the same number.

`log_derivatives` still exists and still builds the matrix, because stochastic reconfiguration
needs `O` itself and not just this one contraction of it. [`local_estimators`](@ref) is the entry
point that returns it.

## Parameter containers

Lux hands parameters back as nested `NamedTuple`s, while every piece of linear algebra
downstream wants a flat vector. [`flatten_parameters`](@ref) provides the bijection and the
function that rebuilds the original structure:

```@example logderiv
using NQSCore

θ = (W=[1.0 2.0; 3.0 4.0], b=[5.0, 6.0])
flat, restore = flatten_parameters(θ)
flat
```

```@example logderiv
restore(flat)
```

It rebuilds a `NamedTuple`, not a `ComponentArray`: the promise downstream is that a gradient
has the same structure as the parameters it belongs to, so that `fmap` can walk the two
together.

## Backends

Differentiation goes through
[DifferentiationInterface.jl](https://github.com/JuliaDiff/DifferentiationInterface.jl), so
`backend` is any of its objects — `AutoForwardDiff()`, `AutoZygote()`, `AutoEnzyme()`. None of
those packages is a dependency here; load the one you want and pass it.
