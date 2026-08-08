"""
    QuantumGeometricTensor(X, diag_scale=0.0)

The quantum geometric tensor `S = XᵀX`, represented by `X` and never formed.

`X` is the stacked, centered, probability-weighted design matrix from `_weighted_design`, of
size `2N × P`. The tensor it defines is `P × P`, and *that* is the matrix a variational run
cannot afford: for a network with ten thousand parameters `S` is a hundred million entries
before anything has been solved, and building it costs `O(N P²)` on top.

Nothing here needs `S` itself. Its action on a vector is two matrix–vector products,

```math
S v = X^T (X v)
```

each `O(N P)`, and an iterative solver asks for nothing else. That is the whole type: an
[`AbstractLinearSolver`](@ref) that only multiplies — [`ConjugateGradientSolver`](@ref) — can
solve against it directly, with peak memory set by `X` rather than by `S`.

# What this is and is not

This is NetKet's `QGTJacobianDense` with a matrix-free apply: the Jacobian is still materialized
(it comes from [`local_estimators`](@ref)), only the geometric tensor is not. NetKet's
`QGTOnTheFly` goes one step further and stores neither, re-differentiating on every product.
That would remove the remaining `2N × P`, at the cost of an automatic-differentiation pass per
solver iteration, and is a separate change.

# Fields
- `X`: the design matrix.
- `diagonal`: `diag(XᵀX)`, precomputed once, for the relative regularization.
- `diag_scale`: the relative shift `ε₁`, applied as `S + ε₁ diag(S)`. The absolute shift is the
  solver's business, since it is what makes the system solvable at all.
"""
struct QuantumGeometricTensor{M<:AbstractMatrix,V<:AbstractVector}
    X::M
    diagonal::V
    diag_scale::Float64
end

function QuantumGeometricTensor(X::AbstractMatrix, diag_scale::Real=0.0)
    diagonal = vec(sum(abs2, X; dims=1))
    return QuantumGeometricTensor(X, diagonal, Float64(diag_scale))
end

Base.size(S::QuantumGeometricTensor) = (size(S.X, 2), size(S.X, 2))
Base.size(S::QuantumGeometricTensor, i::Integer) = size(S)[i]
Base.eltype(S::QuantumGeometricTensor) = eltype(S.X)

function Base.:*(S::QuantumGeometricTensor, v::AbstractVector)
    Sv = transpose(S.X) * (S.X * v)
    iszero(S.diag_scale) && return Sv
    return Sv .+ (S.diag_scale .* S.diagonal) .* v
end

"""
    to_dense(S) -> Matrix

Materialize the geometric tensor. Defeats the purpose of the type, and exists so that a direct
solver, or a test, can have the matrix when it genuinely needs one.
"""
function to_dense(S::QuantumGeometricTensor)
    A = transpose(S.X) * S.X
    iszero(S.diag_scale) && return A
    return A + S.diag_scale * Diagonal(S.diagonal)
end
