"""
    QuantumGeometricTensor(X, diag_scale=0.0)

The quantum geometric tensor `S = XᵀX`, represented by `X` and never formed.

`X` is the stacked, centered, probability-weighted design matrix from `_weighted_design`, of
size `2N × P`. The tensor it defines is `P × P`, which for a network with ten thousand
parameters is a hundred million entries costing `O(N P²)` to build.

Nothing here needs `S` itself. Its action on a vector is `S v = Xᵀ(X v)`, two products of
`O(N P)` each, and an iterative solver asks for nothing else — so
[`ConjugateGradientSolver`](@ref) can solve against this directly with peak memory set by `X`.

This is NetKet's `QGTJacobianDense` with a matrix-free apply: the Jacobian is still materialized
by [`local_estimators`](@ref), only the geometric tensor is not.

# Fields
- `X`: the design matrix.
- `diagonal`: `diag(XᵀX)`, precomputed once, for the relative regularization.
- `diag_scale`: the relative shift `ε₁`, applied as `S + ε₁ diag(S)`. The absolute shift is the
  solver's business.
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

Materialize the geometric tensor. Defeats the purpose of the type, and exists for a direct
solver or a test that genuinely needs the matrix.
"""
function to_dense(S::QuantumGeometricTensor)
    A = transpose(S.X) * S.X
    iszero(S.diag_scale) && return A
    return A + S.diag_scale * Diagonal(S.diagonal)
end
