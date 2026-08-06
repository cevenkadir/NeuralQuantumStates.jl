"""
    AbstractLinearSolver

How to solve the linear system at the heart of a preconditioned update.

Every stochastic-reconfiguration variant reduces to solving `(A + λI) x = b` for a symmetric
positive-semidefinite `A`. Keeping that behind an interface is not ceremony: the geometric
tensor is routinely **singular**, because redundant parameter directions and directions no
sample explores both produce exact zero modes. Which regularization is used is therefore a real
choice with real consequences, not an implementation detail — and it is also the seam where a
batched GPU solver such as BatchSolve.jl would attach.

# Interface
    solve(solver, A, b, shift) -> x
"""
abstract type AbstractLinearSolver end

"""
    solve(solver, A, b, shift) -> x

Solve `(A + shift*I) x = b`.
"""
function solve end

"""
    CholeskySolver() <: AbstractLinearSolver

Dense Cholesky factorization of the shifted system.

The fastest option when the matrix fits in memory and the shift is large enough to make it
positive definite. Falls back to a symmetric indefinite factorization if Cholesky fails, which
happens when the shift is too small to cover the zero modes — silently returning garbage there
would be worse than the extra cost.
"""
struct CholeskySolver <: AbstractLinearSolver end

function solve(::CholeskySolver, A::AbstractMatrix, b::AbstractVector, shift::Real)
    M = Hermitian(A + shift * I)
    F = cholesky(M; check=false)
    issuccess(F) && return F \ b
    # Not positive definite: the shift did not cover the zero modes. Bunch-Kaufman handles the
    # indefinite case rather than returning a meaningless answer.
    return bunchkaufman(M; check=false) \ b
end

"""
    PseudoInverseSolver(; rtol=1e-10) <: AbstractLinearSolver

Moore–Penrose pseudo-inverse via SVD, discarding singular values below `rtol` times the largest.

The most robust option, and the one to reach for when the geometric tensor is genuinely
singular. Rather than inflating the small singular values as a diagonal shift does, this
**projects out** the corresponding directions entirely: the update simply makes no progress
along directions the samples carry no information about, instead of making a wildly large
and arbitrary move along them.

Costs an SVD, so it is for small parameter counts or for diagnosis.
"""
struct PseudoInverseSolver <: AbstractLinearSolver
    rtol::Float64
end
PseudoInverseSolver(; rtol::Real=1e-10) = PseudoInverseSolver(Float64(rtol))

function solve(s::PseudoInverseSolver, A::AbstractMatrix, b::AbstractVector, shift::Real)
    return pinv(Matrix(A + shift * I); rtol=s.rtol) * b
end

"""
    ConjugateGradientSolver(; maxiter=1000, tol=1e-10) <: AbstractLinearSolver

Iterative solve by conjugate gradients, through KrylovKit.

Never forms a factorization, so it scales to parameter counts where a dense solve is
impossible. It needs only the action of the matrix on a vector, which is what makes a
matrix-free geometric tensor usable.
"""
struct ConjugateGradientSolver <: AbstractLinearSolver
    maxiter::Int
    tol::Float64
end
ConjugateGradientSolver(; maxiter::Integer=1000, tol::Real=1e-10) =
    ConjugateGradientSolver(Int(maxiter), Float64(tol))

function solve(s::ConjugateGradientSolver, A, b::AbstractVector, shift::Real)
    op = x -> A * x + shift * x
    x, _ = linsolve(op, b; maxiter=s.maxiter, tol=s.tol, isposdef=true)
    return x
end
