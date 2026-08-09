"""
    AbstractLinearSolver

How to solve the linear system at the heart of a preconditioned update.

Every stochastic-reconfiguration variant reduces to solving `(A + λI) x = b` for a symmetric
positive-semidefinite `A`. The geometric tensor is routinely **singular** — redundant parameter
directions and directions no sample explores both give exact zero modes — so how it is handled
is a real choice, and this is also where a GPU solver attaches.

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
positive definite. Falls back to a symmetric indefinite factorization when Cholesky fails, which
means the shift did not cover the zero modes.
"""
struct CholeskySolver <: AbstractLinearSolver end

function solve(::CholeskySolver, A::AbstractMatrix, b::AbstractVector, shift::Real)
    M = Hermitian(A + shift * I)
    F = cholesky(M; check=false)
    issuccess(F) && return F \ b
    return indefinite_fallback(M, A, b, shift)
end

"""
    indefinite_fallback(M, A, b, shift) -> x

What [`CholeskySolver`](@ref) does when the shift did not cover the tensor's zero modes.

Split out because it is the one step with no accelerator equivalent — cuSOLVER exposes Cholesky
and LU but not Bunch–Kaufman — so a GPU extension can replace exactly this.
"""
indefinite_fallback(M, A, b::AbstractVector, shift::Real) = bunchkaufman(M; check=false) \ b

"""
    PseudoInverseSolver(; rtol=1e-10) <: AbstractLinearSolver

Moore–Penrose pseudo-inverse via SVD, discarding singular values below `rtol` times the largest.

The most robust option. Rather than inflating small singular values as a diagonal shift does,
this projects the corresponding directions out: the update makes no progress along directions
the samples carry no information about, instead of an arbitrarily large move along them.

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

Never forms a factorization, so it scales to parameter counts where a dense solve is impossible.
Needs only the action of the matrix on a vector, which is what makes a matrix-free geometric
tensor usable.
"""
struct ConjugateGradientSolver <: AbstractLinearSolver
    maxiter::Int
    tol::Float64
end
ConjugateGradientSolver(; maxiter::Integer=1000, tol::Real=1e-10) =
    ConjugateGradientSolver(Int(maxiter), Float64(tol))

function solve(s::ConjugateGradientSolver, A, b::AbstractVector, shift::Real)
    op = x -> A * x + shift * x
    # These flags are what select conjugate gradients. KrylovKit cannot see inside a function
    # operator, so without them it falls back to GMRES — same answer, orders of magnitude slower
    # on an ill-conditioned tensor. `isposdef` alone is not enough.
    x, _ = linsolve(
        op, b; maxiter=s.maxiter, tol=s.tol,
        ishermitian=true, issymmetric=eltype(b) <: Real, isposdef=true
    )
    return x
end
