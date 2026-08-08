"""
What stochastic reconfiguration has to do differently when its matrices live on a CUDA device.

Deliberately small. Everything that can be written as a matrix-vector product already works on
a device array without help — which is why `ConjugateGradientSolver` and
[`QuantumGeometricTensor`](@ref) appear nowhere in this file — and the Cholesky path is served
by cuSOLVER through the ordinary `LinearAlgebra` methods. What is left is the two places where
the CPU implementation reaches for something CUDA does not provide.
"""
module NQSOptimisersCUDAExt

using CUDA
using LinearAlgebra
using NQSOptimisers
using NQSOptimisers: ConjugateGradientSolver, PseudoInverseSolver, indefinite_fallback, solve

"""
cuSOLVER has Cholesky and LU but no Bunch–Kaufman, so an indefinite tensor cannot be factorized
on the device. Falling back to conjugate gradients keeps the solve where the data already is;
copying the matrix to the host to factorize it there would cost more than the solve.
"""
NQSOptimisers.indefinite_fallback(
    ::Hermitian{<:Any,<:CUDA.CuMatrix}, A, b::CUDA.CuVector, shift::Real
) = solve(ConjugateGradientSolver(), A, b, shift)

"""
`pinv` has no CUDA method. The generic `LinearAlgebra` fallback would run its singular-value
post-processing by scalar indexing, which on a device array either errors outright or silently
serializes into a host round-trip per element — so this refuses rather than appearing to work.
"""
function NQSOptimisers.solve(
    ::PseudoInverseSolver, ::CUDA.CuMatrix, ::CUDA.CuVector, ::Real
)
    throw(ArgumentError(
        "PseudoInverseSolver is CPU-only: LinearAlgebra.pinv has no CUDA method and its " *
        "fallback indexes scalars. On a device use CholeskySolver, or ConjugateGradientSolver " *
        "— which with mode=:matrixfree never forms the tensor at all."
    ))
end

end # module NQSOptimisersCUDAExt
