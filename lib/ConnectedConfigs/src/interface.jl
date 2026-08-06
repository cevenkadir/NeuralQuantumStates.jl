# The extension seam. Everything the compiler and the kernel need from the outside world goes
# through the generic functions below, and nothing else. SymBasis and OperatorAlgebra supply
# the default methods (see `backends/`), but they are ordinary methods on ordinary generic
# functions: another library plugs in by adding its own, with no package extension to load and
# no capability flag to set.

"""
    expand_terms(operator) -> Vector{<:AbstractVector{<:Pair}}

Flatten `operator` into a plain sum of products of single-site factors.

The return value is one entry per term; each term is a vector of `position => matrix` pairs
listing that term's factors in **matrix-product order** (the leftmost factor of `A * B * C`
first). `position` is an `Integer` digit position — site `i` is digit `i`, counting from the
least significant — and `matrix` is a square `AbstractMatrix` in `mat[out, in]` convention with
any scalar coefficient already folded in. An empty factor list means the identity.

This is called **once**, by [`compile`](@ref), and never on the hot path. An implementation is
therefore free to do arbitrary work here — the default one resolves Jordan–Wigner strings for
fermionic sites, which is the expensive part of the whole pipeline and precisely what
compiling exists to hoist out of the inner loop.

# Implementing a backend
Distinct terms are summed and factors within a term are multiplied, so an implementation must
distribute sums nested inside products rather than leaving them for the caller.

See also [`amplitude_type`](@ref), [`compile`](@ref).
"""
function expand_terms end

"""
    amplitude_type(operator) -> Type

Element type the matrix elements of `operator` should be accumulated in.

Must be wide enough to hold any product of `operator`'s matrix entries. The default promotes to
at least `Float64` so that integer operator matrices — `PAULI_X` and friends are `Int` — do not
truncate the square roots a bosonic ladder operator introduces.
"""
function amplitude_type end

"""
    read_digit(state, position) -> Int

The digit of the packed `state` at `position`, as a **zero-based** local index: digit `d` means
local value `local_values(spec)[d + 1]`.

Called in the innermost loop of the kernel, so an implementation should be allocation-free and
inlineable.
"""
function read_digit end

"""
    write_digit(state, position, digit) -> state

A copy of the packed `state` with `position` set to the zero-based `digit`.

Packed states are immutable values, so this returns a new one rather than mutating. Called in
the innermost loop of the kernel, so an implementation should be allocation-free and
inlineable.
"""
function write_digit end

"""
    fold_state(state, basis) -> (representative, phase)

Map `state` onto the representative of its symmetry orbit in `basis`, together with the phase
picked up by the symmetry operation that gets there.

Only used by the symmetry-reduced path. See also [`state_lookup`](@ref),
[`orbit_norms`](@ref).
"""
function fold_state end

"""
    state_lookup(basis) -> lookup

Build, once, whatever structure [`lookup_index`](@ref) needs to find a state's position in
`basis`. The default is a `Dict` from state to index.

Hoisting this out of the kernel is the entire point: the pre-compilation version of this
package rebuilt an equivalent dictionary over the whole basis on every single call.
"""
function state_lookup end

"""
    lookup_index(lookup, state) -> Int

Position of `state` in the basis `lookup` was built from, or `0` when it is absent — which is
how a configuration outside the symmetry sector is detected and dropped.
"""
function lookup_index end

lookup_index(lookup::AbstractDict, state) = get(lookup, state, 0)

"""
    orbit_norms(basis) -> AbstractVector

Orbit norms of `basis`, indexed as [`lookup_index`](@ref) indexes it. Matrix elements between
representatives are rescaled by `sqrt(orbit_norms[m] / orbit_norms[n])`.
"""
function orbit_norms end
