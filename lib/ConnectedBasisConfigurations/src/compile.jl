"""
    LocalFactor{T}

One single-site factor of a compiled term, stored column-compressed.

Column `d` (zero-based, the incoming local digit) occupies `colptr[d+1]:colptr[d+2]-1` of
`outs` and `vals`. The kernel therefore visits only the non-zero entries of the one column it
needs, instead of scanning the whole matrix for each input configuration.
"""
struct LocalFactor{T}
    position::Int
    colptr::Vector{Int}
    outs::Vector{Int}       # zero-based output digits
    vals::Vector{T}
end

"""
    CompiledTerm{T}

One product of single-site factors, ready to run.

`factors` act on pairwise distinct positions — same-site factors were multiplied together at
compile time — so they commute and the kernel may apply them in any order. `max_branch` bounds
how many configurations the term can produce from one input.
"""
struct CompiledTerm{T}
    factors::Vector{LocalFactor{T}}
    max_branch::Int
end

"""
    CompiledOperator{T}

An operator flattened into the form the kernel consumes, with everything that does not depend
on the configuration hoisted out: Jordan–Wigner strings resolved, sums distributed over
products, same-site factors multiplied out, matrix elements promoted to the single concrete
element type `T`, and each local matrix column-compressed.

Build one with [`compile`](@ref) and hand it to [`connected_padded`](@ref) or
[`connected`](@ref) instead of the operator itself. Doing so is worth it whenever the same
operator is used more than once, which in variational Monte Carlo is always.

Terms are split by whether they can change the configuration at all. Diagonal terms — `σᶻσᶻ`,
`n_i n_j`, an uncancelled Jordan–Wigner tail — are the majority in a typical lattice
Hamiltonian and are evaluated by a branch-free loop that never touches a state.
"""
struct CompiledOperator{T}
    diagonal::Vector{CompiledTerm{T}}
    offdiagonal::Vector{CompiledTerm{T}}
    max_conn::Int
    max_branch::Int
end

Base.eltype(::CompiledOperator{T}) where {T} = T

function Base.show(io::IO, op::CompiledOperator{T}) where {T}
    print(io, "CompiledOperator{", T, "}(",
        length(op.diagonal), " diagonal + ", length(op.offdiagonal), " off-diagonal terms, ",
        "max_conn ", op.max_conn, ")")
end
Base.show(io::IO, ::MIME"text/plain", op::CompiledOperator) = show(io, op)

"""
    max_conn_size(compiled) -> Int

Upper bound on the number of connected configurations any single configuration can have under
`compiled`, known without looking at a single sample.

This is what [`connected_padded!`](@ref) buffers must be sized by. The bound counts one slot
for the diagonal plus the widest branching each off-diagonal term can produce, so it is exact
for the usual lattice Hamiltonian, where every off-diagonal term is a single hop or flip.
"""
max_conn_size(op::CompiledOperator) = op.max_conn

"""
    compile(operator) -> CompiledOperator
    compile(operator, T::Type) -> CompiledOperator{T}
    compile(operator, basis) -> CompiledSector

Flatten `operator` once into the form the kernel runs on, so that the per-sample work contains
no operator-tree traversal, no Jordan–Wigner expansion, and no dynamic dispatch.

`compile` is the entire reason this package is fast. It is also the only place the operator
backend is consulted — see [`expand_terms`](@ref) — so an operator type from any library works
here as long as it implements that interface.

Passing a `basis` compiles the symmetry-reduced path instead, caching the state-to-index
lookup that would otherwise be rebuilt over the whole basis on every call.

# Example
```julia
compiled = compile(H)
for step in 1:1000
    res = connected_padded(compiled, states)   # no setup cost per step
end
```
"""
function compile end

compile(operator) = compile(operator, amplitude_type(operator))

# Idempotent, so that a caller can accept "an operator, compiled or not" and compile it
# unconditionally.
compile(op::CompiledOperator) = op

function compile(operator, ::Type{T}) where {T<:Number}
    diagonal = CompiledTerm{T}[]
    offdiagonal = CompiledTerm{T}[]
    max_branch = 0

    for raw in expand_terms(operator)
        factors = _merge_factors(raw, T)
        factors === nothing && continue          # a zero factor kills the whole term

        term = CompiledTerm{T}(factors, isempty(factors) ? 1 : prod(_branching, factors))
        if all(_is_diagonal, factors)
            push!(diagonal, term)
        else
            push!(offdiagonal, term)
            max_branch = max(max_branch, term.max_branch)
        end
    end

    # One slot for the diagonal — off-diagonal terms can land back on the input configuration
    # for particular inputs, so the slot is needed even with no diagonal term at all.
    max_conn = 1 + sum(t -> t.max_branch, offdiagonal; init=0)
    return CompiledOperator{T}(diagonal, offdiagonal, max_conn, max_branch)
end

"""
Multiply out a term's same-site factors and drop the identities, returning `nothing` if the
term is identically zero.

Factors on distinct sites commute — Jordan–Wigner strings are already explicit by the time
[`expand_terms`](@ref) returns — so gathering same-site factors is sound as long as their
relative order within the product is preserved, which it is.
"""
function _merge_factors(raw, ::Type{T}) where {T}
    merged = Pair{Int,Matrix{T}}[]
    for (position, mat) in raw
        k = findfirst(p -> first(p) == position, merged)
        if k === nothing
            push!(merged, Int(position) => Matrix{T}(mat))
        else
            merged[k] = position => last(merged[k]) * Matrix{T}(mat)
        end
    end

    # An identity factor contributes nothing and only widens the term; cancelling
    # Jordan-Wigner strings show up here as exactly that.
    filter!(p -> !_isidentity(last(p)), merged)
    any(p -> iszero(last(p)), merged) && return nothing

    return [_compress(position, mat) for (position, mat) in merged]
end

_isidentity(mat::AbstractMatrix) =
    size(mat, 1) == size(mat, 2) &&
    all(mat[i, j] == (i == j) for i in axes(mat, 1), j in axes(mat, 2))

function _compress(position::Integer, mat::Matrix{T}) where {T}
    size(mat, 1) == size(mat, 2) || throw(DimensionMismatch(
        "local matrix at site $position is $(size(mat, 1))×$(size(mat, 2)), not square"
    ))
    dim = size(mat, 1)

    colptr = Vector{Int}(undef, dim + 1)
    outs = Int[]
    vals = T[]
    for d in 1:dim
        colptr[d] = length(vals) + 1
        for j in 1:dim
            v = mat[j, d]
            iszero(v) && continue
            push!(outs, j - 1)
            push!(vals, v)
        end
    end
    colptr[dim+1] = length(vals) + 1

    return LocalFactor{T}(Int(position), colptr, outs, vals)
end

"""Whether every non-zero entry of a compressed factor sits on the diagonal."""
function _is_diagonal(f::LocalFactor)
    for d in 1:(length(f.colptr)-1)
        for p in f.colptr[d]:(f.colptr[d+1]-1)
            f.outs[p] == d - 1 || return false
        end
    end
    return true
end

"""How many outputs a factor can produce from one input: its widest column."""
_branching(f::LocalFactor) =
    maximum(d -> f.colptr[d+1] - f.colptr[d], 1:(length(f.colptr)-1); init=0)

"""
    CompiledSector{T,C,L,N}

A [`CompiledOperator`](@ref) together with the symmetry-reduced basis it is to be evaluated
in: the state-to-index lookup and the orbit norms, both built once.

Matrix elements are complex here even when the operator is real, because a symmetry
representative carries the character of the operation that reached it.
"""
struct CompiledSector{T,C<:CompiledOperator{T},B,L,N}
    operator::C
    basis::B
    lookup::L
    norms::N
end

Base.eltype(::CompiledSector{T}) where {T} = T
max_conn_size(sector::CompiledSector) = max_conn_size(sector.operator)

function Base.show(io::IO, sector::CompiledSector{T}) where {T}
    print(io, "CompiledSector{", T, "}(", sector.operator,
        ", ", length(sector.norms), " representatives)")
end
Base.show(io::IO, ::MIME"text/plain", sector::CompiledSector) = show(io, sector)

function compile(operator, basis)
    T = complex(amplitude_type(operator))
    return CompiledSector(compile(operator, T), basis, state_lookup(basis), orbit_norms(basis))
end
