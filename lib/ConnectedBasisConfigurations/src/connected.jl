# The hot loop. Everything here runs once per sample per optimization step, so it allocates
# nothing per sample: the operator arrived pre-flattened from `compile`, and the branching
# intermediates live in a scratch buffer reused across the whole batch.

"""
Ping-pong buffers for one term's intermediate `(configuration, amplitude)` pairs.

A term is applied one factor at a time, each reading the current pairs and writing the next, so
two buffers suffice however many factors there are.
"""
struct Scratch{S,T}
    a_states::Vector{S}
    a_vals::Vector{T}
    b_states::Vector{S}
    b_vals::Vector{T}
end

function Scratch(op::CompiledOperator{T}, ::Type{S}) where {S,T}
    width = max(op.max_branch, 1)
    return Scratch(
        Vector{S}(undef, width), Vector{T}(undef, width),
        Vector{S}(undef, width), Vector{T}(undef, width),
    )
end

"""
Apply one term to a single configuration, returning `(n, states, vals)` — the number of
`(configuration, amplitude)` pairs produced and the scratch buffers holding them.

Which buffer holds the result depends on the number of factors, so it is returned.
"""
function _run_term(term::CompiledTerm{T}, state::S, scratch::Scratch{S,T}) where {S,T}
    states, vals = scratch.a_states, scratch.a_vals
    spare_states, spare_vals = scratch.b_states, scratch.b_vals

    n = 1
    @inbounds states[1] = state
    @inbounds vals[1] = one(T)

    for f in term.factors
        m = 0
        @inbounds for k in 1:n
            s = states[k]
            amp = vals[k]
            d = read_digit(s, f.position)
            for p in f.colptr[d+1]:(f.colptr[d+2]-1)
                m += 1
                spare_states[m] = write_digit(s, f.position, f.outs[p])
                spare_vals[m] = amp * f.vals[p]
            end
        end
        # An empty column annihilates the configuration -- lowering an already-empty bosonic
        # site, say -- and with it the whole term.
        m == 0 && return (0, states, vals)

        n = m
        states, spare_states = spare_states, states
        vals, spare_vals = spare_vals, vals
    end

    return (n, states, vals)
end

"""
Value a purely diagonal term takes on one configuration: the product of the diagonal entries it
selects. Nothing branches, so no configuration is written.
"""
function _diagonal_value(term::CompiledTerm{T}, state) where {T}
    v = one(T)
    @inbounds for f in term.factors
        d = read_digit(state, f.position)
        lo = f.colptr[d+1]
        lo < f.colptr[d+2] || return zero(T)     # a zero on the diagonal kills the term
        v *= f.vals[lo]
    end
    return v
end

"""
Write every connected configuration of `state` into column `b`, returning how many slots were
used. Slots above the returned count are left untouched; padding is a separate step.
"""
function _fill_column!(
    configs::AbstractMatrix{S}, mels::AbstractMatrix{T}, b::Integer,
    op::CompiledOperator{T}, state::S, scratch::Scratch{S,T}
) where {S,T}
    diagonal = zero(T)
    for term in op.diagonal
        diagonal += _diagonal_value(term, state)
    end

    # Slot 1 is held for the diagonal and filled last, once its value is known.
    k = 1
    for term in op.offdiagonal
        n, states, vals = _run_term(term, state, scratch)
        for t in 1:n
            v = @inbounds vals[t]
            iszero(v) && continue
            s′ = @inbounds states[t]
            # A non-diagonal matrix can still map a particular digit to itself, so this is a
            # property of the sample, not of the term.
            if s′ == state
                diagonal += v
            else
                k += 1
                configs[k, b] = s′
                mels[k, b] = v
            end
        end
    end

    if iszero(diagonal)
        # Close the gap rather than leave an inert row: downstream, every padded row costs one
        # full evaluation of the wavefunction.
        for t in 2:k
            configs[t-1, b] = configs[t, b]
            mels[t-1, b] = mels[t, b]
        end
        return k - 1
    end

    configs[1, b] = state
    mels[1, b] = diagonal
    return k
end

"""Fill slots `from:to` of column `b` with the inert padding: the sample, and a zero mel."""
function _pad_column!(
    configs::AbstractMatrix{S}, mels::AbstractMatrix{T},
    b::Integer, from::Integer, to::Integer, state::S
) where {S,T}
    @inbounds for j in from:to
        configs[j, b] = state
        mels[j, b] = zero(T)
    end
    return nothing
end

"""
    connected(operator, state) -> Dict

Every configuration connected to `state`, mapped to its matrix element
``\\langle s' \\vert \\hat{O} \\vert s \\rangle``.

Contributions from different terms reaching the same configuration are summed, and entries
that cancel to exactly zero are dropped.

`operator` may be a [`CompiledOperator`](@ref), which skips the per-call setup entirely — the
right choice when the same operator is queried repeatedly, as a Hamiltonian-driven Metropolis
rule does twice per step.

Use [`connected_padded`](@ref) for a whole batch: it returns arrays rather than a dictionary
and is what a local energy should be built on.
"""
function connected end

connected(operator, state) = connected(compile(operator), state)

function connected(op::CompiledOperator{T}, state::S) where {S,T}
    out = Dict{S,T}()
    scratch = Scratch(op, S)

    diagonal = zero(T)
    for term in op.diagonal
        diagonal += _diagonal_value(term, state)
    end
    iszero(diagonal) || (out[state] = diagonal)

    for term in op.offdiagonal
        n, states, vals = _run_term(term, state, scratch)
        for t in 1:n
            v = @inbounds vals[t]
            iszero(v) && continue
            s′ = @inbounds states[t]
            out[s′] = get(out, s′, zero(T)) + v
        end
    end

    filter!(p -> !iszero(p.second), out)
    return out
end

"""
    connected_padded(operator, states) -> (; configs, mels, counts)

Connected configurations and matrix elements for a **batch** of packed configurations — the
counterpart of NetKet's `get_conn_padded`, and what a local energy is built from.

`states` may have any shape. `configs` and `mels` gain one leading axis for the connections and
keep the batch shape otherwise, so a vector of `B` states gives `(max_conn, B)`; `counts` has
the shape of `states` and says how many entries per sample are real rather than padding. NetKet
puts the connection axis last; here it comes first, so one sample's connections stay contiguous.

`operator` may be a [`CompiledOperator`](@ref), which skips the per-call flattening.

# Padding

Padded slots repeat the sample itself with a matrix element of **zero**. That keeps `mels` a
concrete numeric array rather than a `Union{T,Missing}` one, and makes

```julia
E_loc(s) = sum(mels[:, b] .* exp.(logψ.(configs[:, b]) .- logψ(s)))
```

correct with no masking: a zero matrix element contributes nothing, and a repeated sample is
something the wavefunction can safely be evaluated on.

# Ordering

Entries appear in term order, the diagonal first when non-zero. Two terms reaching the same
configuration produce two entries rather than one summed entry, which keeps the kernel free of
a per-sample hash table; consumers sum over the connection axis, so only the row count differs.
[`connected`](@ref) gives the summed, deduplicated form.
"""
function connected_padded end

connected_padded(operator, states::AbstractArray) =
    connected_padded(compile(operator), states)

function connected_padded(op::CompiledOperator{T}, states::AbstractArray{S}) where {S,T}
    flat = vec(states)
    n = length(flat)
    height = op.max_conn

    configs = Matrix{S}(undef, height, n)
    mels = Matrix{T}(undef, height, n)
    counts = Vector{Int}(undef, n)
    scratch = Scratch(op, S)

    for b in 1:n
        counts[b] = _fill_column!(configs, mels, b, op, @inbounds(flat[b]), scratch)
    end

    max_conn = isempty(counts) ? 0 : maximum(counts)
    for b in 1:n
        _pad_column!(configs, mels, b, counts[b] + 1, max_conn, @inbounds(flat[b]))
    end

    return _reshape_result(configs, mels, counts, max_conn, size(states))
end

"""
    connected_padded!(configs, mels, counts, compiled, states) -> (; configs, mels, counts)

In-place [`connected_padded`](@ref), writing into caller-owned buffers.

`configs` and `mels` must have `max_conn_size(compiled)` rows and the batch shape of `states`
otherwise; `counts` must have the shape of `states`. Unlike the allocating form the leading
axis is **not** trimmed to the batch's actual maximum, so the same buffers can be reused across
calls whose connection counts differ.

This is the form for an optimization loop. The per-sample work allocates nothing, and the one
scratch buffer is sized by the operator's branching rather than by the batch.
"""
function connected_padded!(
    configs::AbstractArray{S}, mels::AbstractArray{T}, counts::AbstractArray{Int},
    op::CompiledOperator{T}, states::AbstractArray{S}
) where {S,T}
    height = op.max_conn
    expected = (height, size(states)...)
    size(configs) == expected || throw(DimensionMismatch(
        "configs is $(size(configs)); for these states and operator it must be $expected"
    ))
    size(mels) == expected || throw(DimensionMismatch(
        "mels is $(size(mels)); for these states and operator it must be $expected"
    ))
    size(counts) == size(states) || throw(DimensionMismatch(
        "counts is $(size(counts)); it must have the shape of states, $(size(states))"
    ))

    flat = vec(states)
    n = length(flat)
    flat_configs = reshape(configs, height, n)
    flat_mels = reshape(mels, height, n)
    flat_counts = vec(counts)
    scratch = Scratch(op, S)

    for b in 1:n
        s = @inbounds flat[b]
        k = _fill_column!(flat_configs, flat_mels, b, op, s, scratch)
        @inbounds flat_counts[b] = k
        _pad_column!(flat_configs, flat_mels, b, k + 1, height, s)
    end

    return (; configs=configs, mels=mels, counts=counts)
end

"""
    connected_padded(operator, states, basis) -> (; configs, mels, counts)

As above, but for configurations living in a **symmetry-reduced** basis.

Each connected configuration is mapped back to its orbit representative and its matrix element
rescaled by the symmetry character and the orbit-norm ratio `sqrt(norm[m] / norm[n])`.
Configurations whose representative is absent from `basis` fall outside the sector and are
dropped.

Contributions landing on the same representative **are** summed here, unlike the unreduced
path: distinct configurations of one orbit are the same basis state, so leaving them separate
would misreport how many basis states the sector connects.

Matrix elements are complex even for a real operator, because the character need not be. Pass
`compile(operator, basis)` to hoist the state-to-index lookup out of the call.
"""
connected_padded(operator, states::AbstractArray, basis) =
    connected_padded(compile(operator, basis), states)

function connected_padded(sector::CompiledSector{T}, states::AbstractArray{S}) where {S,T}
    op = sector.operator
    flat = vec(states)
    n = length(flat)
    height = op.max_conn

    configs = Matrix{S}(undef, height, n)
    mels = Matrix{T}(undef, height, n)
    counts = Vector{Int}(undef, n)
    scratch = Scratch(op, S)

    for b in 1:n
        counts[b] = _fill_sector_column!(configs, mels, b, sector, @inbounds(flat[b]), scratch)
    end

    max_conn = isempty(counts) ? 0 : maximum(counts)
    for b in 1:n
        _pad_column!(configs, mels, b, counts[b] + 1, max_conn, @inbounds(flat[b]))
    end

    return _reshape_result(configs, mels, counts, max_conn, size(states))
end

function _fill_sector_column!(
    configs::AbstractMatrix{S}, mels::AbstractMatrix{T}, b::Integer,
    sector::CompiledSector{T}, state::S, scratch::Scratch{S,T}
) where {S,T}
    op = sector.operator
    norms = sector.norms

    n = lookup_index(sector.lookup, state)
    n == 0 && throw(ArgumentError(
        "sample $b is not a representative state of the given basis"
    ))

    diagonal = zero(T)
    for term in op.diagonal
        diagonal += _diagonal_value(term, state)
    end
    k = _accumulate_folded!(configs, mels, b, 0, sector, n, state, diagonal)

    for term in op.offdiagonal
        m, states, vals = _run_term(term, state, scratch)
        for t in 1:m
            k = _accumulate_folded!(
                configs, mels, b, k, sector, n, @inbounds(states[t]), @inbounds(vals[t])
            )
        end
    end

    # Cancellation between orbit members is common; compact the exact zeros out rather than
    # spend a wavefunction evaluation on each.
    kept = 0
    for j in 1:k
        v = mels[j, b]
        iszero(v) && continue
        kept += 1
        configs[kept, b] = configs[j, b]
        mels[kept, b] = v
    end
    return kept
end

"""
Fold one connected configuration onto its orbit representative and add it to column `b`,
returning the updated number of entries.

Many configurations fold onto one representative, so entries are summed. The search is linear
because a column holds only a handful, where scanning beats hashing.
"""
function _accumulate_folded!(
    configs::AbstractMatrix{S}, mels::AbstractMatrix{T}, b::Integer, k::Int,
    sector::CompiledSector{T}, n::Int, s′::S, v
) where {S,T}
    iszero(v) && return k

    repr, phase = fold_state(s′, sector.basis)
    m = lookup_index(sector.lookup, repr)
    m == 0 && return k                           # outside this symmetry sector

    w = v * phase * sqrt(sector.norms[m] / sector.norms[n])
    for j in 1:k
        if configs[j, b] == repr
            mels[j, b] += w
            return k
        end
    end

    configs[k+1, b] = repr
    mels[k+1, b] = w
    return k + 1
end

"""
Trim the connection axis to what the batch actually used, and restore the batch shape.

A batch that reaches the operator's compile-time bound — the usual case for a lattice
Hamiltonian — reshapes without copying.
"""
function _reshape_result(
    configs::Matrix{S}, mels::Matrix{T}, counts::Vector{Int},
    max_conn::Integer, batch::Dims
) where {S,T}
    trim(a) = size(a, 1) == max_conn ? a : a[1:max_conn, :]
    return (;
        configs=reshape(trim(configs), max_conn, batch...),
        mels=reshape(trim(mels), max_conn, batch...),
        counts=reshape(counts, batch),
    )
end
