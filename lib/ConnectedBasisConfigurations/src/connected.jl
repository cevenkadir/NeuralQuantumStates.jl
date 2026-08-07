# The hot loop. Everything here runs once per sample per optimization step, so it allocates
# nothing per sample: the operator arrived pre-flattened from `compile`, and the branching
# intermediates live in a scratch buffer reused across the whole batch.

"""
Ping-pong buffers for one term's intermediate `(configuration, amplitude)` pairs.

A term is applied one factor at a time, each factor reading the current pairs and writing the
next ones, so two buffers are enough no matter how many factors there are.
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

Which of the two buffers ends up holding the result depends on how many factors the term has,
so it is returned rather than assumed.
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
Value a purely diagonal term takes on one configuration.

Every factor is diagonal, so nothing branches and no configuration is ever written: the term
contributes one number, the product of the diagonal entries the configuration selects.
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
counterpart of NetKet's `get_conn_padded`, and the kernel a variational Monte Carlo local
energy is built from.

`states` may be an array of any shape. `configs` and `mels` gain one leading axis for the
connections and keep the batch shape otherwise, so a vector of `B` states gives
`(max_conn, B)` and an `(A, B)` matrix gives `(max_conn, A, B)`. `counts` has the shape of
`states`. NetKet puts the connection axis last; here it comes first, because Julia is
column-major and this is what keeps one sample's connections contiguous.

`operator` may be a [`CompiledOperator`](@ref) from [`compile`](@ref), which skips the
flattening work; passing the operator itself compiles it on every call.

# Returns
- `configs`: `(max_conn, size(states)...)` connected configurations.
- `mels`: `(max_conn, size(states)...)` matrix elements.
- `counts`: how many entries per sample are real rather than padding.

# Padding

Different configurations have different numbers of connections, so the leading axis is padded
to a common length. Padded slots repeat the sample itself and carry a matrix element of
**zero**, following NetKet. That choice matters twice over: it keeps `mels` a concrete numeric
array rather than a `Union{T,Missing}` one, and it makes the local energy

```julia
E_loc(s) = sum(mels[:, b] .* exp.(logψ.(configs[:, b]) .- logψ(s)))
```

correct with no masking at all, because a zero matrix element contributes nothing — and a
repeated sample is a configuration the wavefunction can safely be evaluated on.

# Ordering, and repeated configurations

Entries appear in term order: the diagonal first when it is non-zero, then each off-diagonal
term's contribution. Two terms reaching the **same** configuration produce two entries rather
than one summed entry, which is what NetKet does and what keeps the kernel free of any
per-sample hash table. Every consumer sums over the connection axis, so the result is
unchanged; only the row count differs. Use [`connected`](@ref) when you want the summed,
deduplicated matrix elements.
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
otherwise; `counts` must have the shape of `states`. Unlike the allocating form, the leading
axis is **not** trimmed to the batch's actual maximum — the buffers keep their full height,
with every slot above a sample's count padded inert — so the same buffers can be reused across
calls whose connection counts differ.

This is the form to use inside an optimization loop, where the batch shape is fixed and
allocating a fresh pair of arrays per step is pure overhead. The per-sample work allocates
nothing; what remains is one small scratch buffer whose size is set by the operator's
branching and **not** by the batch, so the cost per step stops growing with the batch size.
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

Each connected configuration is mapped back to the representative of its symmetry orbit, and
its matrix element is rescaled by the character of the symmetry operation together with the
ratio of orbit norms — the standard factor `sqrt(norm[m] / norm[n])`. Configurations whose
representative is absent from `basis` fall outside the sector and are dropped.

Contributions landing on the same representative **are** summed here, unlike the unreduced
path: distinct configurations of one orbit are the same basis state, so leaving them separate
would not merely be redundant, it would misreport how many basis states the sector connects.

Matrix elements are complex even for a real operator, because the character need not be.

Pass `compile(operator, basis)` instead of `operator` to hoist the state-to-index lookup out of
the call.
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

    # Cancellation between orbit members is real and common, so compact the exact zeros out
    # rather than spend a wavefunction evaluation on each.
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
Fold one connected configuration back onto its orbit representative and add it to column `b`,
returning the updated number of entries.

Folding maps many configurations onto one representative, so entries must be summed. The
search for an existing entry is linear because a column holds a handful of them, where
scanning beats hashing outright.
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

The buffers were allocated at the operator's compile-time bound, which for a lattice
Hamiltonian is usually exactly what the batch reaches — so the common case is a reshape with no
copy at all, and only a batch that falls short of the bound pays for one.
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
