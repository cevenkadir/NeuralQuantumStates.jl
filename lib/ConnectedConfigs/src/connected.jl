"""
    _apply_packed!(out, op, states)

Accumulate `op * states` into `out`, where both map a packed `BaseInt` configuration to an
amplitude.

This is the single-site kernel of the whole package. It mirrors the reference implementation in
OperatorAlgebra's `OperatorAlgebraSymBasisExt`, reimplemented here so it can be driven over a
batch and so the hot loop is under this package's control.

Two input configurations differing only at `op`'s site can be carried onto the same output
configuration whenever a row of `op.mat` has more than one entry, so contributions must be
**accumulated** rather than assigned. That accumulation is what makes the intermediate index of
a product of two operators on the same site get summed over.
"""
function _apply_packed!(out::Dict{S,T}, op::Op, states::Dict{S,T}) where {S,T}
    # `rawsite` strips any fermionic tag, leaving the plain site identifier, which is used
    # directly as the digit position.
    idx = OperatorAlgebra.rawsite(op.site)
    mat = op.mat
    for (s, amp) in states
        d = Int(read(s, idx))                    # 0-based digit currently at this site
        @inbounds for j in axes(mat, 1)
            v = mat[j, d+1]
            iszero(v) && continue
            s′ = write(s, idx, j - 1)
            out[s′] = get(out, s′, zero(T)) + v * amp
        end
    end
    return out
end

function _apply_packed(op::Op, states::Dict{S,T}) where {S,T}
    return _apply_packed!(Dict{S,T}(), op, states)
end

function _apply_packed(chain::OpChain, states::Dict{S,T}) where {S,T}
    # OpChain([A, B]) is the matrix product A*B, so the rightmost factor acts first.
    for op in Iterators.reverse(chain.ops)
        states = _apply_packed(op, states)
    end
    return states
end

function _apply_packed(sum_op::OpSum, states::Dict{S,T}) where {S,T}
    out = Dict{S,T}()
    for op in sum_op.ops
        for (s, v) in _apply_packed(op, states)
            out[s] = get(out, s, zero(T)) + v
        end
    end
    return out
end

"""
    _amplitude_type(H)

Element type of the amplitudes `H` produces, promoted to at least `Float64` so that integer
operator matrices (`PAULI_X` and friends are `Int`) do not truncate square roots.
"""
_amplitude_type(H::AbstractOp) = promote_type(float(eltype(H)), Float64)

"""
    connected(H, state) -> Dict{S,T}

Every configuration connected to `state` by `H`, mapped to its matrix element
``\\langle s' \\vert H \\vert s \\rangle``.

Contributions from different terms of `H` reaching the same configuration are summed, and
entries that cancel to exactly zero are dropped — a vanishing matrix element contributes
nothing to any observable, and keeping it would only waste a slot in the padded output.

Use [`connected_padded`](@ref) for a whole batch of states.
"""
function connected(H::AbstractOp, state::S) where {S}
    T = _amplitude_type(H)
    # Resolve Jordan-Wigner strings for fermionic sites. For an ordinary commuting basis this
    # returns `H` unchanged.
    flat = OperatorAlgebra._jw_expand(H, basis_info(H))
    result = _apply_packed(flat, Dict{S,T}(state => one(T)))
    filter!(p -> !iszero(p.second), result)
    return result
end

"""
    connected_padded(H, states) -> (; configs, mels, counts)

Connected configurations and matrix elements for a **batch** of packed configurations — the
counterpart of NetKet's `get_conn_padded`, and the kernel a variational Monte Carlo local
energy is built from.

# Returns
A named tuple of
- `configs::Matrix{S}`: `(max_connections, batch)` connected configurations.
- `mels::Matrix{T}`: `(max_connections, batch)` matrix elements.
- `counts::Vector{Int}`: how many entries of each column are real rather than padding.

# Padding

Different configurations have different numbers of connections, so the columns are padded to a
common height. Padded slots repeat the sample itself and carry a matrix element of **zero**,
rather than being marked `missing`. That choice matters: it keeps `mels` a concrete numeric
array instead of a `Union{T,Missing}` one, and it makes the local energy

```julia
E_loc(s) = sum(mels[:, b] .* exp.(logψ.(configs[:, b]) .- logψ(s)))
```

correct with no masking at all, because a zero matrix element contributes nothing.

# Ordering

Within a column, entries are sorted by their packed configuration value. Nothing physical
depends on the order, but a deterministic one makes results reproducible and diffable rather
than dependent on hash iteration order.
"""
function connected_padded(H::AbstractOp, states::AbstractVector{S}) where {S}
    T = _amplitude_type(H)
    flat = OperatorAlgebra._jw_expand(H, basis_info(H))

    per_state = Vector{Vector{Pair{S,T}}}(undef, length(states))
    for (b, s) in pairs(states)
        d = _apply_packed(flat, Dict{S,T}(s => one(T)))
        filter!(p -> !iszero(p.second), d)
        entries = collect(d)
        sort!(entries; by=p -> first(p).value)
        per_state[b] = entries
    end

    counts = length.(per_state)
    max_conn = isempty(counts) ? 0 : maximum(counts)

    configs = Matrix{S}(undef, max_conn, length(states))
    mels = zeros(T, max_conn, length(states))
    for (b, entries) in pairs(per_state)
        for (j, (s′, v)) in pairs(entries)
            configs[j, b] = s′
            mels[j, b] = v
        end
        # Pad with the sample itself; the zero matrix element makes the slot inert.
        for j in (length(entries)+1):max_conn
            configs[j, b] = states[b]
        end
    end

    return (; configs=configs, mels=mels, counts=counts)
end

"""
    connected_padded(H, states, basis) -> (; configs, mels, counts)

As above, but for configurations living in a **symmetry-reduced** `SymBasis.Basis`.

Each connected configuration is mapped back to the representative of its symmetry orbit, and
its matrix element is rescaled by the character of the symmetry operation together with the
ratio of orbit norms — the standard factor `sqrt(norm[m] / norm[n])`. Configurations whose
representative is absent from `basis` fall outside the sector and are dropped.

Contributions landing on the same representative are summed after rescaling, which is
essential: distinct configurations in the same orbit are the same basis state here.
"""
function connected_padded(
    H::AbstractOp, states::AbstractVector{S}, basis::SymBasis.Bases.Basis{S}
) where {S}
    T = complex(_amplitude_type(H))
    flat = OperatorAlgebra._jw_expand(H, basis_info(H))

    index_of = Dict(s => i for (i, s) in pairs(basis.states))

    per_state = Vector{Vector{Pair{S,T}}}(undef, length(states))
    for (b, s) in pairs(states)
        n = get(index_of, s, 0)
        n == 0 && throw(ArgumentError(
            "sample $b is not a representative state of the given basis"
        ))

        raw = _apply_packed(flat, Dict{S,T}(s => one(T)))

        folded = Dict{S,T}()
        for (s′, v) in raw
            repr, phase = representative(s′, basis)
            m = get(index_of, repr, 0)
            m == 0 && continue                      # outside this symmetry sector
            factor = sqrt(basis.norms[m] / basis.norms[n])
            folded[repr] = get(folded, repr, zero(T)) + v * phase * factor
        end

        filter!(p -> !iszero(p.second), folded)
        entries = collect(folded)
        sort!(entries; by=p -> first(p).value)
        per_state[b] = entries
    end

    counts = length.(per_state)
    max_conn = isempty(counts) ? 0 : maximum(counts)

    configs = Matrix{S}(undef, max_conn, length(states))
    mels = zeros(T, max_conn, length(states))
    for (b, entries) in pairs(per_state)
        for (j, (s′, v)) in pairs(entries)
            configs[j, b] = s′
            mels[j, b] = v
        end
        for j in (length(entries)+1):max_conn
            configs[j, b] = states[b]
        end
    end

    return (; configs=configs, mels=mels, counts=counts)
end
