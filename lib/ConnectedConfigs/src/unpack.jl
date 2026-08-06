"""
    configurations(spec, states, nsites) -> Array

Unpack `BaseInt` configurations into an array of **physical local values** — magnetic quantum
numbers for a spin, occupation numbers for a boson — as taken from `spec`.

Packed integers are the right representation for the operator kernel: compact, hashable, and
cheap to permute. A neural network wants the opposite, an explicit numeric vector per
configuration, and this is the boundary between the two.

The degree-of-freedom axis comes **first**: a vector of `M` states becomes `(nsites, M)`, and
an `(M, batch)` matrix as returned by [`connected_padded`](@ref) becomes
`(nsites, M, batch)`. Julia is column-major, so this puts each configuration in contiguous
memory, which is the layout a batched forward pass wants.

# Example
```julia
spec = Boson(3)
res = connected_padded(H, states)
x = configurations(spec, res.configs, nsites)   # (nsites, max_conn, batch)
```
"""
function configurations end

function configurations(spec, state::S, nsites::Integer) where {S}
    values = local_values(spec)
    return [values[Int(read(state, i))+1] for i in 1:nsites]
end

function configurations(spec, states::AbstractVector{S}, nsites::Integer) where {S}
    values = local_values(spec)
    out = Matrix{eltype(values)}(undef, nsites, length(states))
    for (j, s) in pairs(states)
        @inbounds for i in 1:nsites
            out[i, j] = values[Int(read(s, i))+1]
        end
    end
    return out
end

function configurations(spec, states::AbstractMatrix{S}, nsites::Integer) where {S}
    values = local_values(spec)
    m, b = size(states)
    out = Array{eltype(values),3}(undef, nsites, m, b)
    for k in 1:b, j in 1:m
        s = states[j, k]
        @inbounds for i in 1:nsites
            out[i, j, k] = values[Int(read(s, i))+1]
        end
    end
    return out
end

"""
    packed(spec, configuration) -> BaseInt

Inverse of [`configurations`](@ref) for a single configuration: pack a vector of physical local
values back into a `BaseInt`.

Throws an `ArgumentError` if a value is not one of `spec`'s local values, which catches the
common mistake of handing over a configuration built with the wrong convention — spins as
`0, 1` rather than `-1//2, 1//2`, say.
"""
function packed(spec, configuration::AbstractVector; T::Type=UInt, Ti::Type=Int)
    B = local_dimension(spec)

    value = zero(T)
    power = one(T)
    for (i, v) in pairs(configuration)
        d = _digit_of(spec, v)
        (d === nothing || d < 0 || d >= B) && throw(ArgumentError(
            "$v at site $i is not a local value of $spec (allowed: $(local_values(spec)))"
        ))
        value += T(d) * power
        power *= T(B)
    end
    return BaseInt{T,Ti,B}(value)
end

"""
    _digit_of(spec, value) -> Union{Int,Nothing}

Zero-based local digit holding `value`, or `nothing` when `spec` has no such local state.

The generic fallback searches `local_values`, which is what any new specification type gets
for free. `Spin` and `Boson` skip the search: their local values are an arithmetic sequence,
so the digit is a subtraction away, and `packed` sits on the boundary every batch of network
inputs crosses.
"""
_digit_of(spec, value) = something(findfirst(==(value), local_values(spec)), 0) - 1

function _digit_of(spec::Spin, value)
    d = value + spec.s
    return isinteger(d) ? Int(d) : nothing
end

function _digit_of(::Boson, value)
    return isinteger(value) ? Int(value) : nothing
end
