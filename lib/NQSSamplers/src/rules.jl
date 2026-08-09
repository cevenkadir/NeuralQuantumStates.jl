"""
    AbstractRule

A Metropolis transition rule: how to propose a move from one configuration to another.

The rule is where the physics knowledge lives: which moves are proposed decides whether the
chain explores the relevant part of configuration space, and whether it stays inside a
conserved-quantity sector.

# Interface
    propose(rule, state, dof, nsites, rng) -> (new_state, log_correction)

`log_correction` is `log[T(s|s') / T(s'|s)]`, the asymmetry of the proposal. It is zero for a
symmetric rule and must be accounted for otherwise, or the chain converges to the wrong
distribution.
"""
abstract type AbstractRule end

"""
    propose(rule, state, dof, nsites, rng) -> (new_state, log_correction)

Propose a move away from `state`.
"""
function propose end

"""
    LocalRule() <: AbstractRule

Change the value at a single randomly chosen site to a different randomly chosen local value.

The workhorse rule for unconstrained systems, and symmetric, so the proposal correction is zero.

It does **not** conserve particle number or magnetization: inside a conserved sector every
proposal leaves the sector, is rejected by the basis check, and the chain stalls. Use
[`ExchangeRule`](@ref) there.
"""
struct LocalRule <: AbstractRule end

function propose(::LocalRule, s, dof, nsites::Integer, rng::AbstractRNG)
    B = local_dimension(dof)
    site = rand(rng, 1:nsites)
    current = Int(read(s, site))
    # Uniform over the other B-1 values, so the move is never a no-op.
    d = rand(rng, 0:(B-2))
    d >= current && (d += 1)
    return write(s, site, d), 0.0
end

"""
    ExchangeRule() <: AbstractRule

Swap the values on two randomly chosen sites.

Conserves any quantity that is a sum over sites — magnetization, particle number — so this is
the rule for a symmetry sector, where [`LocalRule`](@ref) would propose nothing but rejections.
Symmetric, so the proposal correction is zero.

!!! note "It cannot leave its sector"
    That is the point, but it also means the chain only ever explores the sector it starts in.
    The initial configuration therefore determines which sector is sampled.
"""
struct ExchangeRule <: AbstractRule end

function propose(::ExchangeRule, s, dof, nsites::Integer, rng::AbstractRNG)
    i = rand(rng, 1:nsites)
    j = rand(rng, 1:nsites)
    i == j && return s, 0.0
    a, b = Int(read(s, i)), Int(read(s, j))
    a == b && return s, 0.0
    return write(write(s, i, b), j, a), 0.0
end

"""
    HamiltonianRule(operator) <: AbstractRule

Propose moves along the configurations the Hamiltonian actually connects to.

Never proposes a move the Hamiltonian gives zero amplitude to, so it respects whatever the
Hamiltonian conserves without being told what that is, and for a sparse Hamiltonian explores far
better than blind local moves.

The price is that the proposal is **asymmetric** — `s` and `s'` generally have different numbers
of connections — so the correction `log[n_conn(s) / n_conn(s')]` is non-zero. Dropping it biases
the sampled distribution while leaving everything else looking reasonable.

The operator is compiled when the rule is constructed, since every step queries it twice.
"""
struct HamiltonianRule{O} <: AbstractRule
    operator::O

    HamiltonianRule{O}(operator::O) where {O} = new{O}(operator)
end

function HamiltonianRule(operator)
    compiled = ConnectedBasisConfigurations.compile(operator)
    return HamiltonianRule{typeof(compiled)}(compiled)
end

function propose(rule::HamiltonianRule, s, dof, nsites::Integer, rng::AbstractRNG)
    forward = connected(rule.operator, s)
    delete!(forward, s)
    isempty(forward) && return s, 0.0

    candidates = collect(keys(forward))
    s′ = candidates[rand(rng, 1:length(candidates))]

    backward = connected(rule.operator, s′)
    delete!(backward, s′)
    n_back = length(backward)
    n_back == 0 && return s, 0.0

    # T(s'|s) = 1/n_forward and T(s|s') = 1/n_backward.
    return s′, log(length(candidates) / n_back)
end
