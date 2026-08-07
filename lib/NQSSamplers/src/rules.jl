"""
    AbstractRule

A Metropolis transition rule: how to propose a move from one configuration to another.

Separating the rule from the Metropolis machinery is what NetKet does, and the reason is that
the rule is where all the physics knowledge lives. Which moves are proposed determines whether
the chain explores the relevant part of configuration space at all — and, critically, whether
it stays inside a conserved-quantity sector.

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

The workhorse rule for unconstrained systems. It is symmetric — every move has an equally
likely reverse — so the proposal correction is zero.

It does **not** conserve particle number or magnetization. Using it inside a conserved sector
proposes moves that leave the sector; those get rejected by the basis check, and the chain
stalls. Use [`ExchangeRule`](@ref) there instead.
"""
struct LocalRule <: AbstractRule end

function propose(::LocalRule, s, dof, nsites::Integer, rng::AbstractRNG)
    B = local_dimension(dof)
    site = rand(rng, 1:nsites)
    current = Int(read(s, site))
    # Pick uniformly among the other B-1 values, so the move is never a no-op.
    d = rand(rng, 0:(B-2))
    d >= current && (d += 1)
    return write(s, site, d), 0.0
end

"""
    ExchangeRule() <: AbstractRule

Swap the values on two randomly chosen sites.

Conserves any quantity that is a sum over sites — total magnetization, total particle number —
because swapping leaves that sum untouched. This is the rule to use inside a symmetry sector,
where [`LocalRule`](@ref) would propose nothing but rejections.

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

The most physically informed rule available: it never proposes a move the Hamiltonian gives
zero amplitude to, so it automatically respects whatever the Hamiltonian conserves without
being told what that is. For a sparse Hamiltonian it explores far more efficiently than blind
local moves.

The price is that the proposal is **asymmetric** — `s` may have a different number of
connections than `s'` — so the correction `log[n_conn(s) / n_conn(s')]` is non-zero and must be
carried into the acceptance test. Dropping it silently biases the sampled distribution, which
is exactly the kind of error that produces plausible-looking but wrong expectation values.

The operator is compiled once, when the rule is constructed. Every Metropolis step queries it
twice — forwards and backwards — so leaving that work in the step is the difference between
flattening the Hamiltonian a handful of times and flattening it a few million.
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
    delete!(forward, s)                              # a move must go somewhere else
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
