# Default operator backend: OperatorAlgebra's `Op` / `OpChain` / `OpSum` tree.

amplitude_type(operator::AbstractOp) = promote_type(float(eltype(operator)), Float64)

"""Factors of one term: digit position paired with its local matrix."""
const TermFactors = Vector{Pair{Int,AbstractMatrix}}

# `rawsite` strips a fermionic tag, leaving the plain site identifier, which doubles as the
# digit position. Anything that is not an integer site has no digit to address.
function _digit_position(site)
    raw = OperatorAlgebra.rawsite(site)
    raw isa Integer || throw(ArgumentError(
        "site $site is not an integer; ConnectedConfigs addresses digit positions directly, " *
        "so sites must be integers 1:nsites"
    ))
    return Int(raw)
end

_flatten(op::Op) = [TermFactors([_digit_position(op.site) => op.mat])]

# `append!` rather than `reduce(vcat, ...)`: a lattice Hamiltonian is a sum of O(nsites) terms,
# and repeated `vcat` would copy the accumulated list once per term.
function _flatten(os::OpSum)
    terms = TermFactors[]
    for o in os.ops
        append!(terms, _flatten(o))
    end
    return terms
end

function _flatten(oc::OpChain)
    # A chain may contain nested sums, which have to be distributed here rather than left for
    # the kernel: (A + B) * C is two terms, and the kernel only ever sees terms.
    terms = [TermFactors()]
    for factor in oc.ops
        expanded = _flatten(factor)
        # The overwhelmingly common case is a plain product, where each factor contributes
        # exactly one alternative and the term list neither grows nor needs rebuilding.
        if length(expanded) == 1
            only_factors = only(expanded)
            for t in terms
                append!(t, only_factors)
            end
        else
            terms = TermFactors[vcat(t, f) for t in terms for f in expanded]
        end
    end
    return terms
end

"""Whether any site carries statistics that make Jordan-Wigner strings necessary."""
_needs_strings(op::Op) =
    !(OperatorAlgebra.exchange_style(op.site) isa OperatorAlgebra.Commuting)
_needs_strings(os::OpSum) = any(_needs_strings, os.ops)
_needs_strings(oc::OpChain) = any(_needs_strings, oc.ops)

function expand_terms(operator::AbstractOp)
    # Resolve Jordan-Wigner strings for fermionic sites. Deciding whether they are needed by
    # walking the sites is much cheaper than asking `basis_info`, which has to work out every
    # site's local dimension as well -- and for the spin and boson models that dominate use,
    # the answer is no.
    _needs_strings(operator) || return _flatten(operator)
    return _flatten(OperatorAlgebra._jw_expand(operator, basis_info(operator)))
end
