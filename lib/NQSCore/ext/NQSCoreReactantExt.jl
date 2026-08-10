"""
A whole variational step compiled by XLA, behind [`NQSCore.Compiled`](@ref).

`expect_and_grad` is five operations — unpack the samples, evaluate the network, compute the
connected configurations, reduce them into local energies, differentiate a scalar — and Zygote
runs each of them as its own sequence of kernels. Compiled together, the gradient is 7x faster
than Zygote on CUDA on the same card, and the step 1.7x.

# Three regions, not one

The connected-configuration kernel does not raise to StableHLO: its loops are data-dependent, and
Reactant answers `cannot raise op to stablehlo` on the `scf.for`. Raising is what would let XLA
fuse across it — it is needed for differentiation and for fusion, not for running — so the step is
three regions rather than one. Nothing forces them together: only `energy`'s output crosses to
`gradient`, and connected configurations are data that nothing differentiates.

# The parameters stay on the host

They are ordinary arrays and this uploads what a step needs, per step. That is not a compromise
for `MCState`, it is the reason it works: a sampler evaluates the ansatz once per sweep step, and
those calls belong to no compiled region — held as Reactant arrays they do not run slowly, they
fail. Sixty kilobytes per step against three milliseconds of compute is a good trade, and it keeps
the sampler, the optimizer and the parameters themselves untouched.
"""
module NQSCoreReactantExt

using ConnectedBasisConfigurations: ConnectedBasisConfigurations, flatten, max_conn_size
using Enzyme: Enzyme, Const
using NQSCore
using NQSCore: Compiled
using Reactant
using SymBasis.DigitBase: BaseInt

# ------------------------------------------------------------------------------ the regions

# Top-level functions of explicit arguments, which is the form Reactant caches on. A closure over
# the batch would hand it over as a traced constant and recompile on every call.

"""
Connected configurations and both batches, from packed states and an operator's fields.

`op` is a `NamedTuple` and not a `FlatOperator`: that type's parameters demand
`AbstractVector{Int32}`, and a `TracedRArray{Int32,1}` has element type `TracedRNumber{Int32}`, so
the struct cannot cross into a compiled region at all. Its contents can.

The states arrive as the raw integer a `BaseInt` wraps, for the same reason — `BaseInt` is not an
XLA element type, and `to_rarray` passes an array of them through unconverted rather than refusing
it, which is a silence worth knowing about.
"""
function _prepare(op, states, values_conn, values_sample, height::Int, nsites::Int, base::Val)
    n = length(states)
    backend = Reactant.KernelAbstractions.get_backend(states)

    configs = similar(states, height, n)
    mels = similar(op.vals, height, n)
    counts = similar(states, Int, n)
    ConnectedBasisConfigurations.connected_padded!(
        configs, mels, counts, op, states, backend; base=base
    )

    # The connected configurations are unpacked real and the samples in the network's own type,
    # which is what `configurations_of` and `device_connections` do: only the sample batch is
    # followed by a derivative, and only there does a mixed complex-real product cost anything.
    x_conn = similar(values_conn, nsites, height * n)
    ConnectedBasisConfigurations.configurations!(
        x_conn, values_conn, configs, nsites, backend; base=base
    )
    xs = similar(values_sample, nsites, n)
    ConnectedBasisConfigurations.configurations!(
        xs, values_sample, states, nsites, backend; base=base
    )
    return xs, x_conn, mels
end

"""Local energies, weights and the cotangent — everything the gradient needs except `θ`."""
function _energy(a, θ, xs, x_conn, mels, born::Bool)
    logψ_s = log_amplitude(a, θ, xs)
    logψ_sp = reshape(log_amplitude(a, θ, x_conn), size(mels))
    # The column reduces without a mask because the padding is inert: a padded slot repeats the
    # sample with a zero matrix element.
    E = vec(sum(mels .* exp.(logψ_sp .- transpose(logψ_s)); dims=1))
    # `FullSumState` weights by the Born probabilities; Monte Carlo samples already carry them.
    p = born ? NQSCore.born_probabilities(logψ_s) : nothing
    return E, p, NQSCore._gradient_cotangent(E, p)
end

_loss(a, x, c, θ) = NQSCore._gradient_loss(a, θ, x, c)

"""
The gradient with respect to `θ`, everything before it held constant.

`θ` is a `NamedTuple` and stays one: Enzyme returns a structure gradient, and that is measured to
be the same gradient the library's real/imag split produces — Zygote to 0.0, Reactant to
2.206e-15 — so the flattening the host path does is not needed, which is fortunate, because
`ComponentArrays` scalar-indexes and a compiled region will not have it.
"""
_gradient(a, x, c, θ) =
    Enzyme.gradient(Enzyme.Reverse, Const(_loss), Const(a), Const(x), Const(c), θ)[end]

# -------------------------------------------------------------------------------- the cache

"""
Compiled thunks and the uploaded operator, by shape.

Reactant compiles per shape and the regions here take 14-18 seconds each, so a step that compiled
every time would be slower than the path it replaces by four orders of magnitude. A different
shape mints a different key and nothing is ever invalidated: the cache is global, it grows with
the number of distinct shapes a session uses, and [`NQSCore.clear_compiled_cache!`](@ref) empties
it. Setting `persistent_cache_enabled` in `LocalPreferences.toml` makes the first compilation a
once-per-machine cost rather than a once-per-session one.
"""
const CACHE = Dict{Any,Any}()

NQSCore.clear_compiled_cache!() = (empty!(CACHE); nothing)

"""What a compilation depends on, and nothing else."""
_key(a, θ, states, op, born) = (
    typeof(a), size(states), keys(θ), map(size, values(θ)), map(eltype, values(θ)),
    op.max_conn, op.n_terms, op.max_branch, eltype(op.vals), born,
)

"""
Everything uploaded and compiled for this shape, built once.

The operator is uploaded here rather than per call. Seven small transfers is not much beside a
millisecond, but an optimisation run pays them once per step for as many steps as it takes.
"""
function _prepared(vs, operator, states, born::Bool)
    a, θ = NQSCore.ansatz(vs), NQSCore.parameters(vs)
    op = flatten(operator)
    key = _key(a, θ, states, op, born)
    get!(CACHE, key) do
        S = eltype(states)
        V, B = S.parameters[1], S.parameters[3]
        base = Val(B)
        height, n = max_conn_size(op), length(states)
        nsites = NQSCore.n_sites(a)

        # Precision follows the parameters, as it does everywhere else. The connected batch is
        # real because nothing differentiates it; the sample batch takes the network's own type.
        Tψ = NQSCore.input_type(a, θ)
        Tψ === nothing && (Tψ = ComplexF64)
        values_conn = Reactant.to_rarray(collect(real(Tψ), ConnectedBasisConfigurations.local_values(NQSCore.dof(a))))
        values_sample = Reactant.to_rarray(collect(Tψ, ConnectedBasisConfigurations.local_values(NQSCore.dof(a))))
        op_ra = Reactant.to_rarray((;
            op.term_start, op.factor_position, op.factor_col_start,
            op.colptr, op.outs, op.vals, op.n_diagonal, op.n_terms,
            op.max_conn, op.max_branch,
        ))

        probe_states = Reactant.to_rarray(collect(reinterpret(V, vec(states))))
        prep_args = (op_ra, probe_states, values_conn, values_sample, height, nsites, base)
        prep = Reactant.compile(_prepare, prep_args; sync=true)
        xs, x_conn, mels = prep(prep_args...)

        θ_ra = Reactant.to_rarray(θ)
        energy = Reactant.compile(_energy, (a, θ_ra, xs, x_conn, mels, born); sync=true)
        _, _, c = energy(a, θ_ra, xs, x_conn, mels, born)
        grad = Reactant.compile(_gradient, (a, xs, c, θ_ra); sync=true)
        (; prep, energy, grad, op_ra, values_conn, values_sample, height, nsites, base, V, born)
    end
end

"""Run the regions that produce the local energies, and hand back what a caller needs."""
function _run(vs, operator, states, born::Bool)
    a, θ = NQSCore.ansatz(vs), NQSCore.parameters(vs)
    P = _prepared(vs, operator, states, born)
    st = Reactant.to_rarray(collect(reinterpret(P.V, vec(states))))
    θ_ra = Reactant.to_rarray(θ)
    xs, x_conn, mels = P.prep(P.op_ra, st, P.values_conn, P.values_sample,
                              P.height, P.nsites, P.base)
    E, p, c = P.energy(a, θ_ra, xs, x_conn, mels, born)
    return (; a, θ_ra, xs, E=Array(E), p=(p === nothing ? nothing : Array(p)), c, P)
end

# ------------------------------------------------------------------------------- the seams

function NQSCore.compiled_expect(vs, operator, states, ::Compiled)
    born = vs isa NQSCore.FullSumState
    r = _run(vs, operator, states, born)
    return born ? NQSCore.weighted_statistics(r.E, r.p) :
           NQSCore.statistics(reshape(r.E, size(states)))
end

function NQSCore.compiled_expect_and_grad(vs, operator, states, ::Compiled)
    born = vs isa NQSCore.FullSumState
    r = _run(vs, operator, states, born)
    ∇ = map(Array, r.P.grad(r.a, r.xs, r.c, r.θ_ra))
    stats = born ? NQSCore.weighted_statistics(r.E, r.p) :
            NQSCore.statistics(reshape(r.E, size(states)))
    return stats, ∇
end

end # module NQSCoreReactantExt
