```@meta
CurrentModule = NQSCore
```

# Extending the interfaces

Everything downstream of this package is written against three small contracts. Each is
deliberately narrow — that narrowness is what lets `NQSAnsatze`, `NQSSamplers` and
`NQSOptimisers` be developed without depending on one another.

## An ansatz

Subtype [`AbstractAnsatz`](@ref) and define one method:

```julia
struct MyAnsatz{D} <: NQSCore.AbstractAnsatz
    dof::D
    nsites::Int
    # ...whatever the functional form needs
end

function NQSCore.log_amplitude(a::MyAnsatz, θ, x::AbstractMatrix)
    # x is (nsites, batch) physical local values; return a length-batch vector
end
```

`x` holds **physical local values** — magnetic quantum numbers for a spin, occupation numbers
for a boson — as produced by `ConnectedBasisConfigurations.configurations`. Return the
*logarithm* of the amplitude: amplitudes underflow catastrophically for any interesting system
size, whereas the differences of logarithms that local energies actually need stay bounded.

Parameters are not stored in the ansatz. An ansatz describes the functional form; the parameters
live in the variational state that owns it. That separation is what makes it possible to
evaluate the same ansatz at perturbed parameters, which is exactly what computing a Jacobian
does.

Two accessors describe the space the ansatz lives on, and both default to the correspondingly
named field:

```julia
NQSCore.dof(a::MyAnsatz)      # defaults to a.dof
NQSCore.n_sites(a::MyAnsatz)  # defaults to a.nsites
```

Define them if your fields are named differently or the values are computed. Neither is
exported: `n_sites` is a lattice's word as much as an ansatz's — `LatticeSpaceGroups` exports its
own — and a package whose job is to be depended on should not make that choice for everyone
downstream.

Optionally define [`default_basis`](@ref), if there is an unambiguous answer to what a
`FullSumState` on this ansatz should sum over, and [`init_parameters`](@ref) and
[`n_parameters`](@ref).

## A sampler

Subtype [`AbstractSampler`](@ref) and define [`sample`](@ref):

```julia
function NQSCore.sample(s::MySampler, a::AbstractAnsatz, θ, rng::AbstractRNG, state=nothing)
    # ...
    return drawn, new_state
end
```

Return the drawn configurations as **packed** states, together with whatever the sampler needs
to resume. `state` is what a previous call returned, or `nothing` to start cold; a sampler with
nothing to carry returns `nothing` and ignores the argument. Returning a `(steps, chains)` matrix
is supported and preferred where it makes sense — the chain structure survives into
[`statistics`](@ref), which uses it for split-R̂ and for a between-chain error bar.

## A preconditioner

Subtype [`AbstractPreconditioner`](@ref) and define [`precondition`](@ref):

```julia
function NQSCore.precondition(p::MyPreconditioner, vs::AbstractVariationalState, operator)
    est = local_estimators(vs, operator)   # (; E, O, weights)
    # ...
    return stats, update
end
```

`update` must have the same structure as the state's parameters — [`match_parameter_shape`](@ref)
is there for exactly that — so a step is `θ .- η .* update`. Both returned values must come from
a *single* set of samples: [`local_estimators`](@ref) exists to guarantee that, since recomputing
would silently mix two sample sets and show up as an optimizer that mysteriously fails to
descend.
