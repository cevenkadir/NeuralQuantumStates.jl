"""
    VMCLog

The record of an optimization run: one `Stats` per iteration, plus the wall-clock time.

Kept in memory rather than streamed to a file. Serialization is deliberately left to the caller
— writing JSON or HDF5 would mean a dependency the rest of the stack does not need, and the
fields here are plain arrays that any format handles.

# Fields
- `energies::Vector{Stats}`: the energy at every iteration.
- `elapsed::Float64`: seconds taken.
"""
struct VMCLog
    energies::Vector{Stats}
    elapsed::Float64
end

Base.length(log::VMCLog) = length(log.energies)
Base.lastindex(log::VMCLog) = length(log.energies)
Base.getindex(log::VMCLog, i) = log.energies[i]

"""Final energy of a run, or `nothing` if it never got started."""
final_energy(log::VMCLog) = isempty(log.energies) ? nothing : last(log.energies)

function Base.show(io::IO, log::VMCLog)
    print(io, "VMCLog(", length(log), " iterations, ", round(log.elapsed, digits=2), "s")
    isempty(log.energies) || print(io, ", final E = ", final_energy(log))
    print(io, ")")
    return nothing
end

"""
    VMC(state, hamiltonian; preconditioner=Identity(), optimizer=Descent(0.05))

A variational Monte Carlo ground-state optimization.

The counterpart of NetKet's `driver.VMC`. Each iteration computes the energy and its gradient,
passes the gradient through `preconditioner`, and hands the result to `optimizer`.

# Choosing a preconditioner

`Identity()` is plain gradient descent. It is adequate in weakly correlated regimes and
**stalls** in strongly ordered ones, because the gradient carries a factor of the Born
probability that vanishes for exactly the configurations whose amplitude must grow.
`StochasticReconfiguration()` divides that factor out and is the right default for anything
hard; see `NQSOptimisers`.

# Fields
- `state`: an `AbstractVariationalState` — `FullSumState` for exact summation, `MCState` for
  sampling.
- `hamiltonian`: the operator to minimize.
- `preconditioner`: an `AbstractPreconditioner`.
- `optimizer`: any Optimisers.jl rule (`Descent`, `Adam`, ...).
"""
struct VMC{S,H,P,O}
    state::S
    hamiltonian::H
    preconditioner::P
    optimizer::O
end

function VMC(
    state, hamiltonian;
    preconditioner=NQSOptimisers.Identity(), optimizer=Optimisers.Descent(0.05)
)
    return VMC(state, hamiltonian, preconditioner, optimizer)
end

"""
    run!(driver; iterations, callbacks=[], verbose=false) -> VMCLog

Run the optimization.

Each callback is called as `callback(iteration, stats, state)` after the parameters are
updated; returning `false` stops the run. That is how [`EarlyStopping`](@ref) and
[`InvalidLossStopping`](@ref) work, and it is enough for a user callback to do the same.
"""
function run!(
    driver::VMC; iterations::Integer=100, callbacks=(), verbose::Bool=false
)
    vs = driver.state
    opt_state = Optimisers.setup(driver.optimizer, parameters(vs))

    energies = Stats[]
    started = time()

    for it in 1:iterations
        stats, δ = NQSOptimisers.precondition(driver.preconditioner, vs, driver.hamiltonian)
        push!(energies, stats)

        opt_state, θ = Optimisers.update(opt_state, parameters(vs), δ)
        setparameters!(vs, θ)
        # A Monte Carlo state must redraw: its cached samples came from the old parameters.
        vs isa MCState && resample!(vs)

        verbose && println("iter ", lpad(it, 5), "   E = ", stats)

        stop = false
        for cb in callbacks
            cb(it, stats, vs) === false && (stop = true)
        end
        stop && break
    end

    return VMCLog(energies, time() - started)
end

# ------------------------------------------------------------------------------- callbacks

"""
    EarlyStopping(; patience=50, min_delta=1e-6)

Stop when the energy has not improved by more than `min_delta` for `patience` iterations.

Judged on the energy *mean*. With a Monte Carlo state the mean fluctuates, so `min_delta`
should be set comfortably above the error bar — otherwise noise alone will look like an
improvement and the run will never stop.
"""
mutable struct EarlyStopping
    patience::Int
    min_delta::Float64
    best::Float64
    waited::Int
end
EarlyStopping(; patience::Integer=50, min_delta::Real=1e-6) =
    EarlyStopping(Int(patience), Float64(min_delta), Inf, 0)

function (cb::EarlyStopping)(iteration, stats, state)
    current = real(stats.mean)
    if current < cb.best - cb.min_delta
        cb.best = current
        cb.waited = 0
        return true
    end
    cb.waited += 1
    return cb.waited < cb.patience
end

"""
    InvalidLossStopping()

Stop as soon as the energy becomes `NaN` or `Inf`.

Worth having on by default. A diverged run otherwise keeps going for its full iteration count,
propagating `NaN` through every parameter, and the failure is then reported as a finished run
with a meaningless answer rather than as the divergence it was.
"""
struct InvalidLossStopping end

function (::InvalidLossStopping)(iteration, stats, state)
    return isfinite(real(stats.mean))
end
