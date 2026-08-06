```@meta
CurrentModule = NeuralQuantumStates
```

# Basics

Installing `NeuralQuantumStates` gives you the whole stack; one `using` reaches all of it.

```@example 1
using NeuralQuantumStates
```

## A lattice

Lattices are described by a *specification* — a struct saying what you want — which
`build` turns into a lattice. Here is a periodic chain of 8 sites:

```@example 1
lat = build(Hypercube([8], 1.0; periodic=[true]))
```

Its bonds are pairs of site indices, and those indices are what operators act on:

```@example 1
bonds(lat)
```

## A model

Models work the same way: a specification, then `build`.

```@example 1
model = build(TransverseFieldIsing(lat; J=1.0, h_x=1.0, h_z=1.0))
```

A [`Model`](@ref) pairs the Hamiltonian with the space it acts on — the Hamiltonian alone knows
which sites it touches, but not that they are spin-1/2:

```@example 1
model.dof, model.nsites
```

The Hamiltonian is an ordinary OperatorAlgebra `OpSum`, so everything that package offers
applies to it — sparse and dense conversion, commutators, ITensor export:

```@example 1
typeof(model.hamiltonian)
```

## The computational basis

```@example 1
b = basis(model)
length(b.states)
```

States are packed integers rather than vectors, which is what makes symmetry-reduced bases
practical. `configurations` unpacks one into physical local values:

```@example 1
configurations(model.dof, b.states[6], model.nsites)
```

## Connected configurations

The kernel of variational Monte Carlo: for a configuration ``\vert s \rangle``, which
``\vert s' \rangle`` does the Hamiltonian connect it to, and with what matrix element
``\langle s' \vert \hat{H} \vert s \rangle``?

```@example 1
res = connected_padded(model.hamiltonian, b.states[1:3])
res.mels
```

Columns are padded to a common height with a **zero** matrix element, so the local energy

```julia
E_loc(s) = sum(res.mels[:, i] .* exp.(logψ.(res.configs[:, i]) .- logψ(s)))
```

needs no masking — a zero contributes nothing. `res.counts` records how many entries of each
column are real rather than padding.

## Symmetry sectors

A lattice knows its own symmetries, and can hand SymBasis the site permutations it needs:

```@example 1
dofo = dof_object(model.dof)
sector = basis(model, sym(Translational(0, lat), dofo))
length(sector.states)
```

Summing the dimensions of every momentum sector recovers the full space:

```@example 1
sum(length(basis(model, sym(Translational(k, lat), dofo)).states) for k in 0:7)
```

`connected_padded` takes such a basis directly, folding each connected configuration back onto
its orbit representative with the right norm and phase factors:

```@example 1
connected_padded(model.hamiltonian, sector.states[1:3], sector).mels
```

## Optimizing

Everything above is setup; this is the point of it. A [`VMC`](@ref) driver ties a variational
state, a Hamiltonian, a preconditioner, and an optimizer together.

```@example 1
using DifferentiationInterface, ForwardDiff, Random

nsites = 6
lat6 = build(Hypercube([nsites], 1.0; periodic=[true]))
model6 = build(TransverseFieldIsing(lat6; J=1.0, h_x=1.0))
b6 = basis(model6)

# The exact ansatz: one parameter per basis state. Useless for anything large, indispensable
# for checking that the machinery is right, since it can represent any state.
ansatz = LogStateVector(model6.dof, nsites, b6)
state = FullSumState(ansatz, init_parameters(ansatz, Xoshiro(0); scale=0.1), AutoForwardDiff())

driver = VMC(state, model6.hamiltonian;
    preconditioner=StochasticReconfiguration(; diag_shift=1e-2),
    optimizer=Descent(0.05))

log = run!(driver; iterations=2000, callbacks=(InvalidLossStopping(),))
final_energy(log)
```

Compare against exact diagonalization:

```@example 1
using LinearAlgebra
H6 = zeros(ComplexF64, length(b6.states), length(b6.states))
let index = Dict(s => i for (i, s) in pairs(b6.states)),
    res = connected_padded(model6.hamiltonian, b6.states)
    for n in eachindex(b6.states), j in 1:res.counts[n]
        H6[index[res.configs[j, n]], n] += res.mels[j, n]
    end
end
minimum(real(eigvals(Hermitian(H6))))
```

!!! warning "Stochastic reconfiguration wants a smaller step"
    The preconditioned update is routinely one to two orders of magnitude larger than the raw
    gradient. A learning rate that works for plain gradient descent will overshoot, leaving the
    run plateaued above the ground state with a large energy variance. Start near `0.05` with
    `diag_shift = 0.01`.

## A neural ansatz

For anything beyond a toy the ansatz has to compress. An RBM has a number of parameters linear
in the system size rather than exponential:

```@example 1
rbm = LuxAnsatz(RBM(nsites, 2), model6.dof, nsites; rng=Xoshiro(1))
NQSCore.n_parameters(rbm), length(b6.states)
```

## Sampling instead of summing

`FullSumState` enumerates the whole Hilbert space, so it is limited to small systems. Swapping
in an `MCState` changes nothing else about the code:

```@example 1
starts = random_configurations(model6.dof, nsites, 8, Xoshiro(2))
sampler = MetropolisSampler(LocalRule(), starts;
    n_chains=8, n_samples=5_000, burn_in=1_000)
mc = MCState(ansatz, parameters(state), sampler; backend=AutoForwardDiff(), rng=Xoshiro(3))
expect(mc, model6.hamiltonian)
```

The reported error bar is corrected for autocorrelation, and `r_hat` compares the chains
against one another — anything much above `1.01` means the error bar cannot be trusted.
