```@meta
CurrentModule = NeuralQuantumStates
```

# Lattices and symmetries

*Provided by `LatticeSpaceGroups`, which has [its own
site](https://cevenkadir.github.io/NeuralQuantumStates.jl/LatticeSpaceGroups/dev/).*

Symmetry-reduced exact diagonalization needs a **site permutation**: which site does site `i`
become under a translation, a reflection, a rotation? Writing one by hand only ever works for a
chain — `mod1.((1:N) .+ 1, N)` and nothing else. `LatticeSpaceGroups` derives them from the
lattice geometry instead, so a kagome torus is no harder than a chain.

It is re-exported, so everything below is available from `using NeuralQuantumStates`.

```@example lsg
using NeuralQuantumStates
```

## Specifications, then `build`

A lattice is described by a *specification* — a struct saying what you want, validating its
arguments where you write them rather than deep inside the construction — which `build` turns
into a `Lattice`:

```@example lsg
lat = build(Kagome([2, 2], 1.0; periodic=true))
```

The zoo covers `Hypercube` (with `Square` and `Cube` shorthands), `Triclinic`, `Triangular`,
`Honeycomb`, `Kagome`, `BCC`, `FCC`, `Diamond` and `Pyrochlore`. Anything else you build
yourself from a `LatticeBasis` and, where the connectivity is not distance-derived, an explicit
bond list.

Bonds come out as site-index pairs, listed once and sorted — ready to be the site identifiers
of an operator term:

```@example lsg
bonds(build(Hypercube([4]; periodic=true)))
```

## Permutations for a symmetry-reduced basis

The point of the geometry is what it produces: the permutation vectors SymBasis's symmetry
constructors take.

```@example lsg
translation_permutation(lat, 1)
```

`SymBasis`'s `Translational`, `SpatialReflection` and `Rotational` also accept a lattice
directly — a package extension that loads when both packages are present, which in this stack
they always are:

```@example lsg
chain = build(Hypercube([8]; periodic=true))
dofo = dof_object(Spin(1 // 2))
sector = basis(dofo, 8, sym(Translational(0, chain), dofo))
length(sector.states)
```

## Real space groups, not idealised ones

`point_group` and `space_group` return the operations that actually permute *this* finite
lattice's sites under *its* boundary conditions, rather than the infinite-lattice answer. That
includes operations whose symmetry centre is not the coordinate origin — a honeycomb's six-fold
axis passes through a hexagon centre, an open chain's mirror sits at its midpoint — which is
why those lattices report their full point group rather than a subgroup of it.

```@example lsg
length(point_group(lat))     # D₆, found without being told the symmetry centre
```

## Full documentation

`LatticeSpaceGroups` is registered as a package in its own right and has a complete site:

- [Lattices](https://cevenkadir.github.io/NeuralQuantumStates.jl/LatticeSpaceGroups/dev/manual/lattices/)
  — the predefined zoo, bonds and neighbour shells, and building one by hand
- [Symmetries](https://cevenkadir.github.io/NeuralQuantumStates.jl/LatticeSpaceGroups/dev/manual/symmetries/)
  — translations, point groups, space groups, and the bridge to SymBasis
- [Symmetry-reduced bases](https://cevenkadir.github.io/NeuralQuantumStates.jl/LatticeSpaceGroups/dev/examples/symmetry_reduced_bases/)
- [Three-dimensional lattices](https://cevenkadir.github.io/NeuralQuantumStates.jl/LatticeSpaceGroups/dev/examples/three_dimensional_lattices/)
- [Boundary conditions and symmetry](https://cevenkadir.github.io/NeuralQuantumStates.jl/LatticeSpaceGroups/dev/examples/boundaries_and_symmetry/)
- [API reference](https://cevenkadir.github.io/NeuralQuantumStates.jl/LatticeSpaceGroups/dev/api/lattices/)
