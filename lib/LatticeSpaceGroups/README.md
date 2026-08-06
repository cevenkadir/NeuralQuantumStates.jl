# LatticeSpaceGroups.jl

[![Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://cevenkadir.github.io/NeuralQuantumStates.jl/stable/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

Lattice geometry and the space-group machinery that turns it into **site permutations**.

Symmetry-reduced exact diagonalization needs a permutation vector: which site does site `i`
become under a translation, a reflection, a rotation? Writing those by hand only ever works for
a one-dimensional chain — `mod1.((1:N) .+ 1, N)` and nothing else. This package derives them
from the lattice geometry, so a kagome torus is no harder than a chain.

Its only dependency is StaticArrays.

## Installation

```julia
using Pkg
Pkg.add("LatticeSpaceGroups")
```

## Quick start

```julia
using LatticeSpaceGroups

lat = build(Honeycomb([3, 3], 1.0; periodic=true))

n_sites(lat)                     # 18
bonds(lat)                       # nearest-neighbour pairs, as site indices
length(point_group(lat))         # 12 — the D₆ point group
translation_permutation(lat, 1)  # the site permutation of a one-cell shift
```

Bonds come back as site-index pairs in the same numbering as `site_positions`, so they drop
straight into a Hamiltonian:

```julia
H = sum(J * Sz(i) * Sz(j) for (i, j) in bonds(lat))
```

## The lattice zoo

| Spec | Dim | Sites/cell | Coordination |
|---|---|---|---|
| `Hypercube`, and `Square` / `Cube` | any | 1 | `2D` |
| `Triclinic` | 3 | 1 | 6 |
| `Triangular` | 2 | 1 | 6 |
| `Honeycomb` | 2 | 2 | 3 |
| `Kagome` | 2 | 3 | 4 |
| `BCC` | 3 | 1 | 8 |
| `FCC` | 3 | 1 | 12 |
| `Diamond` | 3 | 2 | 4 |
| `Pyrochlore` | 3 | 4 | 6 |

A spec *describes* a lattice; `build` turns it into one. Specs validate on construction, so a
periodic dimension of extent 1 is rejected where you wrote it, not deep inside the build.

```julia
build(Hypercube([8]; periodic=true))                        # a ring
build(Hypercube([4, 4], 1.0; periodic=[true, false])) # a cylinder
build(Square(4; periodic=true); max_order=2)          # nearest and next-nearest bonds
```

## Symmetries

Point groups are found by enumerating the integer matrices that preserve the metric tensor of
the primitive basis, then keeping those that actually permute *this* lattice's sites — which
accounts for its supercell shape, boundary conditions, and multi-site unit cell.

That last part matters. A honeycomb lattice's six-fold axis passes through a hexagon centre,
not through the coordinate origin, so its rotations are symmetries only when paired with a
compensating fractional translation. `compensating_translation` finds it; `space_group` returns
operations complete with theirs.

```julia
lat = build(Kagome([3, 3], 1.0; periodic=true))

point_group(lat)             # orthogonal parts, 12 of them
space_group(lat)             # complete {R|τ}, point group × translations
rotation_permutations(lat)   # proper rotations, as permutations
reflection_permutation(lat)  # mirror — placed where the lattice admits one
```

## Optional integrations

None of these is a dependency; loading one activates an extension.

**[SymBasis.jl](https://github.com/cevenkadir/SymBasis.jl)** — symmetry constructors take a
lattice directly:

```julia
using LatticeSpaceGroups, SymBasis

lat = build(Hypercube([8]; periodic=true))
dofo = dof_object(Spin(1 // 2))
basis(dofo, 8, sym(Translational(0, lat), dofo))    # momentum-zero sector
```

**[Graphs.jl](https://github.com/JuliaGraphs/Graphs.jl)** — `SimpleGraph(lattice)`, plus
`Graphs.nv` and `Graphs.ne` methods.

**[MetaGraphsNext.jl](https://github.com/JuliaGraphs/MetaGraphsNext.jl)** —
`MetaGraph(lattice)`, carrying site labels, Cartesian positions, and neighbour-shell orders.

## Scope

Irreducible representations and character tables at given wave vectors are deliberately out of
scope.

## Part of a larger ecosystem

LatticeSpaceGroups is one of the packages behind
[NeuralQuantumStates.jl](https://github.com/cevenkadir/NeuralQuantumStates.jl), but it depends
on none of them. It is useful to anyone doing symmetry-reduced exact diagonalization, with or
without a neural network anywhere in sight.

## License

MIT
