<div align="center">

# LatticeSpaceGroups.jl

*Lattice geometry, and the space groups it induces, for Julia*

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://cevenkadir.github.io/NeuralQuantumStates.jl/LatticeSpaceGroups/stable/) [![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://cevenkadir.github.io/NeuralQuantumStates.jl/LatticeSpaceGroups/dev/) [![Build Status](https://github.com/cevenkadir/NeuralQuantumStates.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/cevenkadir/NeuralQuantumStates.jl/actions/workflows/CI.yml?query=branch%3Amain) [![Coverage](https://codecov.io/gh/cevenkadir/NeuralQuantumStates.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/cevenkadir/NeuralQuantumStates.jl) [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE) [![Downloads](https://img.shields.io/badge/dynamic/json?url=http%3A%2F%2Fjuliapkgstats.com%2Fapi%2Fv1%2Fmonthly_downloads%2FLatticeSpaceGroups&query=total_requests&suffix=%2Fmonth&label=Downloads)](https://juliapkgstats.com/pkg/LatticeSpaceGroups)
</div>

Symmetry-reduced exact diagonalization needs a **site permutation**: which site does site `i` become under a translation, a reflection, a rotation? Writing one by hand only ever works for a one-dimensional chain — `mod1.((1:N) .+ 1, N)` and nothing else. LatticeSpaceGroups.jl derives them from the lattice geometry, so a kagome torus is no harder than a chain. It finds the point group by solving the integer-matrix condition on the primitive basis and then keeping the operations that genuinely permute *that* lattice's sites under *its* boundary conditions — including the ones whose symmetry centre is not the coordinate origin, which is what a honeycomb, a diamond and an open chain all need. Its only dependency is [StaticArrays.jl](https://github.com/JuliaArrays/StaticArrays.jl).

## Key features
- **A lattice zoo**: hypercubic, triclinic, triangular, honeycomb, kagome, BCC, FCC, diamond and pyrochlore, each described by a spec that validates its arguments where you write them rather than deep inside the build.
- **Site permutations from geometry**: `translation_permutation`, `reflection_permutation` and `rotation_permutations` produce exactly the vectors SymBasis.jl's symmetry constructors take.
- **Real space groups**: `point_group` and `space_group` return the operations that actually act on your finite lattice, with its boundary conditions and its multi-site unit cell — not the idealised infinite-lattice answer.
- **Symmetry centres away from the origin**: a honeycomb's six-fold axis passes through a hexagon centre, and an open chain's mirror sits at its midpoint. Both are found automatically, which is why those lattices report their full point group instead of a subgroup of it.
- **Periodic boundaries done properly**: neighbour shells use the minimum-image convention, so "next-nearest neighbour" means the same thing on a torus as in the bulk.
- **Almost no dependencies**: a plain struct holds the sites and bonds. Graph interop is available, but only if you ask for it.

## Predefined lattices
| Spec | Dimension | Sites per cell | Coordination |
|---|---|---|---|
| `Hypercube` — with `Square` and `Cube` shorthands | any | 1 | `2D` |
| `Triclinic` | 3 | 1 | 6 |
| `Triangular` | 2 | 1 | 6 |
| `Honeycomb` | 2 | 2 | 3 |
| `Kagome` | 2 | 3 | 4 |
| `BCC` | 3 | 1 | 8 |
| `FCC` | 3 | 1 | 12 |
| `Diamond` | 3 | 2 | 4 |
| `Pyrochlore` | 3 | 4 | 6 |

Anything not on that list you can build yourself from a `LatticeBasis` and, if the connectivity is not distance-derived, an explicit bond list.

## Installation
**Requirements**: Julia 1.11 or later.

To install the latest stable version, use the Julia package manager. Either use the Julia REPL package mode (by pressing `]`):
```julia
pkg> add LatticeSpaceGroups
```
or open the Julia REPL and run the following command:
```julia
julia> import Pkg; Pkg.add("LatticeSpaceGroups")
```

## Documentation
For detailed information on using this package, check out the [stable documentation](https://cevenkadir.github.io/NeuralQuantumStates.jl/LatticeSpaceGroups/stable/).

## Manual
- [Lattices](https://cevenkadir.github.io/NeuralQuantumStates.jl/LatticeSpaceGroups/stable/manual/lattices/)
- [Symmetries](https://cevenkadir.github.io/NeuralQuantumStates.jl/LatticeSpaceGroups/stable/manual/symmetries/)

## Examples
- [Symmetry-reduced bases](https://cevenkadir.github.io/NeuralQuantumStates.jl/LatticeSpaceGroups/stable/examples/symmetry_reduced_bases/)
- [Three-dimensional lattices](https://cevenkadir.github.io/NeuralQuantumStates.jl/LatticeSpaceGroups/stable/examples/three_dimensional_lattices/)
- [Boundary conditions and symmetry](https://cevenkadir.github.io/NeuralQuantumStates.jl/LatticeSpaceGroups/stable/examples/boundaries_and_symmetry/)

## Quick example
Build a kagome torus and ask for the permutation that translates it by one unit cell — the vector you would otherwise have to write out by hand:
```julia
julia> using LatticeSpaceGroups

julia> lat = build(Kagome([2, 2], 1.0; periodic=true))
Lattice{Float64,2,3}([2, 2]; periodic=Bool[1, 1], 12 sites, 24 bonds)

julia> translation_permutation(lat, 1)
12-element Vector{Int64}:
  4
  5
  6
  1
  2
  3
 10
 11
 12
  7
  8
  9

# the D₆ point group, found without being told the symmetry centre
julia> length(point_group(lat))
12
```

Bonds come out as site-index pairs, listed once and sorted, ready to be the site identifiers of an operator term:
```julia
julia> bonds(build(Hypercube([4]; periodic=true)))
4-element Vector{Tuple{Int64, Int64}}:
 (1, 2)
 (1, 4)
 (2, 3)
 (3, 4)
```

## Optional integrations
Each is a package extension, loaded only if you load the other package. None is a dependency.
- **[SymBasis.jl](https://github.com/cevenkadir/SymBasis.jl)**: `Translational`, `SpatialReflection` and `Rotational` gain methods taking a lattice directly, so a symmetry-reduced basis is two lines.
- **[Graphs.jl](https://github.com/JuliaGraphs/Graphs.jl)**: `SimpleGraph(lattice)`.
- **[MetaGraphsNext.jl](https://github.com/JuliaGraphs/MetaGraphsNext.jl)**: `MetaGraph(lattice)`, carrying site labels, positions and neighbour-shell orders.

## Scope
Irreducible representations and character tables at given wave vectors are deliberately out of scope. This package produces the permutations; representation theory on top of them belongs elsewhere.

## Part of a larger ecosystem
LatticeSpaceGroups.jl is developed in the [NeuralQuantumStates.jl](https://github.com/cevenkadir/NeuralQuantumStates.jl) monorepo, and is usable entirely on its own — it has no idea the neural-network stack exists, and its test suite asserts as much.

## Important notice
This project is still under active development. While it includes an extensive test suite and is developed with high scientific rigor, you should always benchmark your own code. Please report any issues you encounter via the [GitHub issue tracker](https://github.com/cevenkadir/NeuralQuantumStates.jl/issues/new).

## Citation
If you use this package in your work, we would appreciate the following reference as in [CITATION.bib](https://github.com/cevenkadir/NeuralQuantumStates.jl/blob/main/CITATION.bib).
