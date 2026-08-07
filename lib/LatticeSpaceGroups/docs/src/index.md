```@meta
CurrentModule = LatticeSpaceGroups
```

# LatticeSpaceGroups.jl

*Lattice geometry, and the space groups it induces.*

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://cevenkadir.github.io/NeuralQuantumStates.jl/LatticeSpaceGroups/stable/) [![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://cevenkadir.github.io/NeuralQuantumStates.jl/LatticeSpaceGroups/dev/) [![Build Status](https://github.com/cevenkadir/NeuralQuantumStates.jl/actions/workflows/CI-LatticeSpaceGroups.yml/badge.svg?branch=main)](https://github.com/cevenkadir/NeuralQuantumStates.jl/actions/workflows/CI-LatticeSpaceGroups.yml?query=branch%3Amain) [![Coverage](https://codecov.io/gh/cevenkadir/NeuralQuantumStates.jl/branch/main/graph/badge.svg?flag=LatticeSpaceGroups)](https://codecov.io/gh/cevenkadir/NeuralQuantumStates.jl/flags?flag=LatticeSpaceGroups) [![Downloads](https://img.shields.io/badge/dynamic/json?url=http%3A%2F%2Fjuliapkgstats.com%2Fapi%2Fv1%2Fmonthly_downloads%2FLatticeSpaceGroups&query=total_requests&suffix=%2Fmonth&label=Downloads)](https://juliapkgstats.com/pkg/LatticeSpaceGroups) [![Total Downloads](https://img.shields.io/badge/dynamic/json?url=http%3A%2F%2Fjuliapkgstats.com%2Fapi%2Fv1%2Ftotal_downloads%2FLatticeSpaceGroups&query=total_requests&label=Total%20Downloads)](https://juliapkgstats.com/pkg/LatticeSpaceGroups)

Symmetry-reduced exact diagonalization needs a **site permutation**: which site does site `i`
become under a translation, a reflection, a rotation? Writing one by hand only ever works for a
one-dimensional chain — `mod1.((1:N) .+ 1, N)` and nothing else. This package derives them from
the lattice geometry, so a kagome torus is no harder than a chain.

Its only dependency is [StaticArrays.jl](https://github.com/JuliaArrays/StaticArrays.jl).

## Key features

- **A lattice zoo** — hypercubic (with `Square` and `Cube` shorthands), triclinic, triangular,
  honeycomb, kagome, BCC, FCC, diamond and pyrochlore, each described by a spec that validates
  its arguments where you write them rather than deep inside the build.
- **Site permutations from geometry** — [`translation_permutation`](@ref),
  [`reflection_permutation`](@ref) and [`rotation_permutations`](@ref) produce exactly the
  vectors SymBasis.jl's symmetry constructors take.
- **Real space groups** — [`point_group`](@ref) and [`space_group`](@ref) find the operations by
  solving the integer-matrix condition ``M^{\mathsf{T}} G M = G`` on the primitive basis, then
  keeping those that actually permute *this* lattice's sites under *its* boundary conditions.
- **Symmetry centres that are not the origin** — a honeycomb's six-fold axis passes through a
  hexagon centre, and an open chain's mirror sits at its midpoint. Both are found automatically
  via [`compensating_translation`](@ref), which is why those lattices report their full point
  group rather than a subgroup of it.
- **Periodic boundaries done properly** — neighbour shells are computed under the minimum-image
  convention, so "next-nearest neighbour" means the same thing on a torus as in the bulk.
- **Almost no dependencies** — a plain struct holds the sites and bonds. Graph interop is
  available, but only if you ask for it.

## The lattice zoo

| Spec | Dim | Sites per cell | Coordination |
|---|---|---|---|
| [`Hypercube`](@ref) — with [`Square`](@ref) and [`Cube`](@ref) | any | 1 | `2D` |
| [`Triclinic`](@ref) | 3 | 1 | 6 |
| [`Triangular`](@ref) | 2 | 1 | 6 |
| [`Honeycomb`](@ref) | 2 | 2 | 3 |
| [`Kagome`](@ref) | 2 | 3 | 4 |
| [`BCC`](@ref) | 3 | 1 | 8 |
| [`FCC`](@ref) | 3 | 1 | 12 |
| [`Diamond`](@ref) | 3 | 2 | 4 |
| [`Pyrochlore`](@ref) | 3 | 4 | 6 |

Anything not on that list you can build yourself — see
[Building a lattice by hand](@ref).

## Installation

**Requirements**: Julia 1.11 or later.

From the Julia REPL package mode (press `]`):

```julia
pkg> add LatticeSpaceGroups
```

or equivalently:

```julia
julia> import Pkg; Pkg.add("LatticeSpaceGroups")
```

## Quick example

The site permutation for translating a periodic chain by one cell is the cyclic shift — and the
same call gives you the corresponding permutation on a kagome torus, where writing it by hand
would not be realistic:

```@example index
using LatticeSpaceGroups

chain = build(Hypercube([6]; periodic=true))
translation_permutation(chain, 1)
```

```@example index
kagome = build(Kagome([2, 2], 1.0; periodic=true))
translation_permutation(kagome, 1)
```

Point groups come out at their crystallographic orders, boundary conditions included:

```@example index
length(point_group(build(Triangular([3, 3], 1.0; periodic=true))))   # D6
```

## Optional integrations

Each is a package extension, loaded only if you load the other package — none is a dependency.

- **[SymBasis.jl](https://github.com/cevenkadir/SymBasis.jl)** — `Translational`,
  `SpatialReflection` and `Rotational` gain methods taking a lattice directly, so a
  symmetry-reduced basis is two lines. See [Symmetry-reduced bases](@ref).
- **[Graphs.jl](https://github.com/JuliaGraphs/Graphs.jl)** — `SimpleGraph(lattice)`.
- **[MetaGraphsNext.jl](https://github.com/JuliaGraphs/MetaGraphsNext.jl)** —
  `MetaGraph(lattice)`, carrying site labels, positions and neighbour-shell orders.

## Scope

Irreducible representations and character tables at given wave vectors are deliberately out of
scope. This package produces the permutations; representation theory on top of them belongs
elsewhere.

## Part of a larger ecosystem

LatticeSpaceGroups is developed in the
[NeuralQuantumStates.jl](https://github.com/cevenkadir/NeuralQuantumStates.jl) monorepo and is
usable entirely on its own — it has no idea the neural-network stack exists, and its test suite
asserts as much. If you want the whole stack, see the
[NeuralQuantumStates.jl documentation](https://cevenkadir.github.io/NeuralQuantumStates.jl/).

## Important notice

This project is under active development. While it carries an extensive test suite, the API is
not yet stable. Please report anything you run into via the
[issue tracker](https://github.com/cevenkadir/NeuralQuantumStates.jl/issues/new).

## Module

What `?LatticeSpaceGroups` shows at the REPL:

```@docs
LatticeSpaceGroups
```
