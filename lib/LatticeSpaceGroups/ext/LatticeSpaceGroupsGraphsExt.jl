"""
    LatticeSpaceGroupsGraphsExt

Turns a [`Lattice`](@ref) into a Graphs.jl `SimpleGraph`.

A lattice already knows its sites and bonds, so this is a conversion rather than a
capability — but it is the conversion that gives access to everything Graphs.jl can do:
shortest paths, connectivity, colourings, plotting. It lives in an extension because a package
about lattice symmetry has no business depending on a graph library to describe a bond list.

```julia
using LatticeSpaceGroups, Graphs

g = SimpleGraph(build(Square(4; periodic=true)))
Graphs.nv(g), Graphs.ne(g)     # (16, 32)
```
"""
module LatticeSpaceGroupsGraphsExt

using LatticeSpaceGroups: Lattice, bonds, n_sites

using Graphs: add_edge!
# `import`, not `using`: these are being extended with new methods, and `using` would define
# shadowing functions in this module instead.
import Graphs: SimpleGraph, nv, ne

"""
    SimpleGraph(lattice) -> Graphs.SimpleGraph

The lattice's bond structure as an undirected graph, with vertex `i` the lattice's site `i`.

Site positions and neighbour-shell orders are not carried over — a `SimpleGraph` has nowhere to
put them. Use `MetaGraph(lattice)` from the MetaGraphsNext extension if you need them.
"""
function SimpleGraph(lattice::Lattice)
    g = SimpleGraph(n_sites(lattice))
    for (i, j) in bonds(lattice)
        add_edge!(g, i, j)
    end
    return g
end

"""
    Graphs.nv(lattice) -> Int

Number of sites, for code written against the Graphs.jl interface. This package's own name for
it is [`n_sites`](@ref), which does not collide with `Graphs.nv` on `using`.
"""
nv(lattice::Lattice) = n_sites(lattice)

"""
    Graphs.ne(lattice) -> Int

Number of bonds, for code written against the Graphs.jl interface.
"""
ne(lattice::Lattice) = length(lattice.edges)

end # module LatticeSpaceGroupsGraphsExt
