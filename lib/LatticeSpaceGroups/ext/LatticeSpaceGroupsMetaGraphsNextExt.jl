"""
    LatticeSpaceGroupsMetaGraphsNextExt

Turns a [`Lattice`](@ref) into a MetaGraphsNext.jl `MetaGraph` carrying its geometry.

Unlike a plain `SimpleGraph`, a `MetaGraph` has somewhere to put the things that make a lattice
a lattice: site labels as vertex labels, Cartesian positions as vertex data, and neighbour-shell
orders as edge data. That makes it the right target for plotting and for algorithms that need
to know *where* a site is, not merely what it connects to.

```julia
using LatticeSpaceGroups, MetaGraphsNext

mg = MetaGraph(build(Honeycomb([3, 3], 1.0; periodic=true)))
mg[(1, 1, 1)]                  # position of sublattice 1 in cell (1, 1)
mg[(1, 1, 1), (2, 1, 1)]       # 1, a nearest-neighbour bond
```
"""
module LatticeSpaceGroupsMetaGraphsNextExt

using LatticeSpaceGroups: Lattice, n_sites, site_labels, site_positions
using StaticArrays: SVector

# This extension triggers on Graphs *and* MetaGraphsNext: the underlying `SimpleGraph` a
# `MetaGraph` wraps comes from Graphs, and an extension may only use its declared triggers.
# MetaGraphsNext depends on Graphs anyway, so the pair is always loaded together in practice.
using Graphs: SimpleGraph
# `import`, not `using`: `MetaGraph` is being extended with a new method.
import MetaGraphsNext: MetaGraph

"""
    MetaGraph(lattice) -> MetaGraphsNext.MetaGraph

The lattice as a labelled graph.

Vertices are labelled by `(sublattice, n₁, …, n_D)` — the same labels [`site_labels`](@ref)
returns — and carry their Cartesian position as data. Edges carry their neighbour shell: `1`
for nearest neighbours, `2` for next-nearest, and so on.

Vertex *codes* follow insertion order, which is the lattice's own site numbering, so
`code_for(mg, label)` agrees with the site index that [`site_permutation`](@ref) permutes.
"""
function MetaGraph(lattice::Lattice{T,D,O}) where {T<:Real,D,O}
    labels = site_labels(lattice)
    positions = site_positions(lattice)

    mg = MetaGraph(
        SimpleGraph();
        label_type=eltype(labels),
        vertex_data_type=SVector{D,T},
        edge_data_type=Int
    )

    # Inserted in site order, so vertex codes match the lattice's site indices.
    for (label, position) in zip(labels, positions)
        mg[label] = position
    end
    for ((i, j), order) in zip(lattice.edges, lattice.edge_orders)
        mg[labels[i], labels[j]] = order
    end
    return mg
end

end # module LatticeSpaceGroupsMetaGraphsNextExt
