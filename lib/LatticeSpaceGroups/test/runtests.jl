using LatticeSpaceGroups
using LatticeSpaceGroups: _site_key, _integer_point_candidates, compensating_translation

using Aqua
using LinearAlgebra: det, norm
using Test

"""Names of packages currently loaded in this session."""
loaded_package_names() = Set(m.name for m in keys(Base.loaded_modules))

"""Compose two site permutations: first `a`, then `b`."""
compose(b, a) = [b[a[i]] for i in eachindex(a)]

"""Average number of bonds per site."""
coordination(lat) = 2 * length(bonds(lat)) / n_sites(lat)

"""Multiplicative order of a site permutation."""
function permutation_order(lat, p)
    identity_perm = collect(1:n_sites(lat))
    order, q = 1, p
    while q != identity_perm
        q = compose(p, q)
        order += 1
    end
    return order
end

@testset "LatticeSpaceGroups.jl" begin
    @testset "dependency weight" begin
        # Two properties at once. This package is Tier 1 -- usable for plain exact
        # diagonalization without dragging in the neural-network stack, which is the single
        # property the whole package split exists to buy. And its graph/spatial-index
        # dependencies were deliberately removed in favour of a plain struct, so their absence
        # is asserted rather than left to drift back in via someone's `Pkg.add`.
        #
        # This testset must run FIRST: the extension testsets below load Graphs, MetaGraphsNext
        # and SymBasis, after which these assertions would be meaningless.
        loaded = loaded_package_names()
        @testset "no neural-network stack" begin
            for heavy in ("Lux", "Enzyme", "Reactant", "Zygote", "CUDA", "Metal", "ComponentArrays")
                @test heavy ∉ loaded
            end
        end
        @testset "no graph or spatial-index libraries" begin
            for removed in ("Graphs", "MetaGraphsNext", "NearestNeighbors", "Distances", "SymBasis")
                @test removed ∉ loaded
            end
        end
    end

    @testset "lattice specs" begin
        # A spec describes a lattice and validates it on construction, so an invalid lattice
        # is rejected where it is described rather than deep inside `build`.
        @testset "validation happens at spec construction" begin
            @test_throws ArgumentError Hypercube([0, 4], 1.0)          # non-positive extent
            @test_throws ArgumentError Hypercube([4, 4], -1.0)         # non-positive spacing
            @test_throws ArgumentError Hypercube([1, 4], 1.0; periodic=true)
            @test_throws ArgumentError Hypercube([4, 4], 1.0; periodic=[true])  # wrong length
            @test_throws ArgumentError Triangular([2, 2, 2], 1.0)      # wrong dimension
            @test_throws ArgumentError Triclinic([2, 2], [1.0, 1.0], [90.0, 90.0])
            @test_throws ArgumentError Pyrochlore([2, 2], 1.0)         # wrong dimension
            @test_throws ArgumentError FCC([2, 2, 2], 0.0)             # non-positive spacing
        end

        @testset "scalar periodic applies to every dimension" begin
            @test Hypercube([3, 3], 1.0; periodic=true).periodic == [true, true]
            @test Hypercube([3, 3], 1.0).periodic == [false, false]
        end

        @testset "dimension is carried in the type" begin
            @test Hypercube([4], 1.0) isa AbstractLatticeSpec{1}
            @test Hypercube([4, 4], 1.0) isa AbstractLatticeSpec{2}
            @test Kagome([2, 2], 1.0) isa AbstractLatticeSpec{2}
            @test Diamond([2, 2, 2], 1.0) isa AbstractLatticeSpec{3}
            @test Triclinic([2, 2, 2], [1.0, 1.0, 1.0], [90.0, 90.0, 90.0]) isa
                  AbstractLatticeSpec{3}
        end

        @testset "edge length defaults to 1" begin
            @test Triangular([2, 2]).edge_length == 1.0
            @test FCC([2, 2, 2]).edge_length == 1.0
        end

        @testset "specs are printable" begin
            @test occursin("Hypercube", sprint(show, Hypercube([4, 4], 1.0; periodic=true)))
            @test occursin("Pyrochlore", sprint(show, Pyrochlore([2, 2, 2], 1.0)))
        end

        @testset "Square and Cube are Hypercube shorthands" begin
            @test Square(4).shape == [4, 4]
            @test Square([3, 5]).shape == [3, 5]
            @test Cube(3).shape == [3, 3, 3]
            @test Cube([2, 3, 4]).shape == [2, 3, 4]
            # Identical lattices, not merely identical specs.
            a, b = build(Square(4; periodic=true)), build(Hypercube([4, 4], 1.0; periodic=true))
            @test site_positions(a) == site_positions(b)
            @test bonds(a) == bonds(b)
            @test_throws ArgumentError Square([2, 2, 2])
            @test_throws ArgumentError Cube([2, 2])
            # `Chain` is deliberately not exported: the name collides with `Lux.Chain`, and
            # this package is meant to be loaded alongside Lux.
            @test !isdefined(LatticeSpaceGroups, :Chain)
        end
    end

    @testset "lattice bases" begin
        @testset "a singular set of primitive vectors is rejected" begin
            # It spans fewer than D dimensions, so it describes no D-dimensional lattice, and
            # `_site_key` could not invert it to get fractional coordinates.
            @test_throws ArgumentError LatticeBasis([1.0 2.0; 2.0 4.0])
            @test_throws ArgumentError LatticeBasis([1.0 0.0; 0.0 0.0])
        end

        @testset "equivalent spellings agree" begin
            m = LatticeBasis([1.0 0.5; 0.0 sqrt(0.75)])
            v = LatticeBasis([[1.0, 0.0], [0.5, sqrt(0.75)]])
            @test m.vectors == v.vectors
            @test LatticeBasis(2.0).vectors == LatticeBasis([2.0;;]).vectors
            @test size(LatticeBasis(1.0, [0.0, 0.5]).site_offsets) == (1, 2)
        end
    end

    @testset "geometry" begin
        @test n_sites(build(Hypercube([8]; periodic=true))) == 8
        @test n_sites(build(Square(3))) == 9
        @test n_sites(build(Cube(3))) == 27
        @test n_sites(build(Triangular([3, 3], 1.0))) == 9
        # Multi-site unit cells: honeycomb has 2 sites per cell, kagome 3, pyrochlore 4.
        @test n_sites(build(Honeycomb([3, 3], 1.0))) == 18
        @test n_sites(build(Kagome([3, 3], 1.0))) == 27
        @test n_sites(build(Diamond([2, 2, 2], 1.0))) == 16
        @test n_sites(build(Pyrochlore([2, 2, 2], 1.0))) == 32
        @test n_sites(build(Triclinic([2, 2, 2], [1.0, 1.5, 2.0], [80.0, 70.0, 60.0]))) == 8

        @testset "site_positions is indexed by site number" begin
            lat = build(Hypercube([4]; periodic=true))
            @test site_positions(lat) == [[0.0], [1.0], [2.0], [3.0]]
            @test length(site_labels(lat)) == n_sites(lat)
            # Labels are (sublattice, cell...) and match the position ordering.
            @test site_labels(lat) == [(1, 1), (1, 2), (1, 3), (1, 4)]
            @test site_labels(build(Honeycomb([2, 1], 1.0)))[1:2] == [(1, 1, 1), (2, 1, 1)]
        end

        @testset "bonds are listed once, sorted, in the site indexing" begin
            lat = build(Hypercube([4]; periodic=true))
            @test bonds(lat) == [(1, 2), (1, 4), (2, 3), (3, 4)]
            @test all(i < j for (i, j) in bonds(lat))
            @test issorted(bonds(lat))
            @test length(bonds(build(Hypercube([8]; periodic=false)))) == 7
        end

        @testset "a dimension of extent 1 cannot be periodic" begin
            @test_throws ArgumentError build(Hypercube([1, 4], 1.0; periodic=true))
            # ...but it is fine when that dimension is left open.
            @test n_sites(build(Hypercube([1, 4], 1.0; periodic=[false, true]))) == 4
        end

        @testset "multi-site bases scale with edge_length" begin
            # The honeycomb site offsets used to be written out as literal numbers and so did
            # not scale with `edge_length`: every honeycomb but the unit one had the wrong
            # bond length. They are now derived from the primitive vectors.
            for a in (1.0, 2.0, 3.0)
                pos = site_positions(build(Honeycomb([2, 2], a)))
                @test norm(pos[2] - pos[1]) ≈ a / sqrt(3)
            end
            for a in (1.0, 2.5)
                pos = site_positions(build(Kagome([2, 2], a)))
                @test norm(pos[2] - pos[1]) ≈ a / 2
            end
            for a in (1.0, 2.5)
                pos = site_positions(build(Diamond([2, 2, 2], a)))
                @test norm(pos[2] - pos[1]) ≈ sqrt(3) * a / 4
                pos = site_positions(build(Pyrochlore([2, 2, 2], a)))
                @test norm(pos[2] - pos[1]) ≈ a / (2 * sqrt(2))
            end
        end

        @testset "honeycomb sits at the textbook sublattice positions" begin
            # Sublattices at fractional (1/3,1/3) and (2/3,2/3) of the triangular cell.
            lat = build(Honeycomb([1, 1], 1.0))
            third = _site_key(lat, site_positions(lat)[1])
            twothirds = _site_key(lat, site_positions(lat)[2])
            @test all(third .== round(Int, 10^12 / 3))
            @test all(twothirds .== round(Int, 2 * 10^12 / 3))
        end
    end

    @testset "coordination numbers are the textbook ones" begin
        # The single most informative check on the geometry: a wrong primitive vector, a wrong
        # site offset, or a broken minimum-image metric all show up here as the wrong number of
        # neighbours per site. This is the class of test that would have caught the honeycomb
        # offset bug at the point it was introduced.
        cases = [
            ("chain", build(Hypercube([8]; periodic=true)), 2),
            ("square", build(Square(4; periodic=true)), 4),
            ("cube", build(Cube(3; periodic=true)), 6),
            ("triangular", build(Triangular([4, 4], 1.0; periodic=true)), 6),
            ("honeycomb", build(Honeycomb([3, 3], 1.0; periodic=true)), 3),
            ("kagome", build(Kagome([3, 3], 1.0; periodic=true)), 4),
            ("bcc", build(BCC([3, 3, 3], 1.0; periodic=true)), 8),
            ("fcc", build(FCC([3, 3, 3], 1.0; periodic=true)), 12),
            ("diamond", build(Diamond([3, 3, 3], 1.0; periodic=true)), 4),
            ("pyrochlore", build(Pyrochlore([2, 2, 2], 1.0; periodic=true)), 6),
        ]
        for (name, lat, expected) in cases
            @test coordination(lat) == expected
        end
    end

    @testset "neighbour shells" begin
        @testset "max_order adds further shells" begin
            # On a triangular lattice the second shell sits at sqrt(3) times the first, and
            # there are as many next-nearest neighbours as nearest ones.
            near = build(Triangular([4, 4], 1.0; periodic=true))
            far = build(Triangular([4, 4], 1.0; periodic=true); max_order=2)
            @test bonds(far; order=1) == bonds(near)
            @test length(bonds(far; order=2)) == length(bonds(far; order=1))
            @test length(bonds(far)) == 2 * length(bonds(near))
        end

        @testset "the second shell of a square lattice is the cell diagonal" begin
            # A grid lattice takes the O(N) grid path at max_order == 1 and the general
            # distance search beyond it. The two must agree on the first shell, or the
            # keyword would be quietly changing bonds it has no business touching.
            lat = build(Square(4, 1.0; periodic=true); max_order=2)
            @test bonds(lat; order=1) == bonds(build(Square(4, 1.0; periodic=true)))
            @test length(bonds(lat; order=1)) == 32       # 4 neighbours x 16 sites / 2
            @test length(bonds(lat; order=2)) == 32       # 4 diagonals x 16 sites / 2

            positions = site_positions(lat)
            # A first-shell bond moves along exactly one axis; a second-shell bond moves along
            # both. Displacements of 3 are the same step taken the long way round the torus.
            axial(Δ) = count(x -> abs(x) ≈ 1.0 || abs(x) ≈ 3.0, Δ)
            for (i, j) in bonds(lat; order=1)
                @test axial(positions[j] - positions[i]) == 1
            end
            for (i, j) in bonds(lat; order=2)
                @test axial(positions[j] - positions[i]) == 2
            end
        end

        @testset "asking for more shells than exist is an error" begin
            @test_throws ArgumentError build(Triangular([2, 2], 1.0; periodic=true); max_order=10)
        end
    end

    @testset "explicit edges" begin
        # Connectivity that is not distance-derived: given as site-index pairs.
        basis = LatticeBasis([1.0;;])
        lat = Lattice([4], basis, [(1, 3), (2, 4)])
        @test bonds(lat) == [(1, 3), (2, 4)]
        @test n_sites(lat) == 4

        @test_throws ArgumentError Lattice([4], basis, [(1, 5)])      # out of range
        @test_throws ArgumentError Lattice([4], basis, [(2, 2)])      # self-loop
        @test_throws ArgumentError Lattice([4], basis, [(1, 2)]; orders=[1, 2])
    end

    @testset "translations" begin
        @testset "a periodic chain translates cyclically" begin
            # This is the exact permutation that was previously hand-written as
            # `mod1.((1:N) .+ 1, N)`; reproducing it is the whole point of the package.
            for n in (4, 6, 8)
                lat = build(Hypercube([n]; periodic=true))
                @test translation_permutation(lat, 1) == mod1.((1:n) .+ 1, n)
            end
        end

        @testset "translating by the full extent is the identity" begin
            lat = build(Hypercube([6]; periodic=true))
            @test translation_permutation(lat, 1; cells=6) == collect(1:6)
        end

        @testset "translations compose additively" begin
            lat = build(Hypercube([6]; periodic=true))
            t1 = translation_permutation(lat, 1; cells=1)
            t2 = translation_permutation(lat, 1; cells=2)
            @test compose(t1, t1) == t2
        end

        @testset "open directions admit no translation" begin
            lat = build(Hypercube([6]; periodic=false))
            @test isempty(translation_generators(lat))
            @test_throws ArgumentError translation_permutation(lat, 1)
        end

        @testset "one generator per periodic direction" begin
            lat = build(Hypercube([3, 4], 1.0; periodic=[true, false]))
            @test length(translation_generators(lat)) == 1
            # The group is the product of the cyclic groups along periodic directions.
            @test length(translation_group(lat)) == 3
            @test length(translation_group(build(Hypercube([3, 4], 1.0; periodic=true)))) == 12
        end
    end

    @testset "point groups match the known crystallographic orders" begin
        # Bravais lattices, then lattices with a multi-site basis whose symmetry centre is
        # not the coordinate origin, then an open lattice whose mirror is at its midpoint.
        cases = [
            ("chain PBC", build(Hypercube([8]; periodic=true)), 2),           # C_i
            ("chain open", build(Hypercube([8]; periodic=false)), 2),         # C_i about midpoint
            ("square PBC", build(Square(4; periodic=true)), 8),         # D4
            ("square open", build(Square(4; periodic=false)), 8),       # D4
            ("triangular", build(Triangular([3, 3], 1.0; periodic=true)), 12),  # D6
            ("honeycomb", build(Honeycomb([3, 3], 1.0; periodic=true)), 12),    # D6
            ("kagome", build(Kagome([3, 3], 1.0; periodic=true)), 12),          # D6
            ("cube", build(Cube(3; periodic=true)), 48),                        # O_h
            # The 3-D lattices with non-cubic primitive cells and multi-site bases: these
            # exercise the compensating-translation machinery hardest, since none of their
            # symmetry centres is the coordinate origin.
            ("bcc", build(BCC([3, 3, 3], 1.0; periodic=true)), 48),             # O_h
            ("fcc", build(FCC([3, 3, 3], 1.0; periodic=true)), 48),             # O_h
            ("diamond", build(Diamond([2, 2, 2], 1.0; periodic=true)), 48),     # O_h
            ("pyrochlore", build(Pyrochlore([2, 2, 2], 1.0; periodic=true)), 48),  # O_h
        ]
        for (name, lat, order) in cases
            @test length(point_group(lat)) == order
        end

        @testset "half the elements are proper rotations" begin
            for (name, lat, order) in cases
                pg = point_group(lat)
                @test count(o -> det(o.matrix) ≈ 1, pg) == order ÷ 2
            end
        end
    end

    @testset "the space group's permutations form a group" begin
        # Closure is the property that caught a real bug: an earlier float-keyed
        # implementation admitted C2 and C3 for the honeycomb but rejected C6 = C2 * C3,
        # producing a "group" of order 8 that could not divide 12.
        #
        # The closure statement has to be made about the *space* group, not the point group.
        # Composing {R₁|τ₁} and {R₂|τ₂} gives {R₁R₂ | R₁τ₂ + τ₁}, whose translation part
        # differs in general from the canonical τ picked for R₁R₂ by a lattice vector -- so
        # the point-group permutations alone are genuinely not closed, and only become so
        # once every translation is included.
        for lat in (
            build(Square(3; periodic=true)),
            build(Triangular([3, 3], 1.0; periodic=true)),
            build(Honeycomb([2, 2], 1.0; periodic=true)),
            build(Diamond([2, 2, 2], 1.0; periodic=true)),
        )
            perms = Set(site_permutation(lat, op) for op in space_group(lat))
            identity_perm = collect(1:n_sites(lat))

            @test identity_perm in perms
            @test all(compose(p, q) in perms for p in perms, q in perms)
            @test all(invperm(p) in perms for p in perms)
        end
    end

    @testset "symmetries preserve the bond set" begin
        for lat in (
            build(Hypercube([6]; periodic=true)),
            build(Square(4; periodic=true)),
            build(Triangular([3, 3], 1.0; periodic=true)),
            build(Honeycomb([3, 3], 1.0; periodic=true)),
            build(Pyrochlore([2, 2, 2], 1.0; periodic=true)),
        )
            for t in translation_group(lat)
                @test preserves_edges(lat, site_permutation(lat, t))
            end
            for op in space_group(lat)
                @test preserves_edges(lat, site_permutation(lat, op))
            end
        end
    end

    @testset "space group order is |point| x |translations|" begin
        for lat in (
            build(Square(4; periodic=true)),
            build(Triangular([3, 3], 1.0; periodic=true)),
        )
            @test length(space_group(lat)) ==
                  length(point_group(lat)) * length(translation_group(lat))
        end
    end

    @testset "reflections" begin
        @testset "an open chain reflects about its midpoint" begin
            lat = build(Hypercube([6]; periodic=false))
            @test reflection_permutation(lat, 1) == collect(6:-1:1)
        end

        @testset "a reflection is its own inverse" begin
            for lat in (
                build(Hypercube([6]; periodic=false)),
                build(Hypercube([6]; periodic=true)),
                build(Square(4; periodic=true)),
                build(Honeycomb([3, 3], 1.0; periodic=true)),
            )
                p = reflection_permutation(lat, 1)
                @test compose(p, p) == collect(1:n_sites(lat))
            end
        end
    end

    @testset "rotations" begin
        @testset "a square lattice has a four-fold axis" begin
            lat = build(Square(4; periodic=true))
            perms = rotation_permutations(lat)
            @test !isempty(perms)
            @test 4 in map(p -> permutation_order(lat, p), perms)
        end

        @testset "a triangular lattice has a six-fold axis" begin
            lat = build(Triangular([3, 3], 1.0; periodic=true))
            @test 6 in map(p -> permutation_order(lat, p),
                rotation_permutations(lat))
        end

        @testset "a cubic lattice has three-fold and four-fold axes" begin
            # The body diagonals of a cube carry C3 axes and its faces carry C4 -- the pair
            # that distinguishes O from the lower cubic groups.
            lat = build(Cube(3; periodic=true))
            orders = map(p -> permutation_order(lat, p), rotation_permutations(lat))
            @test 3 in orders
            @test 4 in orders
        end
    end

    @testset "non-symmetries are rejected" begin
        lat = build(Hypercube([8]; periodic=true))
        # A chain of 8 sites has no 90-degree rotation to speak of; asking for a reflection
        # about a non-existent second axis is an argument error, not a silent wrong answer.
        @test_throws ArgumentError reflection_permutation(lat, 2)
        @test_throws ArgumentError translation_permutation(lat, 2)

        # A translation along an open direction moves sites off the lattice.
        open_lat = build(Hypercube([4, 4], 1.0; periodic=[true, false]))
        @test !is_symmetry(open_lat, Translation([0, 1]))
        @test is_symmetry(open_lat, Translation([1, 0]))
        @test_throws ArgumentError site_permutation(open_lat, Translation([0, 1]))
    end

    # ------------------------------------------------------------------------ extensions
    # Everything below loads an optional package. The dependency-weight testset above must
    # already have run.

    @testset "SymBasis extension" begin
        # Loading SymBasis activates LatticeSpaceGroupsSymBasisExt, which lets the symmetry
        # specifications be built from a lattice instead of a hand-written permutation.
        using SymBasis

        lat = build(Hypercube([8]; periodic=true))
        dofo = dof_object(Spin(1 // 2))

        @testset "constructors accept a lattice" begin
            @test Translational(0, lat).perm == mod1.((1:8) .+ 1, 8)
            @test SpatialReflection(1, lat).perm == reflection_permutation(lat, 1)
            # A chain's only proper rotation about a point is the identity, so there is
            # nothing for `Rotational` to be generated from.
            @test_throws ArgumentError Rotational(1, lat)
        end

        @testset "momentum sectors partition the Hilbert space" begin
            # The strongest available check that the permutations are right: if the generator
            # were wrong, the sector dimensions would not sum to the full space.
            full = length(basis(dofo, 8).states)
            @test full == 2^8

            dims = [
                length(basis(dofo, 8, sym(Translational(k, lat), dofo)).states) for k in 0:7
            ]
            @test sum(dims) == full

            # ...and again within a magnetization sector.
            sz = sym(TotalMagnetization(0 // 1, 8), dofo)
            @test length(basis(dofo, 8, sz).states) == binomial(8, 4)
            sector_dims = [
                length(basis(dofo, 8, sz ∘ sym(Translational(k, lat), dofo)).states)
                for k in 0:7
            ]
            @test sum(sector_dims) == binomial(8, 4)
        end

        @testset "a 2-D lattice works the same way" begin
            square = build(Hypercube([2, 3], 1.0; periodic=true))
            n = n_sites(square)
            dims = [
                length(basis(dofo, n, sym(Translational(k, square; axis=2), dofo)).states)
                for k in 0:2
            ]
            @test sum(dims) == 2^n
        end
    end

    @testset "Graphs extension" begin
        using Graphs

        lat = build(Square(4; periodic=true))
        g = SimpleGraph(lat)

        @test Graphs.nv(g) == n_sites(lat) == 16
        @test Graphs.ne(g) == length(bonds(lat)) == 32
        @test Graphs.nv(lat) == n_sites(lat)
        @test Graphs.ne(lat) == length(bonds(lat))
        @test all(has_edge(g, i, j) for (i, j) in bonds(lat))
        # Every site of a periodic square lattice has four neighbours, and the graph agrees.
        @test all(degree(g, v) == 4 for v in Graphs.vertices(g))

        @testset "a multi-site basis converts too" begin
            honey = build(Honeycomb([3, 3], 1.0; periodic=true))
            gh = SimpleGraph(honey)
            @test Graphs.nv(gh) == 18
            @test all(degree(gh, v) == 3 for v in Graphs.vertices(gh))
            @test is_connected(gh)
        end
    end

    @testset "MetaGraphsNext extension" begin
        using MetaGraphsNext

        lat = build(Honeycomb([2, 2], 1.0; periodic=true))
        mg = MetaGraph(lat)

        @test length(collect(MetaGraphsNext.labels(mg))) == n_sites(lat)
        @testset "vertex codes match the lattice site numbering" begin
            # This is what makes the metagraph safe to use alongside `site_permutation`.
            for (i, label) in enumerate(site_labels(lat))
                @test code_for(mg, label) == i
                @test mg[label] == site_positions(lat)[i]
            end
        end
        @testset "edges carry their neighbour shell" begin
            labels = site_labels(lat)
            for (i, j) in bonds(lat)
                @test haskey(mg, labels[i], labels[j])
                @test mg[labels[i], labels[j]] == 1
            end
        end
    end

    @testset "Aqua quality assurance" begin
        Aqua.test_all(LatticeSpaceGroups)
    end
end
