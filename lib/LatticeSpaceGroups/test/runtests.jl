using LatticeSpaceGroups
using LatticeSpaceGroups: _site_key, _integer_point_candidates, compensating_translation

using LinearAlgebra: det, norm
using Test

"""Names of packages currently loaded in this session."""
loaded_package_names() = Set(m.name for m in keys(Base.loaded_modules))

"""Compose two site permutations: first `a`, then `b`."""
compose(b, a) = [b[a[i]] for i in eachindex(a)]

@testset "LatticeSpaceGroups.jl" begin
    @testset "dependency weight" begin
        # This package is Tier 1: it must be usable for plain exact diagonalization without
        # dragging in the neural-network stack. This is the single property the whole package
        # split exists to buy, so it is asserted rather than assumed.
        loaded = loaded_package_names()
        for heavy in ("Lux", "Enzyme", "Reactant", "Zygote", "CUDA", "Metal", "ComponentArrays")
            @test heavy ∉ loaded
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
        end

        @testset "scalar periodic applies to every dimension" begin
            @test Hypercube([3, 3], 1.0; periodic=true).periodic == [true, true]
            @test Hypercube([3, 3], 1.0).periodic == [false, false]
        end

        @testset "dimension is carried in the type" begin
            @test Hypercube([4], 1.0) isa AbstractLatticeSpec{1}
            @test Hypercube([4, 4], 1.0) isa AbstractLatticeSpec{2}
            @test Kagome([2, 2], 1.0) isa AbstractLatticeSpec{2}
            @test Triclinic([2, 2, 2], [1.0, 1.0, 1.0], [90.0, 90.0, 90.0]) isa
                  AbstractLatticeSpec{3}
        end

        @testset "specs are printable" begin
            @test occursin("Hypercube", sprint(show, Hypercube([4, 4], 1.0; periodic=true)))
        end
    end

    @testset "geometry" begin
        @test nv(build(Hypercube([8], 1.0; periodic=[true]))) == 8
        @test nv(build(Hypercube([3, 3], 1.0))) == 9
        @test nv(build(Hypercube([3, 3, 3], 1.0))) == 27
        @test nv(build(Triangular([3, 3], 1.0))) == 9
        # Multi-site unit cells: honeycomb has 2 sites per cell, kagome 3.
        @test nv(build(Honeycomb([3, 3], 1.0))) == 18
        @test nv(build(Kagome([3, 3], 1.0))) == 27
        @test nv(build(Triclinic([2, 2, 2], [1.0, 1.5, 2.0], [80.0, 70.0, 60.0]))) == 8

        @testset "site_positions is indexed by vertex number" begin
            lat = build(Hypercube([4], 1.0; periodic=[true]))
            @test site_positions(lat) == [[0.0], [1.0], [2.0], [3.0]]
            @test length(site_labels(lat)) == nv(lat)
        end

        @testset "a dimension of extent 1 cannot be periodic" begin
            @test_throws ArgumentError build(Hypercube([1, 4], 1.0; periodic=true))
            # ...but it is fine when that dimension is left open.
            @test nv(build(Hypercube([1, 4], 1.0; periodic=[false, true]))) == 4
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

    @testset "translations" begin
        @testset "a periodic chain translates cyclically" begin
            # This is the exact permutation that was previously hand-written as
            # `mod1.((1:N) .+ 1, N)`; reproducing it is the whole point of the package.
            for n in (4, 6, 8)
                lat = build(Hypercube([n], 1.0; periodic=[true]))
                @test translation_permutation(lat, 1) == mod1.((1:n) .+ 1, n)
            end
        end

        @testset "translating by the full extent is the identity" begin
            lat = build(Hypercube([6], 1.0; periodic=[true]))
            @test translation_permutation(lat, 1; cells=6) == collect(1:6)
        end

        @testset "translations compose additively" begin
            lat = build(Hypercube([6], 1.0; periodic=[true]))
            t1 = translation_permutation(lat, 1; cells=1)
            t2 = translation_permutation(lat, 1; cells=2)
            @test compose(t1, t1) == t2
        end

        @testset "open directions admit no translation" begin
            lat = build(Hypercube([6], 1.0; periodic=[false]))
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
            ("chain PBC", build(Hypercube([8], 1.0; periodic=[true])), 2),      # C_i
            ("chain open", build(Hypercube([8], 1.0; periodic=[false])), 2),    # C_i about midpoint
            ("square PBC", build(Hypercube([4, 4], 1.0; periodic=true)), 8),    # D4
            ("square open", build(Hypercube([4, 4], 1.0; periodic=false)), 8),  # D4
            ("triangular", build(Triangular([3, 3], 1.0; periodic=true)), 12),  # D6
            ("honeycomb", build(Honeycomb([3, 3], 1.0; periodic=true)), 12),    # D6
            ("kagome", build(Kagome([3, 3], 1.0; periodic=true)), 12),          # D6
            ("cube", build(Hypercube([3, 3, 3], 1.0; periodic=true)), 48),      # O_h
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
        # differs in general from the canonical τ picked for R₁R₂ by a lattice vector — so
        # the point-group permutations alone are genuinely not closed, and only become so
        # once every translation is included.
        for lat in (
            build(Hypercube([3, 3], 1.0; periodic=true)),
            build(Triangular([3, 3], 1.0; periodic=true)),
            build(Honeycomb([2, 2], 1.0; periodic=true)),
        )
            perms = Set(site_permutation(lat, op) for op in space_group(lat))
            identity_perm = collect(1:nv(lat))

            @test identity_perm in perms
            @test all(compose(p, q) in perms for p in perms, q in perms)
            @test all(invperm(p) in perms for p in perms)
        end
    end

    @testset "symmetries preserve the edge set" begin
        for lat in (
            build(Hypercube([6], 1.0; periodic=[true])),
            build(Hypercube([4, 4], 1.0; periodic=true)),
            build(Triangular([3, 3], 1.0; periodic=true)),
            build(Honeycomb([3, 3], 1.0; periodic=true)),
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
            build(Hypercube([4, 4], 1.0; periodic=true)),
            build(Triangular([3, 3], 1.0; periodic=true)),
        )
            @test length(space_group(lat)) ==
                  length(point_group(lat)) * length(translation_group(lat))
        end
    end

    @testset "reflections" begin
        @testset "an open chain reflects about its midpoint" begin
            lat = build(Hypercube([6], 1.0; periodic=[false]))
            @test reflection_permutation(lat, 1) == collect(6:-1:1)
        end

        @testset "a reflection is its own inverse" begin
            for lat in (
                build(Hypercube([6], 1.0; periodic=[false])),
                build(Hypercube([6], 1.0; periodic=[true])),
                build(Hypercube([4, 4], 1.0; periodic=true)),
            )
                p = reflection_permutation(lat, 1)
                @test compose(p, p) == collect(1:nv(lat))
            end
        end
    end

    @testset "rotations" begin
        @testset "a square lattice has a four-fold axis" begin
            lat = build(Hypercube([4, 4], 1.0; periodic=true))
            perms = rotation_permutations(lat)
            @test !isempty(perms)
            # Some returned rotation must have order exactly 4.
            orders = map(perms) do p
                o, q = 1, p
                while q != collect(1:nv(lat))
                    q = compose(p, q)
                    o += 1
                end
                o
            end
            @test 4 in orders
        end

        @testset "a triangular lattice has a six-fold axis" begin
            lat = build(Triangular([3, 3], 1.0; periodic=true))
            perms = rotation_permutations(lat)
            orders = map(perms) do p
                o, q = 1, p
                while q != collect(1:nv(lat))
                    q = compose(p, q)
                    o += 1
                end
                o
            end
            @test 6 in orders
        end
    end

    @testset "SymBasis extension" begin
        # Loading SymBasis activates LatticeSpaceGroupsSymBasisExt, which lets the symmetry
        # specifications be built from a lattice instead of a hand-written permutation.
        using SymBasis

        lat = build(Hypercube([8], 1.0; periodic=[true]))
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
            n = nv(square)
            dims = [
                length(basis(dofo, n, sym(Translational(k, square; axis=2), dofo)).states)
                for k in 0:2
            ]
            @test sum(dims) == 2^n
        end
    end

    @testset "non-symmetries are rejected" begin
        lat = build(Hypercube([8], 1.0; periodic=[true]))
        # A chain of 8 sites has no 90-degree rotation to speak of; asking for a reflection
        # about a non-existent second axis is a bounds error, not a silent wrong answer.
        @test_throws ArgumentError reflection_permutation(lat, 2)
        @test_throws ArgumentError translation_permutation(lat, 2)

        # A translation along an open direction moves sites off the lattice.
        open_lat = build(Hypercube([4, 4], 1.0; periodic=[true, false]))
        @test !is_symmetry(open_lat, Translation([0, 1]))
        @test is_symmetry(open_lat, Translation([1, 0]))
    end
end
