"""
Tests for the compiled kernel: that compiling changes nothing observable, that the batch
shapes and padding follow the documented contract, that the in-place form allocates nothing,
and that the backend interface really is the only thing the kernel talks to.
"""

using ConnectedBasisConfigurations
using LinearAlgebra
using OperatorAlgebra
using SparseArrays
using SymBasis
using Test

"""
    oa_index(spec, state, nsites) -> Int

Position of a packed configuration in OperatorAlgebra's Kronecker ordering.

The two packages number digits in **opposite directions**: SymBasis — and so this package —
makes site 1 the least significant digit, while OperatorAlgebra's `sparse` and `apply` make it
the most significant. A chain model that is symmetric under site reversal hides the difference
entirely, which is why `exact.jl` can compare the two matrices directly; an operator that
treats sites differently does not, and comparing without this translation silently checks the
wrong matrix element.
"""
function oa_index(spec, state, nsites)
    d = local_dimension(spec)
    idx = 0
    for i in 1:nsites
        idx = idx * d + read_digit(state, i)
    end
    return idx + 1
end

"""Reduce a padded result column to `configuration => summed matrix element`."""
function as_dict(res, b, spec, nsites)
    out = Dict{Vector{Rational{Int}},ComplexF64}()
    for j in 1:res.counts[b]
        key = Rational{Int}.(configurations(spec, res.configs[j, b], nsites))
        out[key] = get(out, key, 0.0im) + res.mels[j, b]
    end
    filter!(p -> !isapprox(p.second, 0; atol=1e-12), out)
    return out
end

@testset "compiled kernel" begin
    spec, nsites = Spin(1 // 2), 6
    states = basis(dof_object(spec), nsites).states
    ops = local_operators(spec)

    bspec, bsites = Boson(3), 5
    bstates = basis(dof_object(bspec), bsites).states

    models = (
        ("transverse-field Ising", spec, nsites, states,
            transverse_field_ising(nsites; J=1.0, h_x=0.7, h_z=0.3)),
        ("extended Bose-Hubbard", bspec, bsites, bstates,
            extended_bose_hubbard(bsites, 3; J=1.0, U=0.5, V=0.25, μ=0.1)),
    )

    @testset "compiling is observationally neutral: $name" for (name, sp, n, sts, H) in models
        compiled = compile(H)
        direct = connected_padded(H, sts)
        viacompiled = connected_padded(compiled, sts)

        @test viacompiled.counts == direct.counts
        @test viacompiled.configs == direct.configs
        @test viacompiled.mels == direct.mels

        for s in sts[1:min(end, 20)]
            @test connected(compiled, s) == connected(H, s)
        end
    end

    @testset "the compile-time bound really bounds: $name" for (name, sp, n, sts, H) in models
        compiled = compile(H)
        res = connected_padded(compiled, sts)
        @test all(res.counts .<= max_conn_size(compiled))
        @test maximum(res.counts) == size(res.configs, 1)
    end

    @testset "batch shapes follow the states" begin
        H = transverse_field_ising(nsites; J=1.0, h_x=0.7, h_z=0.3)
        compiled = compile(H)

        flat = states[1:24]
        asmatrix = reshape(flat, 4, 6)
        ascube = reshape(flat, 2, 3, 4)

        rf = connected_padded(compiled, flat)
        rm = connected_padded(compiled, asmatrix)
        rc = connected_padded(compiled, ascube)

        k = size(rf.configs, 1)
        @test size(rf.configs) == (k, 24)
        @test size(rm.configs) == (k, 4, 6)
        @test size(rc.configs) == (k, 2, 3, 4)
        @test size(rm.mels) == size(rm.configs)
        @test size(rc.counts) == (2, 3, 4)

        # Reshaping the input must reshape the output and change nothing else.
        @test reshape(rm.configs, k, 24) == rf.configs
        @test reshape(rc.mels, k, 24) == rf.mels
        @test vec(rc.counts) == rf.counts
    end

    @testset "padding is inert and repeats the sample" begin
        H = extended_bose_hubbard(bsites, 3; J=1.0, U=0.5, V=0.25, μ=0.1)
        res = connected_padded(H, bstates)
        padded = [(b, j) for b in eachindex(bstates)
                  for j in (res.counts[b]+1):size(res.configs, 1)]

        @test !isempty(padded)          # otherwise the assertions below are vacuous
        # A zero matrix element makes the slot contribute nothing, and repeating the sample
        # keeps it a configuration the wavefunction can safely be evaluated on.
        @test all(iszero(res.mels[j, b]) for (b, j) in padded)
        @test all(res.configs[j, b] == bstates[b] for (b, j) in padded)
    end

    @testset "repeated terms are summed, not deduplicated" begin
        # Two terms reaching the same configuration deliberately produce two rows. What must
        # hold is that their sum is right -- checked against the sparse matrix, which knows
        # nothing about this package's conventions.
        H = OpSum(AbstractOp[
            0.7 * Op(2 .* ops.sx, 1),
            0.3 * Op(2 .* ops.sx, 1),
            1.5 * Op(2 .* ops.sz, 2),
        ])
        res = connected_padded(H, states)
        dense = Matrix(sparse(H, [i => 2 for i in 1:nsites]))

        # This operator is deliberately *not* symmetric under site reversal, so the digit-order
        # difference between the two packages has to be undone explicitly -- see `oa_index`.
        for (b, s) in pairs(states)
            reconstructed = zeros(ComplexF64, length(states))
            for j in 1:res.counts[b]
                reconstructed[oa_index(spec, res.configs[j, b], nsites)] += res.mels[j, b]
            end
            @test reconstructed ≈ dense[:, oa_index(spec, s, nsites)]
        end

        # Two `σˣ₁` terms at the same site both reach the flipped configuration, so the row
        # count exceeds the number of distinct connected configurations.
        b = 1
        distinct = length(unique(res.configs[1:res.counts[b], b]))
        @test res.counts[b] > distinct
        @test length(connected(H, states[b])) == distinct
    end

    @testset "purely diagonal operators use one slot" begin
        H = OpSum(AbstractOp[Op(2 .* ops.sz, i) * Op(2 .* ops.sz, mod1(i + 1, nsites))
                             for i in 1:nsites])
        res = connected_padded(H, states)
        @test all(res.counts .<= 1)
        @test all(res.configs[1, b] == states[b] for b in eachindex(states))
        # Every configuration of this model has non-zero bond energy, so no slot is dropped.
        @test all(res.counts .== 1)
    end

    @testset "an empty batch is handled" begin
        H = transverse_field_ising(nsites)
        res = connected_padded(H, typeof(first(states))[])
        @test isempty(res.counts)
        @test size(res.configs) == (0, 0)
        @test size(res.mels) == (0, 0)
    end

    @testset "in-place form" begin
        H = transverse_field_ising(nsites; J=1.0, h_x=0.7, h_z=0.3)
        compiled = compile(H)
        batch = states[1:32]

        height = max_conn_size(compiled)
        configs = Matrix{eltype(batch)}(undef, height, length(batch))
        mels = Matrix{Float64}(undef, height, length(batch))
        counts = Vector{Int}(undef, length(batch))

        out = connected_padded!(configs, mels, counts, compiled, batch)
        reference = connected_padded(compiled, batch)

        @test out.counts == reference.counts
        # The in-place form keeps the buffer's full height rather than trimming, so compare
        # over the rows the allocating form kept and require the rest to be inert padding.
        k = size(reference.configs, 1)
        @test out.configs[1:k, :] == reference.configs
        @test out.mels[1:k, :] == reference.mels
        for b in eachindex(batch), j in (counts[b]+1):height
            @test iszero(mels[j, b])
            @test configs[j, b] == batch[b]
        end

        @testset "allocation does not grow with the batch" begin
            # The point of the in-place form: the per-sample work allocates nothing, so what
            # is left is one fixed-size scratch buffer whose size depends on the operator's
            # branching, not on how many samples are passed. The allocating form, by
            # contrast, has to allocate the whole output on every call.
            small = states[1:8]
            large = states[1:64]
            small_out = (Matrix{eltype(small)}(undef, height, 8),
                Matrix{Float64}(undef, height, 8), Vector{Int}(undef, 8))
            large_out = (Matrix{eltype(large)}(undef, height, 64),
                Matrix{Float64}(undef, height, 64), Vector{Int}(undef, 64))

            connected_padded!(small_out..., compiled, small)
            connected_padded!(large_out..., compiled, large)
            @test @allocated(connected_padded!(small_out..., compiled, small)) ==
                  @allocated(connected_padded!(large_out..., compiled, large))
            @test @allocated(connected_padded!(large_out..., compiled, large)) <
                  @allocated(connected_padded(compiled, large))
        end

        @testset "mis-sized buffers are rejected" begin
            small = Matrix{eltype(batch)}(undef, height - 1, length(batch))
            @test_throws DimensionMismatch connected_padded!(
                small, mels, counts, compiled, batch)
            @test_throws DimensionMismatch connected_padded!(
                configs, Matrix{Float64}(undef, height, 1), counts, compiled, batch)
            @test_throws DimensionMismatch connected_padded!(
                configs, mels, Vector{Int}(undef, 1), compiled, batch)
        end
    end

    @testset "type stability" begin
        real_H = compile(transverse_field_ising(nsites; J=1.0, h_x=0.7, h_z=0.3))
        complex_H = compile(OpSum(AbstractOp[
            (0.5 + 0.25im) * Op(2 .* ops.sx, i) for i in 1:nsites
        ]))
        @test eltype(real_H) == Float64
        @test eltype(complex_H) == ComplexF64

        batch = states[1:8]
        @test @inferred(connected_padded(real_H, batch)) isa NamedTuple
        @test @inferred(connected_padded(complex_H, batch)) isa NamedTuple
        @test @inferred(connected(real_H, first(batch))) isa Dict
    end

    @testset "fermionic sites" begin
        # Jordan-Wigner strings are resolved once, at compile time. A hopping term far from
        # the chain edge is the case where an unresolved string would blow up, and where a
        # cancelled one must not be left behind.
        n = 6
        bi = [fermion(i) => 2 for i in 1:n]
        hop(i, j) = Op(OperatorAlgebra.RAISE, fermion(i)) * Op(OperatorAlgebra.LOWER, fermion(j))
        H = OpSum(AbstractOp[
            [-1.0 * (hop(i, i + 1) + hop(i + 1, i)) for i in 1:(n-1)]...,
            [0.5 * Op(OperatorAlgebra.OCC_PART, fermion(i)) for i in 1:n]...,
        ])

        compiled = compile(H)
        fstates = basis(dof_object(Spin(1 // 2)), n).states
        res = connected_padded(compiled, fstates)

        dense = Matrix(sparse(H, bi))
        for (b, s) in pairs(fstates)
            column = zeros(ComplexF64, length(fstates))
            for j in 1:res.counts[b]
                column[oa_index(Spin(1 // 2), res.configs[j, b], n)] += res.mels[j, b]
            end
            @test column ≈ dense[:, oa_index(Spin(1 // 2), s, n)]
        end

        # The string of a nearest-neighbour hop cancels, so the term stays two sites wide and
        # the connection bound does not grow with the chain.
        @test max_conn_size(compiled) == 1 + 2 * (n - 1)
    end
end

# --- the backend seam ---------------------------------------------------------------------
#
# A toy operator type that knows nothing about OperatorAlgebra, implementing only the two
# functions the compiler asks of an operator backend. If this works, the kernel really is
# decoupled from the default backend rather than merely appearing to be.

struct ToySum
    terms::Vector{Vector{Pair{Int,Matrix{Float64}}}}
end

ConnectedBasisConfigurations.expand_terms(op::ToySum) = op.terms
ConnectedBasisConfigurations.amplitude_type(::ToySum) = Float64

@testset "operator backend interface" begin
    spec, nsites = Spin(1 // 2), 5
    states = basis(dof_object(spec), nsites).states
    ops = local_operators(spec)
    σx, σz = 2 .* ops.sx, 2 .* ops.sz

    toy = ToySum([
        [[i => Matrix{Float64}(σz), mod1(i + 1, nsites) => Matrix{Float64}(σz)] for i in 1:nsites]...,
        [[i => Matrix{Float64}(σx)] for i in 1:nsites]...,
    ])
    reference = transverse_field_ising(nsites; J=1.0, h_x=1.0, h_z=0.0)

    @test compile(toy) isa CompiledOperator{Float64}

    fromtoy = connected_padded(toy, states)
    fromreference = connected_padded(reference, states)
    for b in eachindex(states)
        @test as_dict(fromtoy, b, spec, nsites) == as_dict(fromreference, b, spec, nsites)
    end

    @testset "the identity term is diagonal with value one" begin
        identity_only = ToySum([Pair{Int,Matrix{Float64}}[]])
        res = connected_padded(identity_only, states)
        @test all(res.counts .== 1)
        @test all(res.mels[1, b] == 1.0 for b in eachindex(states))
        @test all(res.configs[1, b] == states[b] for b in eachindex(states))
    end
end
