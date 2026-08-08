@testset "flattened operators" begin
    """The dense matrix a padded result stands for, which is the property that has to survive."""
    function dense_from(res, states)
        index = Dict(s => i for (i, s) in pairs(states))
        M = zeros(ComplexF64, length(states), length(states))
        for k in eachindex(states), j in 1:res.counts[k]
            M[index[res.configs[j, k]], k] += res.mels[j, k]
        end
        return M
    end

    """Every model worth checking, as (label, operator, dof, nsites)."""
    models = [
        ("TFI, spin-1/2", transverse_field_ising(4; J=1.0, h_x=0.7, h_z=0.2), Spin(1 // 2), 4),
        ("TFI, larger", transverse_field_ising(8; J=1.0, h_x=0.9, h_z=0.1), Spin(1 // 2), 8),
        ("Bose-Hubbard", extended_bose_hubbard(4, 2), Boson(2), 4),
        ("Bose-Hubbard, n_max=3", extended_bose_hubbard(4, 3), Boson(3), 4),
    ]

    @testset "$label agrees with the nested kernel" for (label, H, spec, nsites) in models
        states = basis(dof_object(spec), nsites).states
        nested = connected_padded(compile(H), states)
        flat = connected_padded(flatten(H), states)

        # Not merely the same operator, but the same output: same counts, same entries, in the
        # same order. A flattening that reordered would still be correct and would still be a
        # behaviour change, so it is worth pinning.
        @test flat.counts == nested.counts
        @test size(flat.configs) == size(nested.configs)
        for k in eachindex(states), j in 1:nested.counts[k]
            @test flat.configs[j, k] == nested.configs[j, k]
            @test flat.mels[j, k] == nested.mels[j, k]
        end
        @test dense_from(flat, states) == dense_from(nested, states)
    end

    @testset "an operator with no off-diagonal terms" begin
        # max_conn collapses to one, which is the edge case where the diagonal slot is the only
        # slot and the compaction path has nothing to compact.
        H = transverse_field_ising(6; J=1.0, h_x=0.0, h_z=0.3)
        states = basis(dof_object(Spin(1 // 2)), 6).states
        op = flatten(H)
        @test max_conn_size(op) == 1
        res = connected_padded(op, states)
        @test res.counts == connected_padded(compile(H), states).counts
        @test dense_from(res, states) == dense_from(connected_padded(compile(H), states), states)
    end

    @testset "it reports what it is" begin
        H = transverse_field_ising(4; J=1.0, h_x=0.7, h_z=0.2)
        op = flatten(H)
        @test op isa FlatOperator
        @test eltype(op) === eltype(compile(H))
        @test max_conn_size(op) == max_conn_size(compile(H))
        @test occursin("FlatOperator", sprint(show, op))

        # Flattening something already compiled and something not must agree.
        @test flatten(compile(H)).colptr == op.colptr
    end

    @testset "the layout is flat all the way down" begin
        # The point of the type: every field is a plain array of numbers, so the whole thing can
        # be moved to a device in one step. A nested field would defeat that silently.
        op = flatten(transverse_field_ising(4; J=1.0, h_x=0.7, h_z=0.2))
        for name in fieldnames(typeof(op))
            f = getfield(op, name)
            @test f isa Union{Integer,AbstractVector{<:Number}}
            f isa AbstractVector && @test isbitstype(eltype(f))
        end
    end

    @testset "buffers are checked before they are written" begin
        H = transverse_field_ising(4; J=1.0, h_x=0.7, h_z=0.2)
        op = flatten(H)
        states = basis(dof_object(Spin(1 // 2)), 4).states
        h, n = max_conn_size(op), length(states)

        good_c = Matrix{eltype(states)}(undef, h, n)
        good_m = Matrix{eltype(op)}(undef, h, n)
        good_k = Vector{Int}(undef, n)

        @test_throws DimensionMismatch connected_padded!(
            Matrix{eltype(states)}(undef, h - 1, n), good_m, good_k, op, states)
        @test_throws DimensionMismatch connected_padded!(
            good_c, Matrix{eltype(op)}(undef, h, n - 1), good_k, op, states)
        @test_throws DimensionMismatch connected_padded!(
            good_c, good_m, Vector{Int}(undef, n - 1), op, states)
    end

    @testset "the in-place form does not allocate per sample" begin
        H = transverse_field_ising(6; J=1.0, h_x=0.7, h_z=0.2)
        op = flatten(H)
        states = basis(dof_object(Spin(1 // 2)), 6).states
        h = max_conn_size(op)

        function run(n)
            s = states[1:n]
            c = Matrix{eltype(states)}(undef, h, n)
            m = Matrix{eltype(op)}(undef, h, n)
            k = Vector{Int}(undef, n)
            connected_padded!(c, m, k, op, s)          # warm up
            return @allocated connected_padded!(c, m, k, op, s)
        end

        # Whatever the scratch costs, it is set by the operator's branching and not by the
        # batch, so a batch eight times larger must not allocate more.
        @test run(8) == run(64)
    end
end
