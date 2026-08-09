using KernelAbstractions

@testset "the KernelAbstractions kernel" begin
    # KernelAbstractions runs the same kernel on a CPU backend as on a GPU one, so the kernel
    # that would run on a device can be checked here — against the reference implementation,
    # on real models, without any GPU present. What a device adds beyond this is the compiler
    # and the memory, not the algorithm.
    backend = CPU()

    """Run the kernel and return the padded result, in the shape the allocating form gives."""
    function via_kernel(op, states)
        h, n = max_conn_size(op), length(states)
        configs = Matrix{eltype(states)}(undef, h, n)
        mels = Matrix{eltype(op)}(undef, h, n)
        counts = Vector{Int}(undef, n)
        connected_padded!(configs, mels, counts, op, states, backend)
        return (; configs=configs, mels=mels, counts=counts)
    end

    models = [
        ("TFI, spin-1/2", transverse_field_ising(4; J=1.0, h_x=0.7, h_z=0.2), Spin(1 // 2), 4),
        ("TFI, larger", transverse_field_ising(8; J=1.0, h_x=0.9, h_z=0.1), Spin(1 // 2), 8),
        ("TFI, no off-diagonal", transverse_field_ising(6; J=1.0, h_x=0.0, h_z=0.3), Spin(1 // 2), 6),
        ("Bose-Hubbard", extended_bose_hubbard(4, 2), Boson(2), 4),
        ("Bose-Hubbard, n_max=3", extended_bose_hubbard(4, 3), Boson(3), 4),
    ]

    @testset "$label matches the reference kernel" for (label, H, spec, nsites) in models
        states = basis(dof_object(spec), nsites).states
        op = flatten(H)
        reference = connected_padded!(
            Matrix{eltype(states)}(undef, max_conn_size(op), length(states)),
            Matrix{eltype(op)}(undef, max_conn_size(op), length(states)),
            Vector{Int}(undef, length(states)),
            op, states,
        )
        got = via_kernel(op, states)

        @test got.counts == reference.counts
        # Every slot, padding included: the padded ones are what a consumer reducing over the
        # whole column relies on, so "equal where it counts" is not the property wanted here.
        @test got.configs == reference.configs
        @test got.mels == reference.mels
    end

    @testset "padding is inert" begin
        # The reduction downstream sums the whole column without a mask, which is only correct
        # because a padded slot repeats the sample with a zero matrix element — never a stale
        # configuration, and never something that makes `0 * exp(...)` anything but zero.
        H = transverse_field_ising(6; J=1.0, h_x=0.7, h_z=0.2)
        spec = Spin(1 // 2)
        states = basis(dof_object(spec), 6).states
        op = flatten(H)
        got = via_kernel(op, states)
        for k in eachindex(states), j in (got.counts[k]+1):size(got.mels, 1)
            @test iszero(got.mels[j, k])
            @test got.configs[j, k] == states[k]
        end
    end

    @testset "unchecked digit access agrees with the checked kind" begin
        # The kernel bypasses the bounds checks because they cannot be compiled for a device.
        # That is only safe if it computes the same thing, at every position and every digit.
        ext = Base.get_extension(
            ConnectedBasisConfigurations, :ConnectedBasisConfigurationsKernelAbstractionsExt
        )
        @test ext !== nothing
        for spec in (Spin(1 // 2), Spin(1 // 1), Boson(3))
            nsites = 4
            for state in basis(dof_object(spec), nsites).states
                for pos in 1:nsites
                    @test ext.unchecked_read(state, pos) == read_digit(state, pos)
                    for d in 0:(local_dimension(spec)-1)
                        @test ext.unchecked_write(state, pos, d) == write_digit(state, pos, d)
                    end
                end
            end
        end
    end

    @testset "buffers are checked before the kernel launches" begin
        H = transverse_field_ising(4; J=1.0, h_x=0.7, h_z=0.2)
        op = flatten(H)
        states = basis(dof_object(Spin(1 // 2)), 4).states
        h, n = max_conn_size(op), length(states)
        c = Matrix{eltype(states)}(undef, h, n)
        m = Matrix{eltype(op)}(undef, h, n)
        k = Vector{Int}(undef, n)

        @test_throws DimensionMismatch connected_padded!(
            Matrix{eltype(states)}(undef, h - 1, n), m, k, op, states, backend)
        @test_throws DimensionMismatch connected_padded!(
            c, m, Vector{Int}(undef, n - 1), op, states, backend)
    end

    @testset "an empty batch is not an error" begin
        H = transverse_field_ising(4; J=1.0, h_x=0.7, h_z=0.2)
        op = flatten(H)
        empty = eltype(basis(dof_object(Spin(1 // 2)), 4).states)[]
        res = connected_padded!(
            Matrix{eltype(empty)}(undef, max_conn_size(op), 0),
            Matrix{eltype(op)}(undef, max_conn_size(op), 0),
            Vector{Int}(undef, 0), op, empty, backend,
        )
        @test isempty(res.counts)
    end

    @testset "the unpacking kernel" begin
        # `configurations!` is the second half of keeping a batch on the device: the connected
        # configurations are computed there, and this turns them into the network's input
        # without them ever coming back. It must agree with the host `configurations` — the
        # only intended difference being the element type, since the host returns exact
        # rationals for a spin and no accelerator can hold one.
        @testset "$spec matches the host unpacking" for spec in
                                                        (Spin(1 // 2), Spin(1 // 1), Boson(3))
            nsites = 4
            states = basis(dof_object(spec), nsites).states
            values = collect(Float64, local_values(spec))

            out = Matrix{Float64}(undef, nsites, length(states))
            configurations!(out, values, states, nsites, backend)
            @test out == Float64.(configurations(spec, states, nsites))
        end

        @testset "states of any shape are read in linear order" begin
            # The caller is the local-energy path, whose states are the `(max_conn, batch)`
            # block of connected configurations; its columns have to line up with `vec` of
            # that block, not with some reshaping of it.
            spec, nsites = Spin(1 // 2), 4
            H = transverse_field_ising(nsites; J=1.0, h_x=0.7, h_z=0.2)
            res = connected_padded(H, basis(dof_object(spec), nsites).states)
            values = collect(Float64, local_values(spec))

            out = Matrix{Float64}(undef, nsites, length(res.configs))
            configurations!(out, values, res.configs, nsites, backend)
            @test out == Float64.(configurations(spec, vec(res.configs), nsites))
        end

        @testset "a wrong-sized output is rejected before the kernel launches" begin
            spec, nsites = Spin(1 // 2), 4
            states = basis(dof_object(spec), nsites).states
            values = collect(Float64, local_values(spec))
            @test_throws DimensionMismatch configurations!(
                Matrix{Float64}(undef, nsites + 1, length(states)),
                values, states, nsites, backend)
            @test_throws DimensionMismatch configurations!(
                Matrix{Float64}(undef, nsites, length(states) - 1),
                values, states, nsites, backend)
        end

        @testset "an empty batch is not an error" begin
            spec = Spin(1 // 2)
            empty = eltype(basis(dof_object(spec), 4).states)[]
            out = Matrix{Float64}(undef, 4, 0)
            @test configurations!(out, collect(Float64, local_values(spec)), empty, 4, backend) === out
        end
    end

    @testset "to_backend moves every array" begin
        op = flatten(transverse_field_ising(4; J=1.0, h_x=0.7, h_z=0.2))
        moved = to_backend(op, backend)
        # On a CPU backend the move is a copy rather than a transfer, but the invariant that
        # matters is that the content survives and the metadata is carried across unchanged.
        @test moved.colptr == op.colptr
        @test moved.vals == op.vals
        @test moved.outs == op.outs
        @test moved.n_diagonal == op.n_diagonal
        @test max_conn_size(moved) == max_conn_size(op)

        states = basis(dof_object(Spin(1 // 2)), 4).states
        @test via_kernel(moved, states).counts == via_kernel(op, states).counts
    end
end
