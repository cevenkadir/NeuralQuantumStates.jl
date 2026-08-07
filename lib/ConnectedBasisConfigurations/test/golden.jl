"""
Equivalence with the pre-split `Operators.connected_basis_configs`.

This is the acceptance gate for the whole package: the reference arrays in `test/golden/` (inside this package) were
generated from the deleted implementation, and reproducing them is what licenses deleting it.

# What "equivalent" means here

Not element-wise array equality, because two conventions deliberately changed:

1. **Padding.** The old code padded short columns with `missing`, giving
   `Union{T,Missing}` arrays. This package pads matrix elements with `0`, which keeps the
   arrays concretely typed and makes the local-energy sum correct with no masking.
2. **Axis order.** The old code was internally inconsistent — its Ising batch method put the
   degree-of-freedom axis last and its Bose-Hubbard method put it first (see
   `test/golden/README.md`). This package always puts the DoF axis first.

The physically meaningful object is the map from each connected configuration to its total
matrix element, so that is what is compared: both sides are reduced to a
`Dict(configuration => matrix element)`, with duplicates summed and exact zeros dropped. Two
implementations agree exactly when those dictionaries agree.
"""

using ConnectedBasisConfigurations
using SymBasis
using Test

# Inside the package, not at the monorepo root: a registered package is published as its own
# subdirectory tree, so data reached for across `../../..` simply would not exist for anyone who
# installed it, and `Pkg.test` would fail on a missing file.
const GOLDEN_DIR = joinpath(@__DIR__, "golden")
include(joinpath(GOLDEN_DIR, "tfi_chain8.jl"))
include(joinpath(GOLDEN_DIR, "bhm_chain16.jl"))

"""
    golden_dict(configs, mels; dof_axis) -> Dict{Vector,Float64}

Reduce a pre-split `(configs, mels)` pair to configuration => total matrix element.

`missing` entries are padding and are skipped; duplicate configurations are summed; entries
that end up exactly zero are dropped, matching what this package emits.
"""
function golden_dict(configs::AbstractMatrix, mels::AbstractVector)
    out = Dict{Vector{eltype(configs)},Float64}()
    for c in axes(configs, 2)
        col = configs[:, c]
        any(ismissing, col) && continue
        ismissing(mels[c]) && continue
        key = collect(skipmissing(col))
        out[key] = get(out, key, 0.0) + mels[c]
    end
    filter!(p -> !iszero(p.second), out)
    return out
end

"""Reduce this package's output for batch column `b` to the same dictionary form."""
function packed_dict(spec, res, b::Integer, nsites::Integer)
    out = Dict{Vector{eltype(local_values(spec))},Float64}()
    for j in 1:res.counts[b]
        key = configurations(spec, res.configs[j, b], nsites)
        out[key] = get(out, key, 0.0) + real(res.mels[j, b])
    end
    filter!(p -> !iszero(p.second), out)
    return out
end

"""Compare two configuration => matrix-element maps, reporting the first discrepancy."""
function agree(a::Dict, b::Dict; atol=1e-10)
    if keys(a) != keys(b)
        only_a = setdiff(keys(a), keys(b))
        only_b = setdiff(keys(b), keys(a))
        @info "configuration sets differ" n_only_reference = length(only_a) n_only_new = length(only_b) sample_only_reference = first(only_a, 2) sample_only_new = first(only_b, 2)
        return false
    end
    for (k, v) in a
        if !isapprox(v, b[k]; atol=atol)
            @info "matrix element differs" configuration = k reference = v new = b[k]
            return false
        end
    end
    return true
end

@testset "the spin convention, which the golden data cannot pin down" begin
    # The stored sample has total magnetization zero, so `h_z * Σ σᶻ` vanishes on it and the
    # golden comparison passes whether σᶻ is diag(-1,+1) or diag(+1,-1). Substituting
    # OperatorAlgebra's PAULI_Z for 2Sᶻ is therefore invisible to every test above — verified
    # by mutation, not assumed. These checks use magnetized configurations, where the two
    # conventions differ by the whole field term.
    spec, nsites = Spin(1 // 2), 8
    H = transverse_field_ising(nsites; J=1.0, h_x=1.0, h_z=1.0)

    up = fill(1 // 2, nsites)
    down = fill(-1 // 2, nsites)

    """Diagonal matrix element ⟨s|H|s⟩ of a configuration."""
    function diagonal(H, spec, config, nsites)
        s = packed(spec, config)
        d = connected(H, s)
        return real(get(d, s, 0.0))
    end

    # Fully polarized: every bond contributes J, every site contributes h_z.
    @test diagonal(H, spec, up, nsites) ≈ 1.0 * nsites + 1.0 * nsites
    # Flipping every spin flips the field term but not the bond term.
    @test diagonal(H, spec, down, nsites) ≈ 1.0 * nsites - 1.0 * nsites

    @test diagonal(H, spec, up, nsites) != diagonal(H, spec, down, nsites)

    # And the zero-magnetization golden sample genuinely cannot tell them apart, which is why
    # the checks above are needed.
    @test sum(2 .* tfi_sample_vec) == 0
end

@testset "equivalence with the pre-split implementation" begin
    @testset "transverse-field Ising, 8-site periodic chain" begin
        spec, nsites = Spin(1 // 2), 8
        H = transverse_field_ising(nsites; J=1.0, h_x=1.0, h_z=1.0)

        @testset "single configuration" begin
            state = packed(spec, tfi_sample_vec)
            res = connected_padded(H, [state])

            reference = golden_dict(tfi_single_configs, tfi_single_mels)
            @test agree(reference, packed_dict(spec, res, 1, nsites))

            # The reference has 9 entries: the diagonal plus one flip per site.
            @test length(reference) == 1 + nsites
        end

        @testset "batch" begin
            # The reference batch arrays are (M_max, batch, N): DoF axis LAST.
            samples = [tfi_sample_mat[b, :] for b in axes(tfi_sample_mat, 1)]
            states = [packed(spec, s) for s in samples]
            res = connected_padded(H, states)

            for b in eachindex(samples)
                reference = Dict{Vector{Rational{Int}},Float64}()
                for j in axes(tfi_batch_configs, 1)
                    row = tfi_batch_configs[j, b, :]
                    any(ismissing, row) && continue
                    ismissing(tfi_batch_mels[j, b]) && continue
                    key = collect(skipmissing(row))
                    reference[key] = get(reference, key, 0.0) + tfi_batch_mels[j, b]
                end
                filter!(p -> !iszero(p.second), reference)
                @test agree(reference, packed_dict(spec, res, b, nsites))
            end
        end
    end

    @testset "extended Bose-Hubbard, 16-site periodic chain" begin
        spec, nsites, n_max = Boson(5), 16, 5
        H = extended_bose_hubbard(nsites, n_max; J=1.0, U=1.0, V=1.0, μ=0.0)

        @testset "single configuration" begin
            state = packed(spec, bhm_sample_vec)
            res = connected_padded(H, [state])

            reference = golden_dict(bhm_single_configs, bhm_single_mels)
            @test agree(reference, packed_dict(spec, res, 1, nsites))

            # 1 diagonal + 10 hops for 5 particles on a 16-site ring, none blocked.
            @test length(reference) == 11
        end

        @testset "batch" begin
            # The reference batch arrays are (N, M_max, batch): DoF axis FIRST.
            samples = [bhm_sample_mat[:, b] for b in axes(bhm_sample_mat, 2)]
            states = [packed(spec, s) for s in samples]
            res = connected_padded(H, states)

            for b in eachindex(samples)
                reference = Dict{Vector{Int},Float64}()
                for j in axes(bhm_batch_configs, 2)
                    col = bhm_batch_configs[:, j, b]
                    any(ismissing, col) && continue
                    ismissing(bhm_batch_mels[1, j, b]) && continue
                    key = collect(skipmissing(col))
                    reference[key] = get(reference, key, 0.0) + bhm_batch_mels[1, j, b]
                end
                filter!(p -> !iszero(p.second), reference)
                @test agree(reference, packed_dict(spec, res, b, nsites))
            end
        end
    end
end
