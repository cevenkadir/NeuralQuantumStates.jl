```@meta
CurrentModule = ConnectedBasisConfigurations
```

# Exact diagonalization

The kernel produces a Hamiltonian one column at a time, which is exactly what building a sparse
matrix needs. Doing it this way also exercises the same code path variational Monte Carlo uses,
so agreement with a known spectrum is evidence about the kernel and not only about the model.

## Building the matrix

Each sample contributes one column: `res.configs[j, b]` is a row index once looked up in the
basis, and `res.mels[j, b]` is the entry. Contributions are **added**, not assigned, because
two terms may reach the same configuration.

```@example ed
using ConnectedBasisConfigurations, OperatorAlgebra, SymBasis, LinearAlgebra

function hamiltonian_matrix(H, spec, nsites)
    b = basis(dof_object(spec), nsites)
    index = Dict(s => i for (i, s) in pairs(b.states))

    res = connected_padded(compile(H), b.states)
    M = zeros(ComplexF64, length(b.states), length(b.states))
    for n in eachindex(b.states), j in 1:res.counts[n]
        M[index[res.configs[j, n]], n] += res.mels[j, n]
    end
    return M
end
nothing # hide
```

## A case with a known answer

Non-interacting spins in a transverse field: ``H = h_x \sum_i \sigma^x_i`` has eigenvalues
``h_x (2k - N)`` with multiplicity ``\binom{N}{k}``, because the ``\sigma^x`` eigenvalues are
``\pm 1`` independently on every site.

```@example ed
spec, nsites, h_x = Spin(1 // 2), 6, 0.8
ops = local_operators(spec)

H = OpSum([h_x * Op(2 .* ops.sx, i) for i in 1:nsites])
M = hamiltonian_matrix(H, spec, nsites)

found = sort(real(eigvals(Hermitian(M))))
expected = sort([h_x * (2k - nsites) for k in 0:nsites for _ in 1:binomial(nsites, k)])

found ≈ expected
```

## The transverse-field Ising chain

With the interaction switched on there is no closed form for finite ``N`` at general couplings,
but the matrix is still small enough to diagonalize directly:

```@example ed
H = OpSum(vcat(
    [Op(2 .* ops.sz, i) * Op(2 .* ops.sz, mod1(i + 1, nsites)) for i in 1:nsites],
    [Op(2 .* ops.sx, i) for i in 1:nsites],
))

M = hamiltonian_matrix(H, spec, nsites)
M ≈ M'          # Hermitian, as it must be
```

```@example ed
minimum(real(eigvals(Hermitian(M))))
```

## Checking the sparsity

[`max_conn_size`](@ref) bounds the number of non-zeros per column without looking at a single
configuration, which is what a sparse assembly wants to preallocate by:

```@example ed
compiled = compile(H)
res = connected_padded(compiled, basis(dof_object(spec), nsites).states)

max_conn_size(compiled), maximum(res.counts)
```

For a lattice Hamiltonian whose off-diagonal terms are single hops or flips, the bound is tight
— one slot for the diagonal plus one per off-diagonal term.

## Working in a symmetry sector

Blocking by symmetry gives smaller matrices whose spectra together reproduce the full one. See
[Symmetry sectors](@ref) for the folding rules; the assembly is unchanged apart from indexing
into the reduced basis:

```@example ed
perm = [mod1(i + 1, nsites) for i in 1:nsites]
dofo = dof_object(spec)

spectrum = Float64[]
for k in 0:(nsites-1)
    b = basis(dofo, nsites, sym(Translational(k, perm), dofo))
    isempty(b.states) && continue

    index = Dict(s => i for (i, s) in pairs(b.states))
    res = connected_padded(H, b.states, b)
    block = zeros(ComplexF64, length(b.states), length(b.states))
    for n in eachindex(b.states), j in 1:res.counts[n]
        block[index[res.configs[j, n]], n] += res.mels[j, n]
    end
    append!(spectrum, real(eigvals(Hermitian(block))))
end

sort(spectrum) ≈ sort(real(eigvals(Hermitian(M))))
```
