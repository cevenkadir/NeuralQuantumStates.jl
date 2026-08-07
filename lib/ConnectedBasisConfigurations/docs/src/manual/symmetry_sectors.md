```@meta
CurrentModule = ConnectedBasisConfigurations
```

# Symmetry sectors

A symmetry-reduced basis keeps one representative per symmetry orbit instead of every
configuration. Matrix elements between representatives are not the raw ones: a connected
configuration generally lands somewhere in the middle of an orbit, and has to be folded back
onto its representative, picking up the character of the operation that got there and a ratio
of orbit norms.

Pass the basis and [`connected_padded`](@ref) does all of that:

```@example sectors
using ConnectedBasisConfigurations, OperatorAlgebra, SymBasis

spec, nsites = Spin(1 // 2), 8
dofo = dof_object(spec)
ops = local_operators(spec)

H = OpSum(vcat(
    [Op(2 .* ops.sz, i) * Op(2 .* ops.sz, mod1(i + 1, nsites)) for i in 1:nsites],
    [Op(2 .* ops.sx, i) for i in 1:nsites],
))

perm = [mod1(i + 1, nsites) for i in 1:nsites]         # translate by one site
b = basis(dofo, nsites, sym(Translational(0, perm), dofo))

res = connected_padded(H, b.states, b)
length(b.states), size(res.configs)
```

The folded matrix element is

```math
\langle m \vert \hat{H} \vert n \rangle =
    \chi \, \sqrt{\frac{\mathcal{N}_m}{\mathcal{N}_n}} \,
    \langle s' \vert \hat{H} \vert s \rangle
```

where ``\chi`` is the character of the symmetry operation mapping ``s'`` onto representative
``m``, and ``\mathcal{N}`` are the orbit norms. Configurations whose representative is absent
from the basis fall outside the sector and are dropped.

!!! note "Matrix elements are complex here"
    Even for a real operator. The character need not be real — a momentum sector's is
    ``e^{i k}`` — so `mels` comes back complex on this path where the unreduced one stays real.

## Representatives are merged, unlike the unreduced path

The unreduced kernel deliberately leaves two terms reaching the same configuration as two rows.
The symmetry path does the opposite and sums them, because here it is not a matter of
convenience: distinct configurations of one orbit *are* the same basis state, so leaving them
separate would misreport how many basis states the sector connects to. Entries that cancel to
exactly zero — which happens often between orbit members — are dropped.

## Compiling a sector

[`compile`](@ref) accepts a basis, and caches the state-to-index lookup along with the
flattened operator. Without it, that lookup is rebuilt over the entire basis on every call:

```@example sectors
sector = compile(H, b)
res2 = connected_padded(sector, b.states)

res2.configs == res.configs, res2.mels == res.mels
```

## Samples must be representatives

Only states that are in the basis can be folded, so a sample that is not one is an error rather
than a silently dropped row:

```@example sectors
outsider = first(basis(dofo, nsites).states)
outsider in b.states
```

## Checking the result

Getting the character or the norm factors wrong yields a matrix that is the right size, sparse,
and Hermitian, and has the wrong spectrum. The invariant worth testing is that the sector
spectra **partition** the full one: concatenate the eigenvalues of every sector and they must
reproduce, with multiplicity, the eigenvalues of the unreduced Hamiltonian. The test suite
checks exactly that for momentum sectors, magnetization sectors, and their combination.
