```@meta
CurrentModule = ConnectedBasisConfigurations
```

# The kernel

The kernel itself: [`connected_padded`](@ref) for a batch, [`connected_padded!`](@ref) for a
batch into buffers you own, and [`connected`](@ref) for a single configuration. See
[Connected configurations](@ref) in the manual for the padding and ordering conventions.

## Index

```@index
Pages = ["connected.md"]
```

```@autodocs
Modules = [ConnectedBasisConfigurations]
Pages = ["connected.jl"]
```
