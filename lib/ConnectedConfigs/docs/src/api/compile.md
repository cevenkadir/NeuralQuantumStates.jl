```@meta
CurrentModule = ConnectedConfigs
```

# Compiling operators

[`compile`](@ref) turns an operator into a [`CompiledOperator`](@ref) — or, given a
symmetry-reduced basis, a [`CompiledSector`](@ref) — hoisting every configuration-independent
decision out of the hot loop. [`max_conn_size`](@ref) is the bound that
[`connected_padded!`](@ref) buffers are sized by.

## Index

```@index
Pages = ["compile.md"]
```

```@autodocs
Modules = [ConnectedConfigs]
Pages = ["compile.jl"]
```
