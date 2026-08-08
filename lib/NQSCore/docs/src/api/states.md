```@meta
CurrentModule = NQSCore
```

# Variational states

The two state types and the estimators built on them, together with the connected-configuration
buffers they keep so that an optimization loop neither recompiles the operator nor reallocates
its workspace on every step.

## Index

```@index
Pages = ["states.md"]
```

```@autodocs
Modules = [NQSCore]
Pages = ["kernel.jl", "states.jl"]
```
