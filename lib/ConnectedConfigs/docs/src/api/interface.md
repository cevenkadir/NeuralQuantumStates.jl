```@meta
CurrentModule = ConnectedConfigs
```

# Backend interface

The generic functions the compiler and the kernel talk to, and nothing else. Adding methods to
them is how a basis or operator library other than SymBasis and OperatorAlgebra plugs in — see
[Backends](@ref) in the manual for what each contract requires.

## Index

```@index
Pages = ["interface.md"]
```

```@autodocs
Modules = [ConnectedConfigs]
Pages = ["interface.jl"]
```
