```@meta
CurrentModule = NQSCore
```

# Reference implementations

An exact ansatz and an exact sampler. Neither is meant for production use: they contain no approximation at all, which is what lets the rest of the machinery be tested against exact answers rather than against itself.

## Index

```@index
Pages = ["reference.md"]
```

```@autodocs
Modules = [NQSCore]
Pages = ["ansatz.jl", "sampler.jl"]
```
