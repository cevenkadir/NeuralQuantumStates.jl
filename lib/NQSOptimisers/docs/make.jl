using Documenter

using NQSOptimisers

DocMeta.setdocmeta!(NQSOptimisers, :DocTestSetup, :(using NQSOptimisers); recursive=true)

makedocs(;
    modules=[NQSOptimisers],
    authors="Kadir Çeven",
    repo=Remotes.GitHub("cevenkadir", "NeuralQuantumStates.jl"),
    sitename="NQSOptimisers.jl",
    format=Documenter.HTML(;
        canonical="https://cevenkadir.github.io/NeuralQuantumStates.jl/NQSOptimisers",
        edit_link="main",
        mathengine=Documenter.MathJax3(),
        size_threshold_ignore=["lib/public.md"],
    ),
    pages=[
        "Home" => "index.md",
        "Reference" => Any["Public API"=>"lib/public.md"],
    ],
    checkdocs=:exports,
)
