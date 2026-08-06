using Documenter

using NQSCore

DocMeta.setdocmeta!(NQSCore, :DocTestSetup, :(using NQSCore); recursive=true)

makedocs(;
    modules=[NQSCore],
    authors="Kadir Çeven",
    repo=Remotes.GitHub("cevenkadir", "NeuralQuantumStates.jl"),
    sitename="NQSCore.jl",
    format=Documenter.HTML(;
        canonical="https://cevenkadir.github.io/NeuralQuantumStates.jl/NQSCore",
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
