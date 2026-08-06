using Documenter

using NQSAnsatze

DocMeta.setdocmeta!(NQSAnsatze, :DocTestSetup, :(using NQSAnsatze); recursive=true)

makedocs(;
    modules=[NQSAnsatze],
    authors="Kadir Çeven",
    repo=Remotes.GitHub("cevenkadir", "NeuralQuantumStates.jl"),
    sitename="NQSAnsatze.jl",
    format=Documenter.HTML(;
        canonical="https://cevenkadir.github.io/NeuralQuantumStates.jl/NQSAnsatze",
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
