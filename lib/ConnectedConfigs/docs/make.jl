using Documenter

using ConnectedConfigs

DocMeta.setdocmeta!(ConnectedConfigs, :DocTestSetup, :(using ConnectedConfigs); recursive=true)

makedocs(;
    modules=[ConnectedConfigs],
    authors="Kadir Çeven",
    repo=Remotes.GitHub("cevenkadir", "NeuralQuantumStates.jl"),
    sitename="ConnectedConfigs.jl",
    format=Documenter.HTML(;
        canonical="https://cevenkadir.github.io/NeuralQuantumStates.jl/ConnectedConfigs",
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
