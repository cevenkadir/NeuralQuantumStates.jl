using Documenter

using LatticeSpaceGroups

DocMeta.setdocmeta!(
    LatticeSpaceGroups, :DocTestSetup, :(using LatticeSpaceGroups); recursive=true
)

makedocs(;
    modules=[LatticeSpaceGroups],
    authors="Kadir Çeven",
    repo=Remotes.GitHub("cevenkadir", "NeuralQuantumStates.jl"),
    sitename="LatticeSpaceGroups.jl",
    format=Documenter.HTML(;
        canonical="https://cevenkadir.github.io/NeuralQuantumStates.jl/LatticeSpaceGroups",
        edit_link="main",
        mathengine=Documenter.MathJax3(),
        size_threshold_ignore=["lib/public.md"],
    ),
    pages=[
        "Home" => "index.md",
        "Manual" => Any[
            "Lattices"=>"manual/lattices.md",
            "Symmetries"=>"manual/symmetries.md",
        ],
        "Reference" => Any["Public API"=>"lib/public.md"],
    ],
    checkdocs=:exports,
)
