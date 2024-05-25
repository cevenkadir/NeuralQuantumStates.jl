using NeuralQuantumStates
using Documenter
using DocumenterVitepress

DocMeta.setdocmeta!(NeuralQuantumStates, :DocTestSetup, :(using NeuralQuantumStates); recursive=true)

makedocs(;
    modules=[NeuralQuantumStates],
    repo=Remotes.GitHub("cevenkadir", "NeuralQuantumStates.jl"),
    authors="Kadir Çeven",
    sitename="NeuralQuantumStates.jl",
    format=DocumenterVitepress.MarkdownVitepress(
        repo="https://github.com/cevenkadir/NeuralQuantumStates.jl",
        devurl="dev",
        deploy_url="https://cevenkadir.github.io/NeuralQuantumStates.jl",
    ),
    pages=[
        "Home" => "index.md",
        "Manual" => Any[
            "Lattices"=>"manual/lattices.md",
        ],
        "Reference" => Any[
            "Public API"=>"lib/public.md",
        ],
    ]
)

deploydocs(;
    repo="github.com/cevenkadir/NeuralQuantumStates.jl",
    target="build", # this is where Vitepress stores its output
    devbranch="main",
    branch="gh-pages",
    push_preview=true
)
