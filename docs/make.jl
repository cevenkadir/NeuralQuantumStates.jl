using NeuralQuantumStates
using Documenter
using DocumenterVitepress

DocMeta.setdocmeta!(NeuralQuantumStates, :DocTestSetup, :(using NeuralQuantumStates); recursive=true)

makedocs(;
    modules=[NeuralQuantumStates],
    authors="Kadir Çeven",
    repo="https://github.com/cevenkadir/NeuralQuantumStates.jl",
    sitename="NeuralQuantumStates.jl",
    format=DocumenterVitepress.MarkdownVitepress(
        repo="https://github.com/cevenkadir/NeuralQuantumStates.jl",
        devurl="dev",
        deploy_url="cevenkadir.github.io/NeuralQuantumStates.jl",
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
