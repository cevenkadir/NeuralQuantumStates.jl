using Documenter, DocumenterVitepress

using NeuralQuantumStates

makedocs(;
    modules=[NeuralQuantumStates],
    authors="Kadir Çeven",
    repo="github.com/cevenkadir/NeuralQuantumStates.jl",
    sitename="NeuralQuantumStates.jl",
    format=DocumenterVitepress.MarkdownVitepress(
        repo="github.com/cevenkadir/NeuralQuantumStates.jl",
        devurl="dev"
    ),
    pages=[
        "Home" => "index.md",
        "Basics" => "basics.md",
        "Manual" => Any[
            "Lattices"=>"manual/lattices.md",
            "Hilberts"=>"manual/hilberts.md",
            "Operators"=>"manual/operators.md",
        ],
        "Reference" => Any[
            "Public API"=>"lib/public.md",
        ],
    ],
)

deploydocs(;
    repo="github.com/cevenkadir/NeuralQuantumStates.jl",
    target="build", # this is where Vitepress stores its output
    branch="gh-pages",
    push_preview=true,
)
