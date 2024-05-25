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
        ],
        "Reference" => Any[
            "Public API"=>"lib/public.md",
        ],
    ],
    warnonly=true,
)

deploydocs(;
    repo="github.com/cevenkadir/NeuralQuantumStates.jl",
    push_preview=true,
)
