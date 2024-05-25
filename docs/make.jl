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
    ],
    warnonly=true,
)

deploydocs(;
    repo="github.com/cevenkadir/NeuralQuantumStates.jl",
    push_preview=true,
)
