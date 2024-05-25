using Documenter, DocumenterVitepress

using NeuralQuantumStates

makedocs(;
    modules=[NeuralQuantumStates],
    authors="Kadir Çeven",
    repo="https://github.com/cevenkadir/NeuralQuantumStates.jl",
    sitename="Chairmarks.jl",
    format=DocumenterVitepress.MarkdownVitepress(
        repo="https://github.com/cevenkadir/NeuralQuantumStates.jl",
        devurl="dev",
        deploy_url="cevenkadir.github.io/NeuralQuantumStates.jl",
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
