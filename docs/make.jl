using NeuralQuantumStates
using Documenter

DocMeta.setdocmeta!(NeuralQuantumStates, :DocTestSetup, :(using NeuralQuantumStates); recursive=true)

makedocs(;
    modules=[NeuralQuantumStates],
    authors="Kadir Çeven",
    repo="https://github.com/cevenkadir/NeuralQuantumStates.jl/blob/{commit}{path}#{line}",
    sitename="NeuralQuantumStates.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://kadirceven.com/NeuralQuantumStates.jl",
        edit_link="main",
        assets=String[]
    ),
    pages=[
        "Home" => "index.md",
    ]
)

deploydocs(;
    repo="github.com/cevenkadir/NeuralQuantumStates.jl",
    devbranch="main"
)
