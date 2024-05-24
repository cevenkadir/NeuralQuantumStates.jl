using NeuralQuantumStates
using MarkdownVitepress

DocMeta.setdocmeta!(NeuralQuantumStates, :DocTestSetup, :(using NeuralQuantumStates); recursive=true)

makedocs(;
    modules=[NeuralQuantumStates],
    authors="Kadir Çeven",
    repo="https://github.com/cevenkadir/NeuralQuantumStates.jl/blob/{commit}{path}#{line}",
    sitename="NeuralQuantumStates.jl",
    format=MarkdownVitepress(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://cevenkadir.github.io/NeuralQuantumStates.jl",
        edit_link="main",
        assets=String[]
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
    devbranch="main"
)
