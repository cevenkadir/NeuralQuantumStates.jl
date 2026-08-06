using Documenter

using LatticeSpaceGroups
using NeuralQuantumStates

DocMeta.setdocmeta!(
    NeuralQuantumStates, :DocTestSetup, :(using NeuralQuantumStates); recursive=true
)

makedocs(;
    modules=[NeuralQuantumStates, LatticeSpaceGroups],
    authors="Kadir Çeven",
    repo=Remotes.GitHub("cevenkadir", "NeuralQuantumStates.jl"),
    sitename="NeuralQuantumStates.jl",
    format=Documenter.HTML(;
        canonical="https://cevenkadir.github.io/NeuralQuantumStates.jl",
        assets=["assets/favicon.ico"],
        edit_link="main",
        # The ansatz and operator pages carry a fair amount of LaTeX.
        mathengine=Documenter.MathJax3(),
        # The API reference is one big `@autodocs` dump; it is expected to be large.
        size_threshold_ignore=["lib/public.md", "lib/latticespacegroups.md"],
    ),
    pages=[
        "Home" => "index.md",
        "Basics" => "basics.md",
        "Manual" => Any[
            "Lattices"=>"manual/lattices.md",
            "Symmetries"=>"manual/symmetries.md",
        ],
        "Reference" => Any[
            "Public API"=>"lib/public.md",
            "LatticeSpaceGroups"=>"lib/latticespacegroups.md",
        ],
    ],
    # The split is in progress: fail the build on a broken cross-reference rather than
    # letting the docs quietly rot while modules move between packages.
    checkdocs=:exports,
)

deploydocs(;
    repo="github.com/cevenkadir/NeuralQuantumStates.jl",
    branch="gh-pages",
    devbranch="main",
    push_preview=true,
)
