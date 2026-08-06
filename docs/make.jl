using Documenter

using ConnectedConfigs
using LatticeSpaceGroups
using NQSAnsatze
using NQSCore
using NQSOptimisers
using NQSSamplers
using NeuralQuantumStates

DocMeta.setdocmeta!(
    NeuralQuantumStates, :DocTestSetup, :(using NeuralQuantumStates); recursive=true
)

# One site for the whole ecosystem. Each `lib/` package used to carry its own Documenter
# environment whose output was never deployed and which no CI job built, so a renamed export
# could break it silently. Listing every module here instead means `checkdocs` covers all of
# them in the build that actually runs.
const ECOSYSTEM = [
    NeuralQuantumStates,
    LatticeSpaceGroups,
    ConnectedConfigs,
    NQSCore,
    NQSAnsatze,
    NQSSamplers,
    NQSOptimisers,
]

const REFERENCE_PAGES = [
    "Public API" => "lib/public.md",
    "LatticeSpaceGroups" => "lib/latticespacegroups.md",
    "ConnectedConfigs" => "lib/connectedconfigs.md",
    "NQSCore" => "lib/nqscore.md",
    "NQSAnsatze" => "lib/nqsansatze.md",
    "NQSSamplers" => "lib/nqssamplers.md",
    "NQSOptimisers" => "lib/nqsoptimisers.md",
]

makedocs(;
    modules=ECOSYSTEM,
    authors="Kadir Çeven",
    repo=Remotes.GitHub("cevenkadir", "NeuralQuantumStates.jl"),
    sitename="NeuralQuantumStates.jl",
    format=Documenter.HTML(;
        canonical="https://cevenkadir.github.io/NeuralQuantumStates.jl",
        assets=["assets/favicon.ico"],
        edit_link="main",
        # The ansatz and operator pages carry a fair amount of LaTeX.
        mathengine=Documenter.MathJax3(),
        # Each reference page is one big `@autodocs` dump; they are expected to be large.
        size_threshold_ignore=[last(p) for p in REFERENCE_PAGES],
    ),
    pages=[
        "Home" => "index.md",
        "Basics" => "basics.md",
        "Manual" => Any[
            "Lattices"=>"manual/lattices.md",
            "Symmetries"=>"manual/symmetries.md",
            "Connected configurations"=>"manual/connectedconfigs.md",
            "Variational states"=>"manual/nqscore.md",
            "Ansätze"=>"manual/nqsansatze.md",
            "Samplers"=>"manual/nqssamplers.md",
            "Optimisers"=>"manual/nqsoptimisers.md",
        ],
        "Reference" => Any[REFERENCE_PAGES...],
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
