using Documenter

using NQSAnsatze
using NQSCore
using NQSOptimisers
using NQSSamplers
using NeuralQuantumStates

DocMeta.setdocmeta!(
    NeuralQuantumStates, :DocTestSetup, :(using NeuralQuantumStates); recursive=true
)

# One site for the ecosystem, minus the separately registered packages. LatticeSpaceGroups and
# ConnectedConfigs each carry their own site under their own directory on the same `gh-pages`
# branch; the rest are internal and documented here. Listing every module means `checkdocs`
# covers them all in the build that actually runs.
const ECOSYSTEM = [
    NeuralQuantumStates,
    NQSCore,
    NQSAnsatze,
    NQSSamplers,
    NQSOptimisers,
]

const REFERENCE_PAGES = [
    "Public API" => "lib/public.md",
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
