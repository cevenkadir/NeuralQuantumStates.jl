using Documenter
using NQSCore

# The navbar's repository link, and the subdirectory rather than the repository root: this
# package lives inside a monorepo, so the root would land a reader on the umbrella's README.
# Only the navbar link is affected -- the per-page "Edit on GitHub" links come from `repo`
# below, which must stay the plain repository for Documenter to resolve source paths.
const REPOLINK = "https://github.com/cevenkadir/NeuralQuantumStates.jl/tree/main/lib/NQSCore"
const CANONICAL = "https://cevenkadir.github.io/NeuralQuantumStates.jl/NQSCore/"

DocMeta.setdocmeta!(
    NQSCore,
    :DocTestSetup,
    :(using NQSCore);
    recursive=true,
)

@info "Generating Documenter.jl site"
makedocs(;
    modules=[NQSCore],
    authors="Kadir Çeven",
    repo=Remotes.GitHub("cevenkadir", "NeuralQuantumStates.jl"),
    sitename="NQSCore.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", nothing) == "true",
        repolink=REPOLINK,
        canonical=CANONICAL,
        edit_link="main",
        mathengine=Documenter.KaTeX(),
    ),
    pages=[
        "Home" => "index.md",
        "Manual" => [
            "Variational states" => "manual/states.md",
            "Statistics" => "manual/statistics.md",
            "Log-derivatives and gradients" => "manual/log_derivatives.md",
            "Extending the interfaces" => "manual/interfaces.md",
        ],
        "API Reference" => [
            "Interfaces" => "api/interface.md",
            "Variational states" => "api/states.md",
            "Statistics" => "api/stats.md",
            "Log-derivatives" => "api/log_derivatives.md",
            "Reference implementations" => "api/reference.md",
        ],
    ],
    checkdocs=:exports,
)

# This package lives in a subdirectory of the NeuralQuantumStates.jl monorepo, so its site is
# deployed alongside the umbrella's rather than in place of it.
#
# `dirname` confines the whole deployment -- pages, `versions.js`, the redirect and the version
# symlinks -- to the `NQSCore/` subfolder of `gh-pages`, leaving the umbrella's root-level
# equivalents untouched.
#
# `tag_prefix` is the other half. Registered subdirectory packages are tagged `NQSCore-v0.1.0`,
# not `v0.1.0`; without the prefix this site would take its version from the umbrella's tags and
# publish under the wrong number, and the umbrella would in turn try to release on this
# package's tags.
@info "Deploying to GitHub"
deploydocs(;
    repo="github.com/cevenkadir/NeuralQuantumStates.jl.git",
    dirname="NQSCore",
    tag_prefix="NQSCore-",
    devbranch="main",
    push_preview=true,
)
