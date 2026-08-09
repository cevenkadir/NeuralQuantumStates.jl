using Documenter
using ConnectedBasisConfigurations

# The navbar's repository link, and the subdirectory rather than the repository root: this
# package lives inside a monorepo, so the root would land a reader on the umbrella's README.
# Only the navbar link is affected -- the per-page "Edit on GitHub" links come from `repo`
# below, which must stay the plain repository for Documenter to resolve source paths.
const REPOLINK =
    "https://github.com/cevenkadir/NeuralQuantumStates.jl/tree/main/lib/ConnectedBasisConfigurations"
const CANONICAL = "https://cevenkadir.github.io/NeuralQuantumStates.jl/ConnectedBasisConfigurations/"

DocMeta.setdocmeta!(
    ConnectedBasisConfigurations,
    :DocTestSetup,
    :(using ConnectedBasisConfigurations);
    recursive=true,
)

@info "Generating Documenter.jl site"
makedocs(;
    modules=[ConnectedBasisConfigurations],
    authors="Kadir Çeven",
    repo=Remotes.GitHub("cevenkadir", "NeuralQuantumStates.jl"),
    sitename="ConnectedBasisConfigurations.jl",
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
            "Connected basis configurations" => "manual/connected_states.md",
            "Local operators" => "manual/local_operators.md",
            "Symmetry sectors" => "manual/symmetry_sectors.md",
            "Backends" => "manual/backends.md",
        ],
        "Examples" => [
            "Exact diagonalization" => "examples/exact_diagonalization.md",
        ],
        "API Reference" => [
            "The kernel" => "api/connected.md",
            "Compiling operators" => "api/compile.md",
            "Configurations and local operators" => "api/configurations.md",
            "Running on a device" => "api/device.md",
            "Backend interface" => "api/interface.md",
        ],
    ],
    checkdocs=:exports,
)

# This package lives in a subdirectory of the NeuralQuantumStates.jl monorepo, so its site is
# deployed alongside the umbrella's rather than in place of it.
#
# `dirname` confines the whole deployment -- pages, `versions.js`, the redirect and the version
# symlinks -- to the `ConnectedBasisConfigurations/` subfolder of `gh-pages`, leaving the
# umbrella's root-level equivalents untouched.
#
# `tag_prefix` is the other half. Registered subdirectory packages are tagged
# `ConnectedBasisConfigurations-v0.1.0`, not `v0.1.0`; without the prefix this site would take
# its version from the umbrella's tags and publish under the wrong number, and the umbrella
# would in turn try to release on this package's tags.
@info "Deploying to GitHub"
deploydocs(;
    repo="github.com/cevenkadir/NeuralQuantumStates.jl.git",
    dirname="ConnectedBasisConfigurations",
    tag_prefix="ConnectedBasisConfigurations-",
    devbranch="main",
    push_preview=true,
)
