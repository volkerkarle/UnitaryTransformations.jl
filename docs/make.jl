using UnitaryTransformations
using Documenter

DocMeta.setdocmeta!(
    UnitaryTransformations,
    :DocTestSetup,
    :(using UnitaryTransformations; using QuantumAlgebra);
    recursive = true,
)

makedocs(;
    modules = [UnitaryTransformations],
    authors = "Volker Karle",
    sitename = "UnitaryTransformations.jl",
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://volkerkarle.github.io/UnitaryTransformations.jl",
        assets = String[],
        edit_link = "main",
        repolink = "https://github.com/volkerkarle/UnitaryTransformations.jl",
    ),
    pages = [
        "Home" => "index.md",
        "Tutorial" => "tutorial.md",
        "Examples" => "examples.md",
        "API Reference" => "api.md",
    ],
    warnonly = [:missing_docs, :cross_references],
)

deploydocs(;
    repo = "github.com/volkerkarle/UnitaryTransformations.jl",
    devbranch = "main",
    push_preview = true,
)
