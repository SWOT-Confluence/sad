using Sad
using Documenter

DocMeta.setdocmeta!(Sad, :DocTestSetup, :(using Sad); recursive=true)

makedocs(;
    modules=[Sad],
    authors="Kostas Andreadis <kandread@umass.edu> and contributors",
    repo="https://gitlab.com/kandread/Sad.jl/blob/{commit}{path}#{line}",
    sitename="Sad.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://kandread.gitlab.io/Sad.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)
