using QuantiKite
using Documenter

DocMeta.setdocmeta!(QuantiKite, :DocTestSetup, :(using QuantiKite); recursive=true)

makedocs(;
    modules=[QuantiKite],
    authors="Fernando Penaranda <fernandopenaranda@github.com> and contributors",
    repo="https://github.com/fernandopenaranda/QuantiKite.jl/blob/{commit}{path}#{line}",
    sitename="QuantiKite.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://fernandopenaranda.github.io/QuantiKite.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/fernandopenaranda/QuantiKite.jl",
    devbranch="main",
)
