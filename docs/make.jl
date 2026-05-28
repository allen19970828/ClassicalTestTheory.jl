using Documenter
using CTT

makedocs(
    sitename = "CTT.jl",
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true"
    ),
    modules = [CTT],
    pages = [
        "Home" => "index.md",
        "API Reference" => "api.md"
    ],
)

deploydocs(
    repo = "github.com/allen19970828/CTT.jl.git",
    devbranch = "main"
)
