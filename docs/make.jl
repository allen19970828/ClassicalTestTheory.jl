using Documenter
using ClassicalTestTheory

makedocs(
    sitename = "ClassicalTestTheory.jl",
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true"
    ),
    modules = [ClassicalTestTheory],
    pages = [
        "Home" => "index.md",
        "API Reference" => "api.md"
    ],
)

deploydocs(
    repo = "github.com/allen19970828/ClassicalTestTheory.jl.git",
    devbranch = "main"
)
