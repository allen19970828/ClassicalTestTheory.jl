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
    remotes = nothing
)

# deploydocs(
#     repo = "github.com/USERNAME/CTT.jl.git",
# )
