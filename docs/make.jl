using Documenter, LocalBarycentric

makedocs(
    sitename = "LocalBarycentric.jl",
    modules = [LocalBarycentric],
    pages = [
        "Home" => "index.md",
        "Theory" => "theory.md",
        "Examples" => "examples.md",
    ],
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
    checkdocs = :exports,
)

# Update `repo` when the package moves to its own repository.
deploydocs(repo = "github.com/Canonical111/LocalBarycentric.jl", devbranch = "main")
