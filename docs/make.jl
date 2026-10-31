using Documenter
using LaplacianOpt

makedocs(
    modules = [LaplacianOpt],
    sitename = "LaplacianOpt",
    format = Documenter.HTML(
        mathengine = Documenter.MathJax(),
        assets = ["assets/extra_styles.css", "assets/publications.js"],
        prettyurls = get(ENV, "CI", nothing) == "true",
    ),
    # strict = true,
    authors = "Harsha Nagarajan",
    pages = [
        "Introduction" => "index.md",
        "Quick Start guide" => "quickguide.md",
        "Function References" => "function_references.md",
        "Research & Publications" => "research.md",
    ],
)

deploydocs(repo = "github.com/harshangrjn/LaplacianOpt.jl.git", push_preview = false)
