using Documenter
using LaplacianOpt

makedocs(
    sitename = "LaplacianOpt",
    format = Documenter.HTML(),
    modules = [LaplacianOpt]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
