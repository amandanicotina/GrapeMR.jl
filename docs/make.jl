using Documenter
using GrapeMR

makedocs(
    sitename = "GrapeMR",
    format = Documenter.HTML(),
    modules = [GrapeMR]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
