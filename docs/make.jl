push!(LOAD_PATH,"../src/")

using Documenter, GrapeMR

# Literate.markdown("./src/tutorial.jl", "./src")

pages = [
"Introduction" => "index.md",
# "Tutorial" => "tutorial.md",
 # "API" => "api.md",
]

# compile to HTML:
makedocs(;pages,
    sitename = "GrapeMR",
    format = Documenter.HTML(),
    modules = [GrapeMR]
)
# makedocs(;pages,  modules = [GrapeMR])
#     # format = Documenter.HTML())
#     # sitename = "GrapeMR")
  
# deploydocs(
#     repo = "github.com/amandanicotina/GrapeMR.jl.git"
# )
# # Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
