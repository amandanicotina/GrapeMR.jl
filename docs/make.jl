using Documenter, DocumenterCitations, Literate, GrapeMR

Literate.markdown("./docs/src/tutorial.jl", "./docs/src/generated"; mdstrings=true)

pages = [
  "Introduction" => "index.md",
  "Tutorial" => "generated/tutorial.md",
  "API" => "api.md",
  "References" => "references.md",
]

bib = CitationBibliography(
    joinpath(".", "CITATIONS.bib");
    style=:numeric
)

# compile to HTML:
makedocs(;pages,
    sitename = "GrapeMR",
    format = Documenter.HTML(
        assets = ["assets/grape.ico"]
    ),
    modules = [GrapeMR],
    checkdocs=:exports,
    plugins = [bib]
)
  
deploydocs(
    repo = "github.com/amandanicotina/GrapeMR.jl.git",
    push_preview = true,
)
