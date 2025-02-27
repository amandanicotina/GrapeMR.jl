using Documenter, DocumenterCitations, Literate, GrapeMR

Literate.markdown("./docs/src/tutorial.jl", "./docs/src/generated"; mdstrings=true)

pages = [
  "Home" => "index.md",
  "Theoretical Background" => "pages/theory.md",
  "Tutorial" => "generated/tutorial.md",
  "Developing" => "pages/developing.md",
  "API" => "pages/api.md",
  # TODO: re-enable after adding references to docstrings and docs
  #"References" => "pages/references.md",
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
