using Documenter, DocumenterTools
using MINOTAUR

makedocs(
    sitename = "MINOTAUR",
    format = Documenter.HTML(),
    modules = [MINOTAUR],
    pages = [
    "Home" => "index.md",
    ],
    doctest = true,
)

deploydocs(repo = "github.com/matthewozon/MINOTAUR.git",branch = "master") 

