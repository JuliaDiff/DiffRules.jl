using Documenter, DiffRules

DocMeta.setdocmeta!(
    DiffRules,
    :DocTestSetup,
    :(using DiffRules);
    recursive=true,
)

makedocs(modules=[DiffRules],
         sitename = "DiffRules",
         pages = ["Documentation" => "index.md"],
         format = Documenter.HTML(
                  prettyurls = get(ENV, "CI", nothing) == "true"
         ),
         strict=true,
         checkdocs=:exports,
)

deploydocs(; repo="github.com/JuliaDiff/DiffRules.jl", push_preview=true)
