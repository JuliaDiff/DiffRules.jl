using Documenter, DiffRules

makedocs(modules=[DiffRules],
         doctest = false,
         sitename = "DiffRules",
         pages = ["Documentation" => "index.md"],
         format = Documenter.HTML(
                  prettyurls = get(ENV, "CI", nothing) == "true"
         ),
)

deploydocs(; repo="github.com/JuliaDiff/DiffRules.jl", push_preview=true)
