using Documenter, DiffRules

makedocs(modules=[DiffRules],
         doctest = false,
         format = :html,
         sitename = "DiffRules",
         pages = ["Documentation" => "index.md"])

deploydocs(repo = "github.com/JuliaDiff/DiffRules.jl.git",
           osname = "linux",
           julia = "0.6",
           target = "build",
           deps = nothing,
           make = nothing)
