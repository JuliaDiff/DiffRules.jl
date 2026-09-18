# Rules are methods, so they survive precompilation of the package defining them, and the
# extensions are loaded while a package depending on SpecialFunctions precompiles. Neither
# is observable in-process, so this builds a throwaway package and loads it from a fresh
# Julia.

@testset "precompilation" begin
    mktempdir() do dir
        mkpath(joinpath(dir, "Probe", "src"))
        write(joinpath(dir, "Probe", "Project.toml"), """
            name = "Probe"
            uuid = "2b9d4c69-3a3d-4f1e-9d84-9f9b5a3f0e01"
            version = "0.1.0"

            [deps]
            DiffRules = "b552c78f-8df3-52c6-915a-8e097449b14b"
            SpecialFunctions = "276daf66-3868-5448-9aa4-cd146d93841b"
            """)
        write(joinpath(dir, "Probe", "src", "Probe.jl"), """
            module Probe
            using DiffRules, SpecialFunctions
            cube(x) = x^3
            DiffRules.@define_diffrule Probe.cube(x) = :(3 * \$x^2)
            const RULES = DiffRules.diffrules(; filter_modules=nothing)
            end
            """)
        write(joinpath(dir, "Project.toml"), """
            [deps]
            DiffRules = "b552c78f-8df3-52c6-915a-8e097449b14b"
            Probe = "2b9d4c69-3a3d-4f1e-9d84-9f9b5a3f0e01"

            [sources]
            Probe = { path = "Probe" }
            DiffRules = { path = "$(escape_string(pkgdir(DiffRules)))" }
            """)

        script = """
            using Pkg; Pkg.instantiate(; io=devnull)
            using Probe, DiffRules
            println((:Probe, :cube, 1) in DiffRules.diffrules(; filter_modules=nothing))
            println(any(r -> r[1] === :SpecialFunctions, Probe.RULES))
            """
        out = read(`$(Base.julia_cmd()) --startup-file=no --project=$dir -e $script`, String)
        rule_survives, ext_loaded_at_precompile = split(strip(out), '\n')

        @test rule_survives == "true"
        @test ext_loaded_at_precompile == "true"
    end
end
