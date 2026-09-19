# Which rules exist depends on what is loaded, which the test suite itself cannot vary: it loads
# all providers. So probe throwaway environments from fresh Julia processes.

const LOAD_PATH_SEP = Sys.iswindows() ? ';' : ':'

# `Pkg.test` runs us with `JULIA_LOAD_PATH="@:<pkgdir>/test"`. A child inherits that, losing
# `@stdlib` and gaining our test environment, so pin it to its own project plus the stdlibs.
function probe(dir, script)
    cmd = `$(Base.julia_cmd()) --startup-file=no --project=$dir -e $script`
    return read(
        addenv(
            cmd,
            "JULIA_LOAD_PATH" => "@$(LOAD_PATH_SEP)@stdlib",
            "JULIA_PROJECT" => nothing,
        ),
        String,
    )
end

const DIFFRULES_UUID = "b552c78f-8df3-52c6-915a-8e097449b14b"
const NANMATH_UUID = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
const SPECIALFUNCTIONS_UUID = "276daf66-3868-5448-9aa4-cd146d93841b"
const PROBE_UUID = "2b9d4c69-3a3d-4f1e-9d84-9f9b5a3f0e01"

function write_env(dir, deps; probe_pkg=false)
    mkpath(dir)
    entries = join(("$name = \"$uuid\"" for (name, uuid) in deps), "\n")
    sources = "DiffRules = { path = \"$(escape_string(pkgdir(DiffRules)))\" }"
    if probe_pkg
        entries *= "\nProbe = \"$PROBE_UUID\""
        sources *= "\nProbe = { path = \"$(escape_string(joinpath(dir, "Probe")))\" }"
    end
    write(joinpath(dir, "Project.toml"), """
        [deps]
        $entries

        [sources]
        $sources
        """)
    return dir
end

@testset "precompilation" begin
    mktempdir() do root
        # `Probe.RULES` records what was visible at `Probe`'s own precompile time.
        mkpath(joinpath(root, "both", "Probe", "src"))
        write(joinpath(root, "both", "Probe", "Project.toml"), """
            name = "Probe"
            uuid = "$PROBE_UUID"
            version = "0.1.0"

            [deps]
            DiffRules = "$DIFFRULES_UUID"
            SpecialFunctions = "$SPECIALFUNCTIONS_UUID"
            """)
        write(joinpath(root, "both", "Probe", "src", "Probe.jl"), """
            module Probe
            using DiffRules, SpecialFunctions
            cube(x) = x^3
            DiffRules.@define_diffrule Probe.cube(x) = :(3 * \$x^2)
            const RULES = DiffRules.diffrules(; filter_modules=nothing)
            end
            """)

        both = write_env(
            joinpath(root, "both"),
            ("DiffRules" => DIFFRULES_UUID,
             "NaNMath" => NANMATH_UUID,
             "SpecialFunctions" => SPECIALFUNCTIONS_UUID);
            probe_pkg=true,
        )
        bare = write_env(joinpath(root, "bare"), ("DiffRules" => DIFFRULES_UUID,))
        nanmath = write_env(
            joinpath(root, "nanmath"),
            ("DiffRules" => DIFFRULES_UUID, "NaNMath" => NANMATH_UUID),
        )

        @testset "rules of a dependent package survive its precompilation" begin
            out = probe(both, """
                using Pkg; Pkg.instantiate(; io=devnull)
                using Probe, DiffRules, NaNMath
                println("rule_survives=",
                        (:Probe, :cube, 1) in DiffRules.diffrules(; filter_modules=nothing))
                println("rule_resolves=", DiffRules.diffrule(:Probe, :cube, :x) == :(3 * x^2))
                println("ext_at_precompile=", any(r -> r[1] === :SpecialFunctions, Probe.RULES))
                println("lgamma=", DiffRules.hasdiffrule(:NaNMath, :lgamma, 1))
                """)
            @test occursin("rule_survives=true", out)
            @test occursin("rule_resolves=true", out)
            @test occursin("ext_at_precompile=true", out)
            @test occursin("lgamma=true", out)
        end

        @testset "no provider loaded, no provider rules" begin
            out = probe(bare, """
                using Pkg; Pkg.instantiate(; io=devnull)
                using DiffRules
                println("erf=", DiffRules.hasdiffrule(:SpecialFunctions, :erf, 1))
                println("nanmath_sin=", DiffRules.hasdiffrule(:NaNMath, :sin, 1))
                println("base_only=", all(M === :Base for (M, _, _) in
                                          DiffRules.diffrules(; filter_modules=nothing)))
                """)
            @test occursin("erf=false", out)
            @test occursin("nanmath_sin=false", out)
            @test occursin("base_only=true", out)
        end

        # `NaNMath.lgamma` differentiates to `SpecialFunctions.digamma` (#106).
        @testset "NaNMath alone does not define the lgamma rule" begin
            out = probe(nanmath, """
                using Pkg; Pkg.instantiate(; io=devnull)
                using DiffRules, NaNMath
                println("nanmath_sin=", DiffRules.hasdiffrule(:NaNMath, :sin, 1))
                println("lgamma=", DiffRules.hasdiffrule(:NaNMath, :lgamma, 1))
                """)
            @test occursin("nanmath_sin=true", out)
            @test occursin("lgamma=false", out)
        end
    end
end
