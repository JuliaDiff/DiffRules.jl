include("baseline.jl")

@testset "registry" begin
    @testset "baseline" begin
        rules = sort!(collect(DiffRules.diffrules(; filter_modules=nothing));
                      by = r -> (string(r[1]), string(r[2]), r[3]))
        @test rules == BASELINE
        @testset "$M.$f/$n" for (M, f, n) in BASELINE
            @test DiffRules.hasdiffrule(M, f, n)
        end
    end

    # `parentmodule` reports `Base.Math` for these, which would break `@eval $M.$f`.
    @testset "Base.Math functions report :Base" begin
        @testset "$f" for f in (:sec, :sind, :hypot, :cot, :deg2rad, :mod2pi, :rem2pi)
            @test any(r -> r[1] === :Base && r[2] === f,
                      DiffRules.diffrules(; filter_modules=nothing))
        end
    end

    @testset "lookup by function" begin
        @test DiffRules.diffrule(Base.sin, :x) == :(cos(x))
        @test DiffRules.diffrule(Base.sec, :x) == :(sec(x) * tan(x))
        @test DiffRules.hasdiffrule(Base.sin, 1)
        @test !DiffRules.hasdiffrule(Base.sin, 2)
        # `log` has a rule at both arities
        @test DiffRules.hasdiffrule(Base.log, 1) && DiffRules.hasdiffrule(Base.log, 2)
    end

    @testset "lookup by module" begin
        @test DiffRules.diffrule(Base, :sin, :x) == :(cos(x))
        @test DiffRules.hasdiffrule(Base, :sin, 1)
        @test !DiffRules.hasdiffrule(Base, :sin, 2)
    end

    # Rules splice in the functions they call, so evaluating one does not require the
    # defining package to be in scope.
    @testset "rules evaluate without the package in scope" begin
        mod = Module(:NoImports)
        d = DiffRules.diffrule(:SpecialFunctions, :digamma, :x)
        @test Base.eval(mod, :(let x = 1.5; $d end)) ≈ SpecialFunctions.trigamma(1.5)
    end

    @testset "unknown rules" begin
        @test !DiffRules.hasdiffrule(:Base, :nonexistent, 1)
        @test !DiffRules.hasdiffrule(:NotLoaded, :f, 1)
        @test_throws Exception DiffRules.diffrule(:Base, :nonexistent, :x)
    end

    @testset "filter_modules" begin
        @test all(M === :Base for (M, _, _) in DiffRules.diffrules(; filter_modules=(:Base,)))
        @test all(M in (:Base, :SpecialFunctions, :NaNMath)
                  for (M, _, _) in DiffRules.diffrules())
    end

    # undocumented, but overloaded by IntervalArithmetic and Zygote
    @testset "_abs_deriv hook" begin
        @test DiffRules._abs_deriv(-2.0) == -1.0
        @test DiffRules.diffrule(Base.abs, :x) == :($(DiffRules._abs_deriv)(x))
    end
end
