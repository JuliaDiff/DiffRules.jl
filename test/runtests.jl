using DiffRules
using Test

import SpecialFunctions, NaNMath, LogExpFunctions
import Random
Random.seed!(1)

function finitediff(f, x)
    ϵ = cbrt(eps(typeof(x))) * max(one(typeof(x)), abs(x))
    return (f(x + ϵ) - f(x - ϵ)) / (ϵ + ϵ)
end

@testset "DiffRules" begin
@testset "check rules" begin

non_numeric_arg_functions = [(:Base, :rem2pi, 2), (:Base, :ifelse, 3)]

for (M, f, arity) in DiffRules.diffrules(; filter_modules=nothing)
    (M, f, arity) ∈ non_numeric_arg_functions && continue
    if arity == 1
        @test DiffRules.hasdiffrule(M, f, 1)
        deriv = DiffRules.diffrule(M, f, :goo)
        modifier = if f in (:asec, :acsc, :asecd, :acscd, :acosh, :acoth)
            1.0
        elseif f === :log1mexp
            -1.0
        elseif f === :log2mexp
            -0.5
        else
            0.0
        end
        @eval begin
            let
                goo = rand() + $modifier
                @test isapprox($deriv, finitediff($M.$f, goo), rtol=0.05)
                # test for 2pi functions
                if "mod2pi" == string($M.$f)
                    goo = 4pi + $modifier
                    @test NaN === $deriv
                end
            end
        end
    elseif arity == 2
        @test DiffRules.hasdiffrule(M, f, 2)
        derivs = DiffRules.diffrule(M, f, :foo, :bar)
        @eval begin
            let
                foo, bar = rand(1:10), rand()
                dx, dy = $(derivs[1]), $(derivs[2])
                if !(isnan(dx))
                    @test isapprox(dx, finitediff(z -> $M.$f(z, bar), float(foo)), rtol=0.05)
                end
                if !(isnan(dy))
                    @test isapprox(dy, finitediff(z -> $M.$f(foo, z), bar), rtol=0.05)
                end
            end
        end
    elseif arity == 3
        #=
        @test DiffRules.hasdiffrule(M, f, 3)
        derivs = DiffRules.diffrule(M, f, :foo, :bar, :goo)
        @eval begin
            foo, bar, goo = randn(3)
            dx, dy, dz = $(derivs[1]), $(derivs[2]), $(derivs[3])
            if !(isnan(dx))
                @test isapprox(dx, finitediff(x -> $M.$f(x, bar, goo), foo), rtol=0.05)
            end
            if !(isnan(dy))
                @test isapprox(dy, finitediff(y -> $M.$f(foo, y, goo), bar), rtol=0.05)
            end
            if !(isnan(dz))
                @test isapprox(dz, finitediff(z -> $M.$f(foo, bar, z), goo), rtol=0.05)
            end
        end
        =#
    end
end

# Treat rem2pi separately because of its non-numeric second argument:
derivs = DiffRules.diffrule(:Base, :rem2pi, :x, :y)
for xtype in [:Float64, :BigFloat, :Int64]
    for mode in [:RoundUp, :RoundDown, :RoundToZero, :RoundNearest]
        @eval begin
            let
                x = $xtype(rand(1 : 10))
                y = $mode
                dx, dy = $(derivs[1]), $(derivs[2])
                @test isapprox(dx, finitediff(z -> rem2pi(z, y), float(x)), rtol=0.05)
                @test isnan(dy)
            end
        end
    end
end
end

    @testset "diffrules" begin
        rules = @test_deprecated(DiffRules.diffrules())
        @test Set(M for (M, _, _) in rules) == Set((:Base, :SpecialFunctions, :NaNMath))

        rules = DiffRules.diffrules(; filter_modules=nothing)
        @test Set(M for (M, _, _) in rules) == Set((:Base, :SpecialFunctions, :NaNMath, :LogExpFunctions))

        rules = DiffRules.diffrules(; filter_modules=(:Base, :LogExpFunctions))
        @test Set(M for (M, _, _) in rules) == Set((:Base, :LogExpFunctions))
    end
end

# Test ifelse separately as first argument is boolean
#=
@test DiffRules.hasdiffrule(:Base, :ifelse, 3)
derivs = DiffRules.diffrule(:Base, :ifelse, :foo, :bar, :goo)
for cond in [true, false]
    @eval begin
        foo = $cond
        bar, gee = randn(2)
        dx, dy, dz = $(derivs[1]), $(derivs[2]), $(derivs[3])
        @test isapprox(dy, finitediff(y -> ifelse(foo, y, goo), bar), rtol=0.05)
        @test isapprox(dz, finitediff(z -> ifelse(foo, bar, z), goo), rtol=0.05)
    end
end
=#
