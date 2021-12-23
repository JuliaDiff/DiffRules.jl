using DiffRules
using FiniteDifferences
using Test

import SpecialFunctions, NaNMath, LogExpFunctions
import Random
Random.seed!(1)

const finitediff = central_fdm(5, 1)

@testset "DiffRules" begin
@testset "check rules" begin

non_diffeable_arg_functions = [(:Base, :rem2pi, 2), (:Base, :ldexp, 2), (:Base, :ifelse, 3)]

for (M, f, arity) in DiffRules.diffrules(; filter_modules=nothing)
    (M, f, arity) ∈ non_diffeable_arg_functions && continue
    if arity == 1
        @test DiffRules.hasdiffrule(M, f, 1)
        deriv = DiffRules.diffrule(M, f, :goo)
        @eval begin
            let
                goo = if $(f in (:asec, :acsc, :asecd, :acscd, :acosh, :acoth))
                    # avoid singularities with finite differencing
                    rand() + 1.5
                elseif $(f === :log)
                    # avoid singularities with finite differencing
                    rand() + 0.5
                elseif $(f === :log1mexp)
                    rand() - 1.0
                elseif $(f in (:log2mexp, :asin, :acos, :erfinv))
                    rand() - 0.5
                else
                    rand()
                end
                @test $deriv ≈ finitediff($M.$f, goo) rtol=1e-9 atol=1e-9
                # test for 2pi functions
                if $(f === :mod2pi)
                    goo = 4 * pi
                    @test NaN === $deriv
                end
            end
        end
    elseif arity == 2
        @test DiffRules.hasdiffrule(M, f, 2)
        derivs = DiffRules.diffrule(M, f, :foo, :bar)
        @eval begin
            let
                foo, bar = if $(f === :mod)
                    rand() + 13, rand() + 5 # make sure x/y is not integer
                elseif $(f === :polygamma)
                    rand(1:10), rand() # only supports integers as first arguments
                elseif $(f === :bessely)
                    # avoid singularities with finite differencing
                    rand(), rand() + 0.5
                elseif $(f === :log)
                    # avoid singularities with finite differencing
                    rand() + 1.5, rand()
                elseif $(f === :^)
                    # avoid singularities with finite differencing
                    rand() + 0.5, rand()
                else
                    rand(), rand()
                end
                dx, dy = $(derivs[1]), $(derivs[2])
                if !isnan(dx)
                    @test dx ≈ finitediff(z -> $M.$f(z, bar), float(foo)) rtol=1e-9 atol=1e-9
                end
                if !isnan(dy)
                    @test dy ≈ finitediff(z -> $M.$f(foo, z), bar) rtol=1e-9 atol=1e-9
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
                @test dx ≈ finitediff(z -> rem2pi(z, y), float(x)) rtol=1e-9 atol=1e-9
                @test isnan(dy)
            end
        end
    end
end

# Treat ldexp separately because of its integer second argument:
derivs = DiffRules.diffrule(:Base, :ldexp, :x, :y)
for xtype in [:Float64, :BigFloat]
    for ytype in [:Integer, :UInt64, :Int64]
        @eval begin
            let
                x = rand($xtype)
                y = $ytype(rand(1 : 10))
                dx, dy = $(derivs[1]), $(derivs[2])
                @test dx ≈ finitediff(z -> ldexp(z, y), x) rtol=1e-9 atol=1e-9
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
