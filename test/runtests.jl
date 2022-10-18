using DiffRules
using Test
using FiniteDifferences

using IrrationalConstants: fourπ

import SpecialFunctions, NaNMath, LogExpFunctions
import Random
Random.seed!(1)

# Set `max_range` to avoid domain errors.
const finitediff = central_fdm(5, 1, max_range=1e-3)

@testset "DiffRules" begin
@testset "check rules" begin

non_diffeable_arg_functions = [(:Base, :rem2pi, 2), (:Base, :ldexp, 2), (:Base, :ifelse, 3)]

@testset "($M, $f, $arity)" for (M, f, arity) in DiffRules.diffrules(; filter_modules=nothing)
    for T in [Float32, Float64]
        (M, f, arity) ∈ non_diffeable_arg_functions && continue
        if arity == 1
            @test DiffRules.hasdiffrule(M, f, 1)
            deriv = DiffRules.diffrule(M, f, :goo)
            @eval begin
                let
                    goo = if $(f in (:asec, :acsc, :asecd, :acscd, :acosh, :acoth))
                        # avoid singularities with finite differencing
                        rand($T) + $T(1.5)
                    elseif $(f in (:log, :airyaix, :airyaiprimex, :logmxp1))
                        # avoid singularities with finite differencing
                        rand($T) + $T(0.5)
                    elseif $(f === :log1mexp)
                        rand($T) - one($T)
                    elseif $(f in (:log2mexp, :erfinv))
                        rand($T) - $T(0.5)
                    else
                        rand($T)
                    end
                    # We're happy with types with the correct promotion behavior, e.g.
                    # it's fine to return `1` as a derivative despite input being `Float64`.
                    @test promote_type(typeof($deriv), $T) === $T
                    # In older versions of LogExpFunctions `log1pmx(::Float32)` and `logmxp1(::Float32)` are not defined
                    if hasmethod($M.$f, Tuple{$T})
                        @test $deriv ≈ finitediff($M.$f, goo) rtol=1e-2 atol=1e-3     
                    end
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
                        rand($T) + 13, rand($T) + 5 # make sure x/y is not integer
                    elseif $(f === :polygamma)
                        rand(1:10), rand($T) # only supports integers as first arguments
                    elseif $(f in (:bessely, :besselyx))
                        # avoid singularities with finite differencing
                        rand($T), rand($T) + $T(0.5)
                    elseif $(f === :log)
                        # avoid singularities with finite differencing
                        rand($T) + $T(1.5), rand($T)
                    elseif $(f === :^)
                        # avoid singularities with finite differencing
                        rand($T) + $T(0.5), rand($T)
                    else
                        rand($T), rand($T)
                    end
                    dx, dy = $(derivs[1]), $(derivs[2])
                    if !(isnan(dx))
                        @test dx ≈ finitediff(z -> $M.$f(z, bar), foo) rtol=1e-2 atol=1e-3

                        # Check type, if applicable.
                        @test promote_type(typeof(real(dx)), $T) === $T
                    end
                    if !(isnan(dy))
                        @test dy ≈ finitediff(z -> $M.$f(foo, z), bar) rtol=1e-2 atol=1e-3

                        # Check type, if applicable.
                        @test promote_type(typeof(real(dy)), $T) === $T
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

# Check negative branch for `airybix` and `airybiprimex`
for f in (:airybix, :airybiprimex)
    deriv = DiffRules.diffrule(:SpecialFunctions, f, :goo)
    @eval begin
        let
            goo = -rand()
            @test $deriv ≈ finitediff(SpecialFunctions.$f, goo) rtol=1e-9 atol=1e-9
        end
    end
end

# Test `iszero(x)` branch of `xlogy`
derivs = DiffRules.diffrule(:LogExpFunctions, :xlogy, :x, :y)
for xytype in [:Float32, :Float64, :BigFloat]
    @eval begin
        let
            x = zero($xytype)
            y = rand($xytype)
            dx, dy = $(derivs[1]), $(derivs[2])
            @test iszero(dy)

            y = one($xytype)
            dx, dy = $(derivs[1]), $(derivs[2])
            @test iszero(dy)
        end
    end
end

# Test `iszero(x)` branch of `xlog1py`
derivs = DiffRules.diffrule(:LogExpFunctions, :xlog1py, :x, :y)
for xytype in [:Float32, :Float64, :BigFloat]
    @eval begin
        let
            x = zero($xytype)
            y = rand($xytype)
            dx, dy = $(derivs[1]), $(derivs[2])
            @test iszero(dy)

            y = -one($xytype)
            dx, dy = $(derivs[1]), $(derivs[2])
            @test iszero(dy)
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
