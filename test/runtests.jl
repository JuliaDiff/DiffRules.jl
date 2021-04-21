if VERSION < v"0.7-"
    using Base.Test
    srand(1)
else
    using Test
    import Random
    Random.seed!(1)
end
import SpecialFunctions, NaNMath
using DiffRules


function finitediff(f, x)
    ϵ = cbrt(eps(typeof(x))) * max(one(typeof(x)), abs(x))
    return (f(x + ϵ) - f(x - ϵ)) / (ϵ + ϵ)
end


non_numeric_arg_functions = [(:Base, :rem2pi, 2)]

for T in [Float32, Float64]
    for (M, f, arity) in DiffRules.diffrules()
        @testset "$M.$f $(arity)" begin
            (M, f, arity) ∈ non_numeric_arg_functions && continue
            if arity == 1
                @test DiffRules.hasdiffrule(M, f, 1)
                deriv = DiffRules.diffrule(M, f, :goo)
                modifier = in(f, (:asec, :acsc, :asecd, :acscd, :acosh, :acoth)) ? 1 : 0
                @eval begin
                    goo = $T(rand() + $modifier)
                    @test $deriv isa $T
                    @test isapprox($deriv, finitediff($M.$f, goo), rtol=0.05)
                    # test for 2pi functions
                    if "mod2pi" == string($M.$f)
                        goo = 4pi + $modifier
                        @test NaN === $deriv
                    end
                end
            elseif arity == 2
                # TODO: Add test for types.
                @test DiffRules.hasdiffrule(M, f, 2)
                derivs = DiffRules.diffrule(M, f, :foo, :bar)
                @eval begin
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
        end
    end
end

# Treat rem2pi separately because of its non-numeric second argument:
derivs = DiffRules.diffrule(:Base, :rem2pi, :x, :y)
for xtype in [:Float64, :BigFloat, :Int64]
    for mode in [:RoundUp, :RoundDown, :RoundToZero, :RoundNearest]
        @eval begin
            x = $xtype(rand(1 : 10))
            y = $mode
            dx, dy = $(derivs[1]), $(derivs[2])
            @test isapprox(dx, finitediff(z -> rem2pi(z, y), float(x)), rtol=0.05)
            @test isnan(dy)
        end
    end
end
