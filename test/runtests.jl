using Base.Test
using RealInterface
using SpecialFunctions
using DiffRules

srand(1)

function finitediff(f, x)
    ϵ = cbrt(eps(typeof(x))) * max(one(typeof(x)), abs(x))
    return (f(x + ϵ) - f(x - ϵ)) / (ϵ + ϵ)
end

for f in vcat(RealInterface.UNARY_ARITHMETIC, RealInterface.UNARY_MATH, RealInterface.UNARY_SPECIAL_MATH)
    if DiffRules.hasdiffrule(f, 1)
        deriv = DiffRules.diffrule(f, :x)
        modifier = in(f, (:asec, :acsc, :asecd, :acscd, :acosh, :acoth)) ? 1 : 0
        @eval begin
            x = rand() + $modifier
            @test isapprox($deriv, finitediff($f, x), rtol=0.05)
        end
    end
end

for f in vcat(RealInterface.BINARY_MATH, RealInterface.BINARY_ARITHMETIC, RealInterface.BINARY_SPECIAL_MATH)
    if DiffRules.hasdiffrule(f, 2)
        derivs = DiffRules.diffrule(f, :x, :y)
        @eval begin
            x, y = rand(1:10), rand()
            dx, dy = $(derivs[1]), $(derivs[2])
            if !(isnan(dx))
                @test isapprox(dx, finitediff(z -> $f(z, y), float(x)), rtol=0.05)
            end
            if !(isnan(dy))
                @test isapprox(dy, finitediff(z -> $f(x, z), y), rtol=0.05)
            end
        end
    end
end
