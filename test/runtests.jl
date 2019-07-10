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


epsreal(x) = cbrt(eps(typeof(real(x)))) * max(one(typeof(real(x))), abs(x))
finitediff(f, x, ϵ) = (f(x + ϵ) - f(x - ϵ)) / (ϵ + ϵ)
finitediff(f, x) = finitediff(f, x, epsreal(x))

non_numeric_arg_functions = [(:Base, :rem2pi, 2)]

# just to satisfy Coverall which is not taking into account macros calls 
DiffRules._getkeyrule(Meta.parse("Base.:+(x) = :(1)"))

print("complex diffrules:\n")
for (M, f, arity) in DiffRules.complex_diffrules()
    (M, f, arity) ∈ non_numeric_arg_functions && continue
    print(M, ".", f, "\n")
    if arity == 1
        @test DiffRules.hascomplex_diffrule(M, f, 1)
        deriv = DiffRules.complex_diffrule(M, f, :goo)
        @eval begin
            for goo in [0.5+0im, 0.5+0.5im, 0.5im, -0.5+0.5im, -0.5+0im, -0.5-0.5im, -0.5im, 0.5-0.5im]
                #
                $M.$f == Base.sqrt  && imag(goo) == 0 && real(goo) <= 0 && continue
                #
                $M.$f == Base.log   && imag(goo) == 0 && real(goo) <= 0 && continue
                $M.$f == Base.log10 && imag(goo) == 0 && real(goo) <= 0 && continue
                $M.$f == Base.log2  && imag(goo) == 0 && real(goo) <= 0 && continue
                #
                $M.$f == Base.asin  && imag(goo) == 0 && real(goo) <= -1 && continue
                $M.$f == Base.asin  && imag(goo) == 0 && real(goo) >=  1 && continue
                $M.$f == Base.acos  && imag(goo) == 0 && real(goo) <= -1 && continue
                $M.$f == Base.acos  && imag(goo) == 0 && real(goo) >=  1 && continue
                $M.$f == Base.atan  && real(goo) == 0 && imag(goo) <= -1 && continue
                $M.$f == Base.atan  && real(goo) == 0 && imag(goo) >=  1 && continue
                $M.$f == Base.asec  && imag(goo) == 0 && abs(real(goo)) <= 1 && continue
                $M.$f == Base.acsc  && imag(goo) == 0 && abs(real(goo)) <= 1 && continue
                $M.$f == Base.acot  && real(goo) == 0 && abs(imag(goo)) >= 1 && continue
                $M.$f == Base.asinh && real(goo) == 0 && abs(imag(goo)) <= 1 && continue
                $M.$f == Base.acosh && imag(goo) == 0 && abs(real(goo)) <= 1 && continue
                $M.$f == Base.acoth && imag(goo) == 0 && abs(real(goo)) <= 1 && continue
                if ($M.$f == Base.acoth && real(goo) < 0 && abs(coth(acoth(goo)) - goo) > 1e-10)
                    println("Skipping bogous acoth function: coth(acoth($goo))  = ", coth(acoth(goo)))
                    continue
                end
                $M.$f == Base.acsch && real(goo) == 0 && continue
                $M.$f == Base.asech && imag(goo) == 0 && real(goo) < 0 && continue
                $M.$f == Base.atanh && imag(goo) == 0 && abs(real(goo)) >= 1.0 && continue
                if ($M.$f == Base.atanh && real(goo) < 0 && abs(tanh(atanh(goo)) - goo) > 1e-10)
                    println("Skipping bogous atanh function: tahh(atahh($goo))  = ", tanh(atanh(goo)))
                    continue
                end
                $M.$f == SpecialFunctions.bessely0 && imag(goo) == 0 && real(goo) < 0 && continue
                $M.$f == SpecialFunctions.bessely1 && imag(goo) == 0 && real(goo) < 0 && continue
                $M.$f == SpecialFunctions.erfc     && imag(goo) == 0 && real(goo) < 0 && continue
                $M.$f == SpecialFunctions.lgamma   && imag(goo) == 0 && real(goo) < 0 && continue
                ϵ = 1e-10
                for phi in 0:11
                    #print(goo, " ", phi, " ", $deriv, " ", finitediff($M.$f, goo, ϵ*exp(1im*pi*phi/6.)), "\n")
                    @test isapprox($deriv, finitediff($M.$f, goo, ϵ*exp(1im*pi*phi/6.)), atol=0.001)
                end
            end
        end
    elseif arity == 2
        @test DiffRules.hascomplex_diffrule(M, f, 2)
        derivs = DiffRules.complex_diffrule(M, f, :foo, :bar)
        @eval begin
            for foo in [0.5+0im, 0.5+0.5im, 0.5im, -0.5+0.5im, -0.5+0im, -0.5-0.5im, -0.5im, 0.5-0.5im]
                for bar in [0.5+0im, 0.5+0.5im, 0.5im, -0.5+0.5im, -0.5+0im, -0.5-0.5im, -0.5im, 0.5-0.5im]
                    $M.$f == "Base.^"  && imag(foo) == 0 && real(foo) < 0 && continue
                    dx, dy = $(derivs[1]), $(derivs[2])
                    if !(isnan(dx))
                        ϵ = 1e-10
                        for phi in 0:11
                            # print(foo, " ", bar, " ", phi, " ", dx, " ", finitediff(z -> $M.$f(z, bar), foo, ϵ*exp(1im*pi*phi/6.)), "\n")
                            @test isapprox(dx, finitediff(z-> $M.$f(z, bar), foo, ϵ*exp(1im*pi*phi/6.)), atol=0.001)
                        end
                    end
                    if !(isnan(dy))
                        ϵ = 1e-10
                        for phi in 0:11
                            # print(foo, " ", bar, " ", phi, " ", dy, " ", finitediff(z -> $M.$f(foo, z), bar, ϵ*exp(1im*pi*phi/6.)), "\n")
                            @test isapprox(dy, finitediff(z-> $M.$f(foo, z), bar, ϵ*exp(1im*pi*phi/6.)), atol=0.001)
                        end
                    end
                end
            end
        end
    end
end

print("\nnon complex diffrules:\n")
for (M, f, arity) in DiffRules.diffrules()
    (M, f, arity) ∈ non_numeric_arg_functions && continue
    print(M, ".", f, "\n")
    if arity == 1
        @test DiffRules.hasdiffrule(M, f, 1)
        deriv = DiffRules.diffrule(M, f, :goo)
        modifier = in(f, (:asec, :acsc, :asecd, :acscd, :acosh, :acoth)) ? 1 : 0
        @eval begin
            goo = rand() + $modifier
            @test isapprox($deriv, finitediff($M.$f, goo), rtol=0.05)
        end
    elseif arity == 2
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

