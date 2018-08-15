using Test, DiffRules, Random, FDM
import Random, SpecialFunctions, NaNMath

Random.seed!(1)

function finitediff(f, x)
    ϵ = cbrt(eps(typeof(x))) * max(one(typeof(x)), abs(x))
    return (f(x + ϵ) - f(x - ϵ)) / (ϵ + ϵ)
end

@testset "DiffRules" begin

    @testset "diffrules" begin
        include("diffrules/rules.jl")
    end

    @testset "forward" begin
        include("forward/api.jl")
    end

    @testset "reverse" begin
        include("reverse/api.jl")
        include("reverse/test_util.jl")
        include("reverse/generic.jl")
    end
end
