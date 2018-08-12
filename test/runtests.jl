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

@testset "DiffRules" begin
    # include("rules.jl")
    include("forward.jl")
    # include("reverse.jl")
end
