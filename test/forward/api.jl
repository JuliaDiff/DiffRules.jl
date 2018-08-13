using DiffRules: @forward_rule, DEFINED_FORWARD_RULES, arity, diffrule, _forward_rule

@testset "api" begin

    # Check that various things that should fail, fail.
    @test_throws AssertionError _forward_rule(:(M.f(x)))
    @test_throws AssertionError _forward_rule(:(f(x::T) where T = :(5x)))
    @test_throws AssertionError _forward_rule(:(f(x::Real; y) = :(5x)))
    @test_throws AssertionError _forward_rule(:(f(x::Real...) = :(5x)))
    @test_throws AssertionError _forward_rule(:(f(x::Real=5) = :(5x)))
    @test_throws ErrorException _forward_rule(:(f(x::Real) = :(5x)))

    # Check that a basic rule works.
    foo(x) = 5x
    @forward_rule Main.foo(x, ẋ) = :(5($x)^2 + $ẋ)
    @test DEFINED_FORWARD_RULES[(:Main, :foo, :(Tuple{Any, Any}))](:g, :h) == :(5g^2 + h)
    delete!(DEFINED_FORWARD_RULES, (:Main, :foo, :(Tuple{Any, Any})))

    non_numeric_arg_functions = [(:Base, :rem2pi, 4)]

    # Check that all forward rules agree with basic diff rules.
    for (key, body) in DEFINED_FORWARD_RULES
        M, f, signature = key
        (M, f, arity(key)) ∈ non_numeric_arg_functions && continue
        if arity(key) == 2
            modifier = f ∈ (:asec, :acsc, :asecd, :acscd, :acosh, :acoth) ? 1 : 0
            simple_rule_code = diffrule(M, f, :x)
            forward_rule_code = body(:x, :ẋ)
            @eval manual_rule = (x, ẋ)->ẋ * $simple_rule_code
            @eval forward_rule = (x, ẋ)->$forward_rule_code
            x, ẋ = rand() + modifier, randn()
            @test manual_rule(x, ẋ) ≈ forward_rule(x, ẋ)
        elseif arity(key) == 4
            ∂f∂x, ∂f∂y = diffrule(M, f, :x, :y)
            forward_rule_code = body(:x, :y, :ẋ, :ẏ)
            @eval manual_rule = (x, y, ẋ, ẏ)->ẋ * $∂f∂x + ẏ * $∂f∂y
            @eval forward_rule = (x, y, ẋ, ẏ)->$forward_rule_code
            x, y, ẋ, ẏ = rand(), rand(), rand(), rand()
            manual, fwd = manual_rule(x, y, ẋ, ẏ), forward_rule(x, y, ẋ, ẏ)
            @test isnan(manual) && isnan(fwd) || manual ≈ fwd
        else
            error("argh")
        end
    end

end
