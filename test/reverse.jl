using DiffRules: DEFINED_REVERSE_RULES, arity, diffrule, _reverse_rule

@testset "reverse" begin

    # Check that various things that should fail, fail.
    @test_throws AssertionError _reverse_rule(:(M.f(x)))
    @test_throws AssertionError _reverse_rule(:(f(x::T) where T = :(5x)))
    @test_throws AssertionError _reverse_rule(:(f(x::Real; y) = :(5x)))
    @test_throws AssertionError _reverse_rule(:(f(x::Real...) = :(5x)))
    @test_throws AssertionError _reverse_rule(:(f(x::Real=5) = :(5x)))
    @test_throws ErrorException _reverse_rule(:(f(x::Real) = :(5x)))

    # Check that a basic rule works.
    foo(x) = 5x
    @forward_rule Main.foo(x, ẋ) = :(5($x)^2 + $ẋ)
    @test DEFINED_FORWARD_RULES[(:Main, :foo, :(Tuple{Any, Any}))](:g, :h) == :(5g^2 + h)
    delete!(DEFINED_FORWARD_RULES, (:Main, :foo, :(Tuple{Any, Any})))


    non_numeric_arg_functions = [(:Base, :rem2pi, 4)]

    # Check that all reverse rules agree with basic diff rules.
    for (key, rule) in DEFINED_REVERSE_RULES
        M, f = key[1], key[2]
        (M, f, arity(key)) ∈ non_numeric_arg_functions && continue
        if arity(key) == 3
            modifier = f ∈ (:asec, :acsc, :asecd, :acscd, :acosh, :acoth) ? 1 : 0
            @eval manual_rule = (z, z̄, g)->z̄ * $(diffrule(M, f, :g))
            @eval reverse_rule = (z, z̄, g)->$(rule(:z, :z̄, :g))
            x, z, z̄ = rand() + modifier, randn(), randn()
            @test manual_rule(z, z̄, x) ≈ reverse_rule(z, z̄, x)
        elseif arity(key) == 4
            @eval manual_rule = (z, z̄, g_x, h_z)->z̄ * $(diffrule(M, f, :g_x, :h_z)[key[4][1]])
            @eval reverse_rule = (z, z̄, g_x, h_z)->$(rule(:z, :z̄, :g_x, :h_z))
            x, y, z, z̄ = rand(), rand(), rand(), rand()
            @test manual_rule(z, z̄, x, y) ≈ reverse_rule(z, z̄, x, y)
        else
            @test 1 === 0
        end
    end

end
