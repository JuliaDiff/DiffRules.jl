using DiffRules: diffrules, DEFINED_REVERSE_RULES, arity, diffrule, @reverse_rule,
    _reverse_rule, ReverseRuleKey

@testset "api" begin

    # Check that various things that should fail, fail.
    @test_throws AssertionError _reverse_rule(:z, :z̄, :(M.f(x)))
    @test_throws AssertionError _reverse_rule(:z, :z̄, :(M.f(x::T) where T = :(5x)))
    @test_throws AssertionError _reverse_rule(:z, :z̄, :(M.f(x::Real; y) = :(5x)))
    @test_throws AssertionError _reverse_rule(:z, :z̄, :(M.f(x::Real...) = :(5x)))
    @test_throws AssertionError _reverse_rule(:z, :z̄, :(M.f(x::Real=5) = :(5x)))
    @test_throws ErrorException _reverse_rule(:z, :z̄, :(f(x::Real) = :(5x)))

    # Check that a basic rule works.
    foo(x) = 5x
    @reverse_rule y ȳ::Real Main.foo(wrt(x::Int)) = :($y * $ȳ * $x)
    key = (:Main, :foo, :(Tuple{Any, Real, Int}), (1,))
    @test DEFINED_REVERSE_RULES[key](:g, :h, :y) == :(g * h * y)
    delete!(DEFINED_REVERSE_RULES, key)


    # non_numeric_arg_functions = [(:Base, :rem2pi, 4)]

    # # Check that all reverse rules agree with basic diff rules.
    # for key in diffrules()
    #     M, f, arity = key
    #     key ∈ non_numeric_arg_functions && continue
    #     if arity == 1
    #         rev_key = ReverseRuleKey(M, f, :(Tuple{Real, Real, Real}), (1,))
    #         rev_rule = DEFINED_REVERSE_RULES[rev_key]
    #         modifier = f ∈ (:asec, :acsc, :asecd, :acscd, :acosh, :acoth) ? 1 : 0
    #         @eval manual_rule = (z, z̄, g)->z̄ * $(diffrule(M, f, :g))
    #         @eval reverse_rule = (z, z̄, g)->$(rev_rule(:z, :z̄, :g))
    #         x, z, z̄ = rand() + modifier, randn(), randn()
    #         @test manual_rule(z, z̄, x) ≈ reverse_rule(z, z̄, x)
    #     elseif arity == 2

    #         ∂f∂x, ∂f∂y = diffrule(M, f, :g_x, :h_z)

    #         # Grab the corresponding reverse rules.
    #         typ = :(Tuple{Real, Real, Real, Real})
    #         key1, key2 = (M, f, typ, (1,)), (M, f, typ, (2,))
    #         if key1 ∈ keys(DEFINED_REVERSE_RULES)
    #             x, y, z, z̄ = rand(), rand(), rand(), rand()
    #             rev_rule_1 = DEFINED_REVERSE_RULES[key1](:z, :z̄, :g_x, :h_z)
    #             @eval manual_rule_1 = (z, z̄, g_x, h_z)->z̄ * $∂f∂x
    #             @eval reverse_rule_1 = (z, z̄, g_x, h_z)->$rev_rule_1
    #             @test manual_rule_1(z, z̄, x, y) ≈ reverse_rule_1(z, z̄, x, y)
    #         end
    #         if key2 ∈ keys(DEFINED_REVERSE_RULES)
    #             x, y, z, z̄ = rand(), rand(), rand(), rand()
    #             rev_rule_2 = DEFINED_REVERSE_RULES[(M, f, typ, (2,))](:z, :z̄, :g_x, :h_z)
    #             @eval manual_rule_2 = (z, z̄, g_x, h_z)->z̄ * $∂f∂y
    #             @eval reverse_rule_2 = (z, z̄, g_x, h_z)->$rev_rule_2
    #             @test manual_rule_2(z, z̄, x, y) ≈ reverse_rule_2(z, z̄, x, y)
    #         end
    #     else
    #         @test 1 === 0
    #     end
    # end

end
