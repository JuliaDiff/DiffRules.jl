using DiffRules: DEFINED_REVERSE_RULES, arity

# for (key, rule_dict) in DEFINED_REVERSE_RULES
#     @show key, arity(key)
#     if arity(key) == 4
#         @show rule_dict[(1,)](:z, :z̄, :g, :h)
#         @show rule_dict[(2,)](:z, :z̄, :g, :h)
#     else
#         @show rule_dict[(1,)](:z, :z̄, :g)
#     end
# end

using DiffRules: @forward_rule, DEFINED_REVERSE_RULES, arity, diffrule

non_numeric_arg_functions = [(:Base, :rem2pi, 4)]

# Check that all forward rules agree with basic diff rules.
for (key, rules) in DEFINED_REVERSE_RULES
    M, f = key[1], key[2]
    (M, f, arity(key)) ∈ non_numeric_arg_functions && continue
    if arity(key) == 3
        modifier = f ∈ (:asec, :acsc, :asecd, :acscd, :acosh, :acoth) ? 1 : 0
        @eval manual_rule = (z, z̄, g)->z̄ * $(diffrule(M, f, :g))
        @eval reverse_rule = (z, z̄, g)->$(rules[(1,)](:z, :z̄, :g))
        x, z, z̄ = rand() + modifier, randn(), randn()
        @test manual_rule(z, z̄, x) ≈ reverse_rule(z, z̄, x)
    elseif arity(key) == 4
        @eval manual_rule_1 = (z, z̄, g_x, h_z)->z̄ * $(diffrule(M, f, :g_x, :h_z)[1])
        @eval manual_rule_2 = (z, z̄, g_x, h_z)->z̄ * $(diffrule(M, f, :g_x, :h_z)[2])
        @eval reverse_rule_1 = (z, z̄, g_x, h_z)->$(rules[(1,)](:z, :z̄, :g_x, :h_z))
        @eval reverse_rule_2 = (z, z̄, g_x, h_z)->$(rules[(2,)](:z, :z̄, :g_x, :h_z))
        x, y, z, z̄ = rand(), rand(), rand(), rand()
        @test manual_rule_1(z, z̄, x, y) ≈ reverse_rule_1(z, z̄, x, y)
        @test manual_rule_2(z, z̄, x, y) ≈ reverse_rule_2(z, z̄, x, y)
    else
        @test 1 === 0
    end
end

# if arity(rule) == 3
#     modifier = f ∈ (:asec, :acsc, :asecd, :acscd, :acosh, :acoth) ? 1 : 0
#     simple_rule_code = diffrule(M, f, :g)
#     reverse_rule_code = rule[2](:x, :ẋ)
#     @eval manual_rule = (x, ẋ)->ẋ * $simple_rule_code
#     @eval forward_rule = (x, ẋ)->$forward_rule_code
#     x, ẋ = rand() + modifier, randn()
#     @test manual_rule(x, ẋ) ≈ forward_rule(x, ẋ)
# elseif arity(rule) == 4
#     ∂f∂x, ∂f∂y = diffrule(M, f, :x, :y)
#     forward_rule_code = rule[2](:x, :y, :ẋ, :ẏ)
#     @eval manual_rule = (x, y, ẋ, ẏ)->ẋ * $∂f∂x + ẏ * $∂f∂y
#     @eval forward_rule = (x, y, ẋ, ẏ)->$forward_rule_code
#     x, y, ẋ, ẏ = rand(), rand(), rand(), rand()
#     manual, fwd = manual_rule(x, y, ẋ, ẏ), forward_rule(x, y, ẋ, ẏ)
#     @test isnan(manual) && isnan(fwd) || manual ≈ fwd
# else

