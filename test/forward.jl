using DiffRules: @forward_rule, DEFINED_FORWARD_RULES, arity, diffrule

non_numeric_arg_functions = [(:Base, :rem2pi, 4)]

# Check that all forward rules agree with basic diff rules.
for ((M, f), rules) in DEFINED_FORWARD_RULES
    for rule in rules
        (M, f, arity(rule)) ∈ non_numeric_arg_functions && continue
        if arity(rule) == 2
            modifier = f ∈ (:asec, :acsc, :asecd, :acscd, :acosh, :acoth) ? 1 : 0
            simple_rule_code = diffrule(M, f, :x)
            forward_rule_code = rule[2](:x, :ẋ)
            @eval manual_rule = (x, ẋ)->ẋ * $simple_rule_code
            @eval forward_rule = (x, ẋ)->$forward_rule_code
            x, ẋ = rand() + modifier, randn()
            @test manual_rule(x, ẋ) ≈ forward_rule(x, ẋ)
        elseif arity(rule) == 4
            ∂f∂x, ∂f∂y = diffrule(M, f, :x, :y)
            forward_rule_code = rule[2](:x, :y, :ẋ, :ẏ)
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
