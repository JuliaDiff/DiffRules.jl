# FORMAT of rule function: (y, ȳ, x₁, x₂, ...) where `y` is output, `ȳ` is sensitivity of
# output, `x₁, x₂, ...` are inputs. If sensitivity computation doesn't require one of the
# arguments, it simply won't be used in the resulting expression.

# A (<module_name>, <function_name>, <argument_signature>, <argument_numbers>)-Tuple.
const ReverseRuleKey = Tuple{SymOrExpr, Symbol, Any, Tuple{Vararg{Int}}}

# All of the defined reverse rules. Keys are of the form:
# (Module, Function, type-tuple, argument numbers)
const DEFINED_REVERSE_RULES = Dict{ReverseRuleKey, Any}()

macro reverse_rule(def::Expr)
    def_ = splitdef(def)
    @assert def_[:whereparams] == () "where parameters not currently supported"
    @assert isempty(def_[:kwargs]) "There exists a keyword argument"

    M, f = _split_qualified_name(def_[:name])

    wrts, args′ = process_args(def_[:args])
    args = splitarg.(def_[:args])
    @assert all(arg->arg[3] === false, args) "At least one argument is slurping"
    @assert all(arg->arg[4] === nothing, args) "At least one argument has a default value"
    signature = :(Tuple{$(getindex.(args, 2)...)})

    body = Expr(:->, Expr(:tuple, args′...), def_[:body])
    return esc(:(add_reverse_rule!(($M, $f, $signature, $wrts), $body)))
end

function process_args(args::Array{Any})
    wrts, args′ = Vector{Int}(), Vector{Any}(undef, length(args))
    for (n, arg) in enumerate(args)
        if arg isa Expr && arg.head == :call && arg.args[1] == :wrt
            @assert length(arg.args) == 2
            push!(wrts, n)
            args′[n] = arg.args[2]
        else
            args′[n] = arg
        end
    end
    return (wrts...,), args′
end

function add_reverse_rule!(key::ReverseRuleKey, rule::Any)
    DEFINED_REVERSE_RULES[key] = rule
end

arity(key::ReverseRuleKey) = length(getfield(key[3], :3))

# Create forward rules from all of the existing diff rules.
for ((M, f, nargs), rules) in DEFINED_DIFFRULES
    if nargs == 1
        reverse_rule = (z::Symbol, z̄::Symbol, x::Symbol)->:($z̄ * $(rules(x)))
        add_reverse_rule!((M, f, Tuple{Real, Real, Real}, (1,)), reverse_rule)
    elseif nargs == 2
        ∂f∂x, ∂f∂y = rules(:x, :y)
        (∂f∂x == :NaN || ∂f∂y == :NaN) && continue
        rev_rule_1 = (z::Symbol, z̄::Symbol, x::Symbol, y::Symbol)->:(z̄ * $(rules(x, y)[1]))
        rev_rule_2 = (z::Symbol, z̄::Symbol, x::Symbol, y::Symbol)->:(z̄ * $(rules(x, y)[2]))
        add_reverse_rule!((M, f, Tuple{Vararg{Real, 4}}, (1,)), rev_rule_1)
        add_reverse_rule!((M, f, Tuple{Vararg{Real, 4}}, (2,)), rev_rule_2)
    else
        error("Arrghh")
    end
end
