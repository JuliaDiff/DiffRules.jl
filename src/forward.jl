# Key indicates the function in terms of (module_name, function_name). Each entry contains
# a vector of implementations for different type signatures. First element of entry is
# method signature, second entry is a (Tuple of) expressions giving the forwards-mode
# sensitivities.
const DEFINED_FORWARD_RULES = Dict{Tuple{SymOrExpr, Symbol}, Vector{Tuple{Any, Any}}}()

macro forward_rule(def::Expr)

    # Split up function definition and assert no whereparams or kwargs.
    def_ = splitdef(def)
    @assert def_[:whereparams] == () "where parameters not currently supported"
    @assert isempty(def_[:kwargs]) "There exists a keyword argument"

    # Split up the arguments and assert no is slurps or default values.
    args = splitarg.(def_[:args])
    @assert all(arg->arg[3] === false, args) "At least one argument is slurping"
    @assert all(arg->arg[4] === nothing, args) "At least one argument has a default value"

    # Construct forward rule.
    spec = _split_qualified_name(def_[:name])
    signature = :(Tuple{$(getindex.(args, 2)...)})
    expr = Expr(:->, Expr(:tuple, def_[:args]...), def_[:body])
    return esc(:(DiffRules.add_forward_rule($spec, $signature, $expr)))
end

function add_forward_rule(key::Tuple{SymOrExpr, Symbol}, signature::DataType, body)
    if key ∉ keys(DEFINED_FORWARD_RULES)
        DEFINED_FORWARD_RULES[key] = Vector{Tuple{Any, Any}}()
    end
    push!(DEFINED_FORWARD_RULES[key], (signature, body))
end

# Compute the arity of a particular method from `DEFINED_FORWARD_RULES`.
arity(meth::Tuple{Vararg{Any, 2}}) = length(getfield(meth[1], :3))

# Create forward rules from all of the existing diff rules.
for ((M, f, nargs), rules) in DEFINED_DIFFRULES
    if nargs == 1
        add_forward_rule(
            (M, f),
            Tuple{Vararg{Real, 2}},
            (x::Symbol, ẋ::Symbol)->:($ẋ * $(rules(x))),
        )
    elseif nargs == 2
        ∂f∂x, ∂f∂y = rules(:x, :y)
        (∂f∂x == :NaN || ∂f∂y == :NaN) && continue
        add_forward_rule(
            (M, f),
            Tuple{Vararg{Real, 4}},
            (x::Symbol, y::Symbol, ẋ::Symbol, ẏ::Symbol)->:($ẋ * $∂f∂x + $ẏ * $∂f∂y),
        )
    else
        error("Arrghh")
    end
end
