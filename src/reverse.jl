# FORMAT of rule function: (y, ȳ, x₁, x₂, ...) where `y` is output, `ȳ` is sensitivity of
# output, `x₁, x₂, ...` are inputs. If sensitivity computation doesn't require one of the
# arguments, it simply won't be used in the resulting expression.

#= A particular method will potentially require a number of different reverse-mode rules.
   Each key is a Tuple of Ints, corresponding to the positions in the signature of the
   arguments w.r.t. which the corresponding value defines sensitivities. For example,
   `d[1]` should return an expression for the sensitivity w.r.t. the first argument of the
   method, and `d[2, 3]` an expression providing the sensitivity w.r.t. the 2nd and 3rd
   arguments respectively. Clearly, not all combinations will typically be provided.
=# 
const ReverseRuleDict = Dict{Tuple{Vararg{Int}}, Any}

# All of the defined reverse rules. Keys are of the form (Module, Function, type-tuple)
const DEFINED_REVERSE_RULES = Dict{Tuple{SymOrExpr, Symbol, Any}, ReverseRuleDict}()

# macro reverse_rule(def::Expr)

#     # Split up function definition and assert no whereparams or kwargs.
#     def_ = splitdef(def)
#     @assert def_[:whereparams] == () "where parameters not currently supported"
#     @assert isempty(def_[:kwargs]) "There exists a keyword argument"

#     # Split up the arguments and assert no is slurps or default values.
#     args = splitarg.(def_[:args])
#     @assert all(arg->arg[3] === false, args) "At least one argument is slurping"
#     @assert all(arg->arg[4] === nothing, args) "At least one argument has a default value"

#     # Construct forward rule.
#     spec = _split_qualified_name(def_[:name])
#     signature = :(Tuple{$(getindex.(args, 2)...)})
#     expr = Expr(:->, Expr(:tuple, def_[:args]...), def_[:body])
#     return esc(:(DiffRules.add_reverse_rule($spec, $signature, $expr)))
# end

"""
    add_reverse_rule!(
        M::SymOrExpr,
        f::Symbol,
        signature::DataType,
        positions::Tuple{Vararg{Int}},
        body::Tuple,
    )

Adds a reverse rule for the method of function `f` in module `M` with signature `signature`
given by `body`, which provides sensitivities w.r.t. arguments with positions `positions` in
the signature.
"""
function add_reverse_rule!(
    M::SymOrExpr,
    f::Symbol,
    signature::DataType,
    positions::Tuple{Vararg{Int}},
    body::Any,
)
    key = (M, f, signature)
    if key ∉ keys(DEFINED_REVERSE_RULES)
        DEFINED_REVERSE_RULES[key] = ReverseRuleDict()
    end
    DEFINED_REVERSE_RULES[key][positions] = body
end

arity(meth::Tuple{SymOrExpr, Symbol, DataType}) = length(getfield(meth[3], :3))

# Create forward rules from all of the existing diff rules.
for ((M, f, nargs), rules) in DEFINED_DIFFRULES
    if nargs == 1
        reverse_rule = (z::Symbol, z̄::Symbol, x::Symbol)->:($z̄ * $(rules(x)))
        add_reverse_rule!(M, f, Tuple{Real, Real, Real}, (1,), reverse_rule)
    elseif nargs == 2
        ∂f∂x, ∂f∂y = rules(:x, :y)
        (∂f∂x == :NaN || ∂f∂y == :NaN) && continue
        rev_rule_1 = (z::Symbol, z̄::Symbol, x::Symbol, y::Symbol)->:(z̄ * $(rules(x, y)[1]))
        rev_rule_2 = (z::Symbol, z̄::Symbol, x::Symbol, y::Symbol)->:(z̄ * $(rules(x, y)[2]))
        @show rev_rule_1(:z, :z̄, :g, :h)
        @show rev_rule_2(:z, :z̄, :g, :h)
        add_reverse_rule!(M, f, Tuple{Vararg{Real, 4}}, (1,), rev_rule_1)
        add_reverse_rule!(M, f, Tuple{Vararg{Real, 4}}, (2,), rev_rule_2)
    else
        error("Arrghh")
    end
end
