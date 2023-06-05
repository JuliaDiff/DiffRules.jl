var documenterSearchIndex = {"docs": [

{
    "location": "index.html#",
    "page": "Documentation",
    "title": "Documentation",
    "category": "page",
    "text": ""
},

{
    "location": "index.html#DiffRules.@define_diffrule",
    "page": "Documentation",
    "title": "DiffRules.@define_diffrule",
    "category": "macro",
    "text": "@define_diffrule M.f(x) = :(df_dx($x))\n@define_diffrule M.f(x, y) = :(df_dx($x, $y)), :(df_dy($x, $y))\n⋮\n\nDefine a new differentiation rule for the function M.f and the given arguments, which should be treated as bindings to Julia expressions. Return the defined rule\'s key.\n\nThe LHS should be a function call with a non-splatted argument list, and the RHS should be the derivative expression, or in the n-ary case, an n-tuple of expressions where the ith expression is the derivative of f w.r.t the ith argument. Arguments should be interpolated wherever they are used on the RHS.\n\nNote that differentiation rules are purely symbolic, so no type annotations should be used.\n\nExamples:\n\n@define_diffrule Base.cos(x)          = :(-sin($x))\n@define_diffrule Base.:/(x, y)        = :(inv($y)), :(-$x / ($y^2))\n@define_diffrule Base.polygamma(m, x) = :NaN,       :(polygamma($m + 1, $x))\n\n\n\n"
},

{
    "location": "index.html#DiffRules.diffrule",
    "page": "Documentation",
    "title": "DiffRules.diffrule",
    "category": "function",
    "text": "diffrule(M::Union{Expr,Symbol}, f::Symbol, args...)\n\nReturn the derivative expression for M.f at the given argument(s), with the argument(s) interpolated into the returned expression.\n\nIn the n-ary case, an n-tuple of expressions will be returned where the ith expression is the derivative of f w.r.t the ith argument.\n\nExamples:\n\njulia> DiffRules.diffrule(:Base, :sin, 1)\n:(cos(1))\n\njulia> DiffRules.diffrule(:Base, :sin, :x)\n:(cos(x))\n\njulia> DiffRules.diffrule(:Base, :sin, :(x * y^2))\n:(cos(x * y ^ 2))\n\njulia> DiffRules.diffrule(:Base, :^, :(x + 2), :c)\n(:(c * (x + 2) ^ (c - 1)), :((x + 2) ^ c * log(x + 2)))\n\n\n\n"
},

{
    "location": "index.html#DiffRules.hasdiffrule",
    "page": "Documentation",
    "title": "DiffRules.hasdiffrule",
    "category": "function",
    "text": "hasdiffrule(M::Union{Expr,Symbol}, f::Symbol, arity::Int)\n\nReturn true if a differentiation rule is defined for M.f and arity, or return false otherwise.\n\nHere, arity refers to the number of arguments accepted by f.\n\nExamples:\n\njulia> DiffRules.hasdiffrule(:Base, :sin, 1)\ntrue\n\njulia> DiffRules.hasdiffrule(:Base, :sin, 2)\nfalse\n\njulia> DiffRules.hasdiffrule(:Base, :-, 1)\ntrue\n\njulia> DiffRules.hasdiffrule(:Base, :-, 2)\ntrue\n\njulia> DiffRules.hasdiffrule(:Base, :-, 3)\nfalse\n\n\n\n"
},

{
    "location": "index.html#DiffRules.diffrules",
    "page": "Documentation",
    "title": "DiffRules.diffrules",
    "category": "function",
    "text": "diffrules()\n\nReturn a list of keys that can be used to access all defined differentiation rules.\n\nEach key is of the form (M::Symbol, f::Symbol, arity::Int).\n\nHere, arity refers to the number of arguments accepted by f.\n\nExamples:\n\njulia> first(DiffRules.diffrules())\n(:Base, :asind, 1)\n\n\n\n"
},

{
    "location": "index.html#DiffRules-1",
    "page": "Documentation",
    "title": "DiffRules",
    "category": "section",
    "text": "CurrentModule = DiffRulesMany differentiation methods rely on the notion of \"primitive\" differentiation rules that can be composed via various formulations of the chain rule. Using DiffRules, you can define new differentiation rules, query whether or not a given rule exists, and symbolically apply rules to simple Julia expressions.Note that DiffRules is not a fully-fledged symbolic differentiation tool. It is a (very) simple global database of common derivative definitions, and was developed with the goal of improving derivative coverage in downstream tools.DiffRules.@define_diffrule\nDiffRules.diffrule\nDiffRules.hasdiffrule\nDiffRules.diffrules"
},

]}
