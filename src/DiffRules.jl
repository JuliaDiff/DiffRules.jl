module DiffRules

using MacroTools

const SymOrExpr = Union{Symbol, Expr}

include("api.jl")
include("rules.jl")

include("forward.jl")

end # module
