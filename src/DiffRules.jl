module DiffRules

using MacroTools

const SymOrExpr = Union{Symbol, Expr}

const AVM = AbstractVecOrMat
const AM = AbstractMatrix

# Simple derivative expressions.
include("diffrules/api.jl")
include("diffrules/rules.jl")

# Forwards-mode stuff.
include("forward/api.jl")

# Reverse-mode stuff.
include("reverse/api.jl")
include("reverse/generic.jl")

end # module
