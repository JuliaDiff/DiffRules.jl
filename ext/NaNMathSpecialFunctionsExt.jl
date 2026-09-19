module DiffRulesNaNMathSpecialFunctionsExt

# `NaNMath.lgamma` is the one rule whose derivative reaches into another package,
# so it needs both loaded (JuliaDiff/DiffRules.jl#106).

using DiffRules: @define_diffrule
using NaNMath
using SpecialFunctions

@define_diffrule NaNMath.lgamma(x) = :(  $(SpecialFunctions.digamma)($x)               )
end # module
