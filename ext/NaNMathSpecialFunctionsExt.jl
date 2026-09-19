module NaNMathSpecialFunctionsExt

# `NaNMath.lgamma` is the one rule whose derivative reaches into another package,
# so it needs both loaded (JuliaDiff/DiffRules.jl#106).

using DiffRules: @define_diffrule
using NaNMath: NaNMath
using SpecialFunctions: digamma

@define_diffrule NaNMath.lgamma(x) = :(  $digamma($x)  )

end # module
