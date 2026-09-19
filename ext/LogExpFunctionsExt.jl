module LogExpFunctionsExt

using DiffRules: @define_diffrule
using LogExpFunctions:
    LogExpFunctions, log1mexp, log2mexp, logaddexp, logexpm1, logistic, logsubexp

###################
# LogExpFunctions #
###################

# unary
@define_diffrule LogExpFunctions.xlogx(x) = :(1 + log($x))
@define_diffrule LogExpFunctions.logistic(x) = :(z = $logistic($x); z * (1 - z))
@define_diffrule LogExpFunctions.logit(x) = :(inv($x * (1 - $x)))
@define_diffrule LogExpFunctions.log1psq(x) = :(2 * $x / (1 + $x^2))
@define_diffrule LogExpFunctions.log1pexp(x) = :($logistic($x))
@define_diffrule LogExpFunctions.log1mexp(x) = :(-exp($x - $log1mexp($x)))
@define_diffrule LogExpFunctions.log2mexp(x) = :(-exp($x - $log2mexp($x)))
@define_diffrule LogExpFunctions.logexpm1(x) = :(exp($x - $logexpm1($x)))
@define_diffrule LogExpFunctions.log1pmx(x) = :(-$x / (1 + $x))
@define_diffrule LogExpFunctions.logmxp1(x) = :((1 - $x) / $x)

# binary
@define_diffrule LogExpFunctions.xlogy(x, y) =
    :(log($y)),
    :(z = $x / $y; iszero($x) && !isnan($y) ? zero(z) : z)
@define_diffrule LogExpFunctions.logaddexp(x, y) =
    :(exp($x - $logaddexp($x, $y))), :(exp($y - $logaddexp($x, $y)))
@define_diffrule LogExpFunctions.logsubexp(x, y) =
    :(z = $logsubexp($x, $y); $x > $y ? exp($x - z) : -exp($x - z)),
    :(z = $logsubexp($x, $y); $x > $y ? -exp($y - z) : exp($y - z))
@define_diffrule LogExpFunctions.xlog1py(x, y) =
    :(log1p($y)),
    :(z = $x / (1 + $y); iszero($x) && !isnan($y) ? zero(z) : z)

end # module
