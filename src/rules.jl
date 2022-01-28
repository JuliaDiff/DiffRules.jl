################
# General Math #
################

# unary #
#-------#

@define_diffrule Base.:+(x)                   = :(   1                                 )
@define_diffrule Base.:-(x)                   = :(  -1                                 )
@define_diffrule Base.sqrt(x)                 = :(  inv(2 * sqrt($x))                  )
@define_diffrule Base.cbrt(x)                 = :(  inv(3 * cbrt($x)^2)                )
@define_diffrule Base.abs2(x)                 = :(  $x + $x                            )
@define_diffrule Base.inv(x)                  = :( -abs2(inv($x))                      )
@define_diffrule Base.log(x)                  = :(  inv($x)                            )
@define_diffrule Base.log10(x)                = :(  inv($x) / $logten      )
@define_diffrule Base.log2(x)                 = :(  inv($x) / $logtwo      )
@define_diffrule Base.log1p(x)                = :(  inv($x + 1)                        )
@define_diffrule Base.exp(x)                  = :(  exp($x)                            )
@define_diffrule Base.exp2(x)                 = :(  exp2($x) * $logtwo     )
@define_diffrule Base.exp10(x)                = :(  exp10($x) * $logten    )
@define_diffrule Base.expm1(x)                = :(  exp($x)                            )
@define_diffrule Base.sin(x)                  = :(  cos($x)                            )
@define_diffrule Base.cos(x)                  = :( -sin($x)                            )
@define_diffrule Base.tan(x)                  = :(  1 + tan($x)^2                      )
@define_diffrule Base.sec(x)                  = :(  sec($x) * tan($x)                  )
@define_diffrule Base.csc(x)                  = :( -csc($x) * cot($x)                  )
@define_diffrule Base.cot(x)                  = :( -(1 + cot($x)^2)                    )
@define_diffrule Base.sind(x)                 = :(  deg2rad(cosd($x))  )
@define_diffrule Base.cosd(x)                 = :( - deg2rad(sind($x))   )
@define_diffrule Base.tand(x)                 = :(  deg2rad(1 + tand($x)^2)  )
@define_diffrule Base.secd(x)                 = :(  deg2rad(secd($x) * tand($x))  )
@define_diffrule Base.cscd(x)                 = :( - deg2rad(cscd($x) * cotd($x))  )
@define_diffrule Base.cotd(x)                 = :( - deg2rad(1 + cotd($x)^2)  )
@define_diffrule Base.sinpi(x)                = :(  π * cospi($x)                      )
@define_diffrule Base.cospi(x)                = :( -(π * sinpi($x))                    )
@define_diffrule Base.asin(x)                 = :(  inv(sqrt(1 - $x^2))                )
@define_diffrule Base.acos(x)                 = :( -inv(sqrt(1 - $x^2))                )
@define_diffrule Base.atan(x)                 = :(  inv(1 + $x^2)                      )
@define_diffrule Base.asec(x)                 = :(  inv(abs($x) * sqrt($x^2 - 1))      )
@define_diffrule Base.acsc(x)                 = :( -inv(abs($x) * sqrt($x^2 - 1))      )
@define_diffrule Base.acot(x)                 = :( -inv(1 + $x^2)                      )
@define_diffrule Base.asind(x)                = :(  inv(deg2rad(sqrt(1 - $x^2)))  )
@define_diffrule Base.acosd(x)                = :(  -inv(deg2rad(sqrt(1 - $x^2)))  )
@define_diffrule Base.atand(x)                = :(  inv(deg2rad(1 + $x^2))  )
@define_diffrule Base.asecd(x)                = :(  inv(deg2rad(abs($x) * sqrt($x^2 - 1)))  )
@define_diffrule Base.acscd(x)                = :(  -inv(deg2rad(abs($x) * sqrt($x^2 - 1)))  )
@define_diffrule Base.acotd(x)                = :(  -inv(deg2rad(1 + $x^2))  )
@define_diffrule Base.sinh(x)                 = :(  cosh($x)                           )
@define_diffrule Base.cosh(x)                 = :(  sinh($x)                           )
@define_diffrule Base.tanh(x)                 = :(  1 - tanh($x)^2                     )
@define_diffrule Base.sech(x)                 = :( -tanh($x) * sech($x)                )
@define_diffrule Base.csch(x)                 = :( -coth($x) * csch($x)                )
@define_diffrule Base.coth(x)                 = :( -(csch($x)^2)                       )
@define_diffrule Base.asinh(x)                = :(  inv(sqrt($x^2 + 1))                )
@define_diffrule Base.acosh(x)                = :(  inv(sqrt($x^2 - 1))                )
@define_diffrule Base.atanh(x)                = :(  inv(1 - $x^2)                      )
@define_diffrule Base.asech(x)                = :( -inv($x * sqrt(1 - $x^2))           )
@define_diffrule Base.acsch(x)                = :( -inv(abs($x) * sqrt(1 + $x^2))      )
@define_diffrule Base.acoth(x)                = :(  inv(1 - $x^2)                      )
@define_diffrule Base.sinc(x)                 = :(  cosc($x)                           )
@define_diffrule Base.deg2rad(x)              = :(   deg2rad(one($x))  )
@define_diffrule Base.mod2pi(x)               = :(  isinteger($x / $twoπ) ? oftype(float($x), NaN) : one(float($x)) )
@define_diffrule Base.rad2deg(x)              = :(  rad2deg(one($x))  )
@define_diffrule SpecialFunctions.gamma(x) =
    :(  SpecialFunctions.digamma($x) * SpecialFunctions.gamma($x)  )
@define_diffrule SpecialFunctions.loggamma(x) =
    :(  SpecialFunctions.digamma($x)  )

@define_diffrule Base.abs(x)                  = :( $(_abs_deriv)($x)            )

# We provide this hook for special number types like `Interval`
# that need their own special definition of `abs`.
_abs_deriv(x) = signbit(x) ? -one(x) : one(x)

# binary #
#--------#

@define_diffrule Base.:+(x, y) = :( one($x)            ), :(  one($y)           )
@define_diffrule Base.:-(x, y) = :( one($x)            ), :( -one($y)           )
@define_diffrule Base.:*(x, y) = :( $y                 ), :(  $x                )
@define_diffrule Base.:/(x, y) = :( one($x) / $y       ), :( -($x / $y / $y)    )
@define_diffrule Base.:\(x, y) = :( -($y / $x / $x)    ), :( one($y) / ($x)     )
@define_diffrule Base.:^(x, y) = :( $y * ($x^($y - 1)) ), :( ($x isa Real && $x<=0) ? Base.oftype(float($x), NaN) : ($x^$y)*log($x) )

@define_diffrule Base.atan(x, y)    = :( $y / ($x^2 + $y^2)                                 ), :( -$x / ($x^2 + $y^2)                                                     )
@define_diffrule Base.hypot(x, y)  = :( $x / hypot($x, $y)                                      ), :(  $y / hypot($x, $y)                                                     )
@define_diffrule Base.log(b, x)    = :( log($x) * inv(-log($b)^2 * $b)                          ), :( inv($x) / log($b)                                                       )
@define_diffrule Base.ldexp(x, y)  = :( exp2($y)                                                ), :NaN

@define_diffrule Base.mod(x, y)    = :( z = $x / $y; ifelse(isinteger(z), oftype(float(z), NaN), one(float(z))) ), :(  z = $x / $y; ifelse(isinteger(z), oftype(float(z), NaN), -floor(float(z))) )
@define_diffrule Base.rem(x, y)    = :( z = $x / $y; ifelse(isinteger(z), oftype(float(z), NaN), one(float(z))) ), :(  z = $x / $y; ifelse(isinteger(z), oftype(float(z), NaN), -trunc(float(z))) )
@define_diffrule Base.rem2pi(x, r) = :( 1                                                       ), :NaN
@define_diffrule Base.max(x, y)    = :( $x > $y ? one($x) : zero($x)                            ), :( $x > $y ? zero($y) : one($y)                                            )
@define_diffrule Base.min(x, y)    = :( $x > $y ? zero($x) : one($x)                            ), :( $x > $y ? one($y) : zero($y)                                            )

# trinary #
#---------#

#=

@define_diffrule Base.muladd(x, y, z) = :($y), :($x), :(one($z))
@define_diffrule Base.fma(x, y, z)    = :($y), :($x), :(one($z))

@define_diffrule Base.ifelse(p, x, y) = false, :($p), :(!$p)

=#

####################
# SpecialFunctions #
####################

# unary #
#-------#

@define_diffrule SpecialFunctions.erf(x)         = :(  2 * ($invsqrtπ * exp(-$x^2))       )
@define_diffrule SpecialFunctions.erfinv(x)      =
    :(  ($sqrtπ * exp(SpecialFunctions.erfinv($x)^2)) / 2  )
@define_diffrule SpecialFunctions.erfc(x)        = :( -($invsqrtπ * exp(-$x^2) * 2)       )
@define_diffrule SpecialFunctions.erfcinv(x)     =
    :( -($sqrtπ * exp(SpecialFunctions.erfcinv($x)^2)) / 2  )
@define_diffrule SpecialFunctions.erfi(x)        = :(  $invsqrtπ * exp($x^2) * 2        )
@define_diffrule SpecialFunctions.erfcx(x)       =
    :(  2 * (($x * SpecialFunctions.erfcx($x)) - $invsqrtπ)  )
@define_diffrule SpecialFunctions.logerfcx(x) =
    :(  2 * ($x - inv(SpecialFunctions.erfcx($x) * $sqrtπ))  )
@define_diffrule SpecialFunctions.dawson(x)      =
    :(  1 - (2 * $x * SpecialFunctions.dawson($x))  )
@define_diffrule SpecialFunctions.digamma(x) =
    :(  SpecialFunctions.trigamma($x)  )
@define_diffrule SpecialFunctions.invdigamma(x)  =
    :(  inv(SpecialFunctions.trigamma(SpecialFunctions.invdigamma($x)))  )
@define_diffrule SpecialFunctions.trigamma(x)    =
    :(  SpecialFunctions.polygamma(2, $x)  )
@define_diffrule SpecialFunctions.airyai(x)      =
    :(  SpecialFunctions.airyaiprime($x)  )
@define_diffrule SpecialFunctions.airyaiprime(x) =
    :(  $x * SpecialFunctions.airyai($x)  )
@define_diffrule SpecialFunctions.airybi(x)      =
    :(  SpecialFunctions.airybiprime($x)  )
@define_diffrule SpecialFunctions.airybiprime(x) =
    :(  $x * SpecialFunctions.airybi($x)  )
@define_diffrule SpecialFunctions.besselj0(x)    =
    :( -SpecialFunctions.besselj1($x)  )
@define_diffrule SpecialFunctions.besselj1(x)    =
    :(  (SpecialFunctions.besselj0($x) - SpecialFunctions.besselj(2, $x)) / 2  )
@define_diffrule SpecialFunctions.bessely0(x)    =
    :( -SpecialFunctions.bessely1($x)  )
@define_diffrule SpecialFunctions.bessely1(x)    =
    :(  (SpecialFunctions.bessely0($x) - SpecialFunctions.bessely(2, $x)) / 2 )

# TODO:
#
# eta
# zeta
# airyaix
# airyaiprimex
# airybix
# airybiprimex

# binary #
#--------#

@define_diffrule SpecialFunctions.besselj(ν, x)   =
    :NaN, :(  (SpecialFunctions.besselj($ν - 1, $x) - SpecialFunctions.besselj($ν + 1, $x)) / 2  )
@define_diffrule SpecialFunctions.besseli(ν, x)   =
    :NaN, :(  (SpecialFunctions.besseli($ν - 1, $x) + SpecialFunctions.besseli($ν + 1, $x)) / 2  )
@define_diffrule SpecialFunctions.bessely(ν, x)   =
    :NaN, :(  (SpecialFunctions.bessely($ν - 1, $x) - SpecialFunctions.bessely($ν + 1, $x)) / 2  )
@define_diffrule SpecialFunctions.besselk(ν, x)   =
    :NaN, :( -(SpecialFunctions.besselk($ν - 1, $x) + SpecialFunctions.besselk($ν + 1, $x)) / 2  )
@define_diffrule SpecialFunctions.hankelh1(ν, x)  =
    :NaN, :(  (SpecialFunctions.hankelh1($ν - 1, $x) - SpecialFunctions.hankelh1($ν + 1, $x)) / 2  )
@define_diffrule SpecialFunctions.hankelh2(ν, x)  =
    :NaN, :(  (SpecialFunctions.hankelh2($ν - 1, $x) - SpecialFunctions.hankelh2($ν + 1, $x)) / 2  )
@define_diffrule SpecialFunctions.polygamma(m, x) =
    :NaN, :(  SpecialFunctions.polygamma($m + 1, $x)  )
@define_diffrule SpecialFunctions.beta(a, b)      =
    :( SpecialFunctions.beta($a, $b)*(SpecialFunctions.digamma($a) - SpecialFunctions.digamma($a + $b)) ), :(  SpecialFunctions.beta($a, $b)*(SpecialFunctions.digamma($b) - SpecialFunctions.digamma($a + $b))     )
@define_diffrule SpecialFunctions.logbeta(a, b)     =
    :( SpecialFunctions.digamma($a) - SpecialFunctions.digamma($a + $b)  ), :(  SpecialFunctions.digamma($b) - SpecialFunctions.digamma($a + $b)  )

# TODO:
#
# zeta
# besseljx
# besselyx
# besselix
# besselkx
# besselh
# besselhx
# hankelh1x
# hankelh2
# hankelh2x

# ternary #
#---------#

# TODO:
#
# besselh
# besselhx

###########
# NaNMath #
###########

# unary #
#-------#

@define_diffrule NaNMath.sqrt(x)   = :(  inv(2 * NaNMath.sqrt($x))                  )
@define_diffrule NaNMath.sin(x)    = :(  NaNMath.cos($x)                            )
@define_diffrule NaNMath.cos(x)    = :( -NaNMath.sin($x)                            )
@define_diffrule NaNMath.tan(x)    = :(  1 + NaNMath.pow(NaNMath.tan($x), 2)        )
@define_diffrule NaNMath.asin(x)   = :(  inv(NaNMath.sqrt(1 - NaNMath.pow($x, 2)))  )
@define_diffrule NaNMath.acos(x)   = :(  -inv(NaNMath.sqrt(1 - NaNMath.pow($x, 2))) )
@define_diffrule NaNMath.acosh(x)  = :(  inv(NaNMath.sqrt(NaNMath.pow($x, 2) - 1))  )
@define_diffrule NaNMath.atanh(x)  = :(  inv(1 - NaNMath.pow($x, 2))                )
@define_diffrule NaNMath.log(x)    = :(  inv($x)                                    )
@define_diffrule NaNMath.log2(x)   = :(  inv($logtwo * $x)                          )
@define_diffrule NaNMath.log10(x)  = :(  inv($logten * $x)                          )
@define_diffrule NaNMath.log1p(x)  = :(  inv($x + 1)                                )
@define_diffrule NaNMath.lgamma(x) = :(  SpecialFunctions.digamma($x)               )


# binary #
#--------#

@define_diffrule NaNMath.pow(x, y) = :( $y * NaNMath.pow($x, ($y - 1)) ), :( NaNMath.pow($x, $y) * NaNMath.log($x) )
@define_diffrule NaNMath.max(x, y) = :(ifelse(($y > $x) | (signbit($y) < signbit($x)), ifelse(isnan($y), one($x), zero($x)), ifelse(isnan($x), zero($x), one($x)))),
                                     :(ifelse(($y > $x) | (signbit($y) < signbit($x)), ifelse(isnan($y), zero($y), one($y)), ifelse(isnan($x), one($y), zero($y))))
@define_diffrule NaNMath.min(x, y) = :(ifelse(($y < $x) | (signbit($y) > signbit($x)), ifelse(isnan($y), one($x), zero($x)), ifelse(isnan($x), zero($x), one($x)))),
                                     :(ifelse(($y < $x) | (signbit($y) > signbit($x)), ifelse(isnan($y), zero($y), one($y)), ifelse(isnan($x), one($x), zero($x))))

###################
# LogExpFunctions #
###################

# unary
@define_diffrule LogExpFunctions.xlogx(x) = :(1 + log($x))
@define_diffrule LogExpFunctions.logistic(x) = :(z = LogExpFunctions.logistic($x); z * (1 - z))
@define_diffrule LogExpFunctions.logit(x) = :(inv($x * (1 - $x)))
@define_diffrule LogExpFunctions.log1psq(x) = :(2 * $x / (1 + $x^2))
@define_diffrule LogExpFunctions.log1pexp(x) = :(LogExpFunctions.logistic($x))
@define_diffrule LogExpFunctions.log1mexp(x) = :(-exp($x - LogExpFunctions.log1mexp($x)))
@define_diffrule LogExpFunctions.log2mexp(x) = :(-exp($x - LogExpFunctions.log2mexp($x)))
@define_diffrule LogExpFunctions.logexpm1(x) = :(exp($x - LogExpFunctions.logexpm1($x)))

# binary
@define_diffrule LogExpFunctions.xlogy(x, y) = :(log($y)), :($x / $y)
@define_diffrule LogExpFunctions.logaddexp(x, y) =
    :(exp($x - LogExpFunctions.logaddexp($x, $y))), :(exp($y - LogExpFunctions.logaddexp($x, $y)))
@define_diffrule LogExpFunctions.logsubexp(x, y) =
    :(z = LogExpFunctions.logsubexp($x, $y); $x > $y ? exp($x - z) : -exp($x - z)),
    :(z = LogExpFunctions.logsubexp($x, $y); $x > $y ? -exp($y - z) : exp($y - z))

@define_diffrule LogExpFunctions.xlog1py(x, y) = :(log1p($y)), :($x / (1 + $y))
