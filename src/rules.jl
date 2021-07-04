################
# General Math #
################

# unary #
#-------#

@define_diffrule Base.:+(x)                   = :(  1                                  )
@define_diffrule Base.:-(x)                   = :( -1                                  )
@define_diffrule Base.sqrt(x)                 = :(  inv(2 * sqrt($x))                  )
@define_diffrule Base.cbrt(x)                 = :(  inv(3 * cbrt($x)^2)                )
@define_diffrule Base.abs2(x)                 = :(  $x + $x                            )
@define_diffrule Base.inv(x)                  = :( -abs2(inv($x))                      )
@define_diffrule Base.log(x)                  = :(  inv($x)                            )
@define_diffrule Base.log10(x)                = :(  inv($x) / log(10)                  )
@define_diffrule Base.log2(x)                 = :(  inv($x) / log(2)                   )
@define_diffrule Base.log1p(x)                = :(  inv($x + 1)                        )
@define_diffrule Base.exp(x)                  = :(  exp($x)                            )
@define_diffrule Base.exp2(x)                 = :(  exp2($x) * log(2)                  )
@define_diffrule Base.exp10(x)                = :(  exp10($x) * log(10)                )
@define_diffrule Base.expm1(x)                = :(  exp($x)                            )
@define_diffrule Base.sin(x)                  = :(  cos($x)                            )
@define_diffrule Base.cos(x)                  = :( -sin($x)                            )
@define_diffrule Base.tan(x)                  = :(  1 + tan($x)^2                      )
@define_diffrule Base.sec(x)                  = :(  sec($x) * tan($x)                  )
@define_diffrule Base.csc(x)                  = :( -csc($x) * cot($x)                  )
@define_diffrule Base.cot(x)                  = :( -(1 + cot($x)^2)                    )
@define_diffrule Base.sind(x)                 = :(  (π / 180) * cosd($x)               )
@define_diffrule Base.cosd(x)                 = :( -(π / 180) * sind($x)               )
@define_diffrule Base.tand(x)                 = :(  (π / 180) * (1 + tand($x)^2)       )
@define_diffrule Base.secd(x)                 = :(  (π / 180) * secd($x) * tand($x)    )
@define_diffrule Base.cscd(x)                 = :( -(π / 180) * cscd($x) * cotd($x)    )
@define_diffrule Base.cotd(x)                 = :( -(π / 180) * (1 + cotd($x)^2)       )
@define_diffrule Base.sinpi(x)                = :(  π * cospi($x)                      )
@define_diffrule Base.cospi(x)                = :( -π * sinpi($x)                      )
@define_diffrule Base.asin(x)                 = :(  inv(sqrt(1 - $x^2))                )
@define_diffrule Base.acos(x)                 = :( -inv(sqrt(1 - $x^2))                )
@define_diffrule Base.atan(x)                 = :(  inv(1 + $x^2)                      )
@define_diffrule Base.asec(x)                 = :(  inv(abs($x) * sqrt($x^2 - 1))      )
@define_diffrule Base.acsc(x)                 = :( -inv(abs($x) * sqrt($x^2 - 1))      )
@define_diffrule Base.acot(x)                 = :( -inv(1 + $x^2)                      )
@define_diffrule Base.asind(x)                = :(  180 / π / sqrt(1 - $x^2)           )
@define_diffrule Base.acosd(x)                = :( -180 / π / sqrt(1 - $x^2)           )
@define_diffrule Base.atand(x)                = :(  180 / π / (1 + $x^2)               )
@define_diffrule Base.asecd(x)                = :(  180 / π / abs($x) / sqrt($x^2 - 1) )
@define_diffrule Base.acscd(x)                = :( -180 / π / abs($x) / sqrt($x^2 - 1) )
@define_diffrule Base.acotd(x)                = :( -180 / π / (1 + $x^2)               )
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
@define_diffrule Base.deg2rad(x)              = :(  π / 180                            )
@define_diffrule Base.mod2pi(x)               = :(  isinteger($x / 2pi) ? NaN : 1      )
@define_diffrule Base.rad2deg(x)              = :(  180 / π                            )

@define_diffrule SpecialFunctions.gamma(x) =
    :(  SpecialFunctions.digamma($x) * SpecialFunctions.gamma($x)  )
@define_diffrule SpecialFunctions.loggamma(x) =
    :(  SpecialFunctions.digamma($x)  )

@define_diffrule Base.identity(x)             = :(  1                                  )
@define_diffrule Base.conj(x)                 = :(  1                                  )
@define_diffrule Base.adjoint(x)              = :(  1                                  )
@define_diffrule Base.transpose(x)            = :(  1                                  )
@define_diffrule Base.abs(x)                  = :( DiffRules._abs_deriv($x)            )

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

if VERSION < v"0.7-"
    @define_diffrule Base.atan2(x, y)   = :( $y / ($x^2 + $y^2)                                 ), :( -$x / ($x^2 + $y^2)                                                     )
else
    @define_diffrule Base.atan(x, y)    = :( $y / ($x^2 + $y^2)                                 ), :( -$x / ($x^2 + $y^2)                                                     )
end
@define_diffrule Base.hypot(x, y)  = :( $x / hypot($x, $y)                                      ), :(  $y / hypot($x, $y)                                                     )
@define_diffrule Base.log(b, x)    = :( log($x) * inv(-log($b)^2 * $b)                          ), :( inv($x) / log($b)                                                       )

@define_diffrule Base.mod(x, y)    = :( first(promote(ifelse(isinteger($x / $y), NaN, 1), NaN)) ), :(  z = $x / $y; first(promote(ifelse(isinteger(z), NaN, -floor(z)), NaN)) )
@define_diffrule Base.rem(x, y)    = :( first(promote(ifelse(isinteger($x / $y), NaN, 1), NaN)) ), :(  z = $x / $y; first(promote(ifelse(isinteger(z), NaN, -trunc(z)), NaN)) )
@define_diffrule Base.rem2pi(x, r) = :( 1                                                       ), :NaN
@define_diffrule Base.max(x, y)    = :( $x > $y ? one($x) : zero($x)                            ), :( $x > $y ? zero($y) : one($y)                                            )
@define_diffrule Base.min(x, y)    = :( $x > $y ? zero($x) : one($x)                            ), :( $x > $y ? one($y) : zero($y)                                            )

# trinary #
#---------#

@define_diffrule Base.muladd(x, y, z) = :($y), :($x), :(one($z))
@define_diffrule Base.fma(x, y, z)    = :($y), :($x), :(one($z))

@define_diffrule Base.ifelse(p, x, y) = false, :($p), :(!$p)

####################
# SpecialFunctions #
####################

# unary #
#-------#

@define_diffrule SpecialFunctions.erf(x)         = :(  (2 / sqrt(π)) * exp(-$x * $x)       )
@define_diffrule SpecialFunctions.erfinv(x)      =
    :(  (sqrt(π) / 2) * exp(SpecialFunctions.erfinv($x)^2)  )
@define_diffrule SpecialFunctions.erfc(x)        = :( -(2 / sqrt(π)) * exp(-$x * $x)       )
@define_diffrule SpecialFunctions.erfcinv(x)     =
    :( -(sqrt(π) / 2) * exp(SpecialFunctions.erfcinv($x)^2)  )
@define_diffrule SpecialFunctions.erfi(x)        = :(  (2 / sqrt(π)) * exp($x * $x)        )
@define_diffrule SpecialFunctions.erfcx(x)       =
    :(  (2 * $x * SpecialFunctions.erfcx($x)) - (2 / sqrt(π))  )
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
@define_diffrule NaNMath.log2(x)   = :(  inv($x) / NaNMath.log(2)                   )
@define_diffrule NaNMath.log10(x)  = :(  inv($x) / NaNMath.log(10)                  )
@define_diffrule NaNMath.log1p(x)  = :(  inv($x + 1)                                )
@define_diffrule NaNMath.lgamma(x) = :(  SpecialFunctions.digamma($x)               )


# binary #
#--------#

@define_diffrule NaNMath.pow(x, y) = :( $y * NaNMath.pow($x, ($y - 1)) ), :( NaNMath.pow($x, $y) * NaNMath.log($x) )
@define_diffrule NaNMath.max(x, y) = :(ifelse(($y > $x) | (signbit($y) < signbit($x)), ifelse(isnan($y), one($x), zero($x)), ifelse(isnan($x), zero($x), one($x)))),
                                     :(ifelse(($y > $x) | (signbit($y) < signbit($x)), ifelse(isnan($y), zero($y), one($y)), ifelse(isnan($x), one($y), zero($y))))
@define_diffrule NaNMath.min(x, y) = :(ifelse(($y < $x) | (signbit($y) > signbit($x)), ifelse(isnan($y), one($x), zero($x)), ifelse(isnan($x), zero($x), one($x)))),
                                     :(ifelse(($y < $x) | (signbit($y) > signbit($x)), ifelse(isnan($y), zero($y), one($y)), ifelse(isnan($x), one($x), zero($x))))
