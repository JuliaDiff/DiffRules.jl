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
@define_diffrule Base.ldexp(x, y)  = :( oftype($x, exp2($y))                                                ), :NaN

@define_diffrule Base.mod(x, y)    = :( z = $x / $y; ifelse(isinteger(z), oftype(float(z), NaN), one(float(z))) ), :(  z = $x / $y; ifelse(isinteger(z), oftype(float(z), NaN), -floor(float(z))) )
@define_diffrule Base.rem(x, y)    = :( z = $x / $y; ifelse(isinteger(z), oftype(float(z), NaN), one(float(z))) ), :(  z = $x / $y; ifelse(isinteger(z), oftype(float(z), NaN), -trunc(float(z))) )
@define_diffrule Base.rem2pi(x, r) = :( 1                                                       ), :NaN
@define_diffrule Base.max(x, y)    = :( $x > $y ), :( !($x > $y) )
@define_diffrule Base.min(x, y)    = :( !($x > $y) ), :( $x > $y )

# trinary #
#---------#

#=

@define_diffrule Base.muladd(x, y, z) = :($y), :($x), :(one($z))
@define_diffrule Base.fma(x, y, z)    = :($y), :($x), :(one($z))

@define_diffrule Base.ifelse(p, x, y) = false, :($p), :(!$p)

=#
