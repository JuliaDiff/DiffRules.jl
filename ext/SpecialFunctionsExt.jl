module SpecialFunctionsExt

using DiffRules: @define_diffrule
using IrrationalConstants: sqrtπ, invsqrtπ
using SpecialFunctions:
    SpecialFunctions,
    airyai,
    airyaiprime,
    airyaiprimex,
    airyaix,
    airybi,
    airybiprime,
    airybiprimex,
    airybix,
    besselh,
    besselhx,
    besseli,
    besselix,
    besselj,
    besselj0,
    besselj1,
    besseljx,
    besselk,
    besselkx,
    bessely,
    bessely0,
    bessely1,
    besselyx,
    beta,
    dawson,
    digamma,
    ellipe,
    ellipk,
    erfcinv,
    erfcx,
    erfinv,
    expint,
    gamma,
    hankelh1,
    hankelh1x,
    hankelh2,
    hankelh2x,
    invdigamma,
    logerfc,
    polygamma,
    trigamma,
    zeta

####################
# SpecialFunctions #
####################

# unary #
#-------#

@define_diffrule SpecialFunctions.gamma(x)       = :(  $digamma($x) * $gamma($x)          )
@define_diffrule SpecialFunctions.loggamma(x)    = :(  $digamma($x)                       )

@define_diffrule SpecialFunctions.erf(x)         = :(  2 * ($invsqrtπ * exp(-$x^2))       )
@define_diffrule SpecialFunctions.erfinv(x)      = :(  ($sqrtπ * exp($erfinv($x)^2)) / 2  )
@define_diffrule SpecialFunctions.erfc(x)        = :( -($invsqrtπ * exp(-$x^2) * 2)       )
@define_diffrule SpecialFunctions.logerfc(x)     =
    :( -2 * ($invsqrtπ * exp(-$x^2 - $logerfc($x)))  )

@define_diffrule SpecialFunctions.erfcinv(x)     = :( -($sqrtπ * exp($erfcinv($x)^2)) / 2 )
@define_diffrule SpecialFunctions.erfi(x)        = :(  $invsqrtπ * exp($x^2) * 2          )
@define_diffrule SpecialFunctions.erfcx(x)       = :(  2 * (($x * $erfcx($x)) - $invsqrtπ))
@define_diffrule SpecialFunctions.logerfcx(x)    = :(  2 * ($x - inv($erfcx($x) * $sqrtπ)))

@define_diffrule SpecialFunctions.dawson(x)      = :(  1 - (2 * $x * $dawson($x))         )
@define_diffrule SpecialFunctions.digamma(x)     = :(  $trigamma($x)                      )
@define_diffrule SpecialFunctions.invdigamma(x)  = :(  inv($trigamma($invdigamma($x)))    )
@define_diffrule SpecialFunctions.trigamma(x)    = :(  $polygamma(2, $x)                  )

# derivatives for `airybix` and `airybiprimex` are only correct for real inputs
# `airyaix` and `airyaiprimex` are only defined for positive real inputs
# `airybix` and `airybiprimex` are unscaled for negative real inputs
@define_diffrule SpecialFunctions.airyai(x)       = :(  $airyaiprime($x)  )
@define_diffrule SpecialFunctions.airyaiprime(x)  = :(  $x * $airyai($x)  )
@define_diffrule SpecialFunctions.airyaix(x)      =
    :(  $airyaiprimex($x) + sqrt($x) * $airyaix($x)  )
@define_diffrule SpecialFunctions.airyaiprimex(x) =
    :(  $x * $airyaix($x) + sqrt($x) * $airyaiprimex($x)  )
@define_diffrule SpecialFunctions.airybi(x)       = :(  $airybiprime($x)  )
@define_diffrule SpecialFunctions.airybiprime(x)  = :(  $x * $airybi($x)  )
@define_diffrule SpecialFunctions.airybix(x)      =
    :(  if $x > zero($x)
            $airybiprimex($x) - sqrt($x) * $airybix($x)
        else
            $airybiprimex($x)
        end  )
@define_diffrule SpecialFunctions.airybiprimex(x) =
    :(  if $x > zero($x)
            $x * $airybix($x) - sqrt($x) * $airybiprimex($x)
        else
            $x * $airybix($x)
        end  )

@define_diffrule SpecialFunctions.besselj0(x)    = :( -$besselj1($x)                         )
@define_diffrule SpecialFunctions.besselj1(x)    = :(  ($besselj0($x) - $besselj(2, $x)) / 2 )
@define_diffrule SpecialFunctions.bessely0(x)    = :( -$bessely1($x)                         )
@define_diffrule SpecialFunctions.bessely1(x)    = :(  ($bessely0($x) - $bessely(2, $x)) / 2 )

@define_diffrule SpecialFunctions.sinint(x)      = :(  sinc($x / π)  )
@define_diffrule SpecialFunctions.cosint(x)      = :(  cos($x) / $x  )

@define_diffrule SpecialFunctions.ellipk(m)      =
    :( ($ellipe($m) / (1 - $m) - $ellipk($m)) / (2 * $m) )
@define_diffrule SpecialFunctions.ellipe(m)      =
    :( ($ellipe($m) - $ellipk($m)) / (2 * $m) )

@define_diffrule SpecialFunctions.expint(x)      = :( -exp(-$x) / $x )

# TODO:
#
# eta
# zeta

# binary #
#--------#

@define_diffrule SpecialFunctions.erf(x, y) =
    :(  -2 * ($invsqrtπ * exp(-$x^2))  ), :(  2 * ($invsqrtπ * exp(-$y^2))  )

# derivatives with respect to the order `ν` exist but are not implemented
# (analogously to the ChainRules definitions in SpecialFunctions)

# derivatives for `besselix`, `besseljx` and `besselyx` are only correct for real inputs
# see https://github.com/JuliaMath/SpecialFunctions.jl/blob/master/src/chainrules.jl
# for forward-mode and reverse-mode derivatives for complex inputs

@define_diffrule SpecialFunctions.besselj(ν, x)   =
    :NaN, :(  ($besselj($ν - 1, $x) - $besselj($ν + 1, $x)) / 2  )
@define_diffrule SpecialFunctions.besseljx(ν, x)  =
    :NaN, :(  ($besseljx($ν - 1, $x) - $besseljx($ν + 1, $x)) / 2  )
@define_diffrule SpecialFunctions.besseli(ν, x)   =
    :NaN, :(  ($besseli($ν - 1, $x) + $besseli($ν + 1, $x)) / 2  )
@define_diffrule SpecialFunctions.besselix(ν, x)  =
    :NaN, :(  ($besselix($ν - 1, $x) + $besselix($ν + 1, $x)) / 2 - sign($x) * $besselix($ν, $x)  )
@define_diffrule SpecialFunctions.bessely(ν, x)   =
    :NaN, :(  ($bessely($ν - 1, $x) - $bessely($ν + 1, $x)) / 2  )
@define_diffrule SpecialFunctions.besselyx(ν, x)  =
    :NaN, :(  ($besselyx($ν - 1, $x) - $besselyx($ν + 1, $x)) / 2  )
@define_diffrule SpecialFunctions.besselk(ν, x)   =
    :NaN, :( -($besselk($ν - 1, $x) + $besselk($ν + 1, $x)) / 2  )
@define_diffrule SpecialFunctions.besselkx(ν, x)  =
    :NaN, :( -($besselkx($ν - 1, $x) + $besselkx($ν + 1, $x)) / 2 + $besselkx($ν, $x)  )
@define_diffrule SpecialFunctions.besselh(ν, x)   =
    :NaN, :(  ($besselh($ν - 1, $x) - $besselh($ν + 1, $x)) / 2  )
@define_diffrule SpecialFunctions.besselhx(ν, x)  =
    :NaN, :(  ($besselhx($ν - 1, $x) - $besselhx($ν + 1, $x)) / 2 - im * $besselhx($ν, $x)  )
@define_diffrule SpecialFunctions.hankelh1(ν, x)  =
    :NaN, :(  ($hankelh1($ν - 1, $x) - $hankelh1($ν + 1, $x)) / 2  )
@define_diffrule SpecialFunctions.hankelh1x(ν, x) =
    :NaN, :(  ($hankelh1x($ν - 1, $x) - $hankelh1x($ν + 1, $x)) / 2 - im * $hankelh1x($ν, $x)  )
@define_diffrule SpecialFunctions.hankelh2(ν, x)  =
    :NaN, :(  ($hankelh2($ν - 1, $x) - $hankelh2($ν + 1, $x)) / 2  )
@define_diffrule SpecialFunctions.hankelh2x(ν, x) =
    :NaN, :(  ($hankelh2x($ν - 1, $x) - $hankelh2x($ν + 1, $x)) / 2 + im * $hankelh2x($ν, $x)  )

@define_diffrule SpecialFunctions.polygamma(m, x) =
    :NaN, :(  $polygamma($m + 1, $x)  )

@define_diffrule SpecialFunctions.beta(a, b)      =
    :(  $beta($a, $b) * ($digamma($a) - $digamma($a + $b))  ),
    :(  $beta($a, $b) * ($digamma($b) - $digamma($a + $b))  )
@define_diffrule SpecialFunctions.logbeta(a, b)   =
    :(  $digamma($a) - $digamma($a + $b)  ), :(  $digamma($b) - $digamma($a + $b)  )

# derivative wrt to `ν` is not implemented
@define_diffrule SpecialFunctions.expint(ν, x)    =
    :NaN, :( -$expint($ν - 1, $x) )

# derivative wrt to `s` is not implemented
@define_diffrule SpecialFunctions.zeta(s, z)      =
    :NaN, :( -$s * $zeta($s + 1, $z) )

# ternary #
#---------#

# TODO:
#
# besselh
# besselhx

end # module
