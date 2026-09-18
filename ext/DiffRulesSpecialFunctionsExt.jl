module DiffRulesSpecialFunctionsExt

using DiffRules: @define_diffrule
using IrrationalConstants: sqrtπ, invsqrtπ
using SpecialFunctions

@define_diffrule SpecialFunctions.gamma(x) =
    :(  $(SpecialFunctions.digamma)($x) * $(SpecialFunctions.gamma)($x)  )
@define_diffrule SpecialFunctions.loggamma(x) =
    :(  $(SpecialFunctions.digamma)($x)  )

####################
# SpecialFunctions #
####################

# unary #
#-------#

@define_diffrule SpecialFunctions.erf(x)         = :(  2 * ($invsqrtπ * exp(-$x^2))       )
@define_diffrule SpecialFunctions.erfinv(x)      =
    :(  ($sqrtπ * exp($(SpecialFunctions.erfinv)($x)^2)) / 2  )
@define_diffrule SpecialFunctions.erfc(x)        = :( -($invsqrtπ * exp(-$x^2) * 2)       )
@define_diffrule SpecialFunctions.logerfc(x)     =
    :( - 2 * ($invsqrtπ * exp(- $x^2 - $(SpecialFunctions.logerfc)($x))) )

@define_diffrule SpecialFunctions.erfcinv(x)     =
    :( -($sqrtπ * exp($(SpecialFunctions.erfcinv)($x)^2)) / 2  )
@define_diffrule SpecialFunctions.erfi(x)        = :(  $invsqrtπ * exp($x^2) * 2        )
@define_diffrule SpecialFunctions.erfcx(x)       =
    :(  2 * (($x * $(SpecialFunctions.erfcx)($x)) - $invsqrtπ)  )
@define_diffrule SpecialFunctions.logerfcx(x) =
    :(  2 * ($x - inv($(SpecialFunctions.erfcx)($x) * $sqrtπ))  )

@define_diffrule SpecialFunctions.dawson(x)      =
    :(  1 - (2 * $x * $(SpecialFunctions.dawson)($x))  )
@define_diffrule SpecialFunctions.digamma(x) =
    :(  $(SpecialFunctions.trigamma)($x)  )
@define_diffrule SpecialFunctions.invdigamma(x)  =
    :(  inv($(SpecialFunctions.trigamma)($(SpecialFunctions.invdigamma)($x)))  )
@define_diffrule SpecialFunctions.trigamma(x)    =
    :(  $(SpecialFunctions.polygamma)(2, $x)  )

# derivatives for `airybix` and `airybiprimex` are only correct for real inputs
# `airyaix` and `airyaiprimex` are only defined for positive real inputs
# `airybix` and `airybiprimex` are unscaled for negative real inputs
@define_diffrule SpecialFunctions.airyai(x)      =
    :(  $(SpecialFunctions.airyaiprime)($x)  )
@define_diffrule SpecialFunctions.airyaiprime(x) =
    :(  $x * $(SpecialFunctions.airyai)($x)  )
@define_diffrule SpecialFunctions.airyaix(x)     =
    :(  $(SpecialFunctions.airyaiprimex)($x) + sqrt($x) * $(SpecialFunctions.airyaix)($x)  )
@define_diffrule SpecialFunctions.airyaiprimex(x) =
    :(  $x * $(SpecialFunctions.airyaix)($x) + sqrt($x) * $(SpecialFunctions.airyaiprimex)($x)  )
@define_diffrule SpecialFunctions.airybi(x)      =
    :(  $(SpecialFunctions.airybiprime)($x)  )
@define_diffrule SpecialFunctions.airybiprime(x) =
    :(  $x * $(SpecialFunctions.airybi)($x)  )
@define_diffrule SpecialFunctions.airybix(x)     =
    :(  if $x > zero($x)
            $(SpecialFunctions.airybiprimex)($x) - sqrt($x) * $(SpecialFunctions.airybix)($x)
        else
            $(SpecialFunctions.airybiprimex)($x)
        end  )
@define_diffrule SpecialFunctions.airybiprimex(x) =
    :(  if $x > zero($x)
            $x * $(SpecialFunctions.airybix)($x) - sqrt($x) * $(SpecialFunctions.airybiprimex)($x)
        else
            $x * $(SpecialFunctions.airybix)($x)
        end  )

@define_diffrule SpecialFunctions.besselj0(x)    =
    :( -$(SpecialFunctions.besselj1)($x)  )
@define_diffrule SpecialFunctions.besselj1(x)    =
    :(  ($(SpecialFunctions.besselj0)($x) - $(SpecialFunctions.besselj)(2, $x)) / 2  )
@define_diffrule SpecialFunctions.bessely0(x)    =
    :( -$(SpecialFunctions.bessely1)($x)  )
@define_diffrule SpecialFunctions.bessely1(x)    =
    :(  ($(SpecialFunctions.bessely0)($x) - $(SpecialFunctions.bessely)(2, $x)) / 2 )

@define_diffrule SpecialFunctions.sinint(x) = :( sinc($x / π) )
@define_diffrule SpecialFunctions.cosint(x) = :( cos($x) / $x )

@define_diffrule SpecialFunctions.ellipk(m) =
    :( ($(SpecialFunctions.ellipe)($m) / (1 - $m) - $(SpecialFunctions.ellipk)($m)) / (2 * $m) )
@define_diffrule SpecialFunctions.ellipe(m) =
    :( ($(SpecialFunctions.ellipe)($m) - $(SpecialFunctions.ellipk)($m)) / (2 * $m) )

@define_diffrule SpecialFunctions.expint(x)    = :( -exp(-$x) / $x )

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
    :NaN, :(  ($(SpecialFunctions.besselj)($ν - 1, $x) - $(SpecialFunctions.besselj)($ν + 1, $x)) / 2  )
@define_diffrule SpecialFunctions.besseljx(ν, x)  =
    :NaN, :(  ($(SpecialFunctions.besseljx)($ν - 1, $x) - $(SpecialFunctions.besseljx)($ν + 1, $x)) / 2  )
@define_diffrule SpecialFunctions.besseli(ν, x)   =
    :NaN, :(  ($(SpecialFunctions.besseli)($ν - 1, $x) + $(SpecialFunctions.besseli)($ν + 1, $x)) / 2  )
@define_diffrule SpecialFunctions.besselix(ν, x)  =
    :NaN, :(  ($(SpecialFunctions.besselix)($ν - 1, $x) + $(SpecialFunctions.besselix)($ν + 1, $x)) / 2 - sign($x) * $(SpecialFunctions.besselix)($ν, $x)  )
@define_diffrule SpecialFunctions.bessely(ν, x)   =
    :NaN, :(  ($(SpecialFunctions.bessely)($ν - 1, $x) - $(SpecialFunctions.bessely)($ν + 1, $x)) / 2  )
@define_diffrule SpecialFunctions.besselyx(ν, x)  =
    :NaN, :(  ($(SpecialFunctions.besselyx)($ν - 1, $x) - $(SpecialFunctions.besselyx)($ν + 1, $x)) / 2  )
@define_diffrule SpecialFunctions.besselk(ν, x)   =
    :NaN, :( -($(SpecialFunctions.besselk)($ν - 1, $x) + $(SpecialFunctions.besselk)($ν + 1, $x)) / 2  )
@define_diffrule SpecialFunctions.besselkx(ν, x)  =
    :NaN, :( -($(SpecialFunctions.besselkx)($ν - 1, $x) + $(SpecialFunctions.besselkx)($ν + 1, $x)) / 2 + $(SpecialFunctions.besselkx)($ν, $x)  )
@define_diffrule SpecialFunctions.besselh(ν, x)   =
    :NaN, :(  ($(SpecialFunctions.besselh)($ν - 1, $x) - $(SpecialFunctions.besselh)($ν + 1, $x)) / 2  )
@define_diffrule SpecialFunctions.besselhx(ν, x) =
    :NaN, :(  ($(SpecialFunctions.besselhx)($ν - 1, $x) - $(SpecialFunctions.besselhx)($ν + 1, $x)) / 2 - im * $(SpecialFunctions.besselhx)($ν, $x)  )
@define_diffrule SpecialFunctions.hankelh1(ν, x)  =
    :NaN, :(  ($(SpecialFunctions.hankelh1)($ν - 1, $x) - $(SpecialFunctions.hankelh1)($ν + 1, $x)) / 2  )
@define_diffrule SpecialFunctions.hankelh1x(ν, x) =
    :NaN, :(  ($(SpecialFunctions.hankelh1x)($ν - 1, $x) - $(SpecialFunctions.hankelh1x)($ν + 1, $x)) / 2 - im * $(SpecialFunctions.hankelh1x)($ν, $x)  )
@define_diffrule SpecialFunctions.hankelh2(ν, x)  =
    :NaN, :(  ($(SpecialFunctions.hankelh2)($ν - 1, $x) - $(SpecialFunctions.hankelh2)($ν + 1, $x)) / 2  )
@define_diffrule SpecialFunctions.hankelh2x(ν, x) =
    :NaN, :(  ($(SpecialFunctions.hankelh2x)($ν - 1, $x) - $(SpecialFunctions.hankelh2x)($ν + 1, $x)) / 2 + im * $(SpecialFunctions.hankelh2x)($ν, $x)  )

@define_diffrule SpecialFunctions.polygamma(m, x) =
    :NaN, :(  $(SpecialFunctions.polygamma)($m + 1, $x)  )

@define_diffrule SpecialFunctions.beta(a, b)      =
    :( $(SpecialFunctions.beta)($a, $b)*($(SpecialFunctions.digamma)($a) - $(SpecialFunctions.digamma)($a + $b)) ), :(  $(SpecialFunctions.beta)($a, $b)*($(SpecialFunctions.digamma)($b) - $(SpecialFunctions.digamma)($a + $b))     )
@define_diffrule SpecialFunctions.logbeta(a, b)   =
    :( $(SpecialFunctions.digamma)($a) - $(SpecialFunctions.digamma)($a + $b)  ), :(  $(SpecialFunctions.digamma)($b) - $(SpecialFunctions.digamma)($a + $b)  )

# derivative wrt to `ν` is not implemented
@define_diffrule SpecialFunctions.expint(ν, x)    = 
    :NaN, :( -$(SpecialFunctions.expint)($ν - 1, $x) )

# derivative wrt to `s` is not implemented
@define_diffrule SpecialFunctions.zeta(s, z)      =
    :NaN, :( - $s * $(SpecialFunctions.zeta)($s + 1, $z)  )
end # module
