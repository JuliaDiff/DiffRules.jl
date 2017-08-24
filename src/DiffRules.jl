module DiffRules

const DEFINED_DIFFRULES = Tuple{Symbol,Symbol,Int}[]

struct DiffRule{M,f} end

(::Type{DiffRule{M,f}})(args...) where {M,f} = error("no derivative rule defined for $(M).$(f) with arguments $args")

"""
    @define_diffrule M.f(x) = :(df_dx(\$x))
    @define_diffrule M.f(x, y) = :(df_dx(\$x, \$y)), :(df_dy(\$x, \$y))
    ⋮

Define a new differentiation rule for the function `M.f` and the given arguments, which should
be treated as bindings to Julia expressions.

The RHS should be a function call with a non-splatted argument list, and the LHS should be
the derivative expression, or in the `n`-ary case, an `n`-tuple of expressions where the
`i`th expression is the derivative of `f` w.r.t the `i`th argument. Arguments should be
interpolated wherever they are used on the RHS.

Note that differentiation rules are purely symbolic, so no type annotations should be used.

Examples:

    @define_diffrule Base.cos(x)          = :(-sin(\$x))
    @define_diffrule Base.:/(x, y)        = :(inv(\$y)), :(-\$x / (\$y^2))
    @define_diffrule Base.polygamma(m, x) = :NaN,       :(polygamma(\$m + 1, \$x))

"""
macro define_diffrule(def)
    @assert isa(def, Expr) && def.head == :(=)
    rhs = def.args[1]
    lhs = def.args[2]
    @assert isa(rhs, Expr) && rhs.head == :call
    qualified_f = rhs.args[1]
    @assert isa(qualified_f, Expr) && qualified_f.head == :(.)
    M, quoted_f = qualified_f.args
    @assert isa(quoted_f, Expr) && quoted_f.head == :quote
    f = first(quoted_f.args)
    args = rhs.args[2:end]
    rhs.args[1] = :(::Type{DiffRule{$(Expr(:quote, M)),$(Expr(:quote, f))}})
    key = (M, f, length(args))
    in(DEFINED_DIFFRULES, key) || push!(DEFINED_DIFFRULES, key)
    return esc(def)
end

"""
    diffrule(M::Symbol, f::Symbol, args...)

Return the derivative expression for `M.f` at the given argument(s), with the argument(s)
interpolated into the returned expression.

In the `n`-ary case, an `n`-tuple of expressions will be returned where the `i`th expression
is the derivative of `f` w.r.t the `i`th argument.

Examples:

    julia> DiffBase.diffrule(:Base, :sin, 1)
    :(cos(1))

    julia> DiffBase.diffrule(:Base, :sin, :x)
    :(cos(x))

    julia> DiffBase.diffrule(:Base, :sin, :(x * y^2))
    :(cos(x * y ^ 2))

    julia> DiffBase.diffrule(:Base, :^, :(x + 2), :c)
    (:(c * (x + 2) ^ (c - 1)), :((x + 2) ^ c * log(x + 2)))
"""
diffrule(M::Symbol, f::Symbol, args...) = DiffRule{M,f}(args...)

"""
    hasdiffrule(M::Symbol, f::Symbol, arity::Int)

Return `true` if a differentiation rule is defined for `M.f` and `arity`, or return `false`
otherwise.

Here, `arity` refers to the number of arguments accepted by `f`.

Examples:

    julia> DiffBase.hasdiffrule(:Base, :sin, 1)
    true

    julia> DiffBase.hasdiffrule(:Base, :sin, 2)
    false

    julia> DiffBase.hasdiffrule(:Base, :-, 1)
    true

    julia> DiffBase.hasdiffrule(:Base, :-, 2)
    true

    julia> DiffBase.hasdiffrule(:Base, :-, 3)
    false
"""
hasdiffrule(M::Symbol, f::Symbol, arity::Int) = in((M, f, arity), DEFINED_DIFFRULES)

"""
    diffrules()

Return a list of keys that can be used to access all defined differentiation rules.

Each key is of the form `(M::Symbol, f::Symbol, arity::Int)`.

Here, `arity` refers to the number of arguments accepted by `f`.

Examples:

    julia> first(diffrules())
    true

"""
diffrules() = DEFINED_DIFFRULES

################
# General Math #
################

# unary #
#-------#

@define_diffrule Base.:+(x)      = :(  1                                  )
@define_diffrule Base.:-(x)      = :( -1                                  )
@define_diffrule Base.sqrt(x)    = :(  inv(2 * sqrt($x))                  )
@define_diffrule Base.cbrt(x)    = :(  inv(3 * cbrt($x)^2)                )
@define_diffrule Base.abs2(x)    = :(  $x + $x                            )
@define_diffrule Base.inv(x)     = :( -abs2(inv($x))                      )
@define_diffrule Base.log(x)     = :(  inv($x)                            )
@define_diffrule Base.log10(x)   = :(  inv($x) / log(10)                  )
@define_diffrule Base.log2(x)    = :(  inv($x) / log(2)                   )
@define_diffrule Base.log1p(x)   = :(  inv($x + 1)                        )
@define_diffrule Base.exp(x)     = :(  exp($x)                            )
@define_diffrule Base.exp2(x)    = :(  exp2($x) * log(2)                  )
@define_diffrule Base.expm1(x)   = :(  exp($x)                            )
@define_diffrule Base.sin(x)     = :(  cos($x)                            )
@define_diffrule Base.cos(x)     = :( -sin($x)                            )
@define_diffrule Base.tan(x)     = :(  1 + tan($x)^2                      )
@define_diffrule Base.sec(x)     = :(  sec($x) * tan($x)                  )
@define_diffrule Base.csc(x)     = :( -csc($x) * cot($x)                  )
@define_diffrule Base.cot(x)     = :( -(1 + cot($x)^2)                    )
@define_diffrule Base.sind(x)    = :(  (π / 180) * cosd($x)               )
@define_diffrule Base.cosd(x)    = :( -(π / 180) * sind($x)               )
@define_diffrule Base.tand(x)    = :(  (π / 180) * (1 + tand($x)^2)       )
@define_diffrule Base.secd(x)    = :(  (π / 180) * secd($x) * tand($x)    )
@define_diffrule Base.cscd(x)    = :( -(π / 180) * cscd($x) * cotd($x)    )
@define_diffrule Base.cotd(x)    = :( -(π / 180) * (1 + cotd($x)^2)       )
@define_diffrule Base.asin(x)    = :(  inv(sqrt(1 - $x^2))                )
@define_diffrule Base.acos(x)    = :( -inv(sqrt(1 - $x^2))                )
@define_diffrule Base.atan(x)    = :(  inv(1 + $x^2)                      )
@define_diffrule Base.asec(x)    = :(  inv(abs($x) * sqrt($x^2 - 1))      )
@define_diffrule Base.acsc(x)    = :( -inv(abs($x) * sqrt($x^2 - 1))      )
@define_diffrule Base.acot(x)    = :( -inv(1 + $x^2)                      )
@define_diffrule Base.asind(x)   = :(  180 / π / sqrt(1 - $x^2)           )
@define_diffrule Base.acosd(x)   = :( -180 / π / sqrt(1 - $x^2)           )
@define_diffrule Base.atand(x)   = :(  180 / π / (1 + $x^2)               )
@define_diffrule Base.asecd(x)   = :(  180 / π / abs($x) / sqrt($x^2 - 1) )
@define_diffrule Base.acscd(x)   = :( -180 / π / abs($x) / sqrt($x^2 - 1) )
@define_diffrule Base.acotd(x)   = :( -180 / π / (1 + $x^2)               )
@define_diffrule Base.sinh(x)    = :(  cosh($x)                           )
@define_diffrule Base.cosh(x)    = :(  sinh($x)                           )
@define_diffrule Base.tanh(x)    = :(  sech($x)^2                         )
@define_diffrule Base.sech(x)    = :( -tanh($x) * sech($x)                )
@define_diffrule Base.csch(x)    = :( -coth($x) * csch($x)                )
@define_diffrule Base.coth(x)    = :( -(csch($x)^2)                       )
@define_diffrule Base.asinh(x)   = :(  inv(sqrt($x^2 + 1))                )
@define_diffrule Base.acosh(x)   = :(  inv(sqrt($x^2 - 1))                )
@define_diffrule Base.atanh(x)   = :(  inv(1 - $x^2)                      )
@define_diffrule Base.asech(x)   = :( -inv($x * sqrt(1 - $x^2))           )
@define_diffrule Base.acsch(x)   = :( -inv(abs($x) * sqrt(1 + $x^2))      )
@define_diffrule Base.acoth(x)   = :(  inv(1 - $x^2)                      )
@define_diffrule Base.deg2rad(x) = :(  π / 180                            )
@define_diffrule Base.rad2deg(x) = :(  180 / π                            )
@define_diffrule Base.gamma(x)   = :(  digamma($x) * gamma($x)            )
@define_diffrule Base.lgamma(x)  = :(  digamma($x)                        )

# binary #
#--------#

@define_diffrule Base.:+(x, y) = :( 1                  ), :(  1                 )
@define_diffrule Base.:-(x, y) = :( 1                  ), :( -1                 )
@define_diffrule Base.:*(x, y) = :( $y                 ), :(  $x                )
@define_diffrule Base.:/(x, y) = :( inv($y)            ), :( -$x / ($y^2)       )
@define_diffrule Base.:^(x, y) = :( $y * ($x^($y - 1)) ), :(  ($x^$y) * log($x) )

@define_diffrule Base.atan2(x, y) = :( $y / ($x^2 + $y^2)                                      ), :( -$x / ($x^2 + $y^2)                                                     )
@define_diffrule Base.hypot(x, y) = :( $x / hypot($x, $y)                                      ), :(  $y / hypot($x, $y)                                                     )
@define_diffrule Base.mod(x, y)   = :( first(promote(ifelse(isinteger($x / $y), NaN, 1), NaN)) ), :(  z = $x / $y; first(promote(ifelse(isinteger(z), NaN, -floor(z)), NaN)) )
@define_diffrule Base.rem(x, y)   = :( first(promote(ifelse(isinteger($x / $y), NaN, 1), NaN)) ), :(  z = $x / $y; first(promote(ifelse(isinteger(z), NaN, -trunc(z)), NaN)) )

####################
# SpecialFunctions #
####################

# unary #
#-------#

@define_diffrule SpecialFunctions.erf(x)         = :(  (2 / sqrt(π)) * exp(-$x * $x)       )
@define_diffrule SpecialFunctions.erfinv(x)      = :(  (sqrt(π) / 2) * exp(erfinv($x)^2)   )
@define_diffrule SpecialFunctions.erfc(x)        = :( -(2 / sqrt(π)) * exp(-$x * $x)       )
@define_diffrule SpecialFunctions.erfcinv(x)     = :( -(sqrt(π) / 2) * exp(erfcinv($x)^2)  )
@define_diffrule SpecialFunctions.erfi(x)        = :(  (2 / sqrt(π)) * exp($x * $x)        )
@define_diffrule SpecialFunctions.erfcx(x)       = :(  (2 * x * erfcx($x)) - (2 / sqrt(π)) )
@define_diffrule SpecialFunctions.dawson(x)      = :(  1 - (2 * x * dawson($x))            )
@define_diffrule SpecialFunctions.digamma(x)     = :(  trigamma($x)                        )
@define_diffrule SpecialFunctions.invdigamma(x)  = :(  inv(trigamma(invdigamma($x)))       )
@define_diffrule SpecialFunctions.trigamma(x)    = :(  polygamma(2, $x)                    )
@define_diffrule SpecialFunctions.airyai(x)      = :(  airyaiprime($x)                     )
@define_diffrule SpecialFunctions.airyaiprime(x) = :(  $x * airyai($x)                     )
@define_diffrule SpecialFunctions.airybi(x)      = :(  airybiprime($x)                     )
@define_diffrule SpecialFunctions.airybiprime(x) = :(  $x * airybi($x)                     )
@define_diffrule SpecialFunctions.besselj0(x)    = :( -besselj1($x)                        )
@define_diffrule SpecialFunctions.besselj1(x)    = :(  (besselj0($x) - besselj(2, $x)) / 2 )
@define_diffrule SpecialFunctions.bessely0(x)    = :( -bessely1($x)                        )
@define_diffrule SpecialFunctions.bessely1(x)    = :(  (bessely0($x) - bessely(2, $x)) / 2 )

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

@define_diffrule SpecialFunctions.besselj(ν, x)   = :NaN,                                               :(  (besselj($ν - 1, $x) - besselj($ν + 1, $x)) / 2   )
@define_diffrule SpecialFunctions.besseli(ν, x)   = :NaN,                                               :(  (besseli($ν - 1, $x) + besseli($ν + 1, $x)) / 2   )
@define_diffrule SpecialFunctions.bessely(ν, x)   = :NaN,                                               :(  (bessely($ν - 1, $x) - bessely($ν + 1, $x)) / 2   )
@define_diffrule SpecialFunctions.besselk(ν, x)   = :NaN,                                               :( -(besselk($ν - 1, $x) + besselk($ν + 1, $x)) / 2   )
@define_diffrule SpecialFunctions.hankelh1(ν, x)  = :NaN,                                               :(  (hankelh1($ν - 1, $x) - hankelh1($ν + 1, $x)) / 2 )
@define_diffrule SpecialFunctions.hankelh2(ν, x)  = :NaN,                                               :(  (hankelh2($ν - 1, $x) - hankelh2($ν + 1, $x)) / 2 )
@define_diffrule SpecialFunctions.polygamma(m, x) = :NaN,                                               :(  polygamma($m + 1, $x)                             )
@define_diffrule SpecialFunctions.beta(a, b)      = :( beta($a, $b)*(digamma($a) - digamma($a + $b)) ), :(  beta($a, $b)*(digamma($b) - digamma($a + $b))     )
@define_diffrule SpecialFunctions.lbeta(a, b)     = :( digamma($a) - digamma($a + $b)                ), :(  digamma($b) - digamma($a + $b)                    )

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

end # module
