module DiffRulesNaNMathExt

using DiffRules: @define_diffrule
using IrrationalConstants: logtwo, logten
using NaNMath

###########
# NaNMath #
###########

# unary #
#-------#

@define_diffrule NaNMath.sqrt(x)   = :(  inv(2 * $(NaNMath.sqrt)($x))                  )
@define_diffrule NaNMath.sin(x)    = :(  $(NaNMath.cos)($x)                            )
@define_diffrule NaNMath.cos(x)    = :( -$(NaNMath.sin)($x)                            )
@define_diffrule NaNMath.tan(x)    = :(  1 + $(NaNMath.pow)($(NaNMath.tan)($x), 2)        )
@define_diffrule NaNMath.asin(x)   = :(  inv($(NaNMath.sqrt)(1 - $(NaNMath.pow)($x, 2)))  )
@define_diffrule NaNMath.acos(x)   = :(  -inv($(NaNMath.sqrt)(1 - $(NaNMath.pow)($x, 2))) )
@define_diffrule NaNMath.acosh(x)  = :(  inv($(NaNMath.sqrt)($(NaNMath.pow)($x, 2) - 1))  )
@define_diffrule NaNMath.atanh(x)  = :(  inv(1 - $(NaNMath.pow)($x, 2))                )
@define_diffrule NaNMath.log(x)    = :(  inv($x)                                    )
@define_diffrule NaNMath.log2(x)   = :(  inv($logtwo * $x)                          )
@define_diffrule NaNMath.log10(x)  = :(  inv($logten * $x)                          )
@define_diffrule NaNMath.log1p(x)  = :(  inv($x + 1)                                )


# binary #
#--------#

@define_diffrule NaNMath.pow(x, y) = :( $y * $(NaNMath.pow)($x, ($y - 1)) ), :( $(NaNMath.pow)($x, $y) * $(NaNMath.log)($x) )
@define_diffrule NaNMath.max(x, y) = :(ifelse(($y > $x) | (signbit($y) < signbit($x)), ifelse(isnan($y), one($x), zero($x)), ifelse(isnan($x), zero($x), one($x)))),
                                     :(ifelse(($y > $x) | (signbit($y) < signbit($x)), ifelse(isnan($y), zero($y), one($y)), ifelse(isnan($x), one($y), zero($y))))
@define_diffrule NaNMath.min(x, y) = :(ifelse(($y < $x) | (signbit($y) > signbit($x)), ifelse(isnan($y), one($x), zero($x)), ifelse(isnan($x), zero($x), one($x)))),
                                     :(ifelse(($y < $x) | (signbit($y) > signbit($x)), ifelse(isnan($y), zero($y), one($y)), ifelse(isnan($x), one($x), zero($x))))
end # module
