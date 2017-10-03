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
    rhs.args[1] = :(::Type{$DiffRules.DiffRule{$(Expr(:quote, M)),$(Expr(:quote, f))}})
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

    julia> DiffResults.diffrule(:Base, :sin, 1)
    :(cos(1))

    julia> DiffResults.diffrule(:Base, :sin, :x)
    :(cos(x))

    julia> DiffResults.diffrule(:Base, :sin, :(x * y^2))
    :(cos(x * y ^ 2))

    julia> DiffResults.diffrule(:Base, :^, :(x + 2), :c)
    (:(c * (x + 2) ^ (c - 1)), :((x + 2) ^ c * log(x + 2)))
"""
diffrule(M::Symbol, f::Symbol, args...) = DiffRule{M,f}(args...)

"""
    hasdiffrule(M::Symbol, f::Symbol, arity::Int)

Return `true` if a differentiation rule is defined for `M.f` and `arity`, or return `false`
otherwise.

Here, `arity` refers to the number of arguments accepted by `f`.

Examples:

    julia> DiffResults.hasdiffrule(:Base, :sin, 1)
    true

    julia> DiffResults.hasdiffrule(:Base, :sin, 2)
    false

    julia> DiffResults.hasdiffrule(:Base, :-, 1)
    true

    julia> DiffResults.hasdiffrule(:Base, :-, 2)
    true

    julia> DiffResults.hasdiffrule(:Base, :-, 3)
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
