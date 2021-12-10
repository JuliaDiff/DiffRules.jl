
const DEFINED_DIFFRULES = Dict{Tuple{Union{Expr,Symbol},Symbol,Int},Any}()

"""
    @define_diffrule M.f(x) = :(df_dx(\$x))
    @define_diffrule M.f(x, y) = :(df_dx(\$x, \$y)), :(df_dy(\$x, \$y))
    ⋮

Define a new differentiation rule for the function `M.f` and the given arguments, which should
be treated as bindings to Julia expressions. Return the defined rule's key.

The LHS should be a function call with a non-splatted argument list, and the RHS should be
the derivative expression, or in the `n`-ary case, an `n`-tuple of expressions where the
`i`th expression is the derivative of `f` w.r.t the `i`th argument. Arguments should be
interpolated wherever they are used on the RHS.

Note that differentiation rules are purely symbolic, so no type annotations should be used.

# Examples

```julia
@define_diffrule Base.cos(x)          = :(-sin(\$x))
@define_diffrule Base.:/(x, y)        = :(inv(\$y)), :(-\$x / (\$y^2))
@define_diffrule Base.polygamma(m, x) = :NaN,       :(polygamma(\$m + 1, \$x))
```
"""
macro define_diffrule(def)
    @assert isa(def, Expr) && def.head == :(=) "Diff rule expression does not have a left and right side"
    lhs = def.args[1]
    rhs = def.args[2]
    @assert isa(lhs, Expr) && lhs.head == :call "LHS is not a function call"
    qualified_f = lhs.args[1]
    @assert isa(qualified_f, Expr) && qualified_f.head == :(.) "Function is not qualified by module"
    M = qualified_f.args[1]
    f = _get_quoted_symbol(qualified_f.args[2])
    args = lhs.args[2:end]
    rule = Expr(:->, Expr(:tuple, args...), rhs)
    key = Expr(:tuple, Expr(:quote, M), Expr(:quote, f), length(args))
    return esc(quote
        $DiffRules.DEFINED_DIFFRULES[$key] = $rule
        $key
    end)
end

"""
    diffrule(M::Union{Expr,Symbol}, f::Symbol, args...)

Return the derivative expression for `M.f` at the given argument(s), with the argument(s)
interpolated into the returned expression.

In the `n`-ary case, an `n`-tuple of expressions will be returned where the `i`th expression
is the derivative of `f` w.r.t the `i`th argument.

# Examples

```jldoctest
julia> DiffRules.diffrule(:Base, :sin, 1)
:(cos(1))

julia> DiffRules.diffrule(:Base, :sin, :x)
:(cos(x))

julia> DiffRules.diffrule(:Base, :sin, :(x * y^2))
:(cos(x * y ^ 2))
```
"""
diffrule(M::Union{Expr,Symbol}, f::Symbol, args...) = DEFINED_DIFFRULES[M,f,length(args)](args...)

"""
    hasdiffrule(M::Union{Expr,Symbol}, f::Symbol, arity::Int)

Return `true` if a differentiation rule is defined for `M.f` and `arity`, or return `false`
otherwise.

Here, `arity` refers to the number of arguments accepted by `f`.

# Examples

```jldoctest
julia> DiffRules.hasdiffrule(:Base, :sin, 1)
true

julia> DiffRules.hasdiffrule(:Base, :sin, 2)
false

julia> DiffRules.hasdiffrule(:Base, :-, 1)
true

julia> DiffRules.hasdiffrule(:Base, :-, 2)
true

julia> DiffRules.hasdiffrule(:Base, :-, 3)
false
```
"""
hasdiffrule(M::Union{Expr,Symbol}, f::Symbol, arity::Int) = haskey(DEFINED_DIFFRULES, (M, f, arity))

# show a deprecation warning if `filter_modules` in `diffrules()` is specified implicitly
# we use a custom singleton to figure out if the keyword argument was set explicitly
struct DefaultFilterModules end

function deprecated_modules(modules)
    return if modules isa DefaultFilterModules
        Base.depwarn(
            "the implicit keyword argument " *
            "`filter_modules=(:Base, :SpecialFunctions, :NaNMath)` in `diffrules()` is " *
            "deprecated and will be changed to `filter_modules=nothing` in an upcoming " *
            "breaking release of DiffRules (i.e., `diffrules()` will return all rules " *
            "defined in DiffRules)",
            :diffrules,
        )
        (:Base, :SpecialFunctions, :NaNMath)
    else
        modules
    end
end

"""
    diffrules(; filter_modules=(:Base, :SpecialFunctions, :NaNMath))

Return a list of keys that can be used to access all defined differentiation rules for
modules in `filter_modules`.

Each key is of the form `(M::Symbol, f::Symbol, arity::Int)`.
Here, `arity` refers to the number of arguments accepted by `f` and `M` is one of the
modules in `filter_modules`.

To include all rules, specify `filter_modules = nothing`.

!!! note
    Calling `diffrules()` with the implicit default keyword argument `filter_modules`
    does *not* return all rules defined by this package but rather only rules for the
    packages for which DiffRules 1.0 provided rules. This is done in order to not to
    break downstream packages that assumed this list would never change.
    It is planned to change `diffrules()` to return all rules, i.e., to use the
    default keyword argument `filter_modules=nothing`, in an upcoming breaking release
    of DiffRules.

# Examples

```jldoctest
julia> first(DiffRules.diffrules()) isa Tuple{Symbol,Symbol,Int}
true

julia> (:Base, :log, 1) in DiffRules.diffrules()
true

julia> (:Base, :*, 2) in DiffRules.diffrules()
true
```

If you call `diffrules()`, only rules for Base, SpecialFunctions, and
NaNMath are returned but no rules for LogExpFunctions:
```jldoctest
julia> any(M === :LogExpFunctions for (M, _, _) in DiffRules.diffrules())
false
```

If you set `filter_modules=nothing`, all rules defined in DiffRules are
returned and in particular also rules for LogExpFunctions:
```jldoctest
julia> any(
           M === :LogExpFunctions
           for (M, _, _) in DiffRules.diffrules(; filter_modules=nothing)
       )
true
```

If you set `filter_modules=(:Base,)` only rules for functions in Base are
returned:
```jldoctest
julia> all(M === :Base for (M, _, _) in DiffRules.diffrules(; filter_modules=(:Base,)))
true
```
"""
function diffrules(; filter_modules=DefaultFilterModules())
    modules = deprecated_modules(filter_modules)
    return if modules === nothing
        keys(DEFINED_DIFFRULES)
    else
        Iterators.filter(keys(DEFINED_DIFFRULES)) do (M, _, _)
            return M in modules
        end
    end
end

# For v0.6 and v0.7 compatibility, need to support having the diff rule function enter as a
# `Expr(:quote...)` and a `QuoteNode`. When v0.6 support is dropped, the function will
# always enter in a `QuoteNode` (#23885).
function _get_quoted_symbol(ex::Expr)
    @assert ex.head == :quote
    @assert length(ex.args) == 1 && isa(ex.args[1], Symbol) "Function not a single symbol"
    ex.args[1]
end

function _get_quoted_symbol(ex::QuoteNode)
    @assert isa(ex.value, Symbol) "Function not a single symbol"
    ex.value
end
