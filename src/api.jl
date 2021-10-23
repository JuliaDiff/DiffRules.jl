
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

"""
    diffrules(; modules=(:Base, :SpecialFunctions, :NaNMath))

Return a list of keys that can be used to access all defined differentiation rules for
functions in the `modules`.

Each key is of the form `(M::Symbol, f::Symbol, arity::Int)`.

Here, `arity` refers to the number of arguments accepted by `f` and `M` is one of the
`modules`.

The default `modules` does *not* include all rules defined by this package, but rather, exactly those packages for which `v1.0` provided rules. This is done in order not to break downstream packages, man or which assumed this list would never change. To include all rules, specify `modules = :all`.

# Examples

```jldoctest
julia> modules = Set(M for (M, f, arity) in DiffRules.diffrules());

julia> modules == Set((:Base, :SpecialFunctions, :NaNMath))
true

julia> modules = Set(M for (M, f, arity) in DiffRules.diffrules(; modules=(:Base,)));

julia> modules == Set((:Base,))
true

julia> isempty(DiffRules.diffrules(; modules=(:StatsFuns,)))
true
```
"""
function diffrules(; modules=(:Base, :SpecialFunctions, :NaNMath))
    return Iterators.filter(keys(DEFINED_DIFFRULES)) do (M, _, _)
        return M in modules
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
