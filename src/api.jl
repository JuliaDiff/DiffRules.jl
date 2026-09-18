"""
    @define_diffrule M.f(x) = :(df_dx(\$x))
    @define_diffrule M.f(x, y) = :(df_dx(\$x, \$y)), :(df_dy(\$x, \$y))
    ⋮

Define a new differentiation rule for the function `M.f` and the given arguments, which should
be treated as bindings to Julia expressions.

The LHS should be a function call with a non-splatted argument list, and the RHS should be
the derivative expression, or in the `n`-ary case, an `n`-tuple of expressions where the
`i`th expression is the derivative of `f` w.r.t the `i`th argument. Arguments should be
interpolated wherever they are used on the RHS.

Note that differentiation rules are purely symbolic, so no type annotations should be used.

Each rule is a method of [`diffrule`](@ref), so rules defined in other packages and in
package extensions are precompiled and visible like any other method.

# Examples

```julia
@define_diffrule Base.cos(x)          = :(-sin(\$x))
@define_diffrule Base.:/(x, y)        = :(inv(\$y)), :(-\$x / (\$y^2))
@define_diffrule Base.polygamma(m, x) = :NaN,       :(polygamma(\$m + 1, \$x))
```
"""
macro define_diffrule(def)
    @assert isa(def, Expr) && def.head == :(=) "Diff rule expression does not have a left and right side"
    lhs, rhs = def.args
    @assert isa(lhs, Expr) && lhs.head == :call "LHS is not a function call"
    f = lhs.args[1]
    @assert isa(f, Expr) && f.head == :(.) "Function is not qualified by module"
    args = lhs.args[2:end]
    return esc(:($DiffRules.diffrule(::typeof($f), $(args...)) = $rhs))
end

"""
    diffrule(f, args...)
    diffrule(M::Union{Module,Symbol}, f::Symbol, args...)

Return the derivative expression for `f` at the given argument(s), with the argument(s)
interpolated into the returned expression.

In the `n`-ary case, an `n`-tuple of expressions will be returned where the `i`th expression
is the derivative of `f` w.r.t the `i`th argument.

# Examples

```jldoctest
julia> DiffRules.diffrule(sin, 1)
:(cos(1))

julia> DiffRules.diffrule(Base, :sin, :x)
:(cos(x))

julia> DiffRules.diffrule(:Base, :sin, :(x * y^2))
:(cos(x * y ^ 2))
```
"""
function diffrule end

diffrule(M::Module, f::Symbol, args...) = diffrule(getproperty(M, f), args...)

function diffrule(M::Symbol, f::Symbol, args...)
    fn = _resolve(M, f)
    fn === nothing && throw(KeyError((M, f, length(args))))
    return diffrule(fn, args...)
end

"""
    hasdiffrule(f, arity::Int)
    hasdiffrule(M::Union{Module,Symbol}, f::Symbol, arity::Int)

Return `true` if a differentiation rule is defined for `f` and `arity`, or return `false`
otherwise. Here, `arity` refers to the number of arguments accepted by `f`.

Rules for a package's functions exist only once that package is loaded, so a query for an
unloaded package returns `false`.

# Examples

```jldoctest
julia> DiffRules.hasdiffrule(sin, 1)
true

julia> DiffRules.hasdiffrule(sin, 2)
false

julia> DiffRules.hasdiffrule(:Base, :-, 2)
true
```
"""
hasdiffrule(f, arity::Int) = hasmethod(diffrule, Tuple{typeof(f),Vararg{Any,arity}})

hasdiffrule(M::Module, f::Symbol, arity::Int) =
    isdefined(M, f) && hasdiffrule(getproperty(M, f), arity)

function hasdiffrule(M::Symbol, f::Symbol, arity::Int)
    fn = _resolve(M, f)
    return fn !== nothing && hasdiffrule(fn, arity)
end

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
modules in `filter_modules`. To include all rules, specify `filter_modules = nothing`.

Each key is of the form `(M::Symbol, f::Symbol, arity::Int)`, where `M` is the name of the
package defining `f` and `arity` is the number of arguments accepted by `f`.

Keys are collected from the method table of [`diffrule`](@ref), so rules defined in other
packages are included. Rules for a package's functions exist only once that package is
loaded: querying before `using SpecialFunctions` will not list its rules.

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

julia> all(M === :Base for (M, _, _) in DiffRules.diffrules(; filter_modules=(:Base,)))
true
```
"""
function diffrules(; filter_modules=DefaultFilterModules())
    modules = deprecated_modules(filter_modules)
    rules = [(_pkgname(fn), nameof(fn), arity) for (fn, arity) in _rules()]
    modules === nothing && return rules
    return filter(r -> r[1] in modules, rules)
end

# `parentmodule` reports submodules such as `Base.Math`, which callers cannot splice into
# `M.f`; the root module is the package name they expect.
_pkgname(fn) = nameof(Base.moduleroot(parentmodule(fn)))

# Rules are the methods of `diffrule` whose first parameter is a singleton function type.
# The `Symbol`/`Module` methods above are not, and are skipped.
function _rules()
    rules = Tuple{Function,Int}[]
    for m in methods(diffrule)
        m.isva && continue
        params = Base.unwrap_unionall(m.sig).parameters
        length(params) >= 2 || continue
        T = params[2]
        T isa DataType && T <: Function && isdefined(T, :instance) || continue
        push!(rules, (T.instance, length(params) - 2))
    end
    return rules
end

function _resolve(M::Symbol, f::Symbol)
    for (fn, _) in _rules()
        nameof(fn) === f && _pkgname(fn) === M && return fn
    end
    return nothing
end
