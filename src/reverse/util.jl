"""
    DiffOp

The information associated with a particular differentiable linear algebra operation. `f` is
either a `Symbol` or `Expr` containing the function name, and `T` is an expression
corresponding to the tuple-type of the arguments of the function. `diff_flags` is a vector
containing flags indicating whether each argument of the function is differentiable or not.
"""
struct DiffOp
    f::Union{Symbol, Expr}
    T::Expr
    diff_flags::Vector{Bool}
end
ops = Set{DiffOp}()

"""
    importable(ex::Expr)

Construct an expression of the form `:(import Package.Subpackage.Foo)` from an expression of
the form `:(Package.Subpackage.Foo)`.
"""
function importable(ex::Expr)
    ex.head === :. || error("Expression is not valid as an import: $ex")
    result = importable(ex.args[1])
    push!(result, ex.args[2].value)
    return result
end
importable(sym::Symbol) = Any[sym]

"""
    import_expr(dop::DiffOp)

Generate an expression to import `dop.f` from the appropriate package.
"""
import_expr(dop::DiffOp) =
    VERSION <= VersionNumber("0.6.2") ?
        Expr(:import, importable(dop.f)...) :
        Expr(:import, Expr(Symbol("."), importable(dop.f)...))
