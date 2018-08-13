import LinearAlgebra: det, logdet, LowerTriangular, UpperTriangular
export det, logdet, LowerTriangular, UpperTriangular

for (ctor, ctor_sym, T, T_sym) in zip([:LowerTriangular, :UpperTriangular],
                               [:(:(LinearAlgebra.LowerTriangular)), :(:(LinearAlgebra.UpperTriangular))],
                               [:(LowerTriangular{<:Real}), :(UpperTriangular{<:Real})],
                               [:(:(LowerTriangular{<:Real})), :(:(UpperTriangular{<:Real}))])

    @eval begin

    push!(ops, DiffOp($ctor_sym, :(Tuple{DLA.AM}), [true]))
    ∇(::Type{$ctor}, ::Arg1, p, Y::$T, Ȳ::$T, X::AM) = Matrix(Ȳ)
    ∇(X̄::AM, ::Type{$ctor}, ::Arg1, p, Y::$T, Ȳ::$T, X::AM) = broadcast!(+, X̄, X̄, Ȳ)

    push!(ops, DiffOp(:(LinearAlgebra.det), Expr(:curly, :Tuple, $T_sym), [true]))
    ∇(::typeof(det), ::Arg1, p, y::Real, ȳ::Real, X::$T) =
        Diagonal(ȳ .* y ./ view(X, diagind(X)))

    # Optimisation for in-place updates.
    function ∇(X̄::AM, ::typeof(det), ::Arg1, p, y::Real, ȳ::Real, X::$T)
        X̄_diag = view(X̄, diagind(X̄))
        broadcast!((x̄, x, y, ȳ)->x̄ + ȳ * y / x,
                   X̄_diag, X̄_diag, view(X, diagind(X)), y, ȳ)
        return X̄
    end

    # Optimisation for in-place updates to `Diagonal` sensitivity cache.
    function ∇(X̄::Diagonal, ::typeof(det), ::Arg1, p, y::Real, ȳ::Real, X::$T)
        X̄.diag .+= ȳ .* y ./ view(X, diagind(X))
        return X̄
    end

    push!(ops, DiffOp(:(LinearAlgebra.logdet), Expr(:curly, :Tuple, $T_sym), [true]))
    ∇(::typeof(logdet), ::Arg1, p, y::Real, ȳ::Real, X::$T) =
        Diagonal(ȳ ./ view(X, diagind(X)))

    # Optimisation for in-place updates.
    function ∇(X̄::AM, ::typeof(logdet), ::Arg1, p, y::Real, ȳ::Real, X::$T)
        X̄_diag = view(X̄, diagind(X̄))
        broadcast!((x̄, x, ȳ)->x̄ + ȳ / x, X̄_diag, X̄_diag, view(X, diagind(X)), ȳ)
        return X̄
    end

    # Optimisation for in-place updates to `Diagonal` sensitivity cache.
    function ∇(X̄::Diagonal, ::typeof(logdet), ::Arg1, p, y::Real, ȳ::Real, X::$T)
        X̄.diag .+= ȳ ./ view(X, diagind(X))
        return X̄
    end

    end
end
