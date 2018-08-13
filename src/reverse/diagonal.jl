import LinearAlgebra: det, logdet, diagm, Diagonal, diag
export diag, diagm, Diagonal

push!(ops, DiffOp(:(LinearAlgebra.diag), :(Tuple{DLA.AM}), [true]))
function ∇(::typeof(diag), ::Arg1, p, y::AV, ȳ::AV, X::AM)
    X̄ = fill!(similar(X), zero(eltype(X)))
    X̄[diagind(X̄)] = ȳ
    return X̄
end
function ∇(X̄::AM, ::typeof(diag), ::Arg1, p, y::AV, ȳ::AV, X::AM)
    X̄_diag = view(X̄, diagind(X̄))
    X̄_diag .+= ȳ
    return X̄
end

push!(ops, DiffOp(:(LinearAlgebra.diag), :(Tuple{DLA.AM, Integer}), [true, false]))
function ∇(::typeof(diag), ::Arg1, p, y::AV, ȳ::AV, X::AM, k::Integer)
    X̄ = fill!(similar(X), zero(eltype(X)))
    X̄[diagind(X̄, k)] = ȳ
    return X̄
end
function ∇(X̄::AM, ::typeof(diag), ::Arg1, p, y::AV, ȳ::AV, X::AM, k::Integer)
    X̄_diag = view(X̄, diagind(X̄, k))
    X̄_diag .+= ȳ
    return X̄
end

push!(ops, DiffOp(:(LinearAlgebra.diagm), :(Tuple{Pair{<:Integer, <:AV}}), [true]))
∇(::typeof(diagm), ::Arg1, p, Y::AM, Ȳ::AM, x::Pair{<:Integer, <:AV}) =
    copyto!(similar(x.second), view(Ȳ, diagind(Ȳ, x.first)))
∇(x̄::AV, ::typeof(diagm), ::Arg1, p, Y::AM, Ȳ::AM, x::Pair{<:Integer, <:AV}) =
    broadcast!(+, x̄, x̄, view(Ȳ, diagind(Ȳ, x.first)))

push!(ops, DiffOp(:(LinearAlgebra.Diagonal), :(Tuple{DLA.AV}), [true]))
∇(::Type{Diagonal}, ::Arg1, p, Y::Diagonal{<:Real}, Ȳ::Diagonal{<:Real}, x::AV) =
    copyto!(similar(x), Ȳ.diag)
∇(x̄::AV, ::Type{Diagonal}, ::Arg1, p, Y::Diagonal{<:Real}, Ȳ::Diagonal{<:Real}, x::AV) =
    broadcast!(+, x̄, x̄, Ȳ.diag)

push!(ops, DiffOp(:(LinearAlgebra.Diagonal), :(Tuple{DLA.AM}), [true]))
function ∇(::Type{Diagonal}, ::Arg1, p, Y::Diagonal{<:Real}, Ȳ::Diagonal{<:Real}, X::AM)
    X̄ = zero(X)
    copyto!(view(X̄, diagind(X)), Ȳ.diag)
    return X̄
end
function ∇(X̄::AM, ::Type{Diagonal}, ::Arg1, p, Y::Diagonal{<:Real}, Ȳ::Diagonal{<:Real}, X::AM)
    X̄_diag = view(X̄, diagind(X̄))
    broadcast!(+, X̄_diag, X̄_diag, Ȳ.diag)
    return X̄
end

push!(ops, DiffOp(:(LinearAlgebra.det), :(Tuple{Diagonal{<:Real}}), [true]))
∇(::typeof(det), ::Arg1, p, y::Real, ȳ::Real, X::Diagonal{<:Real}) =
    Diagonal(ȳ .* y ./ X.diag)
function ∇(X̄::Diagonal{<:Real}, ::typeof(det), ::Arg1, p, y::Real, ȳ::Real, X::Diagonal{<:Real})
    broadcast!((x̄, x, y, ȳ)->x̄ + ȳ * y / x, X̄.diag, X̄.diag, X.diag, y, ȳ)
    return X̄
end

push!(ops, DiffOp(:(LinearAlgebra.logdet), :(Tuple{Diagonal{<:Real}}), [true]))
∇(::typeof(logdet), ::Arg1, p, y::Real, ȳ::Real, X::Diagonal{<:Real}) = Diagonal(ȳ ./ X.diag)
function ∇(X̄::Diagonal{<:Real}, ::typeof(logdet), ::Arg1, p, y::Real, ȳ::Real, X::Diagonal{<:Real})
    broadcast!((x̄, x, ȳ)->x̄ + ȳ / x, X̄.diag, X̄.diag, X.diag, ȳ)
    return X̄
end
