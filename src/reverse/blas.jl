import LinearAlgebra: dot
import LinearAlgebra.BLAS: asum, blascopy!, nrm2, scal, scal!, gemm, gemm!, gemv, gemv!,
    syrk, symm, symm!, symv, symv!, trmm, trsm, trmv, trsv, trsv!, ger!

################################## Level 1 ##################################

# Unit-stride `dot`.
push!(ops, DiffOp(:(LinearAlgebra.dot),
    :(Tuple{DLA.SA{<:DLA.BF}, DLA.SA{<:DLA.BF}}),
    [true, true]
))
∇(::typeof(LinearAlgebra.dot), ::Arg1, p, z::BF, z̄::BF, x::SA{<:BF}, y::SA{<:BF}) = z̄ .* y
∇(::typeof(LinearAlgebra.dot), ::Arg2, p, z::BF, z̄::BF, x::SA{<:BF}, y::SA{<:BF}) = z̄ .* x
function ∇(x̄, ::typeof(LinearAlgebra.dot), ::Arg1, p, z::BF, z̄::BF, x::SA{<:BF}, y::SA{<:BF})
    x̄ .= x̄ .+ z̄ .* y
    return x̄
end
function ∇(ȳ, ::typeof(LinearAlgebra.dot), ::Arg2, p, z::BF, z̄::BF, x::SA{<:BF}, y::SA{<:BF})
    ȳ .= ȳ .+ z̄ .* x
    return ȳ
end

# Arbitrary-stride `dot`.
push!(ops, DiffOp(:(LinearAlgebra.BLAS.dot),
    :(Tuple{Int, DLA.SA{<:DLA.BF}, Int, DLA.SA{<:DLA.BF}, Int}),
    [false, true, false, true, false],
))
∇(::typeof(BLAS.dot), ::Arg2, p, z::BF, z̄::BF, n::Int, x::SA{<:BF}, ix::Int, y::SA{<:BF}, iy::Int) =
    scal!(n, z̄, blascopy!(n, y, iy, fill!(similar(x), zero(eltype(x))), ix), ix)
∇(::typeof(BLAS.dot), ::Arg4, p, z::BF, z̄::BF, n::Int, x::SA{<:BF}, ix::Int, y::SA{<:BF}, iy::Int) =
    scal!(n, z̄, blascopy!(n, x, ix, fill!(similar(y), zero(eltype(y))), iy), iy)
function ∇(x̄, ::typeof(BLAS.dot), ::Arg2, p, z::BF, z̄::BF, n::Int, x::SA{<:BF}, ix::Int, y::SA{<:BF}, iy::Int)
    x̄ .= x̄ .+ scal!(n, z̄, blascopy!(n, y, iy, fill!(similar(x), zero(eltype(x))), ix), ix)
    return x̄
end
function ∇(ȳ, ::typeof(BLAS.dot), ::Arg4, p, z::BF, z̄::BF, n::Int, x::SA{<:BF}, ix::Int, y::SA{<:BF}, iy::Int)
    ȳ .= ȳ .+ scal!(n, z̄, blascopy!(n, x, ix, fill!(similar(y), zero(eltype(y))), iy), iy)
    return ȳ
end

# Unit-stride `nrm2`.
push!(ops, DiffOp(:(LinearAlgebra.BLAS.nrm2),
    :(Tuple{DLA.SA{<:DLA.BF}}),
    [true]
))
∇(::typeof(nrm2), ::Arg1, p, y::BF, ȳ::BF, x::SA{<:BF}) = x .* (ȳ / y)
function ∇(x̄::AA, ::typeof(nrm2), ::Arg1, p, y::BF, ȳ::BF, x::SA{<:BF})
    x̄ .= x̄ .+ x .* (ȳ / y)
    return x̄
end

# Arbitrary-stride `nrm2`.
push!(ops, DiffOp(:(LinearAlgebra.BLAS.nrm2),
    :(Tuple{Integer, DLA.SA{<:DLA.BF}, Integer}),
    [false, true, false]
))
∇(::typeof(nrm2), ::Arg2, p, y::BF, ȳ::BF, n::Integer, x::SA{<:BF}, inc::Integer) =
    scal!(n, ȳ / y, blascopy!(n, x, inc, fill!(similar(x), zero(eltype(x))), inc), inc)
function ∇(x̄::SA{<:BF}, ::typeof(nrm2), ::Arg2, p, y::BF, ȳ::BF, n::Integer, x::SA{<:BF}, inc::Integer)
    x̄ .= x̄ .+ scal!(n, ȳ / y, blascopy!(n, x, inc, fill!(similar(x), zero(eltype(x))), inc), inc)
    return x̄
end

# Unit-stride `asum`.
push!(ops, DiffOp(:(LinearAlgebra.BLAS.asum),
    :(Tuple{DLA.SA{<:DLA.BF}}),
    [true]
))
∇(::typeof(asum), ::Arg1, p, y::BF, ȳ::BF, x::SA{<:BF}) = ȳ .* sign.(x)
function ∇(x̄::AA, ::typeof(asum), ::Arg1, p, y::BF, ȳ::BF, x::SA{<:BF})
    x̄ .= x̄ .+ ȳ .* sign.(x)
    return x̄
end

# Arbitrary-stride `asum`.
push!(ops, DiffOp(:(LinearAlgebra.BLAS.asum),
    :(Tuple{Integer, DLA.SA{<:DLA.BF}, Integer}),
    [false, true, false]
))
∇(::typeof(asum), ::Arg2, p, y::BF, ȳ::BF, n::Integer, x, inc::Integer) =
    scal!(n, ȳ, blascopy!(n, sign.(x), inc, fill!(similar(x), zero(eltype(x))), inc), inc)
function ∇(x̄::SA{<:BF}, ::typeof(asum), ::Arg2, p, y::BF, ȳ::BF, n::Integer, x::SA{<:BF}, inc::Integer)
    x̄ .= x̄ .+ scal!(n, ȳ, blascopy!(n, sign.(x), inc, fill!(similar(x), zero(eltype(x))), inc), inc)
    return x̄
end
# Some weird stuff going on that I haven't figured out yet. This is a very old attempt.
# let f = :(scal{T <: AbstractArray, V <: AbstractFloat})
#     ā = :(blascopy!(n, z̄, inc, zeros(X), inc) .* X)
#     X̄ = :(scal!(n, a, z̄, inc))
# end


################################## Level 2 ##################################

# `gemv` sensitivities implementation.
push!(ops, DiffOp(:(LinearAlgebra.BLAS.gemv),
    :(Tuple{Char, T, DLA.SM{T}, DLA.SV{T}} where T<:DLA.BF),
    [false, true, true, true],
))
∇(::typeof(gemv), ::Arg2, _, y::SV{T}, ȳ::SV, tA::Char, α::T, A::SM{T}, x::SV{T}) where T<:BF =
    dot(ȳ, y) / α
∇(::typeof(gemv), ::Arg3, _, y::SV{T}, ȳ::SV, tA::Char, α::T, A::SM{T}, x::SV{T}) where T<:BF =
    uppercase(tA) == 'N' ? α * ȳ * x' : α * x * ȳ'
∇(Ā::SM{T}, ::typeof(gemv), ::Arg3, _, y::SV{T}, ȳ::SV{T}, tA::Char, α::T, A::SM{T}, x::SV{T}) where T<:BF =
    uppercase(tA) == 'N' ? ger!(α, ȳ, x, Ā) : ger!(α, x, ȳ, Ā)
∇(::typeof(gemv), ::Arg4, _, y::SV{T}, ȳ::SV{T}, tA::Char, α::T, A::SM{T}, x::SV{T}) where T<:BF =
    gemv(uppercase(tA) == 'N' ? 'T' : 'N', α, A, ȳ)
∇(x̄::SV{T}, ::typeof(gemv), ::Arg4, _, y::SV{T}, ȳ::SV{T}, tA::Char, α::T, A::SM{T}, x::SV{T}) where T<:BF =
    gemv!(uppercase(tA) == 'N' ? 'T' : 'N', α, A, ȳ, one(T), x̄)

# `gemv` sensitivities implementation with `α = 1`.
push!(ops, DiffOp(:(LinearAlgebra.BLAS.gemv),
    :(Tuple{Char, DLA.SM{T}, DLA.SV{T}} where T<:DLA.BF),
    [false, true, true],
))
∇(::typeof(gemv), ::Arg2, p, y::SV{T}, ȳ::SV{T}, tA::Char, A::SM{T}, x::SV{T}) where T<:BF =
    ∇(gemv, Val{3}, p, y, ȳ, tA, one(T), A, x)
∇(Ā::SM{T}, ::typeof(gemv), ::Arg2, p, y::SV{T}, ȳ::SV{T}, tA::Char, A::SM{T}, x::SV{T}) where T<:BF =
    ∇(Ā, gemv, Val{3}, p, y, ȳ, tA, one(T), A, x)
∇(::typeof(gemv), ::Arg3, p, y::SV{T}, ȳ::SV{T}, tA::Char, A::SM{T}, x::SV{T}) where T<:BF =
    ∇(gemv, Val{4}, p, y, ȳ, tA, one(T), A, x)
∇(x̄::SV{T}, ::typeof(gemv), ::Arg3, p, y::SV{T}, ȳ::SV{T}, tA::Char, A::SM{T}, x::SV{T}) where T<:BF =
    ∇(x̄, gemv, Val{4}, p, y, ȳ, tA, one(T), A, x)

# `symv` sensitivity implementations.
push!(ops, DiffOp(:(LinearAlgebra.BLAS.symv),
    :(Tuple{Char, T, DLA.SM{T}, DLA.SV{T}} where T<:DLA.BF),
    [false, true, true, true],
))
∇(::typeof(symv), ::Arg2, p, y::SV{T}, ȳ::SV{T}, ul::Char, α::T, A::SM{T}, x::SV{T}) where T<:BF =
    dot(ȳ, y) / α
function ∇(::typeof(symv), ::Arg3, p, y::SV{T}, ȳ::SV{T}, ul::Char, α::T, A::SM{T}, x::SV{T}) where T<:BF
    Y, Ȳ, X = reshape(y, length(y), 1), reshape(ȳ, length(ȳ), 1), reshape(x, length(x), 1)
    return ∇(symm, Val{4}, p, Y, Ȳ, 'L', ul, α, A, X)
end
function ∇(Ā::SM{T}, ::typeof(symv), ::Arg3, p, y::SV{T}, ȳ::SV{T}, ul::Char, α::T, A::SM{T}, x::SV{T}) where T<:BF
    Y, Ȳ, X = reshape(y, length(y), 1), reshape(ȳ, length(ȳ), 1), reshape(x, length(x), 1)
    return ∇(Ā, symm, Val{4}, p, Y, Ȳ, 'L', ul, α, A, X)
end
∇(::typeof(symv), ::Arg4, p, y::SV{T}, ȳ::SV{T}, ul::Char, α::T, A::SM{T}, x::SV{T}) where T<:BF =
    symv(ul, α, A, ȳ)
∇(x̄::SV{T}, ::typeof(symv), ::Arg4, p, y::SV{T}, ȳ::SV{T}, ul::Char, α::T, A::SM{T}, x::SV{T}) where T<:BF =
    symv!(ul, α, A, ȳ, one(T), x̄)

# `symv` sensitivity implementations for `α=1`.
push!(ops, DiffOp(:(LinearAlgebra.BLAS.symv),
    :(Tuple{Char, DLA.SM{T}, DLA.SV{T}} where T<:DLA.BF),
    [false, true, true],
))
∇(::typeof(symv), ::Arg2, p, y::SV{T}, ȳ::SV{T}, ul::Char, A::SM{T}, x::SV{T}) where T<:BF =
    ∇(symv, Val{3}, p, y, ȳ, ul, one(T), A, x)
∇(Ā::SM{T}, ::typeof(symv), ::Arg2, p, y::SV{T}, ȳ::SV{T}, ul::Char, A::SM{T}, x::SV{T}) where T<:BF =
    ∇(Ā, symv, Val{3}, p, y::SV{T}, ȳ::SV{T}, ul, one(T), A, x)
∇(::typeof(symv), ::Arg3, p, y::SV{T}, ȳ::SV{T}, ul::Char, A::SM{T}, x::SV{T}) where T<:BF =
    ∇(symv, Val{4}, p, y, ȳ, ul, one(T), A, x)
∇(B̄::SV{T}, ::typeof(symv), ::Arg3, p, y::SV{T}, ȳ::SV{T}, ul::Char, A::SM{T}, x::SV{T}) where T<:BF =
    ∇(B̄, symv, Val{4}, p, y, ȳ, ul, one(T), A, x)

# `trmv` sensitivity implementations.
push!(ops, DiffOp(:(LinearAlgebra.BLAS.trmv),
    :(Tuple{Char, Char, Char, DLA.SM{T}, DLA.SV{T}} where T<:DLA.BF),
    [false, false, false, true, true],
))
function ∇(::typeof(trmv), ::Arg4, p, y::SV{T}, ȳ::SV{T},
    ul::Char, ta::Char, dA::Char,
    A::SM{T},
    b::SV{T},
) where T<:BF
    Ā = (uppercase(ul) == 'L' ? tril! : triu!)(uppercase(ta) == 'N' ? ȳ * b' : b * ȳ')
    dA == 'U' && fill!(view(Ā, diagind(Ā)), zero(T))
    return Ā
end
∇(::typeof(trmv), ::Arg5, p, y::SV{T}, ȳ::SV{T},
    ul::Char, ta::Char, dA::Char,
    A::SM{T},
    b::SV{T},
) where T<:BF = trmv(ul, uppercase(ta) == 'N' ? 'T' : 'N', dA, A, ȳ)

# `trsv` sensitivity implementations.
push!(ops, DiffOp(:(LinearAlgebra.BLAS.trsv),
    :(Tuple{Char, Char, Char, DLA.SM{T}, DLA.SV{T}} where T<:DLA.BF),
    [false, false, false, true, true],
))
function ∇(::typeof(trsv), ::Arg4, p, y::SV{T}, ȳ::SV{T},
    ul::Char, ta::Char, dA::Char,
    A::SM{T},
    x::SV{T},
) where T<:BF
    Y, Ȳ, X = reshape(y, length(y), 1), reshape(ȳ, length(ȳ), 1), reshape(x, length(x), 1)
    Ā = ∇(trsm, Val{6}, p, Y, Ȳ, 'L', ul, ta, dA, one(T), A, X)
    dA == 'U' && fill!(view(Ā, diagind(Ā)), zero(T))
    return Ā
end
∇(::typeof(trsv), ::Arg5, p, y::SV{T}, ȳ::SV{T},
    ul::Char, ta::Char, dA::Char,
    A::SM{T},
    x::SV{T},
) where T<:BF = trsv(ul, uppercase(ta) == 'N' ? 'T' : 'N', dA, A, ȳ)


################################## Level 3 ##################################

# `gemm` sensitivities implementation.
push!(ops, DiffOp(:(LinearAlgebra.BLAS.gemm),
    :(Tuple{Char, Char, T, DLA.SM{T}, DLA.SM{T}} where T<:DLA.BF),
    [false, false, true, true, true],
))
∇(::typeof(gemm), ::Arg3, p, Y::SM{T}, Ȳ::SM{T},
    tA::Char,
    tB::Char,
    α::T,
    A::SM{T},
    B::SM{T},
) where T<:BF = sum(Ȳ .* Y) / α
∇(::typeof(gemm), ::Arg4, p, Y::SM{T}, Ȳ::SM{T},
    tA::Char,
    tB::Char,
    α::T,
    A::SM{T},
    B::SM{T},
) where T<:BF =
    uppercase(tA) == 'N' ?
        uppercase(tB) == 'N' ?
            gemm('N', 'T', α, Ȳ, B) :
            gemm('N', 'N', α, Ȳ, B) :
        uppercase(tB) == 'N' ?
            gemm('N', 'T', α, B, Ȳ) :
            gemm('T', 'T', α, B, Ȳ)
∇(Ā::SM{T}, ::typeof(gemm), ::Arg4, _, Y::SM{T}, Ȳ::SM{T},
    tA::Char,
    tB::Char,
    α::T,
    A::SM{T},
    B::SM{T},
) where T<:BF =
    uppercase(tA) == 'N' ?
        uppercase(tB) == 'N' ?
            gemm!('N', 'T', α, Ȳ, B, one(T), Ā) :
            gemm!('N', 'N', α, Ȳ, B, one(T), Ā) :
        uppercase(tB) == 'N' ?
            gemm!('N', 'T', α, B, Ȳ, one(T), Ā) :
            gemm!('T', 'T', α, B, Ȳ, one(T), Ā)
∇(::typeof(gemm), ::Arg5, p, Y::SM{T}, Ȳ::SM{T},
    tA::Char,
    tB::Char,
    α::T,
    A::SM{T},
    B::SM{T},
) where T<:BF =
    uppercase(tA) == 'N' ?
        uppercase(tB) == 'N' ?
            gemm('T', 'N', α, A, Ȳ) :
            gemm('T', 'N', α, Ȳ, A) :
        uppercase(tB) == 'N' ?
            gemm('N', 'N', α, A, Ȳ) :
            gemm('T', 'T', α, Ȳ, A)
∇(B̄::SM{T}, ::typeof(gemm), ::Arg5, _, Y::SM{T}, Ȳ::SM{T},
    tA::Char,
    tB::Char,
    α::T,
    A::SM{T},
    B::SM{T},
) where T<:BF =
    uppercase(tA) == 'N' ?
        uppercase(tB) == 'N' ?
            gemm!('T', 'N', α, A, Ȳ, one(T), B̄) :
            gemm!('T', 'N', α, Ȳ, A, one(T), B̄) :
        uppercase(tB) == 'N' ?
            gemm!('N', 'N', α, A, Ȳ, one(T), B̄) :
            gemm!('T', 'T', α, Ȳ, A, one(T), B̄)

# `gemm` sensitivities implementation for `α = 1`.
push!(ops, DiffOp(:(LinearAlgebra.BLAS.gemm),
    :(Tuple{Char, Char, DLA.SM{T}, DLA.SM{T}} where T<:DLA.BF),
    [false, false, true, true],
))
∇(::typeof(gemm), ::Arg3, p, Y::SM{T}, Ȳ::SM{T},
    tA::Char,
    tB::Char,
    A::SM{T},
    B::SM{T}
) where T<:BF = ∇(gemm, Val{4}, p, Y, Ȳ, tA, tB, one(T), A, B)
∇(Ā::SM{T}, ::typeof(gemm), ::Arg3, p, Y::SM{T}, Ȳ::SM{T},
    tA::Char,
    tB::Char,
    A::SM{T},
    B::SM{T},
) where T<:BF = ∇(Ā, gemm, Val{4}, p, Y, Ȳ, tA, tB, one(T), A, B)
∇(::typeof(gemm), ::Arg4, p, Y::SM{T}, Ȳ::SM{T},
    tA::Char,
    tB::Char,
    A::SM{T},
    B::SM{T},
) where T<:BF = ∇(gemm, Val{5}, p, Y, Ȳ, tA, tB, one(T), A, B)
∇(B̄::SM{T}, ::typeof(gemm), ::Arg4, p, Y::SM{T}, Ȳ::SM{T},
    tA::Char,
    tB::Char,
    A::SM{T},
    B::SM{T},
) where T<:BF = ∇(B̄, gemm, Val{5}, p, Y, Ȳ, tA, tB, one(T), A, B)

# # `syrk` sensitivity implementations.
# @explicit_intercepts(
#     syrk,
#     Tuple{Char, Char, ∇Scalar, StridedVecOrMat{<:∇Scalar}},
#     [false, false, true, true],
# )
# function ∇(::typeof(syrk), ::Type{Arg{3}}, p, Y, Ȳ,
#     uplo::Char,
#     trans::Char,
#     α::∇Scalar,
#     A::StridedVecOrMat{<:∇Scalar},
# )
#     g! = uppercase(uplo) == 'L' ? tril! : triu!
#     return sum(g!(Ȳ .* Y)) / α
# end
# function ∇(::typeof(syrk), ::Type{Arg{4}}, p, Y, Ȳ,
#     uplo::Char,
#     trans::Char,
#     α::∇Scalar,
#     A::StridedVecOrMat{<:∇Scalar},
# )
#     triȲ = uppercase(uplo) == 'L' ? tril(Ȳ) : triu(Ȳ)
#     out = gemm('N', trans, α, triȲ .+ triȲ', A)
#     return uppercase(trans) == 'N' ? out : out'
# end
# function ∇(Ā::StridedVecOrMat{T}, ::typeof(syrk), ::Type{Arg{4}}, p, Y, Ȳ,
#     uplo::Char,
#     trans::Char,
#     α::∇Scalar,
#     A::StridedVecOrMat{T},
# ) where T<:∇Scalar
#     triȲ = uppercase(uplo) == 'L' ? tril(Ȳ) : triu(Ȳ)
#     out = gemm('N', trans, α, triȲ .+ triȲ', A)
#     return broadcast!((ā, δā)->ā+δā, Ā, Ā, uppercase(trans) == 'N' ? out : out')
# end

# # `syrk` sensitivity implementations for `α=1`.
# @explicit_intercepts(
#     syrk,
#     Tuple{Char, Char, StridedVecOrMat{<:∇Scalar}},
#     [false, false, true],
# )
# ∇(::typeof(syrk), ::Type{Arg{3}}, p, Y, Ȳ,
#     uplo::Char,
#     trans::Char,
#     A::StridedVecOrMat{<:∇Scalar},
# ) = ∇(syrk, Arg{4}, p, Y, Ȳ, uplo, trans, one(eltype(A)), A)
# ∇(Ā::StridedVecOrMat{T}, ::typeof(syrk), ::Type{Arg{4}}, p, Y, Ȳ,
#     uplo::Char,
#     trans::Char,
#     A::StridedVecOrMat{T},
# ) where T<:∇Scalar = ∇(Ā, syrk, Arg{4}, p, Y, Ȳ, uplo, char, one(eltype(A)), A)

# `symm` sensitivity implementations.
push!(ops, DiffOp(:(LinearAlgebra.BLAS.symm),
    :(Tuple{Char, Char, T, DLA.SM{T}, DLA.SM{T}} where T<:DLA.BF),
    [false, false, true, true, true],
))
∇(::typeof(symm), ::Arg3, p, Y::SM{T}, Ȳ::SM{T},
    side::Char,
    ul::Char,
    α::T,
    A::SM{T},
    B::SM{T},
) where T<:BF = sum(Ȳ .* Y) / α
function ∇(::typeof(symm), ::Arg4, p, Y::SM{T}, Ȳ::SM{T},
    side::Char,
    ul::Char,
    α::T,
    A::SM{T},
    B::SM{T},
) where T<:BF
    tmp = uppercase(side) == 'L' ? Ȳ * B' : B'Ȳ
    g! = uppercase(ul) == 'L' ? tril! : triu!
    return α * g!(tmp + tmp' - Diagonal(tmp))
end
function ∇(Ā::SM{T}, ::typeof(symm), ::Arg4, p, Y::SM{T}, Ȳ::SM{T},
    side::Char,
    ul::Char,
    α::T,
    A::SM{T},
    B::SM{T},
) where T<:BF
    tmp = uppercase(side) == 'L' ? Ȳ * B' : B'Ȳ
    g! = uppercase(ul) == 'L' ? tril! : triu!
    return broadcast!((ā, δā)->ā + δā, Ā, Ā, α * g!(tmp + tmp' - Diagonal(tmp)))
end
∇(::typeof(symm), ::Arg5, p, Y::SM{T}, Ȳ::SM{T},
    side::Char,
    ul::Char,
    α::T,
    A::SM{T},
    B::SM{T},
) where T<:BF = symm(side, ul, α, A, Ȳ)
∇(B̄::SM{T}, ::typeof(symm), ::Arg5, p, Y::SM{T}, Ȳ::SM{T},
    side::Char,
    ul::Char,
    α::T,
    A::SM{T},
    B::SM{T},
) where T<:BF = symm!(side, ul, α, A, Ȳ, one(T), B̄)

# `symm` sensitivity implementations for `α=1`.
push!(ops, DiffOp(:(LinearAlgebra.BLAS.symm),
    :(Tuple{Char, Char, DLA.SM{T}, DLA.SM{T}} where T<:DLA.BF),
    [false, false, true, true],
))
∇(::typeof(symm), ::Arg3, p, Y::SM{T}, Ȳ::SM{T},
    side::Char,
    ul::Char,
    A::SM{T},
    B::SM{T},
) where T<:BF = ∇(symm, Val{4}, p, Y, Ȳ, side, ul, one(T), A, B)
∇(Ā::SM{T}, ::typeof(symm), ::Arg3, p, Y::SM{T}, Ȳ::SM{T},
    side::Char,
    ul::Char,
    A::SM{T},
    B::SM{T},
) where T<:BF = ∇(Ā, symm, Val{4}, p, Y, Ȳ, side, ul, one(T), A, B)
∇(::typeof(symm), ::Arg4, p, Y::SM{T}, Ȳ::SM{T},
    side::Char,
    ul::Char,
    A::SM{T},
    B::SM{T},
) where T<:BF = ∇(symm, Val{5}, p, Y, Ȳ, side, ul, one(T), A, B)
∇(B̄::SM{T}, ::typeof(symm), ::Arg4, p, Y::SM{T}, Ȳ::SM{T},
    side::Char,
    ul::Char,
    A::SM{T},
    B::SM{T},
) where T<:BF = ∇(B̄, symm, Val{5}, p, Y, Ȳ, side, ul, one(T), A, B)

# `trmm` sensitivity implementations.
push!(ops, DiffOp(:(LinearAlgebra.BLAS.trmm),
    :(Tuple{Char, Char, Char, Char, T, DLA.SM{T}, DLA.SM{T}} where T<:DLA.BF),
    [false, false, false, false, true, true, true],
))
∇(::typeof(trmm), ::Arg5, p, Y::SM{T}, Ȳ::SM{T},
    side::Char, ul::Char, ta::Char, dA::Char,
    α::T,
    A::SM{T},
    B::SM{T},
) where T<:BF = sum(Ȳ .* Y) / α
function ∇(::typeof(trmm), ::Arg6, p, Y::SM{T}, Ȳ::SM{T},
    side::Char, ul::Char, ta::Char, dA::Char,
    α::T,
    A::SM{T},
    B::SM{T},
) where T<:BF
    Ā_full = uppercase(side) == 'L' ?
        uppercase(ta) == 'N' ?
            gemm('N', 'T', α, Ȳ, B) :
            gemm('N', 'T', α, B, Ȳ) :
        uppercase(ta) == 'N' ?
            gemm('T', 'N', α, B, Ȳ) :
            gemm('T', 'N', α, Ȳ, B)
    dA == 'U' && fill!(view(Ā_full, diagind(Ā_full)), zero(T))
    return (uppercase(ul) == 'L' ? tril! : triu!)(Ā_full)
end
∇(::typeof(trmm), ::Type{Val{7}}, p, Y::SM{T}, Ȳ::SM{T},
    side::Char, ul::Char, ta::Char, dA::Char,
    α::T,
    A::SM{T},
    B::SM{T},
) where T<:BF =
    uppercase(side) == 'L' ?
        uppercase(ta) == 'N' ?
            trmm('L', ul, 'T', dA, α, A, Ȳ) :
            trmm('L', ul, 'N', dA, α, A, Ȳ) :
        uppercase(ta) == 'N' ?
            trmm('R', ul, 'T', dA, α, A, Ȳ) :
            trmm('R', ul, 'N', dA, α, A, Ȳ)


# `trsm` sensitivity implementations.
push!(ops, DiffOp(:(LinearAlgebra.BLAS.trsm),
    :(Tuple{Char, Char, Char, Char, T, DLA.SM{T}, DLA.SM{T}} where T<:DLA.BF),
    [false, false, false, false, true, true, true],
))
∇(::typeof(trsm), ::Arg5, p, Y::SM{T}, Ȳ::SM{T},
    side::Char, ul::Char, ta::Char, dA::Char,
    α::T,
    A::SM{T},
    X::SM{T},
) where T<:BF = sum(Ȳ .* Y) / α
function ∇(::typeof(trsm), ::Arg6, p, Y::SM{T}, Ȳ::SM{T},
    side::Char, ul::Char, ta::Char, dA::Char,
    α::T,
    A::SM{T},
    X::SM{T},
) where T<:BF
    Ā_full = uppercase(side) == 'L' ?
        uppercase(ta) == 'N' ?
            trsm('L', ul, 'T', dA, -one(T), A, Ȳ * Y') :
            trsm('R', ul, 'T', dA, -one(T), A, Y * Ȳ') :
        uppercase(ta) == 'N' ?
            trsm('R', ul, 'T', dA, -one(T), A, Y'Ȳ) :
            trsm('L', ul, 'T', dA, -one(T), A, Ȳ'Y)
    dA == 'U' && fill!(view(Ā_full, diagind(Ā_full)), zero(T))
    return (uppercase(ul) == 'L' ? tril! : triu!)(Ā_full)
end
∇(::typeof(trsm), ::Type{Val{7}}, p, Y::SM{T}, Ȳ::SM{T},
    side::Char, ul::Char, ta::Char, dA::Char,
    α::T,
    A::SM{T},
    X::SM{T},
) where T<:BF =
    uppercase(side) == 'L' ?
        uppercase(ta) == 'N' ?
            trsm('L', ul, 'T', dA, α, A, Ȳ) :
            trsm('L', ul, 'N', dA, α, A, Ȳ) :
        uppercase(ta) == 'N' ?
            trsm('R', ul, 'T', dA, α, A, Ȳ) :
            trsm('R', ul, 'N', dA, α, A, Ȳ)
