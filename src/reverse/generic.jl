import LinearAlgebra: -, tr, inv, det, logdet, transpose, adjoint, norm, kron

# ############################# Unary sensitivities #############################

@reverse_rule(
    Y::AbstractArray{<:Real}, Ȳ::AbstractArray{<:Real},
    Base.:-(wrt(X::AbstractArray{<:Real})) = :(-$Ȳ),
)

@reverse_rule(
    Y::Real, Ȳ::Real,
    LinearAlgebra.tr(wrt(X::AbstractMatrix{<:Real})) = :(Diagonal(fill!(similar($X), $Ȳ)))
)

@reverse_rule(
    Y::AbstractMatrix{<:Real}, Ȳ::AbstractMatrix{<:Real},
    LinearAlgebra.inv(wrt(X::AbstractMatrix{<:Real})) = :(-$Y' * $Ȳ * $Y'),
)

@reverse_rule(
    Y::Real, Ȳ::Real,
    LinearAlgebra.det(wrt(X::AbstractMatrix{<:Real})) = :($Y * $Ȳ * inv($X)'),
)

@reverse_rule(
    Y::Real, Ȳ::Real,
    LinearAlgebra.logdet(wrt(X::AbstractMatrix{<:Real})) = :($Ȳ * inv($X)'),
)

@reverse_rule(
    Y::AbstractVecOrMat{<:Real}, Ȳ::AbstractVecOrMat{<:Real},
    LinearAlgebra.transpose(wrt(X::AbstractVecOrMat{<:Real})) = :($Ȳ'),
)

@reverse_rule(
    Y::AbstractVecOrMat{<:Real}, Ȳ::AbstractVecOrMat{<:Real},
    LinearAlgebra.adjoint(wrt(X::AbstractVecOrMat{<:Real})) = :($Ȳ'),
)

@reverse_rule(
    Y::Real, Ȳ::Real,
    LinearAlgebra.norm(wrt(X::AbstractArray{<:Real})) = :($Ȳ ./ $Y .* abs2.($X) ./ $X),
)

@reverse_rule(
    Y::Real, Ȳ::Real,
    LinearAlgebra.norm(wrt(X::Real)) = :($Ȳ * sign($X)),
)


############################# Binary sensitivities #############################

@reverse_rule(
    Y::AbstractVecOrMat{<:Real}, Ȳ::AbstractVecOrMat{<:Real},
    LinearAlgebra.:*(wrt(A::AbstractVecOrMat{<:Real}), B::AbstractVecOrMat{<:Real}) = :($Ȳ * $B')
)
@reverse_rule(
    Y::AbstractVecOrMat{<:Real}, Ȳ::AbstractVecOrMat{<:Real},
    LinearAlgebra.:*(A::AbstractVecOrMat{<:Real}, wrt(B::AbstractVecOrMat{<:Real})) = :($A' * $Ȳ)
)

@reverse_rule(
    Y::AbstractVecOrMat{<:Real}, Ȳ::AbstractVecOrMat{<:Real},
    LinearAlgebra.:/(wrt(A::AbstractVecOrMat{<:Real}), B::AbstractVecOrMat{<:Real}) = :($Ȳ / $B')
)
@reverse_rule(
    Y::AbstractVecOrMat{<:Real}, Ȳ::AbstractVecOrMat{<:Real},
    LinearAlgebra.:/(A::AbstractVecOrMat{<:Real}, wrt(B::AbstractVecOrMat{<:Real})) = :(-($Y)' * ($Ȳ / $B'))
)

@reverse_rule(
    Y::AbstractVecOrMat{<:Real}, Ȳ::AbstractVecOrMat{<:Real},
    LinearAlgebra.:\(wrt(A::AbstractVecOrMat{<:Real}), B::AbstractVecOrMat{<:Real}) = :(-($A' \ $Ȳ) * $Y')
)
@reverse_rule(
    Y::AbstractVecOrMat{<:Real}, Ȳ::AbstractVecOrMat{<:Real},
    LinearAlgebra.:\(A::AbstractVecOrMat{<:Real}, wrt(B::AbstractVecOrMat{<:Real})) = :($A' \ $Ȳ)
)

@reverse_rule(
    Y::Real, Ȳ::Real,
    LinearAlgebra.norm(wrt(A::AbstractArray{<:Real}), B::Real) =
        :($Ȳ .* $Y^(1 - $B) .* abs.($A).^$B ./ $A)
)
@reverse_rule(
    Y::Real, Ȳ::Real,
    LinearAlgebra.norm(A::AbstractArray{<:Real}, wrt(B::Real)) =
        :($Ȳ * ($Y^(1 - $B) * sum(abs.($A).^$B .* log.(abs.($A))) - $Y * log($Y)) / $B)
)

@reverse_rule(
    Y::Real, Ȳ::Real,
    LinearAlgebra.norm(wrt(A::Real), B::Real) = :($Ȳ * sign($A))
)
@reverse_rule(
    Y::Real, Ȳ::Real,
    LinearAlgebra.norm(A::Real, wrt(B::Real)) = :(0)
)

@reverse_rule(
    Y::AbstractMatrix{<:Real}, Ȳ::AbstractMatrix{<:Real},
    LinearAlgebra.kron(wrt(A::AbstractMatrix{<:Real}), B::AbstractMatrix{<:Real}) =
        :(_kron_rev_kernel_1($Y, $Ȳ, $A, $B)),
)
@reverse_rule(
    Y::AbstractMatrix{<:Real}, Ȳ::AbstractMatrix{<:Real},
    LinearAlgebra.kron(A::AbstractMatrix{<:Real}, wrt(B::AbstractMatrix{<:Real})) =
        :(_kron_rev_kernel_2($Y, $Ȳ, $A, $B)),
)

function _kron_rev_kernel_1(
    Y::AbstractMatrix{<:Real},
    Ȳ::AbstractMatrix{<:Real},
    A::AbstractMatrix{<:Real},
    B::AbstractMatrix{<:Real},
)
    Ā = similar(A)
    (I, J), (K, L), m = size(A), size(B), length(Y)
    @inbounds for j = reverse(1:J), l = reverse(1:L), i = reverse(1:I)
        āij = Ā[i, j]
        for k = reverse(1:K)
            āij += Ȳ[m] * B[k, l]
            m -= 1
        end
        Ā[i, j] = āij
    end
    return Ā
end

function _kron_rev_kernel_2(
    Y::AbstractMatrix{<:Real},
    Ȳ::AbstractMatrix{<:Real},
    A::AbstractMatrix{<:Real},
    B::AbstractMatrix{<:Real},
)
    B̄ = similar(B)
    (I, J), (K, L), m = size(A), size(B), length(Y)
    @inbounds for j = reverse(1:J), l = reverse(1:L), i = reverse(1:I)
        aij = A[i, j]
        for k = reverse(1:K)
            B̄[k, l] += aij * Ȳ[m]
            m -= 1
        end
    end
    return B̄
end

# push!(ops, DiffOp(:(LinearAlgebra.:*), :(Tuple{DLA.AVM, DLA.AVM}), [true, true]))
# ∇(::typeof(*), ::Arg1, p, Y::ASVM, Ȳ::ASVM, A::AVM, B::AVM) = Ȳ * B'
# ∇(::typeof(*), ::Arg2, p, Y::ASVM, Ȳ::ASVM, A::AVM, B::AVM) = A' * Ȳ

# push!(ops, DiffOp(:(LinearAlgebra.:/), :(Tuple{DLA.AVM, DLA.AVM}), [true, true]))
# ∇(::typeof(/), ::Arg1, p, Y::ASVM, Ȳ::ASVM, A::AVM, B::AVM) = Ȳ / B'
# ∇(::typeof(/), ::Arg2, p, Y::ASVM, Ȳ::ASVM, A::AVM, B::AVM) = -Y' * (Ȳ / B')

# push!(ops, DiffOp(:(LinearAlgebra.:\), :(Tuple{DLA.AVM, DLA.AVM}), [true, true]))
# ∇(::typeof(\), ::Arg1, p, Y::ASVM, Ȳ::ASVM, A::AVM, B::AVM) = -(A' \ Ȳ) * Y'
# ∇(::typeof(\), ::Arg2, p, Y::ASVM, Ȳ::ASVM, A::AVM, B::AVM) = A' \ Ȳ

# push!(ops, DiffOp(:(LinearAlgebra.norm), :(Tuple{DLA.AA, Real}), [true, true]))
# ∇(::typeof(norm), ::Arg1, p, Y::Real, Ȳ::Real, A::AA, B::Real) =
#     Ȳ .* Y^(1 - B) .* abs.(A).^B ./ A
# ∇(::typeof(norm), ::Arg2, p, Y::Real, Ȳ::Real, A::AA, B::Real) =
#     Ȳ * (Y^(1 - B) * sum(abs.(A).^B .* log.(abs.(A))) - Y * log(Y)) / B

# push!(ops, DiffOp(:(LinearAlgebra.norm), :(Tuple{Real, Real}), [true, true]))
# ∇(::typeof(norm), ::Arg1, p, Y::Real, Ȳ::Real, A::Real, B::Real) = Ȳ * sign(A)
# ∇(::typeof(norm), ::Arg2, p, Y::Real, Ȳ::Real, A::Real, B::Real) = 0

# push!(ops, DiffOp(:(LinearAlgebra.kron), :(Tuple{AM, AM}), [true, true]))
# ∇(::typeof(kron), ::Type{Val{1}}, p, Y::AM, Ȳ::AM, A::AM, B::AM) =
#     ∇(zero(A), kron, Val{1}, p, Y, Ȳ, A, B)
# ∇(::typeof(kron), ::Type{Val{2}}, p, Y::AM, Ȳ::AM, A::AM, B::AM) =
#     ∇(zero(B), kron, Val{2}, p, Y, Ȳ, A, B)
# function ∇(Ā::AM, ::typeof(kron), ::Type{Val{1}}, p, Y::AM, Ȳ::AM, A::AM, B::AM)
#     @assert size(Ā) == size(A)
#     (I, J), (K, L), m = size(A), size(B), length(Y)
#     @inbounds for j = reverse(1:J), l = reverse(1:L), i = reverse(1:I)
#         āij = Ā[i, j]
#         for k = reverse(1:K)
#             āij += Ȳ[m] * B[k, l]
#             m -= 1
#         end
#         Ā[i, j] = āij
#     end
#     return Ā
# end
# function ∇(B̄::AM, ::typeof(kron), ::Type{Val{2}}, p, Y::AM, Ȳ::AM, A::AM, B::AM)
#     @assert size(B̄) == size(B)
#     (I, J), (K, L), m = size(A), size(B), length(Y)
#     @inbounds for j = reverse(1:J), l = reverse(1:L), i = reverse(1:I)
#         aij = A[i, j]
#         for k = reverse(1:K)
#             B̄[k, l] += aij * Ȳ[m]
#             m -= 1
#         end
#     end
#     return B̄
# end
