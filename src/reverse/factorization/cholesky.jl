import LinearAlgebra.BLAS: gemv, gemv!, gemm!, trsm!, axpy!, ger!
import LinearAlgebra: chol, copytri!

#=
See [1] for implementation details: pages 5-9 in particular. The derivations presented in
[1] assume column-major layout, whereas Julia primarily uses row-major. We therefore
implement both the derivations in [1] and their transpose, which is more appropriate to
Julia.

[2] suggests that implementing these operations at the highest level of abstraction is the
way forward. There is, therefore, some code below to see what happens when we do this.

[1] - "Differentiation of the Cholesky decomposition", Murray 2016
[2] - "Auto-Differentiating Linear Algebra", Seeger et. al 2017.
=#

const UT = UpperTriangular
∇(::typeof(chol), ::Arg1, p, U::UT{T}, Ū::AM{T}, Σ::AM{T}) where T<:BF =
    chol_blocked_rev(Matrix(Ū), Matrix(U), 25, true)

# Experimental code implementing the algebraic sensitivities discussed in [2].
# """
#     ∇(::typeof(chol), ::Arg1, p, U::UT{T}, Ū::AM{T}, Σ::AM{T}) where T<:BF

# Transform Ū into Σ̄ in a non-allocating manner.
# """
# function ∇(::typeof(chol), ::Arg1, p, U::UT{T}, Ū::AM{T}, Σ::AM{T}, ::Symbol) where T<:BF
#     Σ̄ = A_mul_Bt!(Ū, U)
#     Σ̄ = copytri!(Σ̄, 'U')
#     Σ̄ = A_ldiv_B!(U, Σ̄)
#     BLAS.trsm!('R', 'U', 'T', 'N', one(T), U.data, Σ̄)
#     @inbounds for n in diagind(Σ̄)
#         Σ̄[n] *= 0.5
#     end
#     return Σ̄
# end

"""
    level2partition(A::AbstractMatrix, j::Int, upper::Bool)

Returns views to various bits of the lower triangle of `A` according to the
`level2partition` procedure defined in [1] if `upper` is `false`. If `upper` is `true` then
the transposed views are returned from the upper triangle of `A`.
"""
function level2partition(A::AM, j::Int, upper::Bool)

    # Check that A is square and j is a valid index.
    M, N = size(A)
    (0 >= j || j > M) && throw(ArgumentError("j is out of range."))
    M != N && throw(ArgumentError("A is not square."))

    if upper
        r = view(A, 1:j-1, j)
        d = view(A, j, j)
        B = view(A, 1:j-1, j+1:N)
        c = view(A, j, j+1:N)
    else
        r = view(A, j, 1:j-1)
        d = view(A, j, j)
        B = view(A, j+1:N, 1:j-1)
        c = view(A, j+1:N, j)
    end
    return r, d, B, c
end

"""
    level3partition(A::AbstractMatrix, j::Int, k::Int, upper::Bool)

Returns views to various bits of the lower triangle of `A` according to the
`level3partition` procedure defined in [1] if `upper` is `false`. If `upper` is `true` then
the transposed views are returned from the upper triangle of `A`.
"""
function level3partition(A::AM, j::Int, k::Int, upper::Bool)

    # Check that A is square and j is a valid index.
    M, N = size(A)
    (0 >= j || j > M) && throw(ArgumentError("j is out of range."))
    M != N && throw(ArgumentError("A is not square."))

    # Get views into bits of A.
    if upper
        R = view(A, 1:j-1, j:k)
        D = view(A, j:k, j:k)
        B = view(A, 1:j-1, k+1:N)
        C = view(A, j:k, k+1:N)
    else
        R = view(A, j:k, 1:j-1)
        D = view(A, j:k, j:k)
        B = view(A, k+1:N, 1:j-1)
        C = view(A, k+1:N, j:k)
    end
    return R, D, B, C
end

"""
    chol_unblocked_rev!(
        Ā::AbstractMatrix{T},
        L::AbstractMatrix{T},
        upper::Bool
    ) where T<:Real

Compute the reverse-mode sensitivities of the Cholesky factorisation in an unblocked manner.
If `upper` is `false`, then the sensitivites computed from and stored in the lower triangle
of `Ā` and `L` respectively. If `upper` is `true` then they are computed and stored in the
upper triangles. If at input `upper` is `false` and `tril(Ā) = L̄`, at output
`tril(Ā) = tril(Σ̄)`, where `Σ = LLᵀ`. Analogously, if at input `upper` is `true` and
`triu(Ā) = triu(Ū)`, at output `triu(Ā) = triu(Σ̄)` where `Σ = UᵀU`.
"""
function chol_unblocked_rev!(Σ̄::AM{T}, L::AM{T}, upper::Bool) where T<:Real

    # Check that L is square, that Σ̄ is square and that they are the same size.
    M, N = size(Σ̄)
    M != N && throw(ArgumentError("Σ̄ is not square."))

    # Compute the reverse-mode diff.
    j = N
    for ĵ in 1:N
        r, d, B, c = level2partition(L, j, upper)
        r̄, d̄, B̄, c̄ = level2partition(Σ̄, j, upper)

        # d̄ <- d̄ - c'c̄ / d.
        d̄[1] -= dot(c, c̄) / d[1]

        # [d̄ c̄'] <- [d̄ c̄'] / d.
        d̄ ./= d
        c̄ ./= d

        # r̄ <- r̄ - [d̄ c̄'] [r' B']'.
        r̄ = axpy!(-Σ̄[j, j], r, r̄)
        r̄ = gemv!(upper ? 'N' : 'T', -one(T), B, c̄, one(T), r̄)

        # B̄ <- B̄ - c̄ r.
        B̄ = upper ? ger!(-one(T), r, c̄, B̄) : ger!(-one(T), c̄, r, B̄)
        d̄ ./= 2
        j -= 1
    end
    return (upper ? triu! : tril!)(Σ̄)
end
chol_unblocked_rev(Σ̄::AM, L::AM, upper::Bool) = chol_unblocked_rev!(copy(Σ̄), L, upper)

"""
    chol_blocked_rev!(
        Σ̄::AbstractMatrix{T},
        L::AbstractMatrix{T},
        Nb::Int,
        upper::Bool
    ) where T<:BF

Compute the sensitivities of the Cholesky factorisation using a blocked, cache-friendly 
procedure. `Σ̄` are the sensitivities of `L`, and will be transformed into the sensitivities
of `Σ`, where `Σ = LLᵀ`. `Nb` is the block-size to use. If the upper triangle has been used
to represent the factorization, that is `Σ = UᵀU` where `U := Lᵀ`, then this should be
indicated by passing `upper = true`.
"""
function chol_blocked_rev!(Σ̄::AM{T}, L::AM{T}, Nb::Int, upper::Bool) where T<:BF

    # Check that L is square, that Σ̄ is square and that they are the same size.
    M, N = size(Σ̄)
    M != N && throw(ArgumentError("Σ̄ is not square."))

    tmp = Matrix{T}(undef, Nb, Nb)

    # Compute the reverse-mode diff.
    k = N
    if upper
        for k̂ in 1:Nb:N
            j = max(1, k - Nb + 1)
            R, D, B, C = level3partition(L, j, k, true)
            R̄, D̄, B̄, C̄ = level3partition(Σ̄, j, k, true)

            C̄ = trsm!('L', 'U', 'N', 'N', one(T), D, C̄)
            gemm!('N', 'N', -one(T), R, C̄, one(T), B̄)
            gemm!('N', 'T', -one(T), C, C̄, one(T), D̄)
            chol_unblocked_rev!(D̄, D, true)
            gemm!('N', 'T', -one(T), B, C̄, one(T), R̄)
            if size(D̄, 1) == Nb
                tmp = axpy!(one(T), D̄, transpose!(tmp, D̄))
                gemm!('N', 'N', -one(T), R, tmp, one(T), R̄)
            else
                gemm!('N', 'N', -one(T), R, D̄ + D̄', one(T), R̄)
            end

            k -= Nb
        end
        return triu!(Σ̄)
    else
        for k̂ in 1:Nb:N
            j = max(1, k - Nb + 1)
            R, D, B, C = level3partition(L, j, k, false)
            R̄, D̄, B̄, C̄ = level3partition(Σ̄, j, k, false)

            C̄ = trsm!('R', 'L', 'N', 'N', one(T), D, C̄)
            gemm!('N', 'N', -one(T), C̄, R, one(T), B̄)
            gemm!('T', 'N', -one(T), C̄, C, one(T), D̄)
            chol_unblocked_rev!(D̄, D, false)
            gemm!('T', 'N', -one(T), C̄, B, one(T), R̄)
            if size(D̄, 1) == Nb
                tmp = axpy!(one(T), D̄, transpose!(tmp, D̄))
                gemm!('N', 'N', -one(T), tmp, R, one(T), R̄)
            else
                gemm!('N', 'N', -one(T), D̄ + D̄', R, one(T), R̄)
            end

            k -= Nb
        end
        return tril!(Σ̄)
    end
end
chol_blocked_rev(Σ̄::AbstractMatrix, L::AbstractMatrix, Nb::Int, upper::Bool) =
    chol_blocked_rev!(copy(Σ̄), L, Nb, upper)
