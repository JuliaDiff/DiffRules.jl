using LinearAlgebra
using LinearAlgebra: -, tr, inv, det, logdet, transpose, adjoint, norm
using DiffRules: ReverseRuleKey, DEFINED_REVERSE_RULES, make_named_signature

function unary_ȲD(key::ReverseRuleKey)
    f = @eval $(key[1]).$(key[2])
    arg_names = vcat(gensym(), gensym(), [gensym() for _ in 1:arity(key) - 2])
    typed_args = Expr(:tuple, make_named_signature(arg_names, key)...)
    body = DEFINED_REVERSE_RULES[key](arg_names...)
    return f, eval(Expr(Symbol("->"), typed_args, body))
end


@testset "generic" begin

let
    P, Q, rng, N = 4, 3, MersenneTwister(123456), 100

    # Utility for generating square matrices, vectors, non-square matrices and scalars.
    mPP, mQQ = ()->randn(rng, P, P), ()->randn(rng, Q, Q)
    mPQ, mQP = ()->randn(rng, P, Q), ()->randn(rng, Q, P)
    mPQQP = ()->randn(rng, P * Q, P * Q)
    v, sc = ()->randn(rng, P), ()->randn(rng)
    psd = ()->(A = randn(rng, P, P); transpose(A) * A + 1e-3I)

    # sig = :(Tuple{AbstractArray{<:Real}, AbstractArray{<:Real}, AbstractArray{<:Real}})
    # @show unary_ȲD((:Base, :-, sig, (1,)))
    # f, df = unary_ȲD((:Base, :-, sig, (1,)))
    # @show f(mPP())
    # Ȳ = mPP()
    # @show Ȳ
    # @show df(mPP(), Ȳ, mPP())


    sig = :(Tuple{AbstractArray{<:Real}, AbstractArray{<:Real}, AbstractArray{<:Real}})
    @test check_errs(N, unary_ȲD((:Base, :-, sig, (1,)))..., mPQ, mPQ, mPQ)

    sig = :(Tuple{Real, Real, AbstractMatrix{<:Real}})
    @test check_errs(N, unary_ȲD((:LinearAlgebra, :tr, sig, (1,)))..., sc, mPP, mPP)

    sig = :(Tuple{AbstractMatrix{<:Real}, AbstractMatrix{<:Real}, AbstractMatrix{<:Real}})
    @test check_errs(N, unary_ȲD((:LinearAlgebra, :inv, sig, (1,)))..., mPP, mPP, mPP)

    sig = :(Tuple{Real, Real, AbstractMatrix{<:Real}})
    @test check_errs(N, unary_ȲD((:LinearAlgebra, :det, sig, (1,)))..., sc, mPP, mPP)

    sig = :(Tuple{Real, Real, AbstractMatrix{<:Real}})
    @test check_errs(N, unary_ȲD((:LinearAlgebra, :logdet, sig, (1,)))..., sc, psd, psd)

    sig = :(Tuple{AbstractVecOrMat{<:Real}, AbstractVecOrMat{<:Real}, AbstractVecOrMat{<:Real}})
    @test check_errs(N, unary_ȲD((:LinearAlgebra, :transpose, sig, (1,)))..., mQP, mPQ, mPQ)

    sig = :(Tuple{AbstractVecOrMat{<:Real}, AbstractVecOrMat{<:Real}, AbstractVecOrMat{<:Real}})
    @test check_errs(N, unary_ȲD((:LinearAlgebra, :adjoint, sig, (1,)))..., mQP, mPQ, mPQ)

    sig = :(Tuple{Real, Real, AbstractArray{<:Real}})
    @test check_errs(N, unary_ȲD((:LinearAlgebra, :norm, sig, (1,)))..., sc, mPQ, mPQ)

    sig = :(Tuple{Real, Real, Real})
    @test check_errs(N, unary_ȲD((:LinearAlgebra, :norm, sig, (1,)))..., sc, sc, sc)

    # # Test all of the binary sensitivities.
    # @test check_errs(N, binary_ȲD(*, 1, mQP)..., mPP, mPQ, mPQ)
    # @test check_errs(N, binary_ȲD(*, 2, mPQ)..., mPP, mQP, mQP)
    # @test check_errs(N, binary_ȲD(*, 1, mQP)..., mPP, ()->mQP()', ()->mQP()')
    # @test check_errs(N, binary_ȲD(*, 2, ()->mQP()')..., mPP, mQP, mQP)
    # @test check_errs(N, binary_ȲD(*, 1, ()->mPQ()')..., mPP, mPQ, mPQ)
    # @test check_errs(N, binary_ȲD(*, 2, mPQ)..., mPP, ()->mPQ()', ()->mPQ()')
    # @test check_errs(N, binary_ȲD(*, 1, ()->mPQ()')..., mPP, ()->mQP()', ()->mQP()')
    # @test check_errs(N, binary_ȲD(*, 2, ()->mQP()')..., mPP, ()->mPQ()', ()->mPQ()')
    # @test check_errs(N, binary_ȲD(/, 1, mQQ)..., mPQ, mPQ, mPQ)
    # @test check_errs(N, binary_ȲD(/, 2, mPQ)..., mPQ, mQQ, mQQ)
    # @test check_errs(N, binary_ȲD(/, 1, mQQ)..., mPQ, ()->mQP()', ()->mQP()')
    # @test check_errs(N, binary_ȲD(/, 2, ()->mQP()')..., mPQ, mQQ, mQQ)
    # @test check_errs(N, binary_ȲD(/, 1, ()->mQQ()')..., mPQ, mPQ, mPQ)
    # @test check_errs(N, binary_ȲD(/, 2, mPQ)..., mPQ, ()->mQQ()', ()->mQQ()')
    # @test check_errs(N, binary_ȲD(/, 1, ()->mQQ()')..., mPQ, ()->mQP()', ()->mQP()')
    # @test check_errs(N, binary_ȲD(/, 2, ()->mQP()')..., mPQ, ()->mQQ()', ()->mQQ()')
    # @test check_errs(N, binary_ȲD(\, 1, mQP)..., mQP, mQQ, mQQ)
    # @test check_errs(N, binary_ȲD(\, 2, mQQ)..., mQP, mQP, mQP)
    # @test check_errs(N, binary_ȲD(\, 1, mQP)..., mQP, ()->mQQ()', ()->mQQ()')
    # @test check_errs(N, binary_ȲD(\, 2, ()->mQQ()')..., mQP, mQP, mQP)
    # @test check_errs(N, binary_ȲD(\, 1, ()->mPQ()')..., mQP, mQQ, mQQ)
    # @test check_errs(N, binary_ȲD(\, 2, mQQ)..., mQP, ()->mPQ()', ()->mPQ()')
    # @test check_errs(N, binary_ȲD(\, 1, ()->mPQ()')..., mQP, ()->mQQ()', ()->mQQ()')
    # @test check_errs(N, binary_ȲD(\, 2, ()->mQQ()')..., mQP, ()->mPQ()', ()->mPQ()')
    # @test check_errs(N, binary_ȲD(vecnorm, 1, sc)..., sc, mPQ, mPQ)
    # @test check_errs(N, binary_ȲD(vecnorm, 2, mPQ)..., sc, sc, sc)
    # @test check_errs(N, binary_ȲD(vecnorm, 1, sc)..., sc, sc, sc)
    # @test check_errs(N, binary_ȲD(vecnorm, 2, sc)..., sc, sc, sc)
    # @test check_errs(N, binary_ȲD(kron, 1, mQP)..., mPQQP, mPQ, mPQ)
    # @test check_errs(N, binary_ȲD(kron, 2, mPQ)..., mPQQP, mQP, mQP)
    # @test check_errs(N, binary_ȲD_inplace(kron, 1, mQP, mPQ())..., mPQQP, mPQ, mPQ)
    # @test check_errs(N, binary_ȲD_inplace(kron, 2, mPQ, mQP())..., mPQQP, mQP, mQP)
end

end
