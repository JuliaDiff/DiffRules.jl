using LinearAlgebra: -, tr, inv, det, logdet, transpose, adjoint, norm

@testset "generic" begin

let
    P = 4, Q = 3, rng = MersenneTwister(123456), N = 100

    # Utility for generating square matrices, vectors, non-square matrices and scalars.
    mPP, mQQ = ()->randn(rng, P, P), ()->randn(rng, Q, Q)
    mPQ, mQP = ()->randn(rng, P, Q), ()->randn(rng, Q, P)
    mPQQP = ()->randn(rng, P * Q, P * Q)
    v, sc = ()->randn(rng, P), ()->randn(rng)
    psd = ()->(A = randn(rng, P, P); transpose(A) * A + 1e-3I)

    # Test all of the unary sensitivities.
    @test check_errs(N, unary_ȲD(-)..., mPQ, mPQ, mPQ)
    @test check_errs(N, unary_ȲD(tr)..., sc, mPP, mPP)
    @test check_errs(N, unary_ȲD(inv)..., mPP, mPP, mPP)
    @test check_errs(N, unary_ȲD(det)..., sc, mPP, mPP)
    @test check_errs(N, unary_ȲD(logdet)..., sc, psd, psd)
    @test check_errs(N, unary_ȲD(transpose)..., mQP, mPQ, mPQ)
    @test check_errs(N, unary_ȲD(adjoint)..., mQP, mPQ, mPQ)
    @test check_errs(N, unary_ȲD(norm)..., sc, mPQ, mPQ)
    @test check_errs(N, unary_ȲD(norm)..., sc, sc, sc)

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
