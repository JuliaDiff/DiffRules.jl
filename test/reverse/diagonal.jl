@testset "Diagonal" begin

    # diag:
    let P = 10, rng = MersenneTwister(123456), N = 10, k = 2

        # Some rngs.
        vP, mPP, vPk = ()->randn(rng, P), ()->randn(rng, P, P), ()->randn(rng, P - k)

        # Test on-central-diagonal `diag`.
        X0 = randn!(rng, similar(mPP()))
        @test check_errs(N, unary_ȲD(diag)..., vP, mPP, mPP)
        @test check_errs(N, unary_ȲD_inplace(diag, X0)..., vP, mPP, mPP)

        # Test off-central-diagonal `diag`.
        _diag = X->diag(X, k)
        _∇diag = (y, ȳ, X)->∇(diag, Val{1}, (), y, ȳ, X, k)
        _∇diag_inp = (y, ȳ, X)->∇(copy(X0), diag, Val{1}, (), y, ȳ, X, k) - X0
        @test check_errs(N, _diag, _∇diag, vPk, mPP, mPP)
        @test check_errs(N, _diag, _∇diag_inp, vPk, mPP, mPP)
    end

    # diagm:
    let P = 10, rng = MersenneTwister(123456), N = 10, k = 3

        # Some rngs.
        sc, vP = ()->randn(rng), ()->randn(rng, P)
        mPP, mPPk = ()->randn(rng, P, P), ()->randn(rng, P + k, P + k)

        # Test on-central-diagonal `diagm`.
        x0 = randn!(rng, similar(vP()))
        _diagm = x->diagm(0=>x)
        _∇diagm = (Y, Ȳ, x)->∇(diagm, Val{1}, (), Y, Ȳ, 0=>x)
        _∇diagm_inp = (Y, Ȳ, x)->∇(copy(x0), diagm, Val{1}, (), Y, Ȳ, 0=>x) - x0
        @test check_errs(N, _diagm, _∇diagm, mPP, vP, vP)
        @test check_errs(N, _diagm, _∇diagm_inp, mPP, vP, vP)

        # Test off-central-diagonal `diagm`.
        _diagm = x->diagm(k=>x)
        _∇diagm = (Y, Ȳ, x)->∇(diagm, Val{1}, (), Y, Ȳ, k=>x)
        _∇diagm_inp = (Y, Ȳ, x)->∇(copy(x0), diagm, Val{1}, (), Y, Ȳ, k=>x) - x0
        @test check_errs(N, _diagm, _∇diagm, mPPk, vP, vP)
        @test check_errs(N, _diagm, _∇diagm_inp, mPPk, vP, vP)
    end

    # Diagonal:
    let P = 10, rng = MersenneTwister(123456), N = 10
        sc, vP = ()->abs(randn(rng)), ()->abs.(randn(rng, P))
        mPP, dPP = ()->abs.(randn(rng, P, P)), ()->Diagonal(abs.(randn(rng, P)))
        D0 = Diagonal(randn!(rng, similar(vP())))

        # Construction.
        @test check_errs(N, unary_ȲD(Diagonal)..., dPP, vP, vP)
        @test check_errs(N, unary_ȲD_inplace(Diagonal, vP())..., dPP, vP, vP)
        @test check_errs(N, unary_ȲD(Diagonal)..., dPP, mPP, mPP)
        @test check_errs(N, unary_ȲD_inplace(Diagonal, mPP())..., dPP, mPP, mPP)

        # Determinant.
        @test check_errs(N, unary_ȲD(det)..., sc, dPP, dPP)
        @test check_errs(N, unary_ȲD_inplace(det, D0)..., sc, dPP, dPP)

        # Log Determinant.
        @test check_errs(N, unary_ȲD(logdet)..., sc, dPP, dPP)
        @test check_errs(N, unary_ȲD_inplace(logdet, D0)..., sc, dPP, dPP)
    end
end
