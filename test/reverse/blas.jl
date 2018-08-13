@testset "BLAS" begin

import LinearAlgebra.BLAS: nrm2, asum, gemm, gemv, symm, symv, trmm, trmv, trsm, trsv


################################## Level 1 ##################################
let P = 10, Q = 6, rng = MersenneTwister(123456), N = 10

    # Utility random generators.
    sc, vP, vQ = ()->randn(rng), ()->randn(rng, P), ()->randn(rng, Q)

    # Unit-stride dot.
    @test check_errs(N, binary_ȲD(LinearAlgebra.dot, 1, vP)..., sc, vP, vP)
    @test check_errs(N, binary_ȲD(LinearAlgebra.dot, 2, vP)..., sc, vP, vP)
    @test check_errs(N, binary_ȲD_inplace(LinearAlgebra.dot, 1, vP, vP())..., sc, vP, vP)
    @test check_errs(N, binary_ȲD_inplace(LinearAlgebra.dot, 2, vP, vP())..., sc, vP, vP)

    # Strided dot.
    _x, _y = vP(), vQ()
    _dot2, _dot4 = x->BLAS.dot(5, x, 2, _y, 1), y->BLAS.dot(5, _x, 2, y, 1)
    _∇dot2 = (z, z̄, x)->∇(BLAS.dot, Val{2}, (), z, z̄, 5, x, 2, _y, 1)
    _∇dot4 = (z, z̄, y)->∇(BLAS.dot, Val{4}, (), z, z̄, 5, _x, 2, y, 1)
    @test check_errs(N, _dot2, _∇dot2, sc, vP, vP)
    @test check_errs(N, _dot4, _∇dot4, sc, vQ, vQ)

    # In-place strided dot.
    _δx, _δy = vP(), vQ()
    _∇dot2 = (z, z̄, x)->∇(copy(_δx), BLAS.dot, Val{2}, (), z, z̄, 5, x, 2, _y, 1) - _δx
    _∇dot4 = (z, z̄, y)->∇(copy(_δy), BLAS.dot, Val{4}, (), z, z̄, 5, _x, 2, y, 1) - _δy
    @test check_errs(N, _dot2, _∇dot2, sc, vP, vP)
    @test check_errs(N, _dot4, _∇dot4, sc, vQ, vQ)

    # Unit-stride nrm2.
    @test check_errs(N, unary_ȲD(nrm2)..., sc, vP, vP)
    @test check_errs(N, unary_ȲD_inplace(nrm2, vP())..., sc, vP, vP)

    # Arbitrary-stride nrm2.
    _nrm2 = x->nrm2(5, x, 2)
    _∇nrm2 = (y, ȳ, x)->∇(nrm2, Val{2}, (), y, ȳ, 5, x, 2)
    _∇nrm2_in_place = (y, ȳ, x)->∇(copy(_δx), nrm2, Val{2}, (), y, ȳ, 5, x, 2) - _δx
    @test check_errs(N, _nrm2, _∇nrm2, sc, vP, vP)
    @test check_errs(N, _nrm2, _∇nrm2_in_place, sc, vP, vP)

    # Unit-stride `asum`.
    @test check_errs(N, unary_ȲD(asum)..., sc, vP, vP)
    @test check_errs(N, unary_ȲD_inplace(asum, vP())..., sc, vP, vP)

    # Arbitrary-stride `asum`.
    _asum = x->asum(5, x, 2)
    _∇asum = (y, ȳ, x)->∇(asum, Val{2}, (), y, ȳ, 5, x, 2)
    _∇asum_in_place = (y, ȳ, x)->∇(copy(_δx), asum, Val{2}, (), y, ȳ, 5, x, 2) - _δx
    @test check_errs(N, _asum, _∇asum, sc, vP, vP)
    @test check_errs(N, _asum, _∇asum_in_place, sc, vP, vP)
end


################################## Level 2 ##################################
let P = 10, Q = 6, rng = MersenneTwister(123456), N = 10

    # Utility random generators.
    sc, vP, vQ = ()->randn(rng), ()->randn(rng, P), ()->randn(rng, Q)
    mPQ, mQP = ()->randn(rng, P, Q), ()->randn(rng, Q, P)
    mPP = ()->randn(rng, P, P)

    # gemv:
    for tA in ['T', 'N']
        _α, _x = sc(), vQ()
        _A = tA == 'N' ? mPQ() : mQP()
        A = tA == 'N' ? mPQ : mQP
        _δA, _δx = randn!(rng, similar(_A)), randn!(rng, similar(_x))

        # α != 1 tests.
        _gemv1 = α->gemv(tA, α, _A, _x)
        _gemv2 = A->gemv(tA, _α, A, _x)
        _gemv3 = x->gemv(tA, _α, _A, x)
        _∇gemv1 = (y, ȳ, α)->∇(gemv, Val{2}, (), y, ȳ, tA, α, _A, _x)
        _∇gemv2 = (y, ȳ, A)->∇(gemv, Val{3}, (), y, ȳ, tA, _α, A, _x)
        _∇gemv2_inp = (y, ȳ, A)->∇(copy(_δA), gemv, Val{3}, (), y, ȳ, tA, _α, A, _x) - _δA
        _∇gemv3 = (y, ȳ, x)->∇(gemv, Val{4}, (), y, ȳ, tA, _α, _A, x)
        _∇gemv3_inp = (y, ȳ, x)->∇(copy(_δx), gemv, Val{4}, (), y, ȳ, tA, _α, _A, x) - _δx

        @test check_errs(N, _gemv1, _∇gemv1, vP, sc, sc)
        @test check_errs(N, _gemv2, _∇gemv2, vP, A, A)
        @test check_errs(N, _gemv2, _∇gemv2_inp, vP, A, A)
        @test check_errs(N, _gemv3, _∇gemv3, vP, vQ, vQ)
        @test check_errs(N, _gemv3, _∇gemv3_inp, vP, vQ, vQ)

        # α = 1 tests.
        _gemv1, _gemv2 = A->gemv(tA, A, _x), x->gemv(tA, _A, x)
        _∇gemv1 = (y, ȳ, A)->∇(gemv, Val{2}, (), y, ȳ, tA, A, _x)
        _∇gemv1_inp = (y, ȳ, A)->∇(copy(_δA), gemv, Val{2}, (), y, ȳ, tA, A, _x) - _δA
        _∇gemv2 = (y, ȳ, x)->∇(gemv, Val{3}, (), y, ȳ, tA, _A, x)
        _∇gemv2_inp = (y, ȳ, x)->∇(copy(_δx), gemv, Val{3}, (), y, ȳ, tA, _A, x) - _δx

        @test check_errs(N, _gemv1, _∇gemv1, vP, A, A)
        @test check_errs(N, _gemv1, _∇gemv1_inp, vP, A, A)
        @test check_errs(N, _gemv2, _∇gemv2, vP, vQ, vQ)
        @test check_errs(N, _gemv2, _∇gemv2_inp, vP, vQ, vQ)
    end

    # symv:
    for ul in ['L', 'U']
        α_gen, A_gen, x_gen = sc, mPP, vP
        _α, _A, _x = α_gen(), A_gen(), x_gen()
        _δA, _δx = randn!(rng, similar(_A)), randn!(rng, similar(_x))

        # α != 1 tests.
        _symv_α = α->symv(ul, α, _A, _x)
        _symv_A = A->symv(ul, _α, A, _x)
        _symv_x = x->symv(ul, _α, _A, x)
        _∇symv_α = (y, ȳ, α)->∇(symv, Val{2}, (), y, ȳ, ul, α, _A, _x)
        _∇symv_A = (y, ȳ, A)->∇(symv, Val{3}, (), y, ȳ, ul, _α, A, _x)
        _∇symv_A_inp = (y, ȳ, A)->∇(copy(_δA), symv, Val{3}, (), y, ȳ, ul, _α, A, _x) - _δA
        _∇symv_x = (y, ȳ, x)->∇(symv, Val{4}, (), y, ȳ, ul, _α, _A, x)
        _∇symv_x_inp = (y, ȳ, x)->∇(copy(_δx), symv, Val{4}, (), y, ȳ, ul, _α, _A, x) - _δx

        @test check_errs(N, _symv_α, _∇symv_α, vP, sc, sc)
        @test check_errs(N, _symv_A, _∇symv_A, vP, A_gen, A_gen)
        @test check_errs(N, _symv_A, _∇symv_A_inp, vP, A_gen, A_gen)
        @test check_errs(N, _symv_x, _∇symv_x, vP, x_gen, x_gen)
        @test check_errs(N, _symv_x, _∇symv_x_inp, vP, x_gen, x_gen)

        # α = 1 tests.
        _symv_A = A->symv(ul, A, _x)
        _symv_x = x->symv(ul, _A, x)
        _∇symv_A = (y, ȳ, A)->∇(symv, Val{2}, (), y, ȳ, ul, A, _x)
        _∇symv_A_inp = (y, ȳ, A)->∇(copy(_δA), symv, Val{2}, (), y, ȳ, ul, A, _x) - _δA
        _∇symv_x = (y, ȳ, x)->∇(symv, Val{3}, (), y, ȳ, ul, _A, x)
        _∇symv_x_inp = (y, ȳ, x)->∇(copy(_δx), symv, Val{3}, (), y, ȳ, ul, _A, x) - _δx

        @test check_errs(N, _symv_A, _∇symv_A, vP, A_gen, A_gen)
        @test check_errs(N, _symv_A, _∇symv_A_inp, vP, A_gen, A_gen)
        @test check_errs(N, _symv_x, _∇symv_x, vP, x_gen, x_gen)
        @test check_errs(N, _symv_x, _∇symv_x_inp, vP, x_gen, x_gen)
    end

    for f in [trmv, trsv], ul in ['L', 'U'], tA in ['N', 'T'], dA in ['U', 'N']
        A_gen, x_gen = mPP, vP
        _A, _x = A_gen(), x_gen()

        _f_A = A->f(ul, tA, dA, A, _x)
        _f_x = x->f(ul, tA, dA, _A, x)
        _∇f_A = (y, ȳ, A)->∇(f, Val{4}, (), y, ȳ, ul, tA, dA, A, _x)
        _∇f_x = (y, ȳ, x)->∇(f, Val{5}, (), y, ȳ, ul, tA, dA, _A, x)

        @test check_errs(N, _f_A, _∇f_A, vP, A_gen, A_gen)
        @test check_errs(N, _f_x, _∇f_x, vP, x_gen, x_gen)
    end
end


################################## Level 3 ##################################
let P = 5, Q = 3, rng = MersenneTwister(123456), N = 10

    # Utility random generators.
    sc, mPP, mQQ = ()->randn(rng), ()->randn(rng, P, P), ()->randn(rng, Q, Q)
    mPQ, mQP = ()->randn(rng, P, Q), ()->randn(rng, Q, P)

    # gemm:
    for tA in ['T', 'N'], tB in ['T', 'N']

        # Generate conformal test matrices.
        A_gen = mPQ
        B_gen = tA == 'N' ?
                    (tB == 'N' ? mQP : mPQ) :
                    (tB == 'T' ? mQP : mPQ)
        C_gen = tA == 'N' ? mPP : mQQ

        # α != 1 tests.
        _α, _A, _B = sc(), A_gen(), B_gen()
        _δA, _δB = randn!(rng, similar(_A)), randn!(rng, similar(_B))
        _gemm_α = α->gemm(tA, tB, α, _A, _B)
        _gemm_A = A->gemm(tA, tB, _α, A, _B)
        _gemm_B = B->gemm(tA, tB, _α, _A, B)
        _∇gemm_α = (Y, Ȳ, α)->∇(gemm, Val{3}, (), Y, Ȳ, tA, tB, α, _A, _B)
        _∇gemm_A = (Y, Ȳ, A)->∇(gemm, Val{4}, (), Y, Ȳ, tA, tB, _α, A, _B)
        _∇gemm_B = (Y, Ȳ, B)->∇(gemm, Val{5}, (), Y, Ȳ, tA, tB, _α, _A, B)
        _∇gemm_A_inp = (Y, Ȳ, A)->∇(copy(_δA), gemm, Val{4}, (), Y, Ȳ, tA, tB, _α, A, _B) - _δA
        _∇gemm_B_inp = (Y, Ȳ, B)->∇(copy(_δB), gemm, Val{5}, (), Y, Ȳ, tA, tB, _α, _A, B) - _δB

        @test check_errs(N, _gemm_α, _∇gemm_α, C_gen, sc, sc)
        @test check_errs(N, _gemm_A, _∇gemm_A, C_gen, A_gen, A_gen)
        @test check_errs(N, _gemm_A, _∇gemm_A_inp, C_gen, A_gen, A_gen)
        @test check_errs(N, _gemm_B, _∇gemm_B, C_gen, B_gen, B_gen)
        @test check_errs(N, _gemm_B, _∇gemm_B_inp, C_gen, B_gen, B_gen)

        # α = 1 tests.
        _A, _B = A_gen(), B_gen()
        _gemm_A = A->gemm(tA, tB, A, _B)
        _gemm_B = B->gemm(tA, tB, _A, B)
        _∇gemm_A = (Y, Ȳ, A)->∇(gemm, Val{3}, (), Y, Ȳ, tA, tB, A, _B)
        _∇gemm_B = (Y, Ȳ, B)->∇(gemm, Val{4}, (), Y, Ȳ, tA, tB, _A, B)
        _∇gemm_A_inp = (Y, Ȳ, A)->∇(copy(_δA), gemm, Val{3}, (), Y, Ȳ, tA, tB, A, _B) - _δA
        _∇gemm_B_inp = (Y, Ȳ, B)->∇(copy(_δB), gemm, Val{4}, (), Y, Ȳ, tA, tB, _A, B) - _δB

        @test check_errs(N, _gemm_A, _∇gemm_A, C_gen, A_gen, A_gen)
        @test check_errs(N, _gemm_A, _∇gemm_A_inp, C_gen, A_gen, A_gen)
        @test check_errs(N, _gemm_B, _∇gemm_B, C_gen, B_gen, B_gen)
        @test check_errs(N, _gemm_B, _∇gemm_B_inp, C_gen, B_gen, B_gen)
    end

    # symm:
    for side in ['L', 'R'], ul in ['L', 'U']

        # Fixed qtts.
        _α, _A, _B = sc(), mPP(), side == 'L' ? mPQ() : mQP()
        _δA, _δB = randn!(rng, similar(_A)), randn!(rng, similar(_B))

        # α != 1 tests.
        _symm_α = α->symm(side, ul, α, _A, _B)
        _symm_A = A->symm(side, ul, _α, A, _B)
        _symm_B = B->symm(side, ul, _α, _A, B)
        _∇symm_α = (Y, Ȳ, α)->∇(symm, Val{3}, (), Y, Ȳ, side, ul, α, _A, _B)
        _∇symm_A = (Y, Ȳ, A)->∇(symm, Val{4}, (), Y, Ȳ, side, ul, _α, A, _B)
        _∇symm_B = (Y, Ȳ, B)->∇(symm, Val{5}, (), Y, Ȳ, side, ul, _α, _A, B)
        _∇symm_A_inp = (Y, Ȳ, A)->∇(copy(_δA), symm, Val{4}, (), Y, Ȳ, side, ul, _α, A, _B) - _δA
        _∇symm_B_inp = (Y, Ȳ, B)->∇(copy(_δB), symm, Val{5}, (), Y, Ȳ, side, ul, _α, _A, B) - _δB

        B_gen = side == 'L' ? mPQ : mQP
        @test check_errs(N, _symm_α, _∇symm_α, B_gen, sc, sc)
        @test check_errs(N, _symm_A, _∇symm_A, B_gen, mPP, mPP)
        @test check_errs(N, _symm_A, _∇symm_A_inp, B_gen, mPP, mPP)
        @test check_errs(N, _symm_B, _∇symm_B, B_gen, B_gen, B_gen)
        @test check_errs(N, _symm_B, _∇symm_B_inp, B_gen, B_gen, B_gen)

        # α = 1 tests.
        _symm_A = A->symm(side, ul, A, _B)
        _symm_B = B->symm(side, ul, _A, B)
        _∇symm_A = (Y, Ȳ, A)->∇(symm, Val{3}, (), Y, Ȳ, side, ul, A, _B)
        _∇symm_B = (Y, Ȳ, B)->∇(symm, Val{4}, (), Y, Ȳ, side, ul, _A, B)
        _∇symm_A_inp = (Y, Ȳ, A)->∇(copy(_δA), symm, Val{3}, (), Y, Ȳ, side, ul, A, _B) - _δA
        _∇symm_B_inp = (Y, Ȳ, B)->∇(copy(_δB), symm, Val{4}, (), Y, Ȳ, side, ul, _A, B) - _δB

        B_gen = side == 'L' ? mPQ : mQP
        @test check_errs(N, _symm_A, _∇symm_A, B_gen, mPP, mPP)
        @test check_errs(N, _symm_A, _∇symm_A_inp, B_gen, mPP, mPP)
        @test check_errs(N, _symm_B, _∇symm_B, B_gen, B_gen, B_gen)
        @test check_errs(N, _symm_B, _∇symm_B_inp, B_gen, B_gen, B_gen)
    end

    # trmm / trsm:
    for f in [trmm, trsm], side in ['L', 'R'], ul in ['L', 'U'], tA in ['N', 'T'], dA in ['U', 'N']
        A_gen = ()->mPP() + 1e-3I
        B_gen = side == 'L' ? mPQ : mQP
        C_gen = B_gen

        _α, _A, _B = sc(), A_gen(), B_gen()
        _trmm_α = α->f(side, ul, tA, dA, α, _A, _B)
        _trmm_A = A->f(side, ul, tA, dA, _α, A, _B)
        _trmm_B = B->f(side, ul, tA, dA, _α, _A, B)
        _∇trmm_α = (Y, Ȳ, α)->∇(f, Val{5}, (), Y, Ȳ, side, ul, tA, dA, α, _A, _B)
        _∇trmm_A = (Y, Ȳ, A)->∇(f, Val{6}, (), Y, Ȳ, side, ul, tA, dA, _α, A, _B)
        _∇trmm_B = (Y, Ȳ, B)->∇(f, Val{7}, (), Y, Ȳ, side, ul, tA, dA, _α, _A, B)

        @test check_errs(N, _trmm_α, _∇trmm_α, C_gen, sc, sc)
        @test check_errs(N, _trmm_A, _∇trmm_A, C_gen, A_gen, A_gen)
        @test check_errs(N, _trmm_B, _∇trmm_B, C_gen, B_gen, B_gen)
    end
end
end
