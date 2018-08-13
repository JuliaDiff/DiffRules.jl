@testset "Cholesky" begin
    import DiffLinearAlgebra: level2partition, level3partition, chol_unblocked_rev,
        chol_blocked_rev

    let rng = MersenneTwister(123456), N = 5
        A = randn(rng, N, N)
        r, d, B2, c = level2partition(A, 4, false)
        R, D, B3, C = level3partition(A, 4, 4, false)
        @test all(r .== R')
        @test all(d .== D)
        @test B2[1] == B3[1]
        @test all(c .== C)

        # Check that level2partition with 'U' is consistent with 'L'.
        rᵀ, dᵀ, B2ᵀ, cᵀ = level2partition(transpose(A), 4, true)
        @test r == rᵀ
        @test d == dᵀ
        @test B2' == B2ᵀ
        @test c == cᵀ

        # Check that level3partition with 'U' is consistent with 'L'.
        R, D, B3, C = level3partition(A, 2, 4, false)
        Rᵀ, Dᵀ, B3ᵀ, Cᵀ = level3partition(transpose(A), 2, 4, true)
        @test transpose(R) == Rᵀ
        @test transpose(D) == Dᵀ
        @test transpose(B3) == B3ᵀ
        @test transpose(C) == Cᵀ
    end

    let rng = MersenneTwister(123456), N = 10
        A, Ā = Matrix.(LowerTriangular.(randn.(Ref(rng), [N, N], [N, N])))
        B, B̄ = copy.(transpose.([A, Ā]))
        @test chol_unblocked_rev(Ā, A, false) ≈ chol_blocked_rev(Ā, A, 1, false)
        @test chol_unblocked_rev(Ā, A, false) ≈ chol_blocked_rev(Ā, A, 3, false)
        @test chol_unblocked_rev(Ā, A, false) ≈ chol_blocked_rev(Ā, A, 5, false)
        @test chol_unblocked_rev(Ā, A, false) ≈ chol_blocked_rev(Ā, A, 10, false)
        @test chol_unblocked_rev(Ā, A, false) ≈ chol_unblocked_rev(B̄, B, true)'

        @test chol_unblocked_rev(B̄, B, true) ≈ chol_blocked_rev(B̄, B, 1, true)
        @test chol_unblocked_rev(B̄, B, true) ≈ chol_blocked_rev(B̄, B, 5, true)
        @test chol_unblocked_rev(B̄, B, true) ≈ chol_blocked_rev(B̄, B, 10, true)
    end

    # Check sensitivities for lower-triangular version.
    let P = 15, rng = MersenneTwister(123456), N = 10

        Σ_gen = ()->(A = randn(rng, P, P); A'A + 1e-3I)
        S_gen, H_gen = ()->Symmetric(Σ_gen()), ()->Hermitian(Σ_gen())
        U_gen = ()->chol(S_gen())

        _chol = Σ->chol(Symmetric(Σ))
        _∇chol = (Y, Ȳ, Σ)->∇(chol, Val{1}, (), Y, Ȳ, Symmetric(Σ))

        @test check_errs(N, _chol, _∇chol, U_gen, Σ_gen, Σ_gen)
        @test check_errs(N, _chol, _∇chol, U_gen, H_gen, H_gen)
        @test check_errs(N, _chol, _∇chol, U_gen, S_gen, S_gen)
    end
end
