__precompile__(true)

module DiffLinearAlgebra

    using LinearAlgebra

    # Some aliases used repeatedly throughout the package.
    const AV, AM, AVM, AA = AbstractVector, AbstractMatrix, AbstractVecOrMat, AbstractArray
    const SV, SM, SVM, SA = StridedVector, StridedMatrix, StridedVecOrMat, StridedArray
    const AS, ASVM = Union{Real, AA}, Union{Real, AVM}
    const Arg1, Arg2, Arg3 = Type{Val{1}}, Type{Val{2}}, Type{Val{3}}
    const Arg4, Arg5, Arg6 = Type{Val{4}}, Type{Val{5}}, Type{Val{6}}
    const BF = Union{Float32, Float64}
    const DLA = DiffLinearAlgebra

    export ∇, DLA, import_expr

    # Define container and meta-data for defined operations.
    include("util.jl")

    # Define operations, and log them for external use via `ops`.
    include("generic.jl")
    include("blas.jl")
    include("diagonal.jl")
    include("triangular.jl")
    include("uniformscaling.jl")
    include("factorization/cholesky.jl")

end # module
