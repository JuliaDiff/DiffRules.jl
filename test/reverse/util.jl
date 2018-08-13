@testset "util" begin
    import DiffLinearAlgebra: importable
    # @test Expr(:import, Expr(Symbol("."), importable(:Foo)...)) ==
    #     :(import Foo)
    # @test Expr(:import, Expr(Symbol("."), importable(:(Package.Foo))...)) ==
    #     :(import Package.Foo)
    # @test Expr(:import, Expr(Symbol("."), importable(:(Package.Subpackage.Foo))...)) ==
    #     :(import Package.Subpackage.Foo)

    @test import_expr(DLA.DiffOp(:Foo, :(Tuple{Foo}), [true])) ==
        :(import Foo)
    @test import_expr(DLA.DiffOp(:(Package.Foo), :(Tuple{Foo}), [true])) ==
        :(import Package.Foo)
    @test import_expr(DLA.DiffOp(:(Package.Subpackage.Foo), :(Tuple{Foo}), [true])) ==
        :(import Package.Subpackage.Foo)
end
