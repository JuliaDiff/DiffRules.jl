@testset "testimports" begin
    # Check that everything is imported without an error. This is quite a weak criterion,
    # but I'm not sure what else is easily doable.
    for op in DLA.ops
        @eval $(import_expr(op))
    end
end
