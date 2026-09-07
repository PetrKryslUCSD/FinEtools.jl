module testmortar_linLM
using FinEtools
using Test

function test()
        fensA, fesA = T3block(1.0, 1.0, 1, 1, :a)
        fensB, fesB = T3block(1.0, 1.0, 1, 1, :b)


    D, meta = common_refinement(fensA, fesA, fensB, fesB; lam_order = 1, tri_order = 1, dim_u=2)
        @test size(D, 1) == 8
        @test size(D, 2) == 8
        @test abs(sum(meta["M"])-2) < 1e-12
    return true
end
end

using .testmortar_linLM
testmortar_linLM.test()

module testmortar_pwcLM
using FinEtools
using Test

function test()
        fensA, fesA = T3block(1.0, 1.0, 1, 1, :a)
        fensB, fesB = T3block(1.0, 1.0, 1, 1, :b)


    D, meta = common_refinement(fensA, fesA, fensB, fesB; lam_order = 0, tri_order = 1, dim_u=3)
        @test size(D, 1) == 6
        @test size(D, 2) == 12
        @test abs(sum(meta["M"])-3) < 1e-12
    return true
end
end

using .testmortar_pwcLM
testmortar_pwcLM.test()


# module testmortar_t6tri
# using FinEtools
# using Test

# function test()
#         fensA, fesA = T6block(1.0, 1.0, 1, 1, :a)
#         fensB, fesB = T6block(1.0, 1.0, 1, 1, :a)


#     D, meta = common_refinement(fensA, fesA, fensB, fesB; lam_order = 0, tri_order = 2, dim_u=2)
#         @test size(D, 1) == 4
#         @test size(D, 2) == 18
#         @test abs(sum(meta["M"])-2) < 1e-12
#     return true
# end
# end

# using .testmortar_t6tri
# testmortar_t6tri.test()