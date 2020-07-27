using Test
@time @testset "Meshing" begin include("test_meshing.jl") end
@time @testset "Basics" begin include("test_basics.jl") end
@time @testset "Miscellaneous" begin include("test_miscellaneous.jl") end
@time @testset "Matrix multiplication" begin include("test_matrix_multiplication.jl") end
true
