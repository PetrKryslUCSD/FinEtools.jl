using Test
# @time @testset "Debug" begin include("test_debug.jl") end
# @time @testset "Acoustics" begin include("test_acoustics.jl") end
# @time @testset "Heat diffusion" begin include("test_heat.jl") end
# @time @testset "Linear deformation" begin include("test_linear_deformation.jl") end
@time @testset "Meshing" begin include("test_meshing.jl") end
@time @testset "Miscellaneous" begin include("test_miscellaneous.jl") end
@time @testset "Matrix multiplication" begin include("test_matrix_multiplication.jl") end
# @time @testset "Voxel box" begin include("test_voxel_box.jl") end
# @time @testset "Failing tests" begin include("test_failing.jl") end
true
