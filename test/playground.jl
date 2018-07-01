using Test

# module mmsparsem1
# using FinEtools
# using SparseArrays
# using Test
# function test()
#     N = 10
#     a = vec(rand(N))
#     i = convert(Vector{Int64}, round.(rand(N) * N) .+ 1) 
#     j = convert(Vector{Int64}, round.(rand(N) * N) .+ 1) 
#     # println("i = $(i)")
#     # println("j = $(j)")
#     A = sparse(i, j, a, N+1, N+1)
#     @test typeof(A * A) <: AbstractSparseMatrix
#     true
# end
# end
# using .mmsparsem1
# mmsparsem1.test()


@time @testset "Miscellaneous" begin include("test_miscellaneous.jl") end
@time @testset "Acoustics" begin include("test_acoustics.jl") end
@time @testset "Heat diffusion" begin include("test_heat.jl") end
@time @testset "Linear deformation" begin include("test_linear_deformation.jl") end
@time @testset "Meshing" begin include("test_meshing.jl") end
@time @testset "Voxel box" begin include("test_voxel_box.jl") end

