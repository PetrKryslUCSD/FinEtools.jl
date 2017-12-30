using Compat.Test
t0 = time()
@testset "Miscellaneous" begin include("test_miscellaneous.jl") end
println("timing = $(time() - t0)")
t0 = time()
@testset "Acoustics" begin include("test_acoustics.jl") end
println("timing = $(time() - t0)")
t0 = time()
@testset "Heat diffusion" begin include("test_heat.jl") end
println("timing = $(time() - t0)")
t0 = time()
@testset "Linear deformation" begin include("test_linear_deformation.jl") end
println("timing = $(time() - t0)")
t0 = time()
@testset "Meshing" begin include("test_meshing.jl") end
println("timing = $(time() - t0)")
t0 = time()
@testset "Voxel box" begin include("test_voxel_box.jl") end
println("timing = $(time() - t0)")