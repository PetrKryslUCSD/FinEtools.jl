using Base.Test
@testset "Miscellaneous" begin include("test_miscellaneous.jl") end
@testset "Acoustics" begin include("test_acoustics.jl") end
@testset "Heat diffusion" begin include("test_heat.jl") end
@testset "Linear deformation" begin include("test_linear_deformation.jl") end
