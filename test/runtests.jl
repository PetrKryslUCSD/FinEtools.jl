using Test

@testset "Basics" begin
    include("test_basics.jl")
end

@testset "Miscellaneous" begin
    include("test_miscellaneous.jl")
end
@testset "Miscell. 2" begin
    include("test_miscellaneous2.jl")
end

@testset "Matrix ops" begin
    include("test_matrix_multiplication.jl")
end

@testset "Meshing" begin
    include("test_meshing.jl")
end
@testset "Meshing 2" begin
    include("test_meshing_2.jl")
end

@testset "Forms" begin
    include("test_forms.jl")
end

true
