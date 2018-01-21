# Tested 01/18/2017
# $  ~/AppData/Local/Julia-0.6.1/bin/julia.exe
#                _
#    _       _ _(_)_     |  A fresh approach to technical computing
#   (_)     | (_) (_)    |  Documentation: https://docs.julialang.org
#    _ _   _| |_  __ _   |  Type "?help" for help.
#   | | | | | | |/ _` |  |
#   | | |_| | | | (_| |  |  Version 0.6.1 (2017-10-24 22:15 UTC)
#  _/ |\__'_|_|_|\__'_|  |  Official http://julialang.org/ release
# |__/                   |  x86_64-w64-mingw32

# julia> Pkg.test("FinEtools")
# INFO: Testing FinEtools
# WARNING: Method definition Test Summary: | pairs(Pass  AnyTotal)
#  in module FEMMBaseModule at C:\Users\PetrKrysl\.julia\v0.6\FinEtools\src\FEMMBaseModule.jl:9Miscellaneous overwritten at C:\Users\PetrKrysl\.julia\v0.6\FinEtools\src\FEMMBaseModule.jl:34 | .
#  100    100
# timing = 27.314000129699707
# Test Summary: | Pass  Total
# Acoustics     |   20     20
# timing = 76.84099984169006
# Test Summary:  | Pass  Total
# Heat diffusion |   25     25
# timing = 25.91600012779236
# Test Summary:      | Pass  Total
# Linear deformation |  159    159
# timing = 90.73000001907349
# Test Summary: | Pass  Total
# Meshing       |   97     97
# timing = 661.904000043869
# Test Summary: | Pass  Total
# Voxel box     |   27     27
# timing = 5.4039998054504395
# INFO: FinEtools tests passed

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