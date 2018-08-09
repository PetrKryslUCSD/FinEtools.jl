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

# julia> versioninfo()
# Julia Version 0.6.2
# Commit d386e40c17* (2017-12-13 18:08 UTC)
# Platform Info:
#   OS: Windows (x86_64-w64-mingw32)
#   CPU: Intel(R) Core(TM) i7-6650U CPU @ 2.20GHz
#   WORD_SIZE: 64
#   BLAS: libopenblas (USE64BITINT DYNAMIC_ARCH NO_AFFINITY Haswell)
#   LAPACK: libopenblas64_
#   LIBM: libopenlibm
#   LLVM: libLLVM-3.9.1 (ORCJIT, skylake)
# julia> Pkg.test("FinEtools")
# INFO: Testing FinEtools
# Test Summary: | Pass  Total
# Miscellaneous |  107    107
#  38.122776 seconds (34.92 M allocations: 3.518 GiB, 2.47% gc time)
# Test Summary: | Pass  Total
# Acoustics     |   20     20
#  40.384056 seconds (138.85 M allocations: 6.875 GiB, 3.09% gc time)
# Test Summary:  | Pass  Total
# Heat diffusion |   25     25
#  26.885352 seconds (129.07 M allocations: 3.836 GiB, 4.71% gc time)
# Test Summary:      | Pass  Total
# Linear deformation |  159    159
#  92.007878 seconds (331.79 M allocations: 17.178 GiB, 7.56% gc time)
# Test Summary: | Pass  Total
# Meshing       |   97     97
# 220.561762 seconds (2.32 G allocations: 122.789 GiB, 19.69% gc time)
# Test Summary: | Pass  Total
# Voxel box     |   27     27
#   7.464730 seconds (32.47 M allocations: 3.291 GiB, 16.59% gc time)
# INFO: FinEtools tests passed

# julia> versioninfo()
# Julia Version 0.7.0-DEV.3404
# Commit d569a2923c* (2018-01-14 21:52 UTC)
# Platform Info:
#   OS: Windows (x86_64-w64-mingw32)
#   CPU: Intel(R) Core(TM) i7-6650U CPU @ 2.20GHz
#   WORD_SIZE: 64
#   BLAS: libopenblas (USE64BITINT DYNAMIC_ARCH NO_AFFINITY Haswell)
#   LAPACK: libopenblas64_
#   LIBM: libopenlibm
#   LLVM: libLLVM-3.9.1 (ORCJIT, skylake)
# Environment:

# julia> Pkg.test("FinEtools")
# [ Info: Testing FinEtools
# Test Summary: | Pass  Total
# Miscellaneous |  107    107
#  49.633658 seconds (66.62 M allocations: 5.368 GiB, 3.51% gc time)
# Test Summary: | Pass  Total
# Acoustics     |   20     20
#  58.462955 seconds (148.15 M allocations: 8.753 GiB, 3.70% gc time)
# Test Summary:  | Pass  Total
# Heat diffusion |   25     25
#  25.859226 seconds (132.25 M allocations: 4.375 GiB, 5.74% gc time)
# Test Summary:      | Pass  Total
# Linear deformation |  159    159
# 105.097057 seconds (350.96 M allocations: 18.821 GiB, 7.23% gc time)
# Test Summary: | Pass  Total
# Meshing       |   97     97
# 197.330767 seconds (1.42 G allocations: 98.406 GiB, 15.45% gc time)
# Test Summary: | Pass  Total
# Voxel box     |   27     27
#   4.020409 seconds (24.34 M allocations: 1.920 GiB, 13.51% gc time)
# [ Info: FinEtools tests passed

using Test
# @time @testset "Debug" begin include("test_debug.jl") end
@time @testset "Miscellaneous" begin include("test_miscellaneous.jl") end
@time @testset "Acoustics" begin include("test_acoustics.jl") end
@time @testset "Heat diffusion" begin include("test_heat.jl") end
@time @testset "Linear deformation" begin include("test_linear_deformation.jl") end
@time @testset "Meshing" begin include("test_meshing.jl") end
@time @testset "Voxel box" begin include("test_voxel_box.jl") end
# @time @testset "Failing tests" begin include("test_failing.jl") end
true
