module mbuffers13

using FinEtools
import FinEtools.MatrixUtilityModule: add_mggt_ut_only!, complete_lt!
using LinearAlgebra.BLAS
using LinearAlgebra: Transpose, mul!
using StaticArrays
using BenchmarkTools

function test_tda(N, NLOOP)
    gradN = rand(N, 2)
    # Ke = fill(0.0, size(gradN, 1), size(gradN, 1))
    multiplier = 2.0
    t = @elapsed for i = 1:NLOOP
        # Ke .= gradN * (gradN' * multiplier)
        gradN * (gradN' * multiplier)
    end
    return t ./ NLOOP
end 

function test_tsa(N, NLOOP)
    gradN = MMatrix{N, 2, Float64}(rand(N, 2))
    Ke = MMatrix{N, N, Float64}(fill(0.0, size(gradN, 1), size(gradN, 1)))
    multiplier = 2.0
    t = @elapsed for i = 1:NLOOP
        # Ke .= gradN * (gradN' * multiplier)
        gradN * (gradN' * multiplier)
    end
    return t ./ NLOOP
end 

function test(N)
    println("N = $(N)")
    NLOOP = 1
    tda = test_tda(N, NLOOP)
    tsa = test_tsa(N, NLOOP)
    vec([tda tsa])
end

end

using .mbuffers13

NS = [3, 9, 16, 25] # , 225, 900
ts = []
for N in NS
    push!(ts, mbuffers13.test(N))
end 
@show ts
