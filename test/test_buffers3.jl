module mbuffers13

using FinEtools
import FinEtools.MatrixUtilityModule: add_mggt_ut_only!, complete_lt!
using LinearAlgebra.BLAS
using LinearAlgebra: Transpose, mul!
using StaticArrays
using BenchmarkTools

function makebuffers(N, mdim)
    gradN = MMatrix{N, 2, Float64}(rand(N, 2))
    Ke = MMatrix{N, N, Float64}(fill(0.0, size(gradN, 1), size(gradN, 1)))

    # gradN = MMatrix{N, 2, Float64}(rand(N, 2))
    # if N <= 20
    #     Ke = MMatrix{N, N, Float64}(fill(0.0, size(gradN, 1), size(gradN, 1)))
    # else  
    #     Ke = fill(0.0, size(gradN, 1), size(gradN, 1))
    # end

    # gradN = rand(N, 2)
    # Ke = fill(0.0, size(gradN, 1), size(gradN, 1))
    return gradN, Ke
end 

function add_mggt_ut_only_s!(Ke::T1, gradN::T2, mult::FFlt) where {T1, T2}
    Kedim = size(Ke, 1)
    @assert Kedim == size(Ke, 2) # Square matrix?
    nne, mdim = size(gradN)
    @assert nne == Kedim # compatible matrices?
    @inbounds for nx = 1:Kedim # Do: Ce  =  Ce + gradN*((Jac*w[j]))*gradN' ;
        @inbounds for px = 1:mdim
            a = (mult)*gradN[nx, px]
            @inbounds for mx = 1:nx # only the upper triangle
                Ke[mx, nx] +=  gradN[mx, px] * a
            end
        end
    end
    return true
end

function test_tsa(N, NLOOP)
    gradN, Ke = makebuffers(N, 2)
    multiplier = 2.0
    @code_warntype add_mggt_ut_only_s!(Ke, gradN, multiplier)
    t = @elapsed for i = 1:NLOOP
        add_mggt_ut_only_s!(Ke, gradN, multiplier)
    end
    return t ./ NLOOP
end 

function test(N)
    println("N = $(N)")
    NLOOP = 1
    tsa = test_tsa(N, NLOOP)
    vec([tsa])
end

end

using .mbuffers13

NS = [3, 25] # , 225, 900
ts = []
for N in NS
    push!(ts, mbuffers13.test(N))
end 
@show ts
