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

function test_tlt(N, NLOOP)
    gradN = rand(N, 2)
    Ke = fill(0.0, size(gradN, 1), size(gradN, 1))
    multiplier = 2.0
    t = @elapsed for i = 1:NLOOP
        complete_lt!(Ke)
    end
    return t ./ NLOOP
end 

function test_twi(N, NLOOP)
    gradN = rand(N, 2)
    Ke = fill(0.0, size(gradN, 1), size(gradN, 1))
    multiplier = 2.0
    t = @elapsed for i = 1:NLOOP
        add_mggt_ut_only!(Ke, gradN, multiplier)
    end
    return t ./ NLOOP
end 

function test_tgd(N, NLOOP)
    gradN = rand(N, 2)
    Ke = fill(0.0, size(gradN, 1), size(gradN, 1))
    multiplier = 2.0
    t = @elapsed for i = 1:NLOOP
        BLAS.gemm!('N', 'T', multiplier, gradN, gradN, 0.0, Ke)
    end
    return t ./ NLOOP
end 

function test_tsa(N, NLOOP)
    gradN, Ke = makebuffers(N, 2)
    multiplier = 2.0
    # @code_warntype add_mggt_ut_only_s!(Ke, gradN, multiplier)
    gradNt = gradN'
    t = @elapsed for i = 1:NLOOP
        # Ke .= gradN * (gradNt * multiplier)
        add_mggt_ut_only_s!(Ke, gradN, multiplier)
    end
    return t ./ NLOOP
end 

function test(N)
    println("N = $(N)")
    NLOOP = Int64(round(10 / N^2)) + 10
    tlt = test_tlt(N, NLOOP)
    twi = test_twi(N, NLOOP)
    tgd = test_tgd(N, NLOOP)
    tsa = test_tsa(N, NLOOP)
    vec([tlt twi tgd tsa])
end

end

using .mbuffers13
using Gaston
set(axis="loglog", plotstyle="linespoints", linewidth=2, pointsize = 1, color = "black", xlabel = "N", ylabel = "Time [microseconds]", grid="on", title = "")


NS = [3, 9, 16, 25, 36] # , 225, 900
ts = []
for N in NS
    push!(ts, mbuffers13.test(N))
end 
@show ts
f = figure()
# TS = [1.0e6 * t[1] for t in ts] # Time in Microseconds
# plot(NS, TS, legend = "Complete triangle")
TS = [1.0e6 * (t[2] + t[1]) for t in ts] # Time in Microseconds
plot(NS, TS, legend = "Loops inbounds", gpcom = """set terminal wxt font ",6" """, box = "left top")
TS = [1.0e6 * t[3] for t in ts] # Time in Microseconds
plot!(NS, TS, legend = "gemm!" )
TS = [1.0e6 * (t[4] + t[1]) for t in ts] # Time in Microseconds
plot!(NS, TS, legend = "Loops opt/static" )
figure(f)
