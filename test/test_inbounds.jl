module moinbounds13
using FinEtools
using FinEtools.MatrixUtilityModule: add_mggt_ut_only!, complete_lt!
using LinearAlgebra: Transpose, mul!
using BenchmarkTools
function add_mggt_ut_only_wo!(Ke::FFltMat, gradN::FFltMat, mult::FFlt)
    @assert size(Ke, 1) == size(Ke, 2)
    Kedim = size(Ke, 1)
    nne, mdim = size(gradN)
    mx::FInt = 0
    nx::FInt = 0
    px::FInt = 0
    for nx = 1:Kedim # Do: Ce  =  Ce + gradN*((Jac*w[j]))*gradN' ;
        for px = 1:mdim
            for mx = 1:nx # only the upper triangle
                Ke[mx, nx] +=  gradN[mx, px]*(mult)*gradN[nx, px]
            end
        end
    end
    return true
end
function test(N)
    println("N = $(N)")
    gradN = rand(N, 2)
    Ke = fill(0.0, size(gradN, 1), size(gradN, 1))
    tlt = @belapsed complete_lt!($Ke)
    twi = @belapsed add_mggt_ut_only!($Ke, $gradN, 1.0)
    tsn = @belapsed 1.0 * ($gradN*Transpose($gradN))
    tmn = @belapsed 1.0 * mul!($Ke, $gradN, Transpose($gradN))
    two = @belapsed add_mggt_ut_only_wo!($Ke, $gradN, 1.0)
    tsd = @belapsed 1.0 .* ($gradN*Transpose($gradN))
    tmd = @belapsed 1.0 .* mul!($Ke, $gradN, Transpose($gradN))
    tgd = @belapsed BLAS.gemm!('N', 'T', 1.0, $gradN, $gradN, 0.0, $Ke)
    vec([tlt twi tsn tmn two tsd tmd tgd])
end
end
using .moinbounds13
using Gaston
set(axis="loglog", plotstyle="linespoints", linewidth=2, pointsize = 1, color = "black", xlabel = "N", ylabel = "Time [microseconds]", grid="on", title = "")


NS = [3, 9, 16, 25, 36, 81, 225, 900]
ts = []
for N in NS
    push!(ts, moinbounds13.test(N))
end 
@show ts
f = figure()
# TS = [1.0e6 * t[1] for t in ts] # Time in Microseconds
# plot(NS, TS, legend = "Complete triangle")
TS = [1.0e6 * (t[2] + t[1]) for t in ts] # Time in Microseconds
plot!(NS, TS, legend = "Loops inbounds" )
TS = [1.0e6 * t[3] for t in ts] # Time in Microseconds
plot!(NS, TS, legend = "Matrix mult" )
TS = [1.0e6 * t[4] for t in ts] # Time in Microseconds
plot!(NS, TS, legend = "mul!" )
TS = [1.0e6 * (t[5] + t[1]) for t in ts] # Time in Microseconds
plot!(NS, TS, legend = "Loops NO inbounds" )
TS = [1.0e6 * t[6] for t in ts] # Time in Microseconds
plot!(NS, TS, legend = "Matrix mult w/ dot" )
TS = [1.0e6 * t[7] for t in ts] # Time in Microseconds
plot!(NS, TS, legend = "mul! w/ dot" )
TS = [1.0e6 * t[8] for t in ts] # Time in Microseconds
plot!(NS, TS, legend = "gemm!" )
figure(f)