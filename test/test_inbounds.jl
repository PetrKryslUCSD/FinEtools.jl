module moinbounds13
using FinEtools
import FinEtools.MatrixUtilityModule: add_mggt_ut_only!, complete_lt!
using LinearAlgebra: Transpose, mul!, BLAS
using StaticArrays
using BenchmarkTools

struct Vec{N,T} <: AbstractVector{T}
    # data::NTuple{N,Core.VecElement{T}}
    data::NTuple{N,T}
    Vec(data::Vararg{T,N}) where {N,T} = new{N,T}(data)
end

@inline Base.getindex(x::Vec, i) = x.data[i]
@inline Base.length(::Vec{N}) where N = N
@inline Base.size(::Vec{N}) where N = (N,)
@inline Base.eltype(::Vec{N,T}) where {N,T} = T
@inline function vload(::Type{Vec{N,T}}, x::AbstractArray{T}, i) where {N,T}
    unsafe_load(Base.unsafe_convert(Ptr{Vec{N,T}}, pointer(x)) + sizeof(T)*i)
end
@inline function vstore!(x::AbstractArray{T}, v::Vec{N,T}, i) where {N,T}
    unsafe_store!(Base.unsafe_convert(Ptr{Vec{N,T}}, pointer(x)) + sizeof(T)*i, v)
end

function create_quote()
    q = quote
        @inbounds @fastmath begin
            $(Expr(:meta, :inline))
            Vec()
        end
    end
    q, q.args[2].args[3].args[3].args[4].args
end

@generated function Base.fma(x::Vec{N,T}, y::T, z::Vec{N,T}) where {N,T}
    q, qa = create_quote()
    for n ∈ 1:N
        push!(qa, :(x[$n] * y + z[$n]))
    end
    q
end

@generated function add_mggt_ut_only_g!(Ke::FFltMat, gradN::FFltMat, mult::FFlt, ::Val{mdim}) where mdim
    L = 8
    quote
        @assert size(Ke, 1) == size(Ke, 2)
        Kedim = size(Ke, 1)
        # nne, mdim = size(gradN)
        V = Vec{$L,Float64}
        SIMD_loop_max = Kedim ÷ $L * $L -1
        @inbounds for c = 1:Kedim# ÷ $L * $L
            Base.Cartesian.@nexprs $mdim m -> gradNc_m = gradN[c, m] * mult
            for r = 0:$L:min(c-1,SIMD_loop_max)
                Ktemp = vload(V, Ke, r + (c-1)*Kedim)
                Base.Cartesian.@nexprs $mdim m -> (Ktemp = fma(vload(V, gradN, r+(m-1)*Kedim), gradNc_m, Ktemp))
                vstore!(Ke, Ktemp, r + (c-1)*Kedim)
            end
        end
        remainder = Kedim % $L
        if remainder > 0
            bound = Kedim-remainder+1
            @inbounds for c = bound:Kedim
                Base.Cartesian.@nexprs $mdim m -> gradNc_m = gradN[c, m] * mult
                for r = bound:c
                    Ktemp = Ke[r,c]
                    Base.Cartesian.@nexprs $mdim m -> (Ktemp = fma(gradN[r,m], gradNc_m, Ktemp))
                    Ke[r,c] = Ktemp
                end
            end
        end
        return true
    end
end

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

function add_mggt_ut_only_sa!(Ke::T1, gradN::T2, mult::FFlt) where {T1, T2}
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
    tw2 = @belapsed add_mggt_ut_only_g!($Ke, $gradN, 1.0, Val(2))
    gradN = MMatrix{N, 2, Float64}(rand(N, 2))
    Ke = MMatrix{N, N, Float64}(fill(0.0, size(gradN, 1), size(gradN, 1)))
    tsa = @belapsed add_mggt_ut_only_sa!($Ke, $gradN, 1.0)
    vec([tlt twi tsn tmn two tsd tmd tgd tw2 tsa])
end

end
using .moinbounds13
using Gaston
set(axis="loglog", plotstyle="linespoints", linewidth=2, pointsize = 1, color = "black", xlabel = "N", ylabel = "Time [microseconds]", grid="on", title = "")


NS = [3, 9, 16, 25, 36, 81] # , 225, 900
ts = []
for N in NS
    push!(ts, moinbounds13.test(N))
end 
@show ts
f = figure()
# TS = [1.0e6 * t[1] for t in ts] # Time in Microseconds
# plot(NS, TS, legend = "Complete triangle")
TS = [1.0e6 * (t[2] + t[1]) for t in ts] # Time in Microseconds
plot(NS, TS, legend = "Loops inbounds", gpcom = """set terminal wxt font ",6" """, box = "left top")
# TS = [1.0e6 * t[3] for t in ts] # Time in Microseconds
# plot!(NS, TS, legend = "Matrix mult" )
# TS = [1.0e6 * t[4] for t in ts] # Time in Microseconds
# plot!(NS, TS, legend = "mul!" )
TS = [1.0e6 * (t[5] + t[1]) for t in ts] # Time in Microseconds
plot!(NS, TS, legend = "Loops NO inbounds" )
TS = [1.0e6 * t[6] for t in ts] # Time in Microseconds
plot!(NS, TS, legend = "Matrix mult w/ dot" )
TS = [1.0e6 * t[7] for t in ts] # Time in Microseconds
plot!(NS, TS, legend = "mul! w/ dot" )
TS = [1.0e6 * t[8] for t in ts] # Time in Microseconds
plot!(NS, TS, legend = "gemm!" )
TS = [1.0e6 * t[9] for t in ts] # Time in Microseconds
plot!(NS, TS, legend = "Loops generated" )
TS = [1.0e6 * t[10] for t in ts] # Time in Microseconds
plot!(NS, TS, legend = "Loops static" )
figure(f)
