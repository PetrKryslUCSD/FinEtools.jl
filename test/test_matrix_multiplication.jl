module mmmultest1
using Test
# using BenchmarkTools
using LoopVectorization
using LinearAlgebra
using FinEtools.MatrixUtilityModule: mulCAB!

function gemmavx!(C, A, B)
	M, N = size(C); K = size(B,1)
	C .= 0
    @avx for n ∈ 1:N, k ∈ 1:K
    	Bkn = B[k,n]
        for m ∈ 1:M
            C[m,n] += A[m,k] * Bkn
        end
    end
    return C
end

function gemmblas!(C, A, B)
    mul!(C, A, B)
    return C
end

function test(n)
	# println("n = ", n)
	C, A, B = rand(n, n), rand(n, n), rand(n, n)
	@test norm(mulCAB!(C, A, B) - A * B) / norm(C) <= 1.0e-9
	@test norm(gemmblas!(C, A, B) - A * B) / norm(C) <= 1.0e-9
	@test norm(gemmavx!(C, A, B) - A * B) / norm(C) <= 1.0e-9

	BLAS.set_num_threads(1)
	# println("mulCAB!")
	
# @btime mulCAB!($C, $A, $B)
	# println("gemmblas!")
	
# @btime gemmblas!($C, $A, $B)
	# println("gemmavx!")
	
# @btime gemmavx!($C, $A, $B)
end
end
using .mmmultest1

mmmultest1.test(2^6)
mmmultest1.test(2^5)
mmmultest1.test(2^4)
mmmultest1.test(2^3)

mmmultest1.test(3)
mmmultest1.test(5)
mmmultest1.test(7)


module mmmultest2
using Test
# using BenchmarkTools
using LoopVectorization
using LinearAlgebra
using FinEtools.MatrixUtilityModule: mulCAB!

function gemmavx!(C, A, B)
	M, N = size(C); K = size(B,1)
	C .= 0
    @avx  for n ∈ 1:N, k ∈ 1:K # for k ∈ 1:K, n ∈ 1:N #
    	Bkn = B[k,n]
        for m ∈ 1:M
            C[m,n] += A[m,k] * Bkn
        end
    end
    return C
end

function gemmblas!(C, A, B)
    mul!(C, A, B)
    return C
end



function test(n)
	# println("n = ", n)
	C, A, B = rand(n, n), rand(n, n), rand(n, n)
	@test norm(mulCAB!(C, A, B) - A * B) / norm(C) <= 1.0e-9
	@test norm(gemmblas!(C, A, B) - A * B) / norm(C) <= 1.0e-9
	@test norm(gemmavx!(C, A, B) - A * B) / norm(C) <= 1.0e-9

	BLAS.set_num_threads(1)
	# println("mulCAB!")
	
# @btime mulCAB!($C, $A, $B)
	# println("gemmblas!")
	
# @btime gemmblas!($C, $A, $B)
	# println("gemmavx!")
	
# @btime gemmavx!($C, $A, $B)
end
end
using .mmmultest2

mmmultest2.test(2^6)
mmmultest2.test(2^5)
mmmultest2.test(2^4)
mmmultest2.test(2^3)


module mmmultest3
using Test
# using BenchmarkTools
using LoopVectorization
using LinearAlgebra
using FinEtools.MatrixUtilityModule: mulCAtB!


function gemmblas!(C, A, B)
    mul!(C, Transpose(A), B)
    return C
end



function test(n)
	# println("n = ", n)
	C, A, B = rand(n, n), rand(n, n), rand(n, n)
	mulCAtB!(C, A, B)
	@test norm(mulCAtB!(C, A, B) - A' * B) / norm(C) <= 1.0e-9
	@test norm(gemmblas!(C, A, B) - A' * B) / norm(C) <= 1.0e-9

	# println("mulCAtB!")
	
# @btime mulCAtB!($C, $A, $B)
	# println("gemmblas!")
	
# @btime gemmblas!($C, $A, $B)
end
end
using .mmmultest3
mmmultest3.test(2^6)
mmmultest3.test(2^5)
mmmultest3.test(2^4)
mmmultest3.test(2^3)

# mlatest1.test(3)
# mlatest1.test(5)
# mlatest1.test(7)

module mmmultest4
using Test
# using BenchmarkTools
using LoopVectorization
using LinearAlgebra
using FinEtools.MatrixUtilityModule: mulCAtB!

function gemmblas!(C, A, B)
    mul!(C, Transpose(A), B)
    return C
end



function test(M, N, K)
	# println("M, N, K = ", (M, N, K))
	C, A, B = rand(M, N), rand(K, M), rand(K, N)
	mulCAtB!(C, A, B)
	@test norm(mulCAtB!(C, A, B) - A' * B) / norm(C) <= 1.0e-9
	@test norm(gemmblas!(C, A, B) - A' * B) / norm(C) <= 1.0e-9

	# println("mulCAtB!")
	
# @btime mulCAtB!($C, $A, $B)
	# println("gemmblas!")
	
# @btime gemmblas!($C, $A, $B)
end
end
using .mmmultest4
M, N, K = 2^6, 2^6, 2^6
mmmultest4.test(M, N, K)
M, N, K = 2^5, 2^5, 2^5
mmmultest4.test(M, N, K)
M, N, K = 2^4, 2^4, 2^4
mmmultest4.test(M, N, K)
M, N, K = 2^3, 2^3, 2^3
mmmultest4.test(M, N, K)

M, N, K = 3, 4, 5
mmmultest4.test(M, N, K)
M, N, K = 13, 4, 5
mmmultest4.test(M, N, K)
M, N, K = 6, 4, 9
mmmultest4.test(M, N, K)

module mmmultest6
using Test
# using BenchmarkTools
using LoopVectorization
using LinearAlgebra
using FinEtools.MatrixUtilityModule: mulCABt!

function gemmblas!(C, A, B)
    mul!(C, A, Transpose(B))
    return C
end



function test(M, N, K)
	# println("M, N, K = ", (M, N, K))
	C, A, B = rand(M, N), rand(M, K), rand(N, K)
	@test norm(mulCABt!(C, A, B) - A * B') / norm(C) <= 1.0e-9
	@test norm(gemmblas!(C, A, B) - A * B') / norm(C) <= 1.0e-9

	# println("mulCABt!")
	
# @btime mulCABt!($C, $A, $B)
	# println("gemmblas!")
	
# @btime gemmblas!($C, $A, $B)
end
end
using .mmmultest6

M, N, K = 3, 4, 5
mmmultest6.test(M, N, K)
M, N, K = 13, 4, 5
mmmultest6.test(M, N, K)
M, N, K = 6, 4, 9
mmmultest6.test(M, N, K)
M, N, K = 7, 4, 9
mmmultest6.test(M, N, K)
M, N, K = 1, 4, 9
mmmultest6.test(M, N, K)
M, N, K = 7, 1, 9
mmmultest6.test(M, N, K)
M, N, K = 1, 1, 9
mmmultest6.test(M, N, K)
M, N, K = 1, 1, 1
mmmultest6.test(M, N, K)
M, N, K = 7, 4, 1
mmmultest6.test(M, N, K)

M, N, K = 2^6, 2^6, 2^6
mmmultest6.test(M, N, K)
M, N, K = 2^5, 2^5, 2^5
mmmultest6.test(M, N, K)
M, N, K = 2^4, 2^4, 2^4
mmmultest6.test(M, N, K)
M, N, K = 2^3, 2^3, 2^3
mmmultest6.test(M, N, K)

module mmvtest1
using FinEtools
using FinEtools.MatrixUtilityModule: mulCAB!
using LinearAlgebra
using Test
function test()
	M, N = 7, 9
	C, A, B = rand(M), rand(M, N), rand(N)
	@test norm(mulCAB!(C, A, B) - A * B) / norm(C) <= 1.0e-9
	C, A, B = rand(M), rand(M, N), rand(N)
	@test norm(mulCAB!(C, A, B) - A * B) / norm(C) <= 1.0e-9
	C, A, B = rand(M), rand(M, N), rand(N)
	@test norm(mulCAB!(C, A, B) - A * B) / norm(C) <= 1.0e-9
end
end
using .mmvtest1
mmvtest1.test()