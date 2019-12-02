module maddbts1
using FinEtools
using FinEtools.MatrixUtilityModule: add_btsigma!
using LinearAlgebra: norm
using BenchmarkTools
using Test
function test()
	for i in 1:1
		B = rand(6, 30)
		s = rand(6)
		c = rand()
		F1 = fill(0.0, size(B, 2))
		@btime $F1 .+= $B' * ($c * $s)
		F2 = fill(0.0, size(B, 2))
		@btime add_btsigma!($F2, $B, $c, $s)
	end
	true
end
end
using .maddbts1
maddbts1.test()