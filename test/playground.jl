module msymmx
using FinEtools
using FinEtools.MatrixUtilityModule: symmetrize!
using LinearAlgebra
using Test
function test()
    a = rand(6, 6)
    b = deepcopy(a)
    symmetrize!(a)
    @test norm(a - (b + b')/2) / norm(a) <= 1.0e-9
    true
end
end
using .msymmx
msymmx.test()