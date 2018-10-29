
module mmmmtest_trirules
using FinEtools
using Test
import LinearAlgebra: cholesky
function test()
    for NPTS = [1, 3, 4, 6, 7, 9, 12, 13]
        r = TriRule(NPTS)
        # @show NPTS
        @test abs(sum(r.weights) - 0.5) <1.e-5
    end
    true
end
end
using .mmmmtest_trirules
mmmmtest_trirules.test()