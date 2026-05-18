module testmortar
using FinEtools
using Test
function test()
    @test common_refinement() === nothing
    return true
end
end
using .testmortar
testmortar.test()