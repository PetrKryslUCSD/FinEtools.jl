module mbascs06
using FinEtools
using FinEtools.AlgoBaseModule: fieldnorm
using LinearAlgebra
using Test
function test()
    XYZ = reshape([0.2, 0.3, 0.4], 1, 3)
    rfcsmat = reshape([0.37139067635410367; 0.5570860145311555; 0.7427813527082073], 3, 1)
    tangents = reshape([0.2 0.0; 0.0 0.5; 0.4 0.2], 3, 2)
    csys = CSys(3, 2)
    updatecsmat!(csys, XYZ, tangents, 0)
    # @show csys.csmat
    n1 = cross(vec(tangents[:, 1]), vec(tangents[:, 2]))
    n1 = n1 / norm(n1)
    # @show n1
    n2 = cross(vec(csys.csmat[:, 1]), vec(csys.csmat[:, 2]))
    n2 = n2 / norm(n2)
    # @show n2
    @test norm(n1 -  n2) / norm(n2) <= 1.0e-6
    true
end
end
using .mbascs06
mbascs06.test()