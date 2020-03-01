module mbas001
using FinEtools
using Test
function test()
   fens,fes = Q4block(1.3, 3.1,3,2); # Mesh
   @test manifdim(fes) == 2
   bfes = meshboundary(fes)
   @test manifdim(bfes) == 1
   true
end
end
using .mbas001
mbas001.test()

module mbas002
using FinEtools
using Test
function test()
   fens,fes = H8block(1.3, 3.1, 2.7, 3, 2, 4); # Mesh
   @test manifdim(fes) == 3
   bfes = meshboundary(fes)
   @test manifdim(bfes) == 2
   true
end
end
using .mbas002
mbas002.test()

module mbas003
using FinEtools
using Test
function test()
   fens,fes = T4block(1.3, 3.1, 2.7, 3, 2, 4); # Mesh
   @test manifdim(fes) == 3
   bfes = meshboundary(fes)
   @test manifdim(bfes) == 2
   true
end
end
using .mbas003
mbas003.test()


module mbas004
using FinEtools
using Test
function test()
   fens,fes = Q4block(1.3, 3.1,3,2); # Mesh
   @test nodesperelem(fes) == 4
   bfes = meshboundary(fes)
   @test nodesperelem(bfes) == 2
   true
end
end
using .mbas004
mbas004.test()

module mbas005
using FinEtools
using Test
function test()
   fens,fes = H8block(1.3, 3.1, 2.7, 3, 2, 4); # Mesh
   @test nodesperelem(fes) == 8
      bfes = meshboundary(fes)
      @test nodesperelem(bfes) == 4
   true
end
end
using .mbas005
mbas005.test()

module mbas006
using FinEtools
using Test
function test()
   fens,fes = T4block(1.3, 3.1, 2.7, 3, 2, 4); # Mesh
   @test nodesperelem(fes) == 4
      bfes = meshboundary(fes)
      @test nodesperelem(bfes) == 3
   true
end
end
using .mbas006
mbas006.test()


module mbas007
using FinEtools
using Test
function test()
   fens,fes = T10block(1.3, 3.1, 2.7, 3, 2, 4); # Mesh
   @test nodesperelem(fes) == 10
      bfes = meshboundary(fes)
      @test nodesperelem(bfes) == 6
   true
end
end
using .mbas007
mbas007.test()


module mbas008
using FinEtools
using Test
function test()
 fens,fes = L2block(1.37, 3); # Mesh
 @test nodesperelem(fes) == 2
 bfes = meshboundary(fes)
 @test nodesperelem(bfes) == 1
 true
end
end
using .mbas008
mbas008.test()


module mbas009
using FinEtools
using Test
function test()
   fens,fes = Q8block(1.3, 3.1,3,2); # Mesh
   @test nodesperelem(fes) == 8
   bfes = meshboundary(fes)
   @test nodesperelem(bfes) == 3
   true
end
end
using .mbas009
mbas009.test()

module mbas010
using FinEtools
using FinEtools.AlgoBaseModule: fieldnorm
using Test
function test()
   fens,fes = Q8block(1.3, 3.1,13,12); # Mesh
   geom = NodalField(fens.xyz)
   femm = FEMMBase(IntegDomain(fes, GaussRule(2, 3)))
   regions = [FDataDict("femm"=>femm)]
   # x = NodalField(reshape(fens.xyz[:, 1], count(fens), 1))
   # targetfields = [x]
   c = NodalField(fill(1.0, count(fens), 1))
   targetfields = [c]
   modeldata = FDataDict("fens"=>fens, "regions"=>regions, "targetfields"=>targetfields, "geom"=>geom, "elementsize"=>3.1 / 12)
   @test abs(fieldnorm(modeldata) - sqrt(1.3 * 3.1)) <= 1.0e-5
   true
end
end
using .mbas010
mbas010.test()

module mbascs01
using FinEtools
using FinEtools.AlgoBaseModule: fieldnorm
using LinearAlgebra
using Test
function test()
    center = [0.0, 0.0, 0.0]
    ez = [0.0, 0.0, 1.0]
    v = rand(3)
    XYZ = reshape([0.2, 0.3, 0.4], 1, 3)
    tangents = rand(3, 3)
    rfcsmat = [0.5547001962252291 -0.8320502943378437 0.0; 0.8320502943378436 0.5547001962252293 0.0; 0.0 0.0 1.0]
    function compute!(csmatout::FFltMat, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)
        # Cylindrical coordinate system
        xyz = XYZ[:] .- center
        xyz[3] = 0.0
        csmatout[:, 1] = xyz / norm(xyz)
        cross3!(v, ez, vec(csmatout[:, 1]))
        csmatout[:, 2] = v / norm(v)
        csmatout[:, 3] = ez
        return csmatout
    end
    csys = CSys(3, 3, compute!)
    updatecsmat!(csys, XYZ, tangents, 0)
    @test norm(csys.csmat -  rfcsmat) / norm(rfcsmat) <= 1.0e-6
    true
end
end
using .mbascs01
mbascs01.test()

module mbascs02
using FinEtools
using FinEtools.AlgoBaseModule: fieldnorm
using LinearAlgebra
using Test
function test()
    center = [0.0, 0.0, 0.0]
    ez = [0.0, 0.0, 1.0]
    v = rand(3)
    XYZ = reshape([0.2, 0.3, 0.4], 1, 3)
    tangents = rand(3, 3)
    rfcsmat = [0.5547001962252291 -0.8320502943378437 0.0; 0.8320502943378436 0.5547001962252293 0.0; 0.0 0.0 1.0]
    # function compute!(csmatout::FFltMat, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)
    #     # Cylindrical coordinate system
    #     xyz = XYZ[:] .- center
    #     xyz[3] = 0.0
    #     csmatout[:, 1] = xyz / norm(xyz)
    #     cross3!(v, ez, vec(csmatout[:, 1]))
    #     csmatout[:, 2] = v / norm(v)
    #     csmatout[:, 3] = ez
    #     return csmatout
    # end
    csys = CSys(rfcsmat)
    updatecsmat!(csys, XYZ, tangents, 0)
    @test norm(csys.csmat -  rfcsmat) / norm(rfcsmat) <= 1.0e-6
    true
end
end
using .mbascs02
mbascs02.test()

module mbascs03
using FinEtools
using FinEtools.AlgoBaseModule: fieldnorm
using LinearAlgebra
using Test
function test()
    center = [0.0, 0.0, 0.0]
    ez = [0.0, 0.0, 1.0]
    v = rand(3)
    XYZ = reshape([0.2, 0.3, 0.4], 1, 3)
    tangents = rand(3, 3)
    rfcsmat = Matrix(1.0 * I, 3, 3)
    # function compute!(csmatout::FFltMat, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)
    #     # Cylindrical coordinate system
    #     xyz = XYZ[:] .- center
    #     xyz[3] = 0.0
    #     csmatout[:, 1] = xyz / norm(xyz)
    #     cross3!(v, ez, vec(csmatout[:, 1]))
    #     csmatout[:, 2] = v / norm(v)
    #     csmatout[:, 3] = ez
    #     return csmatout
    # end
    csys = CSys(3)
    updatecsmat!(csys, XYZ, tangents, 0)
    @test norm(csys.csmat -  rfcsmat) / norm(rfcsmat) <= 1.0e-6
    true
end
end
using .mbascs03
mbascs03.test()

module mbascs04
using FinEtools
using FinEtools.AlgoBaseModule: fieldnorm
using LinearAlgebra
using Test
function test()
    XYZ = reshape([0.2, 0.3, 0.4], 1, 3)
    rfcsmat = Matrix(1.0 * I, 3, 3)
    tangents = rand(3, 3)
    csys = CSys(3, 3)
    updatecsmat!(csys, XYZ, tangents, 0)
    @test norm(csys.csmat -  rfcsmat) / norm(rfcsmat) <= 1.0e-6
    true
end
end
using .mbascs04
mbascs04.test()

module mbascs05
using FinEtools
using FinEtools.AlgoBaseModule: fieldnorm
using LinearAlgebra
using Test
function test()
    XYZ = reshape([0.2, 0.3, 0.4], 1, 3)
    rfcsmat = reshape([0.37139067635410367; 0.5570860145311555; 0.7427813527082073], 3, 1)
    tangents = reshape([0.2, 0.3, 0.4], 3, 1)
    csys = CSys(3, 1)
    updatecsmat!(csys, XYZ, tangents, 0)
    # @show csys.csmat
    @test norm(csys.csmat -  rfcsmat) / norm(rfcsmat) <= 1.0e-6
    true
end
end
using .mbascs05
mbascs05.test()

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

