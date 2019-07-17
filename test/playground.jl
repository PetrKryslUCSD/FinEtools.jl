module surfacenormaltest3
using FinEtools
using Test
function test()
   XYZ = reshape([0.0, 0.0], 2, 1)
   tangents = reshape([-1.0, 1.0], 2, 1)
   fe_label = 0
   fi = SurfaceNormal(2)
   @show v = updatenormal!(fi, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)
   @test v == [0.7071067811865475, 0.7071067811865475]
end
end
using .surfacenormaltest3
surfacenormaltest3.test()
