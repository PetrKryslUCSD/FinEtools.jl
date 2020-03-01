
module mexpvecv2
using FinEtools
using FinEtools.MeshExportModule: VTK, VTK.vtkexportvectors
using Test
function test()
    rin::FFlt, rex::FFlt, nr::FInt, nc::FInt, Angl::FFlt, orientation::Symbol = 100.0, 200.0, 3, 6, pi/3, :a
    fens, fes = T3annulus(rin::FFlt, rex::FFlt, nr::FInt, nc::FInt, Angl::FFlt, orientation::
      Symbol)
    fens.xyz = xyz3(fens)
    d = (fens.xyz[:, 1].^2 + fens.xyz[:, 2].^2)
    File = "mesh.vtk"
    # Export of multiple scalar fields
    result =  VTK.vtkexportmesh(File, fens, fes; scalars = [("d", d), ("invd", 1 ./ d)])
    @test result == true
    # rm(File)

    points = [vec(fens.xyz[idx, :]) for idx in 1:count(fens)] 
    vectors = [("grad", [vec(fens.xyz[idx, [2,1,3]]) for idx in 1:count(fens)])]
    filename = "vectors.vtk"
    result =  VTK.vtkexportvectors(filename, points, vectors)
    @test result == true
    # rm(filename)
end
end
using .mexpvecv2
mexpvecv2.test()
