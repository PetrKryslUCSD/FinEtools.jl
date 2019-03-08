module mesh_H8cylindern_1
using FinEtools
using FinEtools.MeshExportModule
using Test
function test()
  Radius::FFlt, Length::FFlt, nperradius, nL = 1.0, 2.0, 17, 3
  
  fens, fes = H8cylindern(Radius, Length, nperradius, nL)
  
  # # Postprocessing
  # vtkexportmesh("cylinder.vtk", fes.conn, fens.xyz,  FinEtools.MeshExportModule.H8)

  geom  =  NodalField(fens.xyz)

  femm  =  FEMMBase(IntegDomain(fes, GaussRule(3, 2)))
  V = integratefunction(femm, geom, (x) ->  1.0)

  @test abs(V - pi * Radius^2 * Length) / (pi * Radius^2 * Length) < 1.30e-3
end

end
using .mesh_H8cylindern_1
mesh_H8cylindern_1.test()
