
module mimportexportm1
using FinEtools
using FinEtools.MeshImportModule
using FinEtools.MeshExportModule
using Test
function test()
  output = MeshImportModule.import_NASTRAN(dirname(@__FILE__) * "/" * "cylinder.nas";
    allocationchunk = 13, expectfixedformat = true)
  # show(fes.conn[count(fes), :])
  File = "cylinder.vtk"
  MeshExportModule.vtkexportmesh(File, output["fens"], output["fesets"][1])
#   rm(File)
#   @test output["fesets"][1].conn[count(output["fesets"][1]), :] == NTuple{10,Int64}[(143, 140, 144, 138, 361, 363, 176, 519, 781, 520)]
  # @async run(`"paraview.exe" $File`)
end
end
using .mimportexportm1
 mimportexportm1.test()