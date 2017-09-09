module junk

using FinEtools
using FinEtools.MeshExportModule
using FinEtools.MeshImportModule: import_ABAQUS



output = import_ABAQUS("$(@__DIR__)/" * "NLE10-coarse-T10.inp")
fens, fes = output["fens"], output["fesets"][1]

File = "LE10NAFEMS_T10.vtk"
MeshExportModule.vtkexportmesh(File, fens, fes)
@async run(`"paraview.exe" $File`)



end
