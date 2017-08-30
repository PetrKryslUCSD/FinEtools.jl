module junk

using FinEtools
using FinEtools.MeshExportModule
using FinEtools.MeshImportModule: import_ABAQUS



fens, fesarray = import_ABAQUS("$(@__DIR__)/" * "NLE10-coarse-T10.inp")

File = "LE10NAFEMS_T10.vtk"
MeshExportModule.vtkexportmesh(File, fens, fesarray[1])
@async run(`"paraview.exe" $File`)



end
