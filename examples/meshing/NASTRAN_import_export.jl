module junk

using FinEtools
using FinEtools.MeshImportModule: import_NASTRAN
using FinEtools.MeshExportModule

output = import_NASTRAN("$(@__DIR__)" * "/Slot-coarser.nas")
File = "Slot-coarser.vtk"
MeshExportModule.vtkexportmesh(File, output["fens"], output["fesets"][1])
@async run(`"paraview.exe" $File`)
end
