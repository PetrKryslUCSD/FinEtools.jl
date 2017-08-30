module junk

using FinEtools
using FinEtools.MeshImportModule: import_NASTRAN
using FinEtools.MeshExportModule

fens, fes = import_NASTRAN("$(@__DIR__)" * "/Slot-coarser.nas")
File = "Slot-coarser.vtk"
MeshExportModule.vtkexportmesh(File, fens, fes)
@async run(`"paraview.exe" $File`)
end
