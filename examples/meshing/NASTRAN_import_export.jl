module junk

using FinEtools
using FinEtools.MeshImportModule
using FinEtools.MeshExportModule

fens, fes = MeshImportModule.import_NASTRAN("C:/Users/Petr Krysl/Dropbox/Julia/FinEtools_Examples/meshing/Slot-coarser.nas")
File = "Slot-coarser.vtk"
MeshExportModule.vtkexportmesh(File, fens, fes)
@async run(`"paraview.exe" $File`)
end
