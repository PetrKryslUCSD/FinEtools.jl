module NASTRAN_import_export_examples
using FinEtools
using FinEtools.MeshImportModule: import_NASTRAN
using FinEtools.MeshExportModule

function NASTRAN_import_export()
    output = import_NASTRAN("$(@__DIR__)" * "/Slot-coarser.nas")
    File = "Slot-coarser.vtk"
    MeshExportModule.vtkexportmesh(File, output["fens"], output["fesets"][1])
    @async run(`"paraview.exe" $File`)
end # NASTRAN_import_export

function allrun()
    println("#####################################################") 
    println("# NASTRAN_import_export ")
    NASTRAN_import_export()
    return true
end # function allrun

end # module NASTRAN_import_export_examples
