module mexpoff1
using FinEtools
using FinEtools.MeshImportModule
using FinEtools.MeshExportModule.OFF: header, vertex, facet
using Test
function test()
    
    output = MeshImportModule.import_ABAQUS("failed-tetgen.inp"; allocationchunk = 11)
    fens, fes = output["fens"], output["fesets"][1]

    filename = "mesh.off"
    e = OFFExporter(filename::AbstractString)
    header(e, count(fens), count(fes))
    for i in eachindex(fens)
        vertex(e, fens.xyz[i, :]...)
    end 
    for i in eachindex(fes)
        c = fes.conn[i]
        facet(e, c[1], c[2], c[3],)
    end
    close(e)
    @test filesize(filename) > 0
    rm(filename)
    # @async run(`"paraview.exe" $File`)
end
end
using .mexpoff1
mexpoff1.test()
