using Test

module mt4refine1
using FinEtools
using FinEtools.MeshExportModule
using Test
function test()
    xs = collect(linearspace(0.0, pi / 2, 5))
    ys = collect(linearspace(0.0, 1.0, 6).^2)
    zs = collect(linearspace(0.0, 1.0, 3))
    fens, fes = T4blockx(xs, ys, zs, :a)
    # for i = 1:count(fens)
    #     a, y, z = fens.xyz[i,:]
    #     fens.xyz[i,1] = sin(a) * (y + 0.5)
    #     fens.xyz[i,2] = cos(a) * (y + 0.5)
    #     fens.xyz[i,3] = z
    # end
    @test count(fes) == 240
    bfes = meshboundary(fes)
    @test count(bfes) == 2*2*(4*5 + 5*2 + 4*2)
    fens, fes = T4refine(fens, fes)
    for i = 1:count(fens)
        a, y, z = fens.xyz[i,:]
        fens.xyz[i,1] = sin(a) * (y + 0.5)
        fens.xyz[i,2] = cos(a) * (y + 0.5)
        fens.xyz[i,3] = z
    end
    @test count(fes) == 240*8
    bfes = meshboundary(fes)
    @test count(bfes) == 4*2*2*(4*5 + 5*2 + 4*2)

    geom  =  NodalField(fens.xyz)
    femm  =  FEMMBase(IntegData(fes, SimplexRule(3, 4)))
    V = integratefunction(femm, geom, (x) ->  1.0)
    
    File = "Refine-T4-a.vtk"
    MeshExportModule.vtkexportmesh(File, fens, bfes)
    rm(File)
    # @async run(`"paraview.exe" $File`)
end
end
using .mt4refine1
mt4refine1.test()


module mt4refine2
using FinEtools
using FinEtools.MeshExportModule
using Test
function test()
    xs = collect(linearspace(0.0, pi / 2, 5))
    ys = collect(linearspace(0.0, 1.0, 6).^2)
    zs = collect(linearspace(0.0, 1.0, 3))
    fens, fes = T4blockx(xs, ys, zs, :a)
    fens, fes = T4toT10(fens, fes)
    # for i = 1:count(fens)
    #     a, y, z = fens.xyz[i,:]
    #     fens.xyz[i,1] = sin(a) * (y + 0.5)
    #     fens.xyz[i,2] = cos(a) * (y + 0.5)
    #     fens.xyz[i,3] = z
    # end
    @test count(fes) == 240
    bfes = meshboundary(fes)
    @test count(bfes) == 2*2*(4*5 + 5*2 + 4*2)
    fens, fes = T10refine(fens, fes)
    for i = 1:count(fens)
        a, y, z = fens.xyz[i,:]
        fens.xyz[i,1] = sin(a) * (y + 0.5)
        fens.xyz[i,2] = cos(a) * (y + 0.5)
        fens.xyz[i,3] = z
    end
    @test count(fes) == 240*8
    bfes = meshboundary(fes)
    @test count(bfes) == 4*2*2*(4*5 + 5*2 + 4*2)

    geom  =  NodalField(fens.xyz)
    femm  =  FEMMBase(IntegData(fes, SimplexRule(3, 4)))
    V = integratefunction(femm, geom, (x) ->  1.0)
    
    File = "Refine-T10-a.vtk"
    MeshExportModule.vtkexportmesh(File, fens, bfes)
    rm(File)
    # @async run(`"paraview.exe" $File`)
end
end
using .mt4refine2
mt4refine2.test()