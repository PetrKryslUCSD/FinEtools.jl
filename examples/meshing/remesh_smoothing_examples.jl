module remesh_smoothing_examples
using FinEtools
using FinEtools.MeshExportModule
using FinEtools.TetRemeshingModule

function remesh1()
    L= 0.3; 
    W = 0.3;
    a = 0.15;
    nL=46; nW=46; na=36;
    
    fens,fes = T4block(a,L,W,nL,nW,na,:a);
    t = deepcopy(connasarray(fes));
    v = deepcopy(fens.xyz);
    tmid = ones(Int, size(t,1));
    
    desired_ts =a;
    bfes = meshboundary(fes);
    f = connectednodes(bfes);
    bv = zeros(Bool, size(v,1));
    bv[f] .= true;
    
    println("Mesh size: initial = $(size(t,1))")
    t0 = time()
    
    t, v, tmid = TetRemeshingModule.coarsen(t, v, tmid; bv = bv, desired_ts = desired_ts);
    
    println("Mesh size: final = $(size(t,1)) [$(time() - t0) sec]")
    
    fens.xyz = deepcopy(v)
    fes = fromarray!(fes, t)
    setlabel!(fes, tmid)
    geom  =  NodalField(fens.xyz)
    
    femm  =  FEMMBase(IntegData(fes, SimplexRule(3, 1)))
    V = integratefunction(femm, geom, (x) ->  1.0)
    println("V = $(V) compared to $(L * W * a)")
    
    # File = "test1.vtk"
    # MeshExportModule.vtkexportmesh(File, t, v, MeshExportModule.T4)
    # @async run(`"paraview.exe" $File`)
end # remesh1

function mesh_smoothing()
    println("""
    Meshing, deforming  and smoothing
    """)
    
    A = 100. # strip width
    N = 16
    tolerance = A / N / 1.0e5
    fens, fes = T3block(A, A, N, N)
    
    bnl = connectednodes(meshboundary(fes))
    for ixxxx = 1:length(bnl)
        x, y = fens.xyz[bnl[ixxxx], :]
        fens.xyz[bnl[ixxxx], 1] += A / N * sin(2 * pi * y / A)
        fens.xyz[bnl[ixxxx], 2] += -A / N * sin(2 * pi * x / A)
    end
    
    File =  "mesh_smoothing_before.vtk"
    vtkexportmesh(File, fens, fes);
    @async run(`"paraview.exe" $File`)
    println("$(fens.xyz[Int(N^2 / 2), :] )")
    
    fixedv = falses(count(fens))
    fixedv[bnl] .= true
    fens = meshsmoothing(fens, fes; fixedv=fixedv, method=:taubin, npass=100)
    
    println("$(fens.xyz[Int(N^2 / 2), :] )")
    
    geom = NodalField(fens.xyz)
    
    File =  "mesh_smoothing_after.vtk"
    vtkexportmesh(File, fes.conn, geom.values, FinEtools.MeshExportModule.T3);
    @async run(`"paraview.exe" $File`)
    
    println("Done")
    true
    
end # mesh_smoothing

function allrun()
    println("#####################################################") 
    println("# remesh1 ")
    remesh1()
    println("#####################################################") 
    println("# mesh_smoothing ")
    mesh_smoothing()
    return true
end # function allrun

end # module remesh_smoothing_examples
