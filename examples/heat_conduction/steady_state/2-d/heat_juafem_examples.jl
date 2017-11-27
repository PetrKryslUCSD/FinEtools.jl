module heat_juafem_examples
using FinEtools
using FinEtools.MeshExportModule
using Compat.Test

function heat_juafem_example()
    #println("""Heat conduction example from JuAFEM.""")
    # t0 = time()
    
    A = 2.0
    thermal_conductivity =  [i==j ? one(FFlt) : zero(FFlt) for i=1:2, j=1:2]; # conductivity matrix
    function getsource!(forceout::FFltVec, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)
        forceout[1] = 1.0; #heat source
    end
    N = 1000;
    
    #println("Mesh generation")
    fens,fes = Q4block(A, A, N, N)
    fens.xyz[:,1] .-= A/2
    fens.xyz[:,2] .-= A/2
    
    geom = NodalField(fens.xyz)
    Temp = NodalField(zeros(size(fens.xyz,1),1))
    
    #println("Searching nodes  for BC")
    l1 = selectnode(fens; box=[-A/2 -A/2 -A/2 A/2], inflate = 1.0/N/100.0)
    l2 = selectnode(fens; box=[A/2 A/2 -A/2 A/2], inflate = 1.0/N/100.0)
    l3 = selectnode(fens; box=[-A/2 A/2 -A/2 -A/2], inflate = 1.0/N/100.0)
    l4 = selectnode(fens; box=[-A/2 A/2 A/2 A/2], inflate = 1.0/N/100.0)
    List = vcat(l1, l2, l3, l4);
    setebc!(Temp, List, true, 1, 0.0)
    applyebc!(Temp)
    
    numberdofs!(Temp)
    
    # t1 = time()
    
    m = MatHeatDiff(thermal_conductivity)
    femm = FEMMHeatDiff(IntegData(fes, GaussRule(2, 2)), m)
    
    #println("Conductivity")
    K=conductivity(femm, geom, Temp)
    
    #println("Internal heat generation")
    fi = ForceIntensity(FFlt, 1, getsource!);
    F1 = distribloads(femm, geom, Temp, fi, 3);
    
    #println("Factorization")
    K = cholfact(K)
    #println("Solution of the factorized system")
    U =  K\(F1)
    scattersysvec!(Temp, U[:])
    
    
    # #println("Total time elapsed = $(time() - t0) [s]")
    # #println("Solution time elapsed = $(time() - t1) [s]")
    
    #println("Maximum temperature = $(maximum(Temp.values)) ")
    # using MeshExportModule
    
    # File =  "a.vtk"
    # MeshExportModule.vtkexportmesh(File, fes.conn, hcat(geom.values,  Temp.values), MeshExportModule.Q4; scalars=[("Temperature", Temp.values)])
    # @async run(`"paraview.exe" $File`)
end # heat_juafem_example

function allrun()
    #println("#####################################################") 
    println("# heat_juafem_example ")
    heat_juafem_example()
    return true
end # function allrun

end # module heat_juafem_examples
