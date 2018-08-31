module T129b_l2_examples
using FinEtools
using FinEtools.AlgoHeatDiffModule
using LinearAlgebra: cholesky
using PGFPlotsX

function T129b_l2_uq()
    L = 6.0;
    kappa = reshape([4.0], 1, 1);
    Q  = 0.1;
    n = 7; # number of elements
    crosssection = 1.0
    
    fens,fes = L2block(L, n)
    
    # Define boundary conditions
    l1  = selectnode(fens; box=[0. 0.], inflate = L/n/100.0)
    l2  = selectnode(fens; box=[L L], inflate = L/n/100.0)
    
    essential1 = FDataDict("node_list"=>vcat(l1, l2), "temperature"=> 0.0);
    material = MatHeatDiff(kappa)
    femm = FEMMHeatDiff(IntegData(fes, GaussRule(1, 2), crosssection), material)
    region1 = FDataDict("femm"=>femm, "Q"=>Q)
    # Make model data
    modeldata= FDataDict("fens"=> fens,    "regions"=>[region1],    "essential_bcs"=>[essential1]);
    
    geom = NodalField(fens.xyz)
    Temp = NodalField(zeros(size(fens.xyz,1),1))
    setebc!(Temp, vcat(l1, l2), true, 1, 0.0)
    applyebc!(Temp)
    numberdofs!(Temp)
    
    K = conductivity(femm, geom, Temp)
    
    F2 = nzebcloadsconductivity(femm, geom, Temp);
    
    fi = ForceIntensity(FFlt[Q]);
    F1 = distribloads(femm, geom, Temp, fi, 3);
    
    K = cholesky(K)
    U = K\(F1+F2)
    scattersysvec!(Temp,U[:])
    
    println("maximum(U)-0.1102 = $(maximum(U)-0.1102)")
    
    @pgf a = Axis({
        xlabel = "Location",
        ylabel = "T",
        title = "Temperature plot"
    },
    Plot(Table([:x => vec(geom.values), :y => vec(Temp.values)])))
    display(a)

    # function errfh(loc,val)
    #     x = loc[1]
    #     exact = (kappa[1,1]/Q/L^2)*Q/(2*kappa[1,1])*x*(L-x)
    #     return ((exact-val)*exact)[1]
    # end
    #
    # femm.integdata.integration_rule = GaussRule(1, 4)
    # E = integratefieldfunction(femm, geom, Temp, errfh, 0.0, m=3)
    # println("Error=$E")
    
    # Postprocessing
    # geom=modeldata["geom"]
    # Temp=modeldata["temp"]
    # MeshExportModule.vtkexportmesh ("a.vtk", fes.conn, [geom.values Temp.values], MeshExportModule.T3; scalars=Temp.values, scalars_name ="Temperature")
    
end # T129b_l2_uq


function T129b_l2_uq_algo()
    L = 6.0;
    kappa = reshape([4.0], 1, 1);
    Q  = 0.1;
    n = 7; # number of elements
    crosssection = 1.0
    
    fens,fes = L2block(L, n)
    
    # Define boundary conditions
    l1  = selectnode(fens; box=[0. 0.], inflate = L/n/100.0)
    l2  = selectnode(fens; box=[L L], inflate = L/n/100.0)
    
    essential1 = FDataDict("node_list"=>vcat(l1, l2), "temperature"=> 0.0);
    material = MatHeatDiff(kappa)
    femm = FEMMHeatDiff(IntegData(fes, GaussRule(1, 2), crosssection), material)
    region1 = FDataDict("femm"=>femm, "Q"=>Q)
    # Make model data
    modeldata= FDataDict("fens"=> fens,    "regions"=>[region1],    "essential_bcs"=>[essential1]);
    
    
    # Call the solver
    modeldata = AlgoHeatDiffModule.steadystate(modeldata)
    
    geom=modeldata["geom"]
    Temp=modeldata["temp"]
    femm=modeldata["regions"][1]["femm"]
    
    # Evaluate error
    function errfh(loc,val)
        x = loc[1]
        exact = (kappa[1,1]/Q/L^2)*Q/(2*kappa[1,1])*x*(L-x)
        return ((exact-val[1])*exact)[1]
    end
    
    femm.integdata.integration_rule = GaussRule(1, 4)
    E = integratefieldfunction(femm, geom, Temp, errfh, 0.0, m=3)
    println("Error=$E")
    
    @pgf a = Axis({
        xlabel = "Location",
        ylabel = "T",
        title = "Temperature plot"
    },
    Plot(Table([:x => vec(geom.values), :y => vec(Temp.values)])))
    display(a)
    
end # T129b_l2_uq_algo

function allrun()
    println("#####################################################") 
    println("# T129b_l2_uq ")
    T129b_l2_uq()
    println("#####################################################") 
    println("# T129b_l2_uq_algo ")
    T129b_l2_uq_algo()
end # function allrun

end # module T129b_l2_examples
