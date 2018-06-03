module Poisson_examples
using FinEtools

function Poisson_FE_example()
    println("""
    
    Heat conduction example described by Amuthan A. Ramabathiran
    http://www.codeproject.com/Articles/579983/Finite-Element-programming-in-Julia:
    Unit square, with known temperature distribution along the boundary,
    and uniform heat generation rate inside.  Mesh of regular linear TRIANGLES,
    in a grid of 1000 x 1000 edges (2M triangles, 1M degrees of freedom).
    Version: 05/29/2017
    """
    )
    t0 = time()
    
    A = 1.0 # dimension of the domain (length of the side of the square)
    thermal_conductivity =  [i==j ? one(FFlt) : zero(FFlt) for i=1:2, j=1:2]; # conductivity matrix
    Q = -6.0; # internal heat generation rate
    tempf(x) = (1.0 .+ x[:,1].^2 .+ 2.0 .* x[:,2].^2);#the exact distribution of temperature
    N = 1000;# number of subdivisions along the sides of the square domain
    
    
    println("Mesh generation")
    @time fens,fes =T3block(A, A, N, N)
    
    geom = NodalField(fens.xyz)
    Temp = NodalField(zeros(size(fens.xyz,1),1))
    
    println("Searching nodes  for BC")
    @time l1 = selectnode(fens; box=[0. 0. 0. A], inflate = 1.0/N/100.0)
    @time l2 = selectnode(fens; box=[A A 0. A], inflate = 1.0/N/100.0)
    @time l3 = selectnode(fens; box=[0. A 0. 0.], inflate = 1.0/N/100.0)
    @time l4 = selectnode(fens; box=[0. A A A], inflate = 1.0/N/100.0)
    List = vcat(l1, l2, l3, l4)
    @time setebc!(Temp, List, true, 1, tempf(geom.values[List,:])[:])
    @time applyebc!(Temp)
    @time numberdofs!(Temp)
    
    t1 = time()
    
    material = MatHeatDiff(thermal_conductivity)
    
    femm = FEMMHeatDiff(IntegData(fes, TriRule(1)), material)
    
    
    println("Conductivity")
    @time K = conductivity(femm, geom, Temp)
    println("Nonzero EBC")
    @time F2 = nzebcloadsconductivity(femm, geom, Temp);
    println("Internal heat generation")
    # function getsource!(forceout::FFltVec, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)
    #   forceout[1] = Q; #heat source
    # end
    # fi = ForceIntensity(FFlt, getsource!);# alternative  specification
    fi = ForceIntensity(FFlt[Q]);
    @time F1 = distribloads(FEMMBase(IntegData(fes, TriRule(1))), geom, Temp, fi, 3);
    
    println("Factorization")
    @time K = cholesky(K)
    println("Solution of the factorized system")
    @time U = K\(F1+F2)
    scattersysvec!(Temp,U[:])
    
    println("Total time elapsed = $(time() - t0) [s]")
    println("Solution time elapsed = $(time() - t1) [s]")
    
    Error= 0.0
    for k=1:size(fens.xyz,1)
        Error = Error + abs.(Temp.values[k,1] - tempf(reshape(fens.xyz[k,:], (1,2)))[1])
    end
    println("Error =$Error")
    
    
    # File =  "a.vtk"
    # MeshExportModule.vtkexportmesh (File, fes.conn, [geom.values Temp.values], MeshExportModule.T3; scalars=Temp.values, scalars_name ="Temperature")
    
    true
end

function Poisson_FE_example_algo()
    A= 1.0
    thermal_conductivity =  [i==j ? one(FFlt) : zero(FFlt) for i=1:2, j=1:2]; # conductivity matrix
    magn = -6.0; #heat source
    truetempf(x)=1.0 .+ x[1].^2 .+ 2.0 .* x[2].^2;
    N=20;
    
    println("""
    Heat conduction example described by Amuthan A. Ramabathiran
    http://www.codeproject.com/Articles/579983/Finite-Element-programming-in-Julia:
    Unit square, with known temperature distribution along the boundary,
    and uniform heat generation rate inside.  Mesh of regular TRIANGLES,
    in a grid of $N x $N edges.
    This version uses the FinEtools algorithm module.
    """
    )
    t0 = time()
    
    fens,fes = T3block(A, A, N, N)
    
    
    # Define boundary conditions
    l1  = selectnode(fens; box=[0. 0. 0. A], inflate = 1.0/N/100.0)
    l2  = selectnode(fens; box=[A A 0. A], inflate = 1.0/N/100.0)
    l3  = selectnode(fens; box=[0. A 0. 0.], inflate = 1.0/N/100.0)
    l4  = selectnode(fens; box=[0. A A A], inflate = 1.0/N/100.0)
    
    essential1 = FDataDict("node_list"=>vcat(l1, l2, l3, l4),
    "temperature"=>truetempf);
    material = MatHeatDiff(thermal_conductivity)
    femm = FEMMHeatDiff(IntegData(fes, TriRule(1)), material)
    region1 = FDataDict("femm"=>femm, "Q"=>magn)
    # Make model data
    modeldata= FDataDict("fens"=> fens,
    "regions"=>[region1],
    "essential_bcs"=>[essential1]);
    
    
    # Call the solver
    modeldata = FinEtools.AlgoHeatDiffModule.steadystate(modeldata)
    
    println("Total time elapsed = ",time() - t0,"s")
    
    geom=modeldata["geom"]
    Temp=modeldata["temp"]
    femm=modeldata["regions"][1]["femm"]
    
    function errfh(loc,val)
        exact = truetempf(loc)
        return ((exact[1]-val[1])*exact)
    end
    E = integratefieldfunction(femm, geom, Temp, errfh, 0.0, m=3)
    println("Error=$E")
    
end # Poisson_FE_example_algo


function Poisson_FE_example_csys_1()
    println("""
    Heat conduction example described by Amuthan A. Ramabathiran
    http://www.codeproject.com/Articles/579983/Finite-Element-programming-in-Julia:
    Unit square, with known temperature distribution along the boundary,
    and uniform heat generation rate inside.  Mesh of regular linear TRIANGLES,
    in a grid of 1000 x 1000 edges (2M triangles, 1M degrees of freedom).
    The material response is defined in a local coordinate system.
    Version: 05/29/2017
    """
    )
    t0 = time()
    
    A = 1.0 # dimension of the domain (length of the side of the square)
    thermal_conductivity =  [i==j ? one(FFlt) : zero(FFlt) for i=1:2, j=1:2]; # conductivity matrix
    Q = -6.0; # internal heat generation rate
    function getsource!(forceout::FFltVec, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)
        forceout[1] = Q; #heat source
    end
    tempf(x) = (1.0 .+ x[:,1].^2 .+ 2.0 .* x[:,2].^2);#the exact distribution of temperature
    N = 1000;# number of subdivisions along the sides of the square domain
    Rm=[-0.9917568452513019 -0.12813414805267656
    -0.12813414805267656 0.9917568452513019]
    Rm=[-0.8020689950104449 -0.5972313850116512
    -0.5972313850116512 0.8020689950104447]
    
    println("Mesh generation")
    @time fens,fes =T3block(A, A, N, N)
    
    geom = NodalField(fens.xyz)
    Temp = NodalField(zeros(size(fens.xyz,1),1))
    
    println("Searching nodes  for BC")
    @time l1 = selectnode(fens; box=[0. 0. 0. A], inflate = 1.0/N/100.0)
    @time l2 = selectnode(fens; box=[A A 0. A], inflate = 1.0/N/100.0)
    @time l3 = selectnode(fens; box=[0. A 0. 0.], inflate = 1.0/N/100.0)
    @time l4 = selectnode(fens; box=[0. A A A], inflate = 1.0/N/100.0)
    List = vcat(l1, l2, l3, l4)
    @time setebc!(Temp, List, true, 1, tempf(geom.values[List,:])[:])
    @time applyebc!(Temp)
    @time numberdofs!(Temp)
    
    t1 = time()
    
    material = MatHeatDiff(thermal_conductivity)
    
    femm = FEMMHeatDiff(IntegData(fes, TriRule(1)), CSys(Rm), material)
    
    println("Conductivity")
    @time K = conductivity(femm, geom, Temp)
    println("Nonzero EBC")
    @time F2 = nzebcloadsconductivity(femm, geom, Temp);
    println("Internal heat generation")
    fi = ForceIntensity(FFlt[Q]);
    @time F1 = distribloads(femm, geom, Temp, fi, 3);
    
    println("Factorization")
    @time K = cholesky(K)
    println("Solution of the factorized system")
    @time U = K\(F1+F2)
    scattersysvec!(Temp,U[:])
    
    println("Total time elapsed = $(time() - t0) [s]")
    println("Solution time elapsed = $(time() - t1) [s]")
    
    Error= 0.0
    for k=1:size(fens.xyz,1)
        Error = Error + abs.(Temp.values[k,1] - tempf(reshape(fens.xyz[k,:], (1,2)))[1])
    end
    println("Error =$Error")
    
    
    # File =  "a.vtk"
    # MeshExportModule.vtkexportmesh (File, fes.conn, [geom.values Temp.values], MeshExportModule.T3; scalars=Temp.values, scalars_name ="Temperature")
    
    true
    
end # Poisson_FE_example_csys_1


function Poisson_FE_Q4_example()
    println("""
    
    Heat conduction example described by Amuthan A. Ramabathiran
    http://www.codeproject.com/Articles/579983/Finite-Element-programming-in-Julia:
    Unit square, with known temperature distribution along the boundary,
    and uniform heat generation rate inside.  Mesh of regular four-node QUADRILATERALS,
    in a grid of 1000 x 1000 edges (1M quads, 1M degrees of freedom).
    Version: 05/29/2017
    """
    )
    t0 = time()
    
    A = 1.0
    thermal_conductivity =  [i==j ? one(FFlt) : zero(FFlt) for i=1:2, j=1:2]; # conductivity matrix
    function getsource!(forceout::FFltVec, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)
        forceout[1] = -6.0; #heat source
    end
    tempf(x) = (1.0 .+ x[:,1].^2 .+ 2.0 .* x[:,2].^2);#the exact distribution of temperature
    N = 1000;
    
    println("Mesh generation")
    @time fens,fes = Q4block(A, A, N, N)
    
    geom = NodalField(fens.xyz)
    Temp = NodalField(zeros(size(fens.xyz,1),1))
    
    
    println("Searching nodes  for BC")
    @time l1 = selectnode(fens; box=[0. 0. 0. A], inflate = 1.0/N/100.0)
    @time l2 = selectnode(fens; box=[A A 0. A], inflate = 1.0/N/100.0)
    @time l3 = selectnode(fens; box=[0. A 0. 0.], inflate = 1.0/N/100.0)
    @time l4 = selectnode(fens; box=[0. A A A], inflate = 1.0/N/100.0)
    List = vcat(l1, l2, l3, l4);
    setebc!(Temp, List, true, 1, tempf(geom.values[List,:])[:])
    applyebc!(Temp)
    
    numberdofs!(Temp)
    
    t1 = time()
    
    m = MatHeatDiff(thermal_conductivity)
    femm = FEMMHeatDiff(IntegData(fes, GaussRule(2, 2)), m)
    
    println("Conductivity")
    @time K=conductivity(femm, geom, Temp)
    #Profile.print()
    
    println("Nonzero EBC")
    @time F2 = nzebcloadsconductivity(femm, geom, Temp);
    println("Internal heat generation")
    fi = ForceIntensity(FFlt, 1, getsource!);
    @time F1 = distribloads(femm, geom, Temp, fi, 3);
    
    println("Factorization")
    @time K = cholesky(K)
    println("Solution of the factorized system")
    @time U=  K\(F1+F2)
    scattersysvec!(Temp, U[:])
    
    
    println("Total time elapsed = $(time() - t0) [s]")
    println("Solution time elapsed = $(time() - t1) [s]")
    
    # using MeshExportModule
    
    # File =  "a.vtk"
    # MeshExportModule.vtkexportmesh (File, fes.conn, [geom.values Temp.values], MeshExportModule.Q4; scalars=Temp.values, scalars_name ="Temperature")
    
    Error= 0.0
    for k=1:size(fens.xyz,1)
        Error = Error + abs.(Temp.values[k,1] - tempf(reshape(fens.xyz[k,:], (1,2)))[1])
    end
    println("Error =$Error")
    
    
    true
end # Poisson_FE_Q4_example

function allrun()
    println("#####################################################") 
    println("# Poisson_FE_example ")
    Poisson_FE_example()
    println("#####################################################") 
    println("# Poisson_FE_example_algo ")
    Poisson_FE_example_algo()
    println("#####################################################") 
    println("# Poisson_FE_example_csys_1 ")
    Poisson_FE_example_csys_1()
    println("#####################################################") 
    println("# Poisson_FE_Q4_example ")
    Poisson_FE_Q4_example()
end # function allrun

end # module Poisson_examples
