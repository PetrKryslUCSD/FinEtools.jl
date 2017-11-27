module Poisson_examples
using FinEtools

function Poisson_FE_H20_example()
    println("""
    
    Heat conduction example described by Amuthan A. Ramabathiran
    http://www.codeproject.com/Articles/579983/Finite-Element-programming-in-Julia:
    Unit cube, with known temperature distribution along the boundary,
    and uniform heat generation rate inside.  Mesh of regular quadratic HEXAHEDRA,
    in a grid of 30 x 30 x 30 edges (100079 degrees of freedom).
    Version: 06/03/2017
    """
    )
    t0 = time()
    
    A = 1.0 # dimension of the domain (length of the side of the square)
    thermal_conductivity = [i==j ? one(FFlt) : zero(FFlt) for i=1:3, j=1:3]; # conductivity matrix
    Q = -6.0; # internal heat generation rate
    function getsource!(forceout::FFltVec, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)
        forceout[1] = Q; #heat source
    end
    tempf(x) = (1.0 .+ x[:,1].^2 + 2.0 .* x[:,2].^2);#the exact distribution of temperature
    N = 30;# number of subdivisions along the sides of the square domain
    
    
    println("Mesh generation")
    @time fens,fes = H20block(A, A, A, N, N, N)
    
    geom = NodalField(fens.xyz)
    Temp = NodalField(zeros(size(fens.xyz,1),1))
    
    println("Searching nodes  for BC")
    Tolerance = 1.0/N/100.0
    @time l1 = selectnode(fens; box=[0. 0. 0. A 0. A], inflate = Tolerance)
    @time l2 = selectnode(fens; box=[A A 0. A 0. A], inflate = Tolerance)
    @time l3 = selectnode(fens; box=[0. A 0. 0. 0. A], inflate = Tolerance)
    @time l4 = selectnode(fens; box=[0. A A A 0. A], inflate = Tolerance)
    @time l5 = selectnode(fens; box=[0. A 0. A 0. 0.], inflate = Tolerance)
    @time l6 = selectnode(fens; box=[0. A 0. A A A], inflate = Tolerance)
    List = vcat(l1, l2, l3, l4, l5, l6)
    @time setebc!(Temp, List, true, 1, tempf(geom.values[List,:])[:])
    @time applyebc!(Temp)
    @time numberdofs!(Temp)
    
    println( "Number of free degrees of freedom: $(Temp.nfreedofs)")
    t1 = time()
    
    material = MatHeatDiff(thermal_conductivity)
    
    femm = FEMMHeatDiff(IntegData(fes, GaussRule(3, 3)), material)
    
    
    println("Conductivity")
    @time K = conductivity(femm, geom, Temp)
    println("Nonzero EBC")
    @time F2 = nzebcloadsconductivity(femm, geom, Temp);
    println("Internal heat generation")
    # fi = ForceIntensity(FFlt, getsource!);# alternative  specification
    fi = ForceIntensity(FFlt[Q]);
    @time F1 = distribloads(femm, geom, Temp, fi, 3);
    
    println("Solution of the system")
    @time U = K\(F1+F2)
    scattersysvec!(Temp,U[:])
    
    println("Total time elapsed = $(time() - t0) [s]")
    println("Solution time elapsed = $(time() - t1) [s]")
    
    Error= 0.0
    for k=1:size(fens.xyz,1)
        Error = Error + abs.(Temp.values[k,1] - tempf(reshape(fens.xyz[k,:], (1,3)))[1])
    end
    println("Error =$Error")
    
    
    # File =  "a.vtk"
    # MeshExportModule.vtkexportmesh (File, fes.conn, [geom.values Temp.values], MeshExportModule.T3; scalars=Temp.values, scalars_name ="Temperature")
    
    true
    
end # Poisson_FE_H20_example


function Poisson_FE_T10_example()
    println("""
    
    Heat conduction example described by Amuthan A. Ramabathiran
    http://www.codeproject.com/Articles/579983/Finite-Element-programming-in-Julia:
    Unit cube, with known temperature distribution along the boundary,
    and uniform heat generation rate inside.  Mesh of regular quadratic HEXAHEDRA,
    in a grid of 30 x 30 x 30 edges (100079 degrees of freedom).
    Version: 06/03/2017
    """
    )
    t0 = time()
    
    A = 1.0 # dimension of the domain (length of the side of the square)
    thermal_conductivity = [i==j ? one(FFlt) : zero(FFlt) for i=1:3, j=1:3]; # conductivity matrix
    Q = -6.0; # internal heat generation rate
    function getsource!(forceout::FFltVec, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)
        forceout[1] = Q; #heat source
    end
    tempf(x) = (1.0 .+ x[:,1].^2 + 2.0 .* x[:,2].^2);#the exact distribution of temperature
    N = 30;# number of subdivisions along the sides of the square domain
    
    
    println("Mesh generation")
    @time fens,fes = H20block(A, A, A, N, N, N)
    
    geom = NodalField(fens.xyz)
    Temp = NodalField(zeros(size(fens.xyz,1),1))
    
    println("Searching nodes  for BC")
    Tolerance = 1.0/N/100.0
    @time l1 = selectnode(fens; box=[0. 0. 0. A 0. A], inflate = Tolerance)
    @time l2 = selectnode(fens; box=[A A 0. A 0. A], inflate = Tolerance)
    @time l3 = selectnode(fens; box=[0. A 0. 0. 0. A], inflate = Tolerance)
    @time l4 = selectnode(fens; box=[0. A A A 0. A], inflate = Tolerance)
    @time l5 = selectnode(fens; box=[0. A 0. A 0. 0.], inflate = Tolerance)
    @time l6 = selectnode(fens; box=[0. A 0. A A A], inflate = Tolerance)
    List = vcat(l1, l2, l3, l4, l5, l6)
    @time setebc!(Temp, List, true, 1, tempf(geom.values[List,:])[:])
    @time applyebc!(Temp)
    @time numberdofs!(Temp)
    
    println( "Number of free degrees of freedom: $(Temp.nfreedofs)")
    t1 = time()
    
    material = MatHeatDiff(thermal_conductivity)
    
    femm = FEMMHeatDiff(IntegData(fes, GaussRule(3, 3)), material)
    
    
    println("Conductivity")
    @time K = conductivity(femm, geom, Temp)
    println("Nonzero EBC")
    @time F2 = nzebcloadsconductivity(femm, geom, Temp);
    println("Internal heat generation")
    # fi = ForceIntensity(FFlt, getsource!);# alternative  specification
    fi = ForceIntensity(FFlt[Q]);
    @time F1 = distribloads(femm, geom, Temp, fi, 3);
    
    println("Solution of the system")
    @time U = K\(F1+F2)
    scattersysvec!(Temp,U[:])
    
    println("Total time elapsed = $(time() - t0) [s]")
    println("Solution time elapsed = $(time() - t1) [s]")
    
    Error= 0.0
    for k=1:size(fens.xyz,1)
        Error = Error + abs.(Temp.values[k,1] - tempf(reshape(fens.xyz[k,:], (1,3)))[1])
    end
    println("Error =$Error")
    
    
    # File =  "a.vtk"
    # MeshExportModule.vtkexportmesh (File, fes.conn, [geom.values Temp.values], MeshExportModule.T3; scalars=Temp.values, scalars_name ="Temperature")
    
    true
    
end # Poisson_FE_T10_example


function Poisson_FE_T4_example()
    t0 = time()
    
    A = 1.0 # dimension of the domain (length of the side of the square)
    thermal_conductivity = [i==j ? one(FFlt) : zero(FFlt) for i=1:3, j=1:3]; # conductivity matrix
    Q = -6.0; # internal heat generation rate
    function getsource!(forceout::FFltVec, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)
        forceout[1] = Q; #heat source
    end
    tempf(x) = (1.0 .+ x[:,1].^2 + 2.0 .* x[:,2].^2);#the exact distribution of temperature
    N = 30;# number of subdivisions along the sides of the square domain
    
    println("Mesh generation")
    @time fens,fes = T4block(A, A, A, N, N, N)
    
    println("""
    Heat conduction example described by Amuthan A. Ramabathiran
    http://www.codeproject.com/Articles/579983/Finite-Element-programming-in-Julia:
    Unit cube, with known temperature distribution along the boundary,
    and uniform heat generation rate inside.  Mesh of regular quadratic TETRAHEDRA,
    in a grid of $(N) x $(N) x $(N) edges ($(count(fens)) degrees of freedom).
    Version: 07/03/2017
    """
    )
    
    geom = NodalField(fens.xyz)
    Temp = NodalField(zeros(size(fens.xyz,1),1))
    
    println("Searching nodes  for BC")
    Tolerance = 1.0/N/100.0
    @time l1 = selectnode(fens; box=[0. 0. 0. A 0. A], inflate = Tolerance)
    @time l2 = selectnode(fens; box=[A A 0. A 0. A], inflate = Tolerance)
    @time l3 = selectnode(fens; box=[0. A 0. 0. 0. A], inflate = Tolerance)
    @time l4 = selectnode(fens; box=[0. A A A 0. A], inflate = Tolerance)
    @time l5 = selectnode(fens; box=[0. A 0. A 0. 0.], inflate = Tolerance)
    @time l6 = selectnode(fens; box=[0. A 0. A A A], inflate = Tolerance)
    List = vcat(l1, l2, l3, l4, l5, l6)
    @time setebc!(Temp, List, true, 1, tempf(geom.values[List,:])[:])
    @time applyebc!(Temp)
    @time numberdofs!(Temp)
    
    println( "Number of free degrees of freedom: $(Temp.nfreedofs)")
    t1 = time()
    
    material = MatHeatDiff(thermal_conductivity)
    
    femm = FEMMHeatDiff(IntegData(fes, TetRule(1)), material)
    
    
    println("Conductivity")
    @time K = conductivity(femm, geom, Temp)
    println("Nonzero EBC")
    @time F2 = nzebcloadsconductivity(femm, geom, Temp);
    println("Internal heat generation")
    # fi = ForceIntensity(FFlt, getsource!);# alternative  specification
    fi = ForceIntensity(FFlt[Q]);
    @time F1 = distribloads(femm, geom, Temp, fi, 3);
    
    println("Solution of the system")
    @time U = K\(F1+F2)
    scattersysvec!(Temp,U[:])
    
    println("Total time elapsed = $(time() - t0) [s]")
    println("Solution time elapsed = $(time() - t1) [s]")
    
    Error= 0.0
    for k=1:size(fens.xyz,1)
        Error = Error + abs.(Temp.values[k,1] - tempf(reshape(fens.xyz[k,:], (1,3)))[1])
    end
    println("Error =$Error")
    
    
    # File =  "a.vtk"
    # MeshExportModule.vtkexportmesh (File, fes.conn, [geom.values Temp.values], MeshExportModule.T3; scalars=Temp.values, scalars_name ="Temperature")
    
    true
    
end # Poisson_FE_T4_example

function allrun()
    println("#####################################################") 
    println("# Poisson_FE_H20_example ")
    Poisson_FE_H20_example()
    println("#####################################################") 
    println("# Poisson_FE_T10_example ")
    Poisson_FE_T10_example()
    println("#####################################################") 
    println("# Poisson_FE_T4_example ")
    Poisson_FE_T4_example()
end # function allrun

end # module Poisson_examples
