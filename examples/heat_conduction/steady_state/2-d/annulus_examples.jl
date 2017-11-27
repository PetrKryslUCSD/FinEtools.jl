module annulus_examples
using FinEtools
using FinEtools.AlgoHeatDiffModule

function annulus_Q4_example()
    println("""
    Annular region,  ingoing and outgoing flux. Minimum/maximum temperature ~(+/-)0.591106.
    Mesh of linear quadrilaterals.
    """)
    
    t0 = time()
    
    kappa = 0.2*[1.0 0; 0 1.0]; # conductivity matrix
    magn = 0.06;# heat flux along the boundary
    rin =  1.0;#internal radius
    rex =  2.0;#external radius
    nr = 10; nc = 80;
    Angle = 2*pi;
    thickness =  1.0;
    tolerance = min(rin/nr,  rin/nc/2/pi)/10000;
    
    
    fens, fes  =  Q4annulus(rin, rex, nr, nc, Angle)
    fens, fes  =  mergenodes(fens,  fes,  tolerance);
    edge_fes  =  meshboundary(fes);
    
    geom = NodalField(fens.xyz)
    Temp = NodalField(zeros(size(fens.xyz, 1), 1))
    
    
    l1  = selectnode(fens; box=[0.0 0.0 -rex -rex],  inflate = tolerance)
    setebc!(Temp, l1, 1; val=zero(FFlt))
    applyebc!(Temp)
    
    numberdofs!(Temp)
    
    
    material = MatHeatDiff(kappa)
    femm = FEMMHeatDiff(IntegData(fes,  GaussRule(2, 2)),  material)
    
    @time K = conductivity(femm,  geom,  Temp)
    
    l1 = selectelem(fens, edge_fes, box=[-1.1*rex -0.9*rex -0.5*rex 0.5*rex]);
    el1femm = FEMMBase(IntegData(subset(edge_fes, l1),  GaussRule(1, 2)))
    fi = ForceIntensity(FFlt[-magn]);#entering the domain
    @time F1 = (-1.0)* distribloads(el1femm,  geom,  Temp,  fi,  2);
    
    l1 = selectelem(fens, edge_fes, box=[0.9*rex 1.1*rex -0.5*rex 0.5*rex]);
    el1femm =  FEMMBase(IntegData(subset(edge_fes, l1),  GaussRule(1, 2)))
    fi = ForceIntensity(FFlt[+magn]);#leaving the domain
    @time F2 = (-1.0)* distribloads(el1femm,  geom,  Temp,  fi,  2);
    
    @time F3 = nzebcloadsconductivity(femm,  geom,  Temp);
    
    
    @time K = cholfact(K)
    @time U = K\(F1+F2+F3)
    @time scattersysvec!(Temp, U[:])
    
    println("Total time elapsed = ", time() - t0, "s")
    
    File =  "annulus.vtk"
    vtkexportmesh(File,  connasarray(fes),  [geom.values Temp.values], FinEtools.MeshExportModule.Q4; scalars=[("Temperature", Temp.values)])
    @async run(`"paraview.exe" $File`)
    
    println("Minimum/maximum temperature= $(minimum(Temp.values))/$(maximum(Temp.values)))")
    
    true
    
end # annulus_Q4_example


function annulus_Q4_example_algo()
    println("""
    Annular region, ingoing and outgoing flux. Temperature at one node prescribed.
    Minimum/maximum temperature ~(+/-)0.500.
    Mesh of serendipity quadrilaterals.
    This version uses the FinEtools algorithm module.
    Version: 05/29/2017
    """)
    
    t0 = time()
    
    kappa = 0.2*[1.0 0; 0 1.0]; # conductivity matrix
    magn = 0.06;# heat flux along the boundary
    rin =  1.0;#internal radius
    
    rex =  2.0; #external radius
    nr = 3; nc = 40;
    Angle = 2*pi;
    thickness =  1.0;
    tolerance = min(rin/nr,  rin/nc/2/pi)/10000;
    
    fens, fes = Q4annulus(rin, rex, nr, nc, Angle)
    fens, fes = mergenodes(fens,  fes,  tolerance);
    edge_fes = meshboundary(fes);
    
    # At a single point apply an essential boundary condition (pin down the temperature)
    l1  = selectnode(fens; box=[0.0 0.0 -rex -rex],  inflate = tolerance)
    essential1 = FDataDict("node_list"=>l1, "temperature"=>0.0)
    
    # The flux boundary condition is applied at two pieces of surface
    # Side 1
    l1 = selectelem(fens, edge_fes, box=[-1.1*rex -0.9*rex -0.5*rex 0.5*rex]);
    el1femm = FEMMBase(IntegData(subset(edge_fes, l1),  GaussRule(1, 2)))
    fi = ForceIntensity(FFlt[-magn]);#entering the domain
    flux1 = FDataDict("femm"=>el1femm, "normal_flux"=>-magn) # entering the domain
    # Side 2
    l2=selectelem(fens,edge_fes,box=[0.9*rex 1.1*rex -0.5*rex 0.5*rex]);
    el2femm = FEMMBase(IntegData(subset(edge_fes, l2),  GaussRule(1, 2)))
    flux2 = FDataDict("femm"=>el2femm, "normal_flux"=>+magn) # leaving the domain
    
    material = MatHeatDiff(kappa)
    femm = FEMMHeatDiff(IntegData(fes,  GaussRule(2, 2)),  material)
    region1 = FDataDict("femm"=>femm)
    
    # Make model data
    modeldata = FDataDict("fens"=>fens,
    "regions"=>[region1], "essential_bcs"=>[essential1],
    "flux_bcs"=>[flux1, flux2]);
    
    # Call the solver
    modeldata = FinEtools.AlgoHeatDiffModule.steadystate(modeldata)
    geom=modeldata["geom"]
    Temp=modeldata["temp"]
    println("Minimum/maximum temperature= $(minimum(Temp.values))/$(maximum(Temp.values)))")
    
    println("Total time elapsed = ",time() - t0,"s")
    
    # Postprocessing
    vtkexportmesh("annulusmod.vtk", connasarray(fes), [geom.values Temp.values], FinEtools.MeshExportModule.Q4; scalars=[("Temperature", Temp.values)])
    
end # annulus_Q4_example_algo


function annulus_Q8_example()
    println("""
    Annular region,  ingoing and outgoing flux. Minimum/maximum temperature ~(+/-)0.501.
    Mesh of quadratic serendipity quadrilaterals.
    Version: 05/29/2017
    """)
    
    t0 = time()
    
    kappa = 0.2*[1.0 0; 0 1.0]; # conductivity matrix
    magn = 0.06;# heat flux along the boundary
    rin =  1.0;#internal radius
    
    rex =  2.0; #external radius
    nr = 5; nc = 40;
    Angle = 2*pi;
    thickness =  1.0;
    tolerance = min(rin/nr,  rin/nc/2/pi)/10000;
    
    
    fens, fes  =  Q8annulus(rin, rex, nr, nc, Angle)
    fens, fes  =  mergenodes(fens,  fes,  tolerance);
    edge_fes  =  meshboundary(fes);
    
    geom = NodalField(fens.xyz)
    Temp = NodalField(zeros(size(fens.xyz, 1), 1))
    
    
    l1  = selectnode(fens; box=[0.0 0.0 -rex -rex],  inflate = tolerance)
    setebc!(Temp, l1, 1; val=zero(FFlt))
    applyebc!(Temp)
    
    numberdofs!(Temp)
    
    
    material = MatHeatDiff(kappa)
    femm = FEMMHeatDiff(IntegData(fes,  GaussRule(2, 3)),  material)
    
    @time K = conductivity(femm,  geom,  Temp)
    
    l1 = selectelem(fens, edge_fes, box=[-1.1*rex -0.9*rex -0.5*rex 0.5*rex]);
    el1femm = FEMMBase(IntegData(subset(edge_fes, l1),  GaussRule(1, 2)))
    fi = ForceIntensity(FFlt[-magn]);#entering the domain
    @time F1 = (-1.0)* distribloads(el1femm,  geom,  Temp,  fi,  2);
    
    l1 = selectelem(fens, edge_fes, box=[0.9*rex 1.1*rex -0.5*rex 0.5*rex]);
    el1femm =  FEMMBase(IntegData(subset(edge_fes, l1),  GaussRule(1, 2)))
    fi = ForceIntensity(FFlt[+magn]);#leaving the domain
    @time F2 = (-1.0)* distribloads(el1femm,  geom,  Temp,  fi,  2);
    
    @time F3 = nzebcloadsconductivity(femm,  geom,  Temp);
    
    
    @time K = cholfact(K)
    @time U = K\(F1+F2+F3)
    @time scattersysvec!(Temp, U[:])
    
    println("Total time elapsed = ", time() - t0, "s")
    
    File =  "annulusq8.vtk"
    vtkexportmesh(File, connasarray(fes),  [geom.values Temp.values], FinEtools.MeshExportModule.Q8; scalars=[("Temperature", Temp.values)])
    
    println("Minimum/maximum temperature= $(minimum(Temp.values))/$(maximum(Temp.values)))")
    
    true
    
end # annulus_Q8_example


function ebc_annulus_Q4_algo()
    println("""
    Annular region, ingoing and outgoing flux. Temperature at one node prescribed.
    Minimum/maximum temperature ~(+/-)0.500.
    Mesh of serendipity quadrilaterals.
    This version uses the FinEtools algorithm module.
    Version: 05/29/2017
    """)
    
    t0 = time()
    
    kappa = 0.2*[1.0 0; 0 1.0]; # conductivity matrix
    te1 = -0.5
    te2 = 0.6
    hconv1, hconv2 = 1000.0, 1000.0
    rin =  1.0;#internal radius
    rex =  2.0; #external radius
    nr = 7; nc = 90;
    Angle = 2*pi;
    thickness =  1.0;
    tolerance = min(rin/nr,  rin/nc/2/pi)/10000;
    
    fens, fes = Q4annulus(rin, rex, nr, nc, Angle)
    fens, fes = mergenodes(fens,  fes,  tolerance);
    edge_fes = meshboundary(fes);
    
    # The convection boundary condition is applied at two pieces of surface
    # Side 1
    l1 = selectelem(fens, edge_fes, box=[-1.1*rex -0.9*rex -0.5*rex 0.5*rex]);
    el1femm = FEMMHeatDiffSurf(IntegData(subset(edge_fes, l1),  GaussRule(1, 2)), hconv1)
    cbc1 = FDataDict("femm"=>el1femm, "ambient_temperature"=>te1)
    # Side 2
    l2=selectelem(fens,edge_fes,box=[0.9*rex 1.1*rex -0.5*rex 0.5*rex]);
    el2femm = FEMMHeatDiffSurf(IntegData(subset(edge_fes, l2),  GaussRule(1, 2)), hconv2)
    cbc2 = FDataDict("femm"=>el2femm, "ambient_temperature"=>te2)
    
    material = MatHeatDiff(kappa)
    femm = FEMMHeatDiff(IntegData(fes,  GaussRule(2, 2)),  material)
    region1 = FDataDict("femm"=>femm)
    
    # Make model data
    modeldata = FDataDict("fens"=>fens,
    "regions"=>[region1],
    "convection_bcs"=>[cbc1, cbc2]);
    
    # Call the solver
    modeldata = FinEtools.AlgoHeatDiffModule.steadystate(modeldata)
    geom=modeldata["geom"]
    Temp=modeldata["temp"]
    println("Minimum/maximum temperature= $(minimum(Temp.values))/$(maximum(Temp.values)))")
    
    println("Total time elapsed = ",time() - t0,"s")
    
    # Postprocessing
    File = "annulusmod.vtk"
    vtkexportmesh(File, connasarray(fes), [geom.values Temp.values],  FinEtools.MeshExportModule.Q4; scalars=[("Temperature", Temp.values)])
    @async run(`"paraview.exe" $File`)
    
end # ebc_annulus_Q4_algo

function allrun()
    println("#####################################################") 
    println("# annulus_Q4_example ")
    annulus_Q4_example()
    println("#####################################################") 
    println("# annulus_Q4_example_algo ")
    annulus_Q4_example_algo()
    println("#####################################################") 
    println("# annulus_Q8_example ")
    annulus_Q8_example()
    println("#####################################################") 
    println("# ebc_annulus_Q4_algo ")
    ebc_annulus_Q4_algo()
end # function allrun

end # module annulus_examples
