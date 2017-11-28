module tank_examples
using FinEtools
using FinEtools.AlgoAcoustModule

function tank_piston_platten()
    println("""
    Axially symmetric geometry  of the tank plus piston plus platten
    Version: 08/26/2017
    """)
    
    t0 = time()
    
    rho = 1001*phun("kg/m^3");# mass density
    c  = 1500.0*phun("m/s");# sound speed
    bulk =  c^2*rho;
    omega =  200*phun("rev/s");      # frequency of the piston
    a_piston = 1.0*phun("m/s")
    
    axisymmetric = true
    Platten_radius = 12*phun("mm")
    Piston_radius = 33*phun("mm")
    @assert Piston_radius > Platten_radius
    Tank_radius = 50*phun("mm")
    @assert Tank_radius > Piston_radius
    Platten_height = 100*phun("mm")
    Piston_height = 33*phun("mm")
    Tank_height = Platten_height + Piston_height + 75*phun("mm")
    Gap_height = Tank_height-(Platten_height + Piston_height)
    Tolerance = Platten_radius/1.0e6
    
    
    xs = vec([0.0 Platten_radius Piston_radius Tank_radius])
    ys = vec(vcat(vec(linspace(0.0, Platten_height, 4)),
    vec(linspace(Platten_height+Gap_height/4, Platten_height+Gap_height, 3)),
    vec(linspace(Platten_height+Gap_height+Piston_height/4, Tank_height, 3))))
    
    fens, fes = Q4blockx(xs, ys)
    
    l1 = selectelem(fens, fes, box=[0.0 Piston_radius-10*Tolerance Tank_height-Piston_height+10*Tolerance Tank_height],
    inflate=Tolerance, allin=false);
    l2 = selectelem(fens, fes, box=[0.0 Platten_radius-10*Tolerance 0.0 Platten_height-10*Tolerance],
    inflate=Tolerance, allin=false);
    fes = subset(fes, setdiff(1:count(fes), vcat(l1, l2)))
    for ixxxx in 1:3
        fens, fes = Q4refine(fens, fes)
    end
    connected = findunconnnodes(fens, fes);
    fens, new_numbering = compactnodes(fens, connected);
    fes = renumberconn!(fes, new_numbering);
    
    # Postprocessing
    # File = "tank_mesh.vtk"
    # vtkexportmesh(File, fes.conn, fens.xyz, FinEtools.MeshExportModule.Q4)
    # @async run(`"paraview.exe" $File`)
    
    
    edge_fes = meshboundary(fes);
    
    material = MatAcoustFluid(bulk, rho)
    # The pressure boundary condition
    l1 = selectelem(fens, edge_fes, box=[0.0 Piston_radius Tank_height-Piston_height Tank_height-Piston_height]);
    flux1  =  FDataDict("femm"=>FEMMAcoustSurf(IntegData(subset(edge_fes, l1), GaussRule(1, 2), axisymmetric),
    material),  "normal_flux"=> -rho*a_piston+0.0im);
    l2 = selectelem(fens, edge_fes, box=[Piston_radius Tank_radius Tank_height Tank_height]);
    ebc2 = FDataDict("node_list"=>connectednodes(subset(edge_fes, l2)),
    "pressure"=>x -> 0.0) # entering the domain
    
    femm = FEMMAcoust(IntegData(fes,  GaussRule(2, 2), axisymmetric),  material)
    region1 = FDataDict("femm"=>femm)
    
    # Make model data
    modeldata = FDataDict("fens"=>fens,"omega"=>omega,
    "regions"=>[region1], "flux_bcs"=>[flux1], "essential_bcs"=>[ebc2]);
    
    # Call the solver
    modeldata = FinEtools.AlgoAcoustModule.steadystate(modeldata)
    geom = modeldata["geom"]
    P = modeldata["P"]
    println("Minimum/maximum pressure, real= $(minimum(real(P.values)))/$(maximum(real(P.values))))")
    println("Minimum/maximum pressure, imag= $(minimum(imag(P.values)))/$(maximum(imag(P.values))))")
    
    println("Total time elapsed = ",time() - t0,"s")
    
    # Postprocessing
    File = "tank_piston_platten.vtk"
    vtkexportmesh(File, fes.conn, geom.values, FinEtools.MeshExportModule.Q4; scalars=[("Pre", real(P.values)), ("Pim", imag(P.values))])
    @async run(`"paraview.exe" $File`)
    
end # tank_piston_platten


function tank_piston_platten_pressure()
   println("""
    Axially symmetric geometry  of the tank plus piston plus platten
    Version: 08/26/2017
    """)
    
    t0 = time()
    
    rho = 1001*phun("kg/m^3");# mass density
    c  = 1500.0*phun("m/s");# sound speed
    bulk =  c^2*rho;
    omega =  200*phun("rev/s");      # frequency of the piston
    
    axisymmetric = true
    Platten_radius = 15*phun("mm")
    Piston_radius = 35*phun("mm")
    @assert Piston_radius > Platten_radius
    Tank_radius = 50*phun("mm")
    @assert Tank_radius > Piston_radius
    Platten_height = 75*phun("mm")
    Piston_height = 75*phun("mm")
    Tank_height = Platten_height + Piston_height + 75*phun("mm")
    Gap_height = Tank_height-(Platten_height + Piston_height)
    Tolerance = Platten_radius/1.0e6
    
    xs = vec([0.0 Platten_radius Piston_radius Tank_radius])
    ys = vec(vcat(vec(linspace(0.0, Platten_height, 4)),
    vec(linspace(Platten_height+Gap_height/4, Platten_height+Gap_height, 3)),
    vec(linspace(Platten_height+Gap_height+Piston_height/4, Tank_height, 3))))
    
    fens, fes = Q4blockx(xs, ys)
    
    l1 = selectelem(fens, fes, box=[0.0 Piston_radius-10*Tolerance Tank_height-Piston_height+10*Tolerance Tank_height],
    inflate=Tolerance, allin=false);
    l2 = selectelem(fens, fes, box=[0.0 Platten_radius-10*Tolerance 0.0 Platten_height-10*Tolerance],
    inflate=Tolerance, allin=false);
    fes = subset(fes, setdiff(1:count(fes), vcat(l1, l2)))
    for ixxxx in 1:3
        fens, fes = Q4refine(fens, fes)
    end
    connected = findunconnnodes(fens, fes);
    fens, new_numbering = compactnodes(fens, connected);
    fes = renumberconn!(fes, new_numbering);
    
    # Postprocessing
    # File = "tank_mesh.vtk"
    # vtkexportmesh(File, fes.conn, fens.xyz, FinEtools.MeshExportModule.Q4)
    # @async run(`"paraview.exe" $File`)
    
    
    edge_fes = meshboundary(fes);
    
    # The pressure boundary condition
    l1 = selectelem(fens, edge_fes, box=[0.0 Piston_radius Tank_height-Piston_height Tank_height-Piston_height]);
    ebc1 = FDataDict("node_list"=>connectednodes(subset(edge_fes, l1)),
    "pressure"=>x -> 1.0) # entering the domain
    l2 = selectelem(fens, edge_fes, box=[Piston_radius Tank_radius Tank_height Tank_height]);
    ebc2 = FDataDict("node_list"=>connectednodes(subset(edge_fes, l2)),
    "pressure"=>x -> 0.0) # entering the domain
    material = MatAcoustFluid(bulk, rho)
    femm = FEMMAcoust(IntegData(fes,  GaussRule(2, 2), axisymmetric),  material)
    region1 = FDataDict("femm"=>femm)
    
    # Make model data
    modeldata = FDataDict("fens"=>fens,"omega"=>omega,
    "regions"=>[region1], "essential_bcs"=>[ebc1 ebc2]);
    
    # Call the solver
    modeldata = FinEtools.AlgoAcoustModule.steadystate(modeldata)
    geom = modeldata["geom"]
    P = modeldata["P"]
    println("Minimum/maximum pressure, real= $(minimum(real(P.values)))/$(maximum(real(P.values))))")
    println("Minimum/maximum pressure, imag= $(minimum(imag(P.values)))/$(maximum(imag(P.values))))")
    
    println("Total time elapsed = ",time() - t0,"s")
    
    # Postprocessing
    File = "tank_piston_platten_pressure.vtk"
    vtkexportmesh(File, fes.conn, geom.values, FinEtools.MeshExportModule.Q4; scalars=[("Pre", real(P.values)), ("Pim", imag(P.values))])
    @async run(`"paraview.exe" $File`)
    
end # tank_piston_platten_pressure

function allrun()
    println("#####################################################") 
    println("# tank_piston_platten ")
    tank_piston_platten()
    println("#####################################################") 
    println("# tank_piston_platten_pressure ")
    tank_piston_platten_pressure()
    return true
end # function allrun

end # module tank_examples
