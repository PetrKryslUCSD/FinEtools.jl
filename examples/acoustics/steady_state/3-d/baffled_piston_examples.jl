module baffled_piston_examples
using FinEtools


function baffled_piston_H8_ABC_example_algo()
    rho = 1.21*phun("kg/m^3");# mass density
    c  = 343.0*phun("m/s");# sound speed
    bulk =  c^2*rho;
    omega =  7500*phun("rev/s");      # frequency of the piston
    a_piston =  -1.0*phun("mm/s")     # amplitude of the piston acceleration
    R = 50.0*phun("mm");# radius of the piston
    Ro = 150.0*phun("mm"); # radius of the enclosure
    nref = 4;#number of refinements of the sphere around the piston
    nlayers = 35;                     # number of layers of elements surrounding the piston
    tolerance = R/(2^nref)/100
    
    println("""
    
    Baffled piston in a half-sphere domain with ABC.
    
    Hexahedral mesh. Algorithm version.
    """)
    
    t0  =  time()
    
    # Hexahedral mesh
    fens,fes  =  H8sphere(R,nref);
    bfes  =  meshboundary(fes)
    File  =   "baffledabc_boundary.vtk"
    vtkexportmesh(File, connasarray(bfes), fens.xyz, FinEtools.MeshExportModule.Q4)
    @async run(`"paraview.exe" $File`)
    
    l = selectelem(fens,bfes,facing = true,direction = [1.0 1.0  1.0], dotmin= 0.001)
    ex(xyz, layer) = (R+layer/nlayers*(Ro-R))*xyz/norm(xyz)
    fens1,fes1  =  H8extrudeQ4(fens, subset(bfes,l), nlayers, ex);
    fens,newfes1,fes2 =  mergemeshes(fens1, fes1, fens, fes, tolerance)
    fes = cat(newfes1,fes2)
    
    # Piston surface mesh
    bfes  =  meshboundary(fes)
    l1 = selectelem(fens, bfes, facing = true, direction = [-1.0 0.0 0.0])
    l2 = selectelem(fens, bfes, distance = R, from = [0.0 0.0 0.0], inflate = tolerance)
    piston_fes = subset(bfes,intersect(l1,l2));
    
    # Outer spherical boundary
    louter = selectelem(fens, bfes, facing = true, direction = [1.0 1.0  1.0], dotmin= 0.001)
    outer_fes = subset(bfes,louter);
    
    println("Pre-processing time elapsed  =  ",time() - t0,"s")
    
    t1  =  time()
    
    material = MatAcoustFluid(bulk, rho)
    # Region of the fluid
    region1 =  FDataDict("femm"=>FEMMAcoust(IntegData(fes, GaussRule(3, 2)), material))
    
    # Surface for the ABC
    abc1  =  FDataDict("femm"=>FEMMAcoustSurf(IntegData(outer_fes, GaussRule(2, 2)),
    material))
    
    # Surface of the piston
    flux1  =  FDataDict("femm"=>FEMMAcoustSurf(IntegData(piston_fes, GaussRule(2, 2)),
    material),  "normal_flux"=> -rho*a_piston+0.0im);
    
    # Make model data
    modeldata =  FDataDict("fens"=>  fens,
    "omega"=>omega,
    "regions"=>[region1],
    "flux_bcs"=>[flux1], "ABCs"=>[abc1])
    
    # Call the solver
    modeldata = FinEtools.AlgoAcoustModule.steadystate(modeldata)
    
    println("Computing time elapsed  =  ",time() - t1,"s")
    println("Total time elapsed  =  ",time() - t0,"s")
    
    geom = modeldata["geom"]
    P = modeldata["P"]
    
    File  =   "baffledabc.vtk"
    vtkexportmesh(File, connasarray(fes), geom.values, FinEtools.MeshExportModule.H8;
    scalars = [("absP", abs.(P.values))])
    @async run(`"paraview.exe" $File`)
end # baffled_piston_H8_ABC_example_algo

function allrun()
    println("#####################################################") 
    println("# baffled_piston_H8_ABC_example_algo ")
    baffled_piston_H8_ABC_example_algo()
    return true
end # function allrun

end # module baffled_piston_examples
