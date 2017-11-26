using FinEtools
using FinEtools.AlgoHeatDiffModule

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
