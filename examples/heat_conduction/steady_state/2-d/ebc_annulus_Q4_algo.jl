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
el1femm = FEMMHeatDiffSurf(GeoD(subset(edge_fes, l1),  GaussRule(1, 2)), hconv1)
cbc1 = FDataDict("femm"=>el1femm, "ambient_temperature"=>te1)
# Side 2
l2=selectelem(fens,edge_fes,box=[0.9*rex 1.1*rex -0.5*rex 0.5*rex]);
el2femm = FEMMHeatDiffSurf(GeoD(subset(edge_fes, l2),  GaussRule(1, 2)), hconv2)
cbc2 = FDataDict("femm"=>el2femm, "ambient_temperature"=>te2)

material = MatHeatDiff(kappa)
femm = FEMMHeatDiff(GeoD(fes,  GaussRule(2, 2)),  material)
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
vtkexportmesh(File, fes.conn, [geom.values Temp.values],
FinEtools.MeshExportModule.Q4; scalars=[("Temperature", Temp.values)])
@async run(`"paraview.exe" $File`)
