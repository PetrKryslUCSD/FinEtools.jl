using FinEtools
using FinEtools.AlgoAcoustModule

println("""
Annular region, pressure BC + rigid wall.
This version uses the FinEtools algorithm module.
Version: 08/21/2017
""")

t0 = time()

rho = 1001*phun("kg/m^3");# mass density
c  = 1500.0*phun("m/s");# sound speed
bulk =  c^2*rho;
omega =  2000*phun("rev/s");      # frequency of the piston
rin =  1.0*phun("m");#internal radius

rex =  2.0*phun("m"); #external radius
nr = 20; nc = 120;
Angle = 2*pi;
thickness =  1.0*phun("m/s");
tolerance = min(rin/nr,  rin/nc/2/pi)/10000;

fens, fes = Q4annulus(rin, rex, nr, nc, Angle)
fens, fes = mergenodes(fens,  fes,  tolerance);
edge_fes = meshboundary(fes);

# The pressure boundary condition
l1 = selectelem(fens, edge_fes, box=[-1.1*rex -0.9*rex -0.5*rex 0.5*rex]);
ebc1 = FDataDict("node_list"=>connectednodes(subset(edge_fes, l1)),
 "pressure"=>x -> cos(2*pi*x[2]/rin)+1im*sin(2*pi*x[2]/rin)) # entering the domain

material = MatAcoustFluid(bulk, rho)
femm = FEMMAcoust(IntegData(fes,  GaussRule(2, 2)),  material)
region1 = FDataDict("femm"=>femm)

# Make model data
modeldata = FDataDict("fens"=>fens,
"omega"=>omega,
"regions"=>[region1], "essential_bcs"=>[ebc1]);

# Call the solver
modeldata = FinEtools.AlgoAcoustModule.steadystate(modeldata)
geom=modeldata["geom"]
P=modeldata["P"]
println("Minimum/maximum pressure, real= $(minimum(real(P.values)))/$(maximum(real(P.values))))")
println("Minimum/maximum pressure, imag= $(minimum(imag(P.values)))/$(maximum(imag(P.values))))")

println("Total time elapsed = ",time() - t0,"s")

# Postprocessing
File = "acou_annulusmod.vtk"
vtkexportmesh(File, fes.conn, geom.values,
FinEtools.MeshExportModule.Q4; scalars=[("Pre", real(P.values)), ("Pim", imag(P.values))])
@async run(`"paraview.exe" $File`)
