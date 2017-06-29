using FinEtools
using FinEtools.AlgoDeforLinearModule
using FinEtools.MeshExportModule

E = 1.0;
nu = 1.0/3;
width = 48.0; height = 44.0; thickness  = 1.0;
free_height  = 16.0;
Mid_edge  = [48.0, 52.0];# Location of tracked  deflection
magn = 1./free_height;# Magnitude of applied load
convutip = 23.97;
n = 20;#*int(round(sqrt(170.)/2.)); # number of elements per side
tolerance = minimum([width, height])/n/1000.;#Geometrical tolerance

fens,fes = T3block(width, height, n, n)

# Reshape into a trapezoidal panel
for i = 1:count(fens)
    fens.xyz[i,2] = fens.xyz[i,2]+(fens.xyz[i,1]/width)*(height -fens.xyz[i,2]/height*(height-free_height));
end

# Clamped edge of the membrane
l1 = selectnode(fens; box=[0.,0.,-Inf, Inf], inflate = tolerance)
ess1 = FDataDict("displacement"=>  0.0, "component"=> 1, "node_list"=>l1)
ess2 = FDataDict("displacement"=>  0.0, "component"=> 2, "node_list"=>l1)

# Traction on the opposite edge
boundaryfes =  meshboundary(fes);
Toplist  = selectelem(fens, boundaryfes, box= [width, width, -Inf, Inf ], inflate=  tolerance);
el1femm = FEMMBase(GeoD(subset(boundaryfes, Toplist), GaussRule(1, 2)))
flux1 = FDataDict("traction_vector"=>[0.0,+magn],
    "femm"=>el1femm
    )

# Make the region
MR = DeforModelRed2DStress
material = MatDeforElastIso(MR,  0.0, E, nu, 0.0)
region1 = FDataDict("femm"=>FEMMDeforLinear(MR,
    GeoD(fes, TriRule(1)), material))

modeldata = FDataDict("fens"=>fens,
 "regions"=>[region1],
 "essential_bcs"=>[ess1, ess2],
 "traction_bcs"=>[flux1]
 )

# Call the solver
modeldata = AlgoDeforLinearModule.linearstatics(modeldata)

u = modeldata["u"]
geom = modeldata["geom"]

# Extract the solution
nl = selectnode(fens, box=[Mid_edge[1],Mid_edge[1],Mid_edge[2],Mid_edge[2]],
          inflate=tolerance);
theutip = u.values[nl,:]
println("displacement =$(theutip[2]) as compared to converged $convutip")

modeldata["postprocessing"] = FDataDict("file"=>"cookstress",
   "quantity"=>:Cauchy, "component"=>:xy)
modeldata = AlgoDeforLinearModule.exportstress(modeldata)

true
