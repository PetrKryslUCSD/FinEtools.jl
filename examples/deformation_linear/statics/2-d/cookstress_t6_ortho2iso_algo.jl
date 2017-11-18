using FinEtools
using FinEtools.AlgoDeforLinearModule
using FinEtools.MeshExportModule

println("Cook plane stress, with quadratic triangles. With orthotropic  material model.")
E = 1.0;
nu = 1.0/3;
E1 = E2 = E3 = E;
nu12 = nu13 = nu23 = nu;
G12 = G13 = G23 = E/2.0/(1+nu);
width = 48.0; height = 44.0; thickness  = 1.0;
free_height  = 16.0;
Mid_edge  = [48.0, 52.0];# Location of tracked  deflection
magn = 1.0/free_height;# Magnitude of applied load
convutip = 23.97;
n = 10;#*int(round(sqrt(170.)/2.)); # number of elements per side
tolerance = minimum([width, height])/n/1000.;#Geometrical tolerance

fens,fes = T6block(width, height, n, n)

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
el1femm = FEMMBase(IntegData(subset(boundaryfes, Toplist), GaussRule(1, 3)))
flux1 = FDataDict("traction_vector"=>[0.0,+magn],
    "femm"=>el1femm
    )

# Make the region
MR = DeforModelRed2DStress
# This material model is orthotropic,  but the input parameters correspond to an
# isotropiic material  model..
material=MatDeforElastOrtho(MR, E1,E2,E3,nu12,nu13,nu23,G12,G13,G23)
region1 = FDataDict("femm"=>FEMMDeforLinear(MR,
    IntegData(fes, TriRule(3)), material))

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
File = modeldata["postprocessing"]["exported"][1]["file"]
@async run(`"paraview.exe" $File`)
fld = modeldata["postprocessing"]["exported"][1]["field"]
println("$(minimum(fld.values)) $(maximum(fld.values))")
true
