using FinEtools
using FinEtools.AlgoDeforLinearModule
using FinEtools.MeshExportModule

mu = 756000
nu = 0.49546;
E = 2*(1.0 + nu) * mu;
width = 48.0; height = 44.0; thickness  = 10.0;
free_height  = 16.0;
Mid_edge  = [48.0 52.0 0.0];# Location of tracked  deflection
magn = 1.0/(free_height*thickness);# Density of applied load
nt = 1
n = 16*nt;# number of elements per side
tolerance = minimum([width, height])/n/1000.;#Geometrical tolerance

fens,fes = T10block(width, height, thickness, n, n, nt)

# Reshape into a trapezoidal panel
for i = 1:count(fens)
    fens.xyz[i,2] = fens.xyz[i,2]+(fens.xyz[i,1]/width)*(height -fens.xyz[i,2]/height*(height-free_height));
end

# Clamped edge of the membrane
l1 = selectnode(fens; box=[0.,0.,-Inf, Inf,-Inf, Inf], inflate = tolerance)
ess1 = FDataDict("displacement"=>  0.0, "component"=> 1, "node_list"=>l1)
ess2 = FDataDict("displacement"=>  0.0, "component"=> 2, "node_list"=>l1)
ess3 = FDataDict("displacement"=>  0.0, "component"=> 3, "node_list"=>l1)

# Traction on the opposite edge
boundaryfes =  meshboundary(fes);
Toplist  = selectelem(fens, boundaryfes, box= [width, width, -Inf, Inf, -Inf, Inf], inflate=  tolerance);
el1femm = FEMMBase(IntegData(subset(boundaryfes, Toplist), SimplexRule(2, 3)))
flux1 = FDataDict("traction_vector"=>[0.0, 0.0, +magn],
    "femm"=>el1femm
    )

# Make the region
MR = DeforModelRed3D
material = MatDeforElastIso(MR,  0.0, E, nu, 0.0)
region1 = FDataDict("femm"=>FEMMDeforLinearMST10(MR,
    IntegData(fes, SimplexRule(3, 2)), material))

modeldata = FDataDict("fens"=>fens,
 "regions"=>[region1],
 "essential_bcs"=>[ess1, ess2, ess3],
 "traction_bcs"=>[flux1]
 )

# Call the solver
modeldata = AlgoDeforLinearModule.linearstatics(modeldata)

u = modeldata["u"]
geom = modeldata["geom"]

# Extract the solution
nl = selectnode(fens, box=[Mid_edge[1],Mid_edge[1],Mid_edge[2],Mid_edge[2], -Inf, Inf],
          inflate=tolerance);
theutip = mean(u.values[nl,:])
println("displacement =$(theutip[2]) as compared to converged $convutip")

modeldata["postprocessing"] = FDataDict("file"=>"cook_like-ew",
   "quantity"=>:Cauchy, "component"=>:xy)
modeldata = AlgoDeforLinearModule.exportstresselementwise(modeldata)
fld = modeldata["postprocessing"]["exported"][1]["field"]
println("range of Cauchy_xy = $((minimum(fld.values), maximum(fld.values)))")
File = modeldata["postprocessing"]["exported"][1]["file"]
@async run(`"paraview.exe" $File`)

modeldata["postprocessing"] = FDataDict("file"=>"cook_like-ew-vm",
   "quantity"=>:vm, "component"=>1)
modeldata = AlgoDeforLinearModule.exportstresselementwise(modeldata)
fld = modeldata["postprocessing"]["exported"][1]["field"]
println("range of vm = $((minimum(fld.values), maximum(fld.values)))")
File = modeldata["postprocessing"]["exported"][1]["file"]
@async run(`"paraview.exe" $File`)

modeldata["postprocessing"] = FDataDict("file"=>"cook_like-ew-pressure",
   "quantity"=>:pressure, "component"=>1)
modeldata = AlgoDeforLinearModule.exportstresselementwise(modeldata)
fld = modeldata["postprocessing"]["exported"][1]["field"]
println("range of pressure = $((minimum(fld.values), maximum(fld.values)))")
File = modeldata["postprocessing"]["exported"][1]["file"]
@async run(`"paraview.exe" $File`)

modeldata["postprocessing"] = FDataDict("file"=>"cook_like-ew-princ1",
   "quantity"=>:princCauchy, "component"=>1)
modeldata = AlgoDeforLinearModule.exportstresselementwise(modeldata)
fld = modeldata["postprocessing"]["exported"][1]["field"]
println("range of princCauchy Max = $((minimum(fld.values), maximum(fld.values)))")
File = modeldata["postprocessing"]["exported"][1]["file"]
@async run(`"paraview.exe" $File`)

modeldata["postprocessing"] = FDataDict("file"=>"cook_like-ew-princ3",
   "quantity"=>:princCauchy, "component"=>3)
modeldata = AlgoDeforLinearModule.exportstresselementwise(modeldata)
fld = modeldata["postprocessing"]["exported"][1]["field"]
println("range of princCauchy Min = $((minimum(fld.values), maximum(fld.values)))")
File = modeldata["postprocessing"]["exported"][1]["file"]
@async run(`"paraview.exe" $File`)

true
