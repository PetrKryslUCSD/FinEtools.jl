using FinEtools

println("""
% Vibration modes of unit cube  of almost incompressible material.
% Mean-strain hexahedron.
% Reference: Puso MA, Solberg J (2006) A stabilized nodally integrated
% tetrahedral. International Journal for Numerical Methods in
% Engineering 67: 841-867.""")
t0 = time()


E = 1*phun("PA");
nu = 0.499;
rho= 1*phun("KG/M^3");
a=1*phun("M"); b=a; h= a;
n1=8 # How many element edges per side?
na= n1; nb= n1; nh =n1;
neigvs=20                   # how many eigenvalues
omega_shift=(0.1*2*pi)^2;

fens,fes = H8block(a,b,h, na,nb,nh)

# Make the region
MR = DeforModelRed3D
material = MatDeforElastIso(MR, rho, E, nu, 0.0)
region1 = FDataDict("femm"=>FEMMDeforLinearMSH8(MR, GeoD(fes, GaussRule(3,2)),
  material), "femm_mass"=>FEMMDeforLinearMSH8(MR, GeoD(fes, GaussRule(3,3)),
  material))

# Make model data
modeldata =  FDataDict(
  "fens"=> fens, "regions"=>  [region1],
  "omega_shift"=>omega_shift, "neigvs"=>neigvs)

# Solve
modeldata = FinEtools.AlgoDeforLinearModule.modal(modeldata)

fs = modeldata["omega"]/(2*pi)
println("Eigenvalues: $fs [Hz]")

modeldata["postprocessing"] = FDataDict("file"=>"unit_cube_mode",
  "mode"=>10)
modeldata=FinEtools.AlgoDeforLinearModule.exportmode(modeldata)
@async run(`"paraview.exe" $(modeldata["postprocessing"]["file"]*"1.vtk")`)

true
