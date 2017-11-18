using FinEtools

println("""
FV32: Cantilevered tapered membrane
This is a test recommended by the National Agency for Finite Element Methods and
Standards (U.K.): Test FV32 from NAFEMS publication TNSB, Rev. 3, “The Standard
NAFEMS Benchmarks,” October 1990.

Reference solution: 44.623	130.03	162.70	246.05	379.90	391.44 for the first
six modes.
""")

t0 = time()


E = 200*phun("GPA");
nu = 0.3;
rho= 8000*phun("KG/M^3");
L = 10*phun("M");
W0 = 5*phun("M");
WL = 1*phun("M");
H = 0.05*phun("M");
nL, nW, nH = 8, 4, 1;# How many element edges per side?
neigvs=20                   # how many eigenvalues
Reffs = [44.623	130.03	162.70	246.05	379.90	391.44]

fens,fes =H20block(1.0, 2.0, 1.0, nL, nW, nH)
for i = 1:count(fens)
    xi, eta, theta = fens.xyz[i,:];
    eta = eta - 1.0
    fens.xyz[i,:] = [xi*L eta*(1.0 - 0.8*xi)*W0/2 theta*H/2];
end
# File =  "mesh.vtk"
# vtkexportmesh(File, fens, fes)
# @async run(`"paraview.exe" $File`)

# Make the region
MR = DeforModelRed3D
material = MatDeforElastIso(MR, rho, E, nu, 0.0)
region1 = FDataDict("femm"=>FEMMDeforLinear(MR, IntegData(fes, GaussRule(3,2)),
  material), "femm_mass"=>FEMMDeforLinear(MR, IntegData(fes, GaussRule(3,3)),
  material))

nl1 = selectnode(fens; plane=[1.0 0.0 0.0 0.0], thickness=H/1.0e4)
ebc1 = FDataDict("node_list"=>nl1, "component"=>1, "displacement"=>0.0)
ebc2 = FDataDict("node_list"=>nl1, "component"=>2, "displacement"=>0.0)
ebc3 = FDataDict("node_list"=>nl1, "component"=>3, "displacement"=>0.0)

nl4 = selectnode(fens; plane=[0.0 0.0 1.0 0.0], thickness=H/1.0e4)
ebc4 = FDataDict("node_list"=>nl4, "component"=>3, "displacement"=>0.0)

# Make model data
modeldata =  FDataDict(
  "fens"=> fens, "regions"=>  [region1], "essential_bcs"=>[ebc1 ebc2 ebc3 ebc4],
  "neigvs"=>neigvs)

# Solve
modeldata = FinEtools.AlgoDeforLinearModule.modal(modeldata)

fs = modeldata["omega"]/(2*pi)
println("Eigenvalues: $fs [Hz]")
println("Percentage frequency errors: $((vec(fs[1:6]) - vec(Reffs))./vec(Reffs)*100)")

modeldata["postprocessing"] = FDataDict("file"=>"FV32-modes", "mode"=>1:10)
modeldata=FinEtools.AlgoDeforLinearModule.exportmode(modeldata)
@async run(`"paraview.exe" $(modeldata["postprocessing"]["file"]*"1.vtk")`)

true
