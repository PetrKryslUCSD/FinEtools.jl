using FinEtools
using FinEtools.AlgoDeforLinearModule

println("""
R0031(2): Wrapped thick cylinder under pressure and thermal loading.
""")

t0 = time()
pu = ustring -> phun(ustring; system_of_units = :SIMM)
# Orthotropic material parameters of the outer cylinder
E1s = 130.0*pu("GPa")
E2s = 5.0*pu("GPa")
E3s = E2s
nu12s = nu13s = 0.25
nu23s = 0.0
G12s = 10.0*pu("GPa")
G13s = G12s
G23s = 5.0*pu("GPa")
CTE1 = 3.0e-6
CTE2 = 2.0e-5
CTE3 = 2.0e-5
# Isotropic material parameters of the inner cylinder
E = 2.1e5*pu("MPa")
nu = 0.3
CTE = 2.0e-5

L = 200.0*pu("mm") # length of the cylinder
ri = 23.0*pu("mm") # inner radius of the cylinder
ti = 2.0*pu("mm") # thickness of the inner cylinder
te = 2.0*pu("mm") # thickness of the outer cylinder
q0 = 200.0*pu("MPa") # inner pressure
dT = 130*pu("K") # temperature rise

tolerance = 0.0001*ti

# Generate mesh
nL = 16 # number of elements lengthwise
nc = 8 # number of elements circumferentially
xs = collect(linspace(0.0, L/2, nL+1))
ys = collect(linspace(0.0, pi/2, nc+1))
ts = [ti; te];# layer thicknesses
nts= 6*ones(Int, length(ts));# number of elements per layer
fens,fes = H8compositeplatex(xs, ys, ts, nts)
fens,fes = H8toH20(fens,fes)
bfes = meshboundary(fes)
# inner surface  for the pressure loading
intl = selectelem(fens, bfes; facing=true, direction = [0.0 0.0 -1.0])
# Shape into a cylinder
for i = 1:count(fens)
  z = fens.xyz[i,1]; a = fens.xyz[i,2]; t = fens.xyz[i,3];
  fens.xyz[i,:] = [(ri+t)*cos(pi/2-a) (ri+t)*sin(pi/2-a) z];
end

MR = DeforModelRed3D
outermaterial = MatDeforElastOrtho(MR,
  0.0, E1s, E2s, E3s,
  nu12s, nu13s, nu23s,
  G12s, G13s, G23s,
  CTE1, CTE2, CTE3)
innermaterial = MatDeforElastIso(MR,
  0.0, E, nu, CTE)

function cylcs!(csmatout::FFltMat, XYZ::FFltMat)
  csmatout[:, 2] = [0.0 0.0 1.0]
  csmatout[:, 3] = XYZ
  csmatout[3, 3] = 0.0
  csmatout[:, 3] = csmatout[:, 3]/norm(csmatout[:, 3])
  csmatout[:, 1] = cross(csmatout[:, 2], csmatout[:, 3])
end

function updatecs!(csmatout::FFltMat, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)
  cylcs!(csmatout, XYZ)
end

gr = GaussRule(3, 3)

rli = selectelem(fens, fes, label=1)
innerregion = FDataDict("femm"=>FEMMDeforLinear(MR,
    IntegData(subset(fes, rli), gr), innermaterial))
rle = selectelem(fens, fes, label=2)
outerregion = FDataDict("femm"=>FEMMDeforLinear(MR,
    IntegData(subset(fes, rle), gr, CSys(3, 3, updatecs!)), outermaterial))

lx0 = selectnode(fens, box=[0.0 0.0 -Inf Inf -Inf Inf], inflate=tolerance)
ly0 = selectnode(fens, box=[-Inf Inf 0.0 0.0 -Inf Inf], inflate=tolerance)
lz0 = selectnode(fens, box=[-Inf Inf -Inf Inf 0.0 0.0], inflate=tolerance)

ex0 = FDataDict( "displacement"=>  0.0, "component"=> 1, "node_list"=>lx0 )
ey0 = FDataDict( "displacement"=>  0.0, "component"=> 2, "node_list"=>ly0 )
ez0 = FDataDict( "displacement"=>  0.0, "component"=> 3, "node_list"=>lz0 )

function getpr!(forceout::FFltVec, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)
  csmatout = zeros(3, 3)
  cylcs!(csmatout, XYZ)
  copy!(forceout, q0*csmatout[:, 3])
end

Trac = FDataDict("traction_vector"=>getpr!,
    "femm"=>FEMMBase(IntegData(subset(bfes, intl), GaussRule(2, 3)))
    )

modeldata = FDataDict("fens"=>fens,
 "regions"=>[innerregion, outerregion],
 "essential_bcs"=>[ex0, ey0, ez0],
 "traction_bcs"=>[Trac],
 "temperature_change"=>FDataDict("temperature"=>dT)
 )
modeldata = AlgoDeforLinearModule.linearstatics(modeldata)

u = modeldata["u"]
geom = modeldata["geom"]
# File =  "NAFEMS-R0031-2-plate.vtk"
# vtkexportmesh(File, fes.conn, geom.values, FinEtools.MeshExportModule.H20;
#     scalars = [("Layer", fes.label)], vectors = [("displacement", u.values)])
# @async run(`"paraview.exe" $File`)

modeldata["postprocessing"] = FDataDict("file"=>"NAFEMS-R0031-2",
  "outputcsys"=>CSys(3, 3, updatecs!), "quantity"=>:Cauchy, "component"=>6)
modeldata = AlgoDeforLinearModule.exportstress(modeldata)
modeldata["postprocessing"] = FDataDict("file"=>"NAFEMS-R0031-2-elem",
  "outputcsys"=>CSys(3, 3, updatecs!), "quantity"=>:Cauchy, "component"=>6)
modeldata = AlgoDeforLinearModule.exportstresselementwise(modeldata)

println("Done")
true
