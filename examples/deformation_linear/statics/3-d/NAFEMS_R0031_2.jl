using FinEtools
using FinEtools.AlgoDeforLinearModule

println("""
R0031(2): Wrapped thick cylinder under pressure and thermal loading.
""")

t0 = time()
# Orthotropic material parameters of the external cylinder
E1s = 130.0*phun("GPa")
E2s = 5.0*phun("GPa")
E3s = E2s
nu12s = nu13s = 0.25
nu23s = 0.0
G12s = 10.0*phun("GPa")
G13s = G12s
G23s = 5.0*phun("GPa")
CTE1 = 3.0e-6
CTE2 = 2.0e-5
CTE3 = 2.0e-5
# Isotropic material parameters of the internal cylinder
E = 2.1e5*phun("MPa")
nu = 0.3
CTE = 2.0e-5

L = 200.0*phun("mm") # length of the cylinder
ri = 23.0*phun("mm") # internal radius of the cylinder
ti = 2.0*phun("mm") # thickness of the internal cylinder
te = 2.0*phun("mm") # thickness of the external cylinder
q0 = 200.0*phun("MPa") # internal pressure
dT = 130*phun("K") # temperature rise

tolerance = 0.0001*ti

# Generate mesh
nL = 8 # number of elements lengthwise
nc = 8 # number of elements circumferentially
xs = collect(linspace(0.0, L/2, nL+1))
ys = collect(linspace(0.0, pi/2, nc+1))
ts = [ti; te];# layer thicknesses
nts= 3*ones(Int, length(ts));# number of elements per layer
fens,fes = H8layeredplatex(xs, ys, ts, nts)
fens,fes = H8toH20(fens,fes)
bfes = meshboundary(fes)
# internal surface  for the pressure loading
intl = selectelem(fens, bfes; facing=true, direction = [0.0 0.0 -1.0])
# Shape into a cylinder
for i = 1:count(fens)
  z = fens.xyz[i,1]; a = fens.xyz[i,2]; t = fens.xyz[i,3];
  fens.xyz[i,:] = [(ri+t)*cos(pi/2-a) (ri+t)*sin(pi/2-a) z];
end

MR = DeforModelRed3D
externalmaterial = MatDeforElastOrtho(MR,
  0.0, E1s, E2s, E3s,
  nu12s, nu13s, nu23s,
  G12s, G13s, G23s,
  CTE1, CTE2, CTE3)
internalmaterial = MatDeforElastIso(MR,
  0.0, E, nu, CTE)

function updatecs!(csmatout::FFltMat, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)
  rotmat3!(csmatout, angles[fe_label]/180.*pi* [0.;0.;1.]);
end

gr = GaussRule(3, 3)

internalregion = FDataDict("femm"=>FEMMDeforLinear(MR,
    FEMMBase(fes, gr), internalmaterial))
externalregion = FDataDict("femm"=>FEMMDeforLinear(MR,
    FEMMBase(fes, gr, CSys(updatecs!)), externalmaterial))

lx0 = selectnode(fens, box=[0.0 0.0 -Inf Inf -Inf Inf], inflate=tolerance)
ly0 = selectnode(fens, box=[-Inf Inf 0.0 0.0 -Inf Inf], inflate=tolerance)
lz0 = selectnode(fens, box=[-Inf Inf -Inf Inf 0.0 0.0], inflate=tolerance)

ex0 = FDataDict( "displacement"=>  0.0, "component"=> 1, "node_list"=>lx0 )
ey0 = FDataDict( "displacement"=>  0.0, "component"=> 2, "node_list"=>ly0 )
ez0 = FDataDict( "displacement"=>  0.0, "component"=> 3, "node_list"=>lz0 )


Trac = FDataDict("traction_vector"=>[0.0; 0.0; -q0],
    "femm"=>FEMMBase(subset(bfes, ttopl), GaussRule(2, 3)))

modeldata = FDataDict("fens"=>fens,
 "regions"=>[region],
 "essential_bcs"=>[ex02, ex03, exa2, exa3, ey01, ey03, eyb1, eyb3],
 "traction_bcs"=> [Trac]
 )
modeldata = AlgoDeforLinearModule.linearstatics(modeldata)

u = modeldata["u"]
geom = modeldata["geom"]
lcenter = selectnode(fens, box=[a/2 a/2  b/2 b/2 -Inf Inf], inflate=tolerance)
cdis = abs(mean(u.values[lcenter, 3]))
println("")
println("Normalized Center deflection: $(cdis/wc_analytical)")

File =  "NAFEMS-R0031-2-plate.vtk"
vtkexportmesh(File, fes.conn, geom.values, FinEtools.MeshExportModule.H20)
@async run(`"paraview.exe" $File`)

println("Done")
true
