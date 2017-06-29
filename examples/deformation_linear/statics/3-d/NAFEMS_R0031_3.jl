using FinEtools
using FinEtools.AlgoDeforLinearModule

println("""
NAFEMS publication R0031/3 Composite plate test.
Simply supported on all four edges.  Uniform transverse  loading.
The modeled part is one quarter of the full plate here.
""")

t0 = time()
# Skin (face) material parameters
E1s = 1.0e7*phun("psi")
E2s = 0.4e7*phun("psi")
E3s = 0.4e7*phun("psi")
nu12s = 0.3
nu13s = 0.3
nu23s = 0.3
G12s = 0.1875e7*phun("psi")
G13s = 0.1875e7*phun("psi")
G23s = 0.1875e7*phun("psi")
# Core material parameters
E1c = 10.*phun("psi")
E2c = 10.*phun("psi")
E3c = 10e4.*phun("psi")
nu12c = 0.
nu13c = 0.
nu23c = 0.
G12c = 10.*phun("psi")
G13c = 3.0e4*phun("psi")
G23c = 1.2e4*phun("psi")
L = 10.0*phun("in") # side of the square plate
nL = 8 # number of elements along the side of the plate
tolerance = 0.0001*phun("in")
xs = collect(linspace(0.0, L/2, nL+1))
ys = collect(linspace(0.0, L/2, nL+1))
ts = [0.028; 0.75; 0.028]*phun("in")
nts = [2; 3;  2; ] # number of elements through the thickness
tmag = 100*phun("psi")

# Generate mesh
fens,fes = H8compositeplatex(xs, ys, ts, nts)
fens,fes = H8toH20(fens,fes)

MR = DeforModelRed3D
skinmaterial = MatDeforElastOrtho(MR,
  0.0, E1s, E2s, E3s,
  nu12s, nu13s, nu23s,
  G12s, G13s, G23s,
  0.0, 0.0, 0.0)
corematerial = MatDeforElastOrtho(MR,
  0.0, E1c, E2c, E3c,
  nu12c, nu13c, nu23c,
  G12c, G13c, G23c,
  0.0, 0.0, 0.0)

gr = GaussRule(3, 3)

rl1 = selectelem(fens, fes, label=1)
skinbot = FDataDict("femm"=>FEMMDeforLinear(MR,
    FEMMBase(subset(fes, rl1), gr), skinmaterial))

rl3 = selectelem(fens, fes, label=3)
skintop = FDataDict("femm"=>FEMMDeforLinear(MR,
    FEMMBase(subset(fes, rl3), gr), skinmaterial))

rl2 = selectelem(fens, fes, label=2)
core = FDataDict("femm"=>FEMMDeforLinear(MR,
    FEMMBase(subset(fes, rl2), gr), corematerial))

lx0 = selectnode(fens, box=[0.0 0.0 -Inf Inf -Inf Inf], inflate=tolerance)
lxL2 = selectnode(fens, box=[L/2 L/2 -Inf Inf -Inf Inf], inflate=tolerance)
ly0 = selectnode(fens, box=[-Inf Inf 0.0 0.0 -Inf Inf], inflate=tolerance)
lyL2 = selectnode(fens, box=[-Inf Inf L/2 L/2 -Inf Inf], inflate=tolerance)

ex0 = FDataDict( "displacement"=>  0.0, "component"=> 3, "node_list"=>lx0 )
exL2 = FDataDict( "displacement"=>  0.0, "component"=> 1, "node_list"=>lxL2 )
ey0 = FDataDict( "displacement"=>  0.0, "component"=> 3, "node_list"=>ly0 )
eyL2 = FDataDict( "displacement"=>  0.0, "component"=> 2, "node_list"=>lyL2 )

bfes = meshboundary(fes)
ttopl = selectelem(fens, bfes; facing=true, direction = [0.0 0.0 1.0])
Trac = FDataDict("traction_vector"=>[0.0; 0.0; -tmag],
    "femm"=>FEMMBase(subset(bfes, ttopl), GaussRule(2, 3)))

modeldata = FDataDict("fens"=>fens,
 "regions"=>[skinbot, core, skintop], "essential_bcs"=>[ex0, exL2, ey0, eyL2],
 "traction_bcs"=> [Trac]
 )
modeldata = AlgoDeforLinearModule.linearstatics(modeldata)

u = modeldata["u"]
geom = modeldata["geom"]
lcenter = selectnode(fens, box=[L/2 L/2  L/2 L/2 -Inf Inf], inflate=tolerance)
cdis = mean(u.values[lcenter, 3])/phun("in")
println("Center node displacements $(cdis) [in]")
println("")

File =  "NAFEMS-R0031-3-plate.vtk"
vtkexportmesh(File, fes.conn, geom.values, FinEtools.MeshExportModule.H20;
    scalars = [("Layer", fes.label)], vectors = [("displacement", u.values)])
@async run(`"paraview.exe" $File`)

println("Done")
true
