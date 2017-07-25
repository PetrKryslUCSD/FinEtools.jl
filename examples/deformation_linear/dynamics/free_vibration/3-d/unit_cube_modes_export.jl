using FinEtools
using FinEtools.AbaqusExportModule

println("""
Vibration modes of unit cube  of almost incompressible material.

This example EXPORTS the model to Abaqus.

Reference: Puso MA, Solberg J (2006) A stabilized nodally integrated
tetrahedral. International Journal for Numerical Methods in
Engineering 67: 841-867.
""")
t0 = time()


E = 1*phun("PA");
nu = 0.499;
rho = 1*phun("KG/M^3");
a = 1*phun("M"); b = a; h =  a;
n1 = 5;# How many element edges per side?
na =  n1; nb =  n1; nh  = n1;
neigvs = 20                   # how many eigenvalues
OmegaShift = (0.01*2*pi)^2;

MR = DeforModelRed3D
fens,fes  = H20block(a,b,h, na,nb,nh)

geom = NodalField(fens.xyz)
u = NodalField(zeros(size(fens.xyz,1),3)) # displacement field

numberdofs!(u)

material=MatDeforElastIso(MR, rho, E, nu, 0.0)

femm = FEMMDeforLinear(MR, GeoD(fes, GaussRule(3,2)), material)

K =stiffness(femm, geom, u)
femm = FEMMDeforLinear(MR, GeoD(fes, GaussRule(3,3)), material)
M =mass(femm, geom, u)
d,v,nev,nconv = eigs(K+OmegaShift*M, M; nev=neigvs, which=:SM)
d = d - OmegaShift;
fs = real(sqrt.(complex(d)))/(2*pi)
println("Eigenvalues: $fs [Hz]")

mode = 7
scattersysvec!(u, v[:,mode])
File =  "unit_cube_modes.vtk"
vtkexportmesh(File, fens, fes; vectors=[("mode$mode", u.values)])


AE = AbaqusExporter("unit_cube_modes_h20");
# AE.ios = STDOUT;
HEADING(AE, "Vibration modes of unit cube  of almost incompressible material.");
COMMENT(AE, "The  first six frequencies are rigid body modes.");
COMMENT(AE, "The  first nonzero frequency (7) should be around 0.26 Hz");
PART(AE, "part1");
END_PART(AE);
ASSEMBLY(AE, "ASSEM1");
INSTANCE(AE, "INSTNC1", "PART1");
NODE(AE, fens.xyz);
COMMENT(AE, "The hybrid form of the serendipity hexahedron is chosen because");
COMMENT(AE, "the material is  nearly incompressible.");
ELEMENT(AE, "c3d20rh", "AllElements", 1, fes.conn)
ORIENTATION(AE, "GlobalOrientation", vec([1. 0 0]), vec([0 1. 0]));
SOLID_SECTION(AE, "elasticity", "GlobalOrientation", "AllElements");
END_INSTANCE(AE);
END_ASSEMBLY(AE);
MATERIAL(AE, "elasticity")
ELASTIC(AE, E, nu)
DENSITY(AE, rho)
STEP_FREQUENCY(AE, neigvs)
END_STEP(AE)
close(AE)

true
