using FinEtools

println("Cook membrane problem,  plane stress."        )
t0 = time()

E = 1.0;
nu = 1.0/3;
width = 48.0; height = 44.0; thickness  = 1.0;
free_height  = 16.0;
Mid_edge  = [48.0, 52.0];# Location of tracked  deflection
magn = 1.0/free_height;# Magnitude of applied load
convutip = 23.97;
n = 32;#*int(round(sqrt(170.)/2.)); # number of elements per side
tolerance = minimum([width, height])/n/1000.;#Geometrical tolerance

fens, fes = T3block(width, height,  n,  n)

# Reshape into a trapezoidal panel
for i=1:count(fens)
  fens.xyz[i, 2]=fens.xyz[i, 2]+(fens.xyz[i, 1]/width)*(height -fens.xyz[i, 2]/height*(height-free_height));
end

geom = NodalField(fens.xyz)
u = NodalField(zeros(size(fens.xyz, 1), 2)) # displacement field

l1 = selectnode(fens; box=[0, 0, -Inf,  Inf],  inflate = tolerance)
setebc!(u, l1, 1, val=0.0)
setebc!(u, l1, 2, val=0.0)
applyebc!(u)
numberdofs!(u)

@time boundaryfes =  meshboundary(fes);
Toplist = selectelem(fens, boundaryfes,  box= [width,  width,  -Inf,  Inf ],  inflate=  tolerance);
el1femm =  FEMMBase(IntegData(subset(boundaryfes, Toplist),  GaussRule(1, 2)))
fi = ForceIntensity([0.0, +magn]);
F2 = distribloads(el1femm,  geom,  u,  fi,  2);


MR = DeforModelRed2DStress
material = MatDeforElastIso(MR,  0.0, E, nu, 0.0)

femm = FEMMDeforLinear(MR, IntegData(fes,  TriRule(1)),  material)

K = stiffness(femm,  geom,  u)
K = cholfact(K)
U=  K\(F2)
scattersysvec!(u, U[:])

nl = selectnode(fens,  box=[Mid_edge[1], Mid_edge[1], Mid_edge[2], Mid_edge[2]], inflate=tolerance);
theutip = zeros(FFlt, 1, 2)
gathervalues_asmat!(u, theutip, nl);
println("$(time()-t0) [s];  displacement =$(theutip[2]) as compared to converged $convutip")

using FinEtools.MeshExportModule

File =  "a.vtk"
vtkexportmesh(File,  fes.conn,  geom.values+u.values,
               FinEtools.MeshExportModule.T3; vectors=[("u", u.values)])

true
