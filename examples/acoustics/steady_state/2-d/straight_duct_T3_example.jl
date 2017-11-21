using FinEtools

t0  =  time()

rho = 1.21*phun("kg/m^3");# mass density
c  = 343.0*phun("m/s");# sound speed
bulk =  c^2*rho;
omega =  54.5901*phun("rev/s")
vn0 =  -1.0*phun("m/s")
Lx = 10.0*phun("m");# length of the box, millimeters
Ly = 1.0*phun("m"); # length of the box, millimeters
n = 20;#number of elements along the length

println("""

Straight duct with anechoic termination.
Example from Boundary element acoustics: Fundamentals and computer codes, TW Wu, page 44.
Both real and imaginary components of the pressure should have amplitude of
rho*c = $(rho*c).

Triangle mesh.
""")

fens,fes  =  T3block(Lx,Ly,n,2); # Mesh
bfes  =  meshboundary(fes)
L0 = selectelem(fens,bfes,facing = true, direction = [-1.0 0.0])
L10 = selectelem(fens,bfes,facing = true, direction = [+1.0 0.0])
nLx = selectnode(fens,box = [0.0 Lx  0.0 0.0], inflate = Lx/1.0e5)

geom  =  NodalField(fens.xyz)
P  =  NodalField(zeros(Complex128,size(fens.xyz,1),1))

numberdofs!(P)


material = MatAcoustFluid(bulk,rho)
femm  =  FEMMAcoust(IntegData(fes, TriRule(1)), material)

S  =  acousticstiffness(femm, geom, P);
C  =  acousticmass(femm, geom, P);


E10femm  =  FEMMAcoustSurf(IntegData(subset(bfes,L10),GaussRule(1, 2)), material)
D  =  acousticABC(E10femm, geom, P);

E0femm  =  FEMMBase(IntegData(subset(bfes,L0), GaussRule(1,  2)))
fi  =  ForceIntensity(-1.0im*omega*rho*vn0);
F  =  distribloads(E0femm, geom, P, fi, 2);

p = (-omega^2*S +omega*1.0im*D + C)\F
scattersysvec!(P, p[:])

println("Pressure amplitude bounds: ")
println("  real $(minimum(real(P.values)))/$(maximum(real(P.values)))")
println("  imag $(minimum(imag(P.values)))/$(maximum(imag(P.values)))")

println("Total time elapsed  =  ",time() - t0,"s")

File  =   "straight_duct.vtk"
scalars = real(P.values);
vtkexportmesh(File, fes.conn, geom.values, FinEtools.MeshExportModule.T3;
scalars = [("Pressure", scalars)])
@async run(`"paraview.exe" $File`)


using PyCall
@pyimport matplotlib.pyplot as plt
plt.style[:use]("seaborn-whitegrid")
fig = plt.figure() 
ax = plt.axes()
ix = sortperm(geom.values[nLx,1])
ax[:plot](geom.values[nLx,1][ix], real(P.values)[nLx][ix], marker=:o, color = :blue, label = "real")
ax[:plot](geom.values[nLx,1][ix], imag(P.values)[nLx][ix], marker=:d,  color = :red, label  =  "imag")
plt.legend()
plt.show()

true
