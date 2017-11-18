using FinEtools

println("""
Example from Sound and Structural Vibration, Second Edition: Radiation, Transmission and Response [Paperback]
Frank J. Fahy, Paolo Gardonio, page 483.

Hexahedral mesh.
""")

t0 = time()

rho=1.21*1e-9;# mass density
c =343.0*1000;# millimeters per second
bulk= c^2*rho;
L=500.0;# length of the box, millimeters
A=200.0; # cross-sectional area of the box
n=40;#
neigvs=8;
OmegaShift=10.0;

fens,fes = H8block(L,sqrt(A),sqrt(A),n,1,1); # Mesh
fens,fes = H8toH20(fens,fes)

geom = NodalField(fens.xyz)
P = NodalField(zeros(size(fens.xyz,1),1))

numberdofs!(P)

femm = FEMMAcoust(IntegD(fes, GaussRule(3, 3)), MatAcoustFluid(bulk, rho))


S = acousticstiffness(femm, geom, P);
C = acousticmass(femm, geom, P);

d,v,nev,nconv = eigs(C+OmegaShift*S, S; nev=neigvs, which=:SM)
d = d - OmegaShift;
fs = real(sqrt.(complex(d)))/(2*pi)
println("Eigenvalues: $fs [Hz]")


println("Total time elapsed = ",time() - t0,"s")

File =  "fahy_H20.vtk"
en = 5;
vtkexportmesh(File, fes.conn, geom.values, FinEtools.MeshExportModule.H20;
          scalars=[("Pressure_mode_$en", v[:,en])])
@async run(`"paraview.exe" $File`)
println("Done")
true
