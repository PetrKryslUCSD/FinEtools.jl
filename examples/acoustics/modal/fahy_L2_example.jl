using FinEtools
println("""
Example from Sound and Structural Vibration, Second Edition: Radiation, Transmission and Response [Paperback]
Frank J. Fahy, Paolo Gardonio, page 483.

1D mesh.
""")

t0  =  time()

rho = 1.21*1e-9;# mass density
c  = 343.0*1000;# millimeters per second
bulk =  c^2*rho;
L = 500.0;# length of the box, millimeters
A = 200.0; # cross-sectional area of the box
graphics =  true;# plot the solution as it is computed?
n = 40;#
neigvs = 8;
OmegaShift = 10.0;

fens,fes  =  L2block(L,n); # Mesh

geom = NodalField(fens.xyz)
P = NodalField(zeros(size(fens.xyz,1),1))

numberdofs!(P)

femm = FEMMAcoust(GeoD(fes, GaussRule(1, 2)), MatAcoustFluid(bulk, rho))

S  =  acousticstiffness(femm, geom, P);
C  =  acousticmass(femm, geom, P);

d,v,nev,nconv  = eigs(C+OmegaShift*S, S; nev = neigvs, which = :SM)
d  =  d - OmegaShift;
fs = real(sqrt.(complex(d)))/(2*pi)
println("Eigenvalues: $fs [Hz]")


println("Total time elapsed  =  ",time() - t0,"s")

using Plots
plotly()
en = 2
ix = sortperm(geom.values[:])
plot(geom.values[:][ix], v[:,en][ix], color = "blue",
title = "Fahy example, mode $en" , xlabel = "x", ylabel = "P")
gui()
# pl  =  FramedPlot(title =)


true
