using FinEtools
using FinEtools.MeshExportModule

rho = 1.21*phun("kg/m^3");# mass density
c  = 343.0*phun("m/s");# sound speed
bulk =  c^2*rho;
R = 50.0*phun("mm");# radius of the piston
Ro = 300.0*phun("mm"); # radius of the enclosure
nPerradius=24;# number of elements along the radius of the scatterer
nlayers=78;                     # number of layers of elements surrounding the piston
tolerance = R/(nPerradius)/100
kR =  1.0;
k = kR/R*[1.0 0.0 0.0]; k = vec(k);
omega = norm(k)*c;    # angular frequency of the incident wave
pincampl = 1.0*phun("Pa");
h = maximum([(Ro-R)/nlayers pi*Ro/(nPerradius)]);

println("""

        Rigid sphere scattering.

        Frequency: $(omega/2/pi). Wave vector: $k. kR = $(norm(k)*R).
        Number of elements per wavelength: $(2*pi/(norm(k)*h)).

        Hexahedral mesh.
        """)

t0  =  time()

# Hexahedral mesh
fens, fes  =  H8spheren(R, nPerradius);
bfes  =  meshboundary(fes)
l = selectelem(fens, bfes, facing = true, direction = [1.0 1.0  1.0], dotmin= 0.001)
ex(xyz, layer) = (R+layer/nlayers*(Ro-R))*xyz/norm(xyz)
fens,fes  =  H8extrudeQ4(fens, subset(bfes,l), nlayers, ex);
fens1,fes1 = mirrormesh(fens, fes, [-1.0, 0.0, 0.0], [0.0, 0.0, 0.0]; renumb =  r(c) = c[[1, 4, 3, 2, 5, 8, 7, 6]])
fens,newfes1,fes2 =  mergemeshes(fens1, fes1, fens, fes, tolerance)
fes = cat(newfes1,fes2)

# Derive the finite element sets for the boundary
bfes  =  meshboundary(fes)
# Outer spherical boundary
function dout(xyz)
    return xyz/norm(xyz)
end

louter = selectelem(fens, bfes, facing = true, direction = dout)
outer_fes = subset(bfes,louter);

# Inner spherical boundary
function din(xyz)
    return -xyz/norm(xyz)
end
linner = selectelem(fens, bfes, facing = true, direction = din)
inner_fes = subset(bfes,linner);

# File  =   "Sphere.vtk"
# vtkexportmesh(File, inner_fes.conn, fens.xyz, MeshExportModule.Q4)
# @async run(`"paraview.exe" $File`)
println("Pre-processing time elapsed  =  ",time() - t0,"s")

t1  =  time()

geom  =  NodalField(fens.xyz)
P  =  NodalField(zeros(FCplxFlt,size(fens.xyz,1),1))

numberdofs!(P)

pinc = deepcopy(P);

material = MatAcoustFluid(bulk,rho)
femm  =  FEMMAcoust(IntegData(fes, GaussRule(3, 3)), material)

S  =  acousticstiffness(femm, geom, P);
C  =  acousticmass(femm, geom, P);

abcfemm  =  FEMMAcoustSurf(IntegData(outer_fes, GaussRule(2, 3)), material)
D  =  acousticABC(abcfemm, geom, P);

# Incident pressure loading
for j = 1:size(geom.values,1)
    xyz = vec(geom.values[j,:]);
    pinc.values[j] = pincampl*exp((-1.0im*(k'*xyz))[1])
end
vpinc = gathersysvec(pinc)

F  =  - (-omega^2*S + C) * vpinc;

#pfemm = FEMMBase(inner_fes, GaussRule(order = 2,dim = 2))
function dpincdn(dpdn, xyz, J, label)
    xyz = vec(xyz);
    n = cross(J[:,1],J[:,2]);
    n = vec(n/norm(n));
    dpdn[1] = pincampl*(-1.0im)*dot(vec(k),n)*exp(-1.0im*dot(vec(k),vec(xyz)))
end

fi  =  ForceIntensity(FCplxFlt, 1, dpincdn);
F  =  F + distribloads(abcfemm, geom, P, fi, 2);

K = ((-omega^2*S + omega*1.0im*D + C))
p = K\F

scattersysvec!(P,p[:])

println("Computing time elapsed  =  ",time() - t1,"s")
println("Total time elapsed  =  ",time() - t0,"s")

ptot = deepcopy(P)
ptot.values = P.values+pinc.values

# File  =   "Sphereptot.vtk"
# vtkexportmesh (File, fes.conn, geom.values, MeshExportModule.H8; scalars = abs(ptot.values), scalars_name  = "ptot")
# @async run(`"C:/Program Files (x86)/ParaView 4.2.0/bin/paraview.exe" $File`)

# File  =   "Spherepinc.vtk"
# vtkexportmesh (File, fes.conn, geom.values, MeshExportModule.H8; scalars = real(pinc.values), scalars_name  = "pinc")
# @async run(`"C:/Program Files (x86)/ParaView 4.2.0/bin/paraview.exe" $File`)

File  =   "SphereP.vtk"
vtkexportmesh(File, connasarray(fes), geom.values, FinEtools.MeshExportModule.H8; scalars = [("P", abs.(P.values))])
@async run(`"paraview.exe" $File`)

# File  =   "SphereP.vtk"
# vtkexportmesh (File, cat(outer_fes,inner_fes).conn, geom.values, MeshExportModule.Q4; scalars = abs(P.values), scalars_name  = "P")
# @async run(`"C:/Program Files (x86)/ParaView 4.2.0/bin/paraview.exe" $File`)

# using Winston
# pl  =  FramedPlot(title = "Matrix",xlabel = "x",ylabel = "Re P, Im P")
# setattr(pl.frame, draw_grid = true)
# add(pl, Curve([1:length(C[:])],vec(C[:]), color = "blue"))

# # pl = plot(geom.values[nLx,1][ix],scalars[nLx][ix])
# # xlabel("x")
# # ylabel("Pressure")
# display(pl)

true
