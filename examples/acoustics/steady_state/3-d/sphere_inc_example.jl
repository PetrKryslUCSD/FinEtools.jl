using FinEtools


rho = 1.21*phun("kg/m^3");# mass density
c  = 343.0*phun("m/s");# sound speed
bulk =  c^2*rho;
omega = 2000*phun("rev/s");      # frequency of the incident wave
Ro = 300.0*phun("mm"); # radius of the enclosure
nPerradius = 18;#number of elements on the radius of the sphere
tolerance = Ro/(2^nPerradius)/100
pincampl = 1.0*phun("kilo*Pa");
k = omega/c*[1.0 0.0 0.0];

println("""
        Spherical domain filled with acoustic medium, no scatterer.
        Incident wave.
Hexahedral mesh.
        """)

t0  =  time()

# Hexahedral mesh
fens,fes  =  H8spheren(Ro, nPerradius);
fens1,fes1 = mirrormesh(fens, fes, [-1.0, 0.0, 0.0], [0.0, 0.0, 0.0];
  renumb =  r(c) = c[[1, 4, 3, 2, 5, 8, 7, 6]])
fens,newfes1,fes2 =  mergemeshes(fens1, fes1, fens, fes, tolerance)
fes = cat(newfes1,fes2)

# Derive the finite element sets for the boundary
bfes  =  meshboundary(fes)
# Outer spherical boundary
function dout(xyz)
    return xyz/norm(xyz)
end

louter = selectelem(fens, bfes, facing = true, direction = dout)
outer_fes = subset(bfes, louter);

# File  =   "Sphere.vtk"
# vtkexportmesh (File, fes.conn, fens.xyz, MeshExportModule.H8)
# @async run(`"C:/Program Files (x86)/ParaView 4.2.0/bin/paraview.exe" $File`)

println("Pre-processing time elapsed  =  ",time() - t0,"s")

t1  =  time()

geom  =  NodalField(fens.xyz)
P  =  NodalField(zeros(FCplxFlt,size(fens.xyz,1),1))

numberdofs!(P)

pinc = deepcopy(P);

material = MatAcoustFluid(bulk,rho)
femm  =  FEMMAcoust(IntegData(fes, GaussRule(3, 2)), material)

@time S  =  acousticstiffness(femm, geom, P);
@time C  =  acousticmass(femm, geom, P);

abcfemm  =  FEMMAcoustSurf(IntegData(outer_fes, GaussRule(2, 2)), material)
@time D  =  acousticABC(abcfemm, geom, P);

pincf(loc)  =  pincampl*exp((-1.0im*(vec(k)'*vec(loc)))[1])

# Incident pressure loading
for j = 1:size(geom.values,1)
    pinc.values[j] = pincampl[1,1]*exp((-1.0im*(vec(k)'*vec(geom.values[j,:])))[1])
end
vpinc = gathersysvec(pinc)

F  =  0.0 * vpinc;                # Start with a zero vector

#pfemm = FEMMBase(inner_fes, GaussRule(order = 2,dim = 2))
function dpincdn(dpdn, xyz, J, label)
    n = cross(J[:,1],J[:,2]);
    n = vec(n/norm(n));
    #println("$( (-1.0im)*dot(vec(k),n)*exp(-1.0im*dot(vec(k),vec(xyz))) )")
    dpdn[1] =  pincampl*(-1.0im)*dot(vec(k),n)*exp(-1.0im*dot(vec(k),vec(xyz)))
end

fi  =  ForceIntensity(FCplxFlt, 1, dpincdn);
#@time F  =  F - distribloads(pfemm, nothing, geom, P, fi, 2);
@time F  =  F + distribloads(abcfemm, geom, P, fi, 2);

@time K = lufact((1.0+0.0im)*(-omega^2*S + C)) # We fake a complex matrix here
@time p = K\F  #+omega*1.0im*D

scattersysvec!(P,p[:])

println("Computing time elapsed  =  ",time() - t1,"s")
println("Total time elapsed  =  ",time() - t0,"s")

ptot = deepcopy(P)
ptot.values = P.values+pinc.values

# File  =   "Sphereptot.vtk"
# vtkexportmesh (File, fes.conn, geom.values, MeshExportModule.H8; scalars = abs(ptot.values), scalars_name  = "ptot")
# @async run(`"C:/Program Files (x86)/ParaView 4.2.0/bin/paraview.exe" $File`)

File  =   "Spherepinc.vtk"
vtkexportmesh(File, connasarray(fes), geom.values, FinEtools.MeshExportModule.H8;
scalars = [("realpinc", real(pinc.values))])
@async run(`"paraview.exe" $File`)

File  =   "SphereP.vtk"
vtkexportmesh(File, connasarray(fes), geom.values, FinEtools.MeshExportModule.H8;
 scalars = [( "realP", real(P.values))])
@async run(`"paraview.exe" $File`)

true
