using FinEtools

println("The interior sphere accelerates in the positive x-direction, generating
positive pressure ahead of it, negative pressure behind
")
rho = 1.21*phun("kg/m^3");# mass density
c  = 343.0*phun("m/s");# sound speed
bulk =  c^2*rho;
a_amplitude=1.*phun("mm/s^2");# amplitude of the  acceleration of the sphere
# omega = 2000*phun("rev/s");      # frequency of the incident wave
R = 5.0*phun("mm"); # radius of the interior sphere
Ro = 4*R # radius of the external sphere
P_amplitude = R*rho*a_amplitude; # pressure amplitude
nref=2;
nlayers=40;
tolerance = Ro/(nlayers)/100
frequency=2000; # Hz
omega=2*pi*frequency;
k = omega/c*[1.0 0.0 0.0];

println("""
        Spherical domain filled with acoustic medium, no scatterer.
Hexahedral mesh.
        """)

t0  =  time()

# Hexahedral mesh
fens,fes  =  H8sphere(R,nref);
fens1,fes1  =  mirrormesh(fens, fes, [-1.0, 0.0, 0.0], [0.0, 0.0, 0.0],
renumb =  r(c) = c[[1, 4, 3, 2, 5, 8, 7, 6]]);
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

fens,fes  = H8extrudeQ4(fens, outer_fes,
  nlayers, (xyz, layer)->(R+layer/nlayers*(Ro-R))*xyz/norm(xyz));
connected = findunconnnodes(fens, fes);
fens, new_numbering = compactnodes(fens, connected);
fess = renumberconn!(fes, new_numbering);

geom  =  NodalField(fens.xyz)
P  =  NodalField(zeros(FCplxFlt,size(fens.xyz,1),1))

bfes  =  meshboundary(fes)

# File  =   "Sphere.vtk"
# vtkexportmesh(File, bfes.conn, geom.values, FinEtools.MeshExportModule.Q4)
# @async run(`"paraview.exe" $File`)

# File  =   "Sphere.vtk"
# vtkexportmesh (File, fes.conn, fens.xyz, MeshExportModule.H8)
# @async run(`"C:/Program Files (x86)/ParaView 4.2.0/bin/paraview.exe" $File`)

linner = selectelem(fens, bfes, distance = R, from = [0.0 0.0 0.0],
  inflate = tolerance)
louter = selectelem(fens, bfes, facing = true, direction = dout)

println("Pre-processing time elapsed  =  ",time() - t0,"s")

t1  =  time()


numberdofs!(P)

material = MatAcoustFluid(bulk,rho)
femm  =  FEMMAcoust(IntegData(fes, GaussRule(3, 2)), material)

@time S  =  acousticstiffness(femm, geom, P);
@time C  =  acousticmass(femm, geom, P);

abcfemm  =  FEMMAcoustSurf(IntegData(subset(bfes, louter), GaussRule(2, 2)), material)
@time D  =  acousticABC(abcfemm, geom, P);

# Inner sphere pressure loading
function dipole(dpdn, xyz, J, label)
    n = cross(J[:,1],J[:,2]);
    n = vec(n/norm(n));
    #println("$( (-1.0im)*dot(vec(k),n)*exp(-1.0im*dot(vec(k),vec(xyz))) )")
    dpdn[1] = -rho*a_amplitude*n[1]
end

fi  =  ForceIntensity(FCplxFlt, 1, dipole);
dipfemm  =  FEMMAcoustSurf(IntegData(subset(bfes, linner), GaussRule(2, 2)), material)
@time F  = distribloads(dipfemm, geom, P, fi, 2);

@time K = lufact((1.0+0.0im)*(-omega^2*S + omega*1.0im*D + C)) # We fake a complex matrix here
@time p = K\F  #

scattersysvec!(P, p[:])

println(" Minimum/maximum pressure= $(minimum(real(p)))/$(maximum(real(p)))")

println("Computing time elapsed  =  ",time() - t1,"s")
println("Total time elapsed  =  ",time() - t0,"s")

File  =   "sphere_dipole_1.vtk"
vtkexportmesh(File, connasarray(fes), geom.values, FinEtools.MeshExportModule.H8;
 scalars = [( "realP", real(P.values))])
@async run(`"paraview.exe" $File`)

true
