
module mmmmmmmmmmrrigid
using FinEtools
using Base.Test
function test()

# println("""
#
# Example from Boundary element acoustics: Fundamentals and computer codes, TW Wu, page 123.
# Internal resonance problem. Reference frequencies: 90.7895, 181.579, 215.625, 233.959, 272.368, 281.895
# Quadrilateral mesh.
# """)

t0 = time()

rho=1.21*1e-9;# mass density
c =345.0*1000;# millimeters per second
bulk= c^2*rho;
Lx=1900.0;# length of the box, millimeters
Ly=800.0; # length of the box, millimeters
n=14;#
neigvs=18;
OmegaShift=10.0;

fens,fes = Q4block(Lx,Ly,n,n); # Mesh

geom = NodalField(fens.xyz)
P = NodalField(zeros(size(fens.xyz,1),1))

numberdofs!(P)

femm = FEMMAcoust(GeoD(fes, GaussRule(2, 2)),
                     MatAcoustFluid(bulk,rho))

S = acousticstiffness(femm, geom, P);
C = acousticmass(femm, geom, P);

d,v,nev,nconv =eigs(C+OmegaShift*S, S; nev=neigvs, which=:SM)
d = d - OmegaShift;
fs=real(sqrt.(complex(d)))/(2*pi)
# println("Eigenvalues: $fs [Hz]")
#
#
# println("Total time elapsed = ",time() - t0,"s")
# println("$fs[2]-90.7895")
@test abs(fs[2]-90.9801) < 1e-4

File =  "rigid_box.vtk"
scalarllist = Any[]
for n  = 1:15
  push!(scalarllist, ("Pressure_mode_$n", v[:,n]));
end
vtkexportmesh(File, fes.conn, geom.values, FinEtools.MeshExportModule.Q4;
          scalars=scalarllist)
sleep(1.0)
rm(File)
# @async run(`"paraview.exe" $File`)

true
end
end
using mmmmmmmmmmrrigid
mmmmmmmmmmrrigid.test()

module fahyL2example

using FinEtools
using Base.Test

function test()
# println("""
# Example from Sound and Structural Vibration, Second Edition: Radiation, Transmission and Response [Paperback]
# Frank J. Fahy, Paolo Gardonio, page 483.
#
# 1D mesh.
# """)

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
# println("Eigenvalues: $fs [Hz]")


# println("Total time elapsed  =  ",time() - t0,"s")

# using Plots
# plotly()
# en = 2
# ix = sortperm(geom.values[:])
# plot(geom.values[:][ix], v[:,en][ix], color = "blue",
# title = "Fahy example, mode $en" , xlabel = "x", ylabel = "P")
# gui()

  @test (fs[2]-343.08)/343.08 < 0.001

  #println("Total time elapsed = ",time() - t0,"s")

  # using Winston
  # en=2
  # pl = FramedPlot(title="Fahy example, mode $en",xlabel="x",ylabel="P")
  # setattr(pl.frame, draw_grid=true)
  # ix=sortperm(geom.values[:])
  # add(pl, Curve(geom.values[:][ix],v[:,en][ix], color="blue"))

  # display(pl)

  true
end
end
using fahyL2example
fahyL2example.test()

module mmfahyH8example
using FinEtools
using Base.Test
function test()

# println("""
# Example from Sound and Structural Vibration, Second Edition: Radiation, Transmission and Response [Paperback]
# Frank J. Fahy, Paolo Gardonio, page 483.
#
# Hexahedral mesh.
# """)

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

geom = NodalField(fens.xyz)
P = NodalField(zeros(size(fens.xyz,1),1))

numberdofs!(P)

femm = FEMMAcoust(GeoD(fes, GaussRule(3, 2)), MatAcoustFluid(bulk, rho))


S = acousticstiffness(femm, geom, P);
C = acousticmass(femm, geom, P);

d,v,nev,nconv = eigs(C+OmegaShift*S, S; nev=neigvs, which=:SM)
d = d - OmegaShift;
fs = real(sqrt.(complex(d)))/(2*pi)
# println("Eigenvalues: $fs [Hz]")
#
#
# println("Total time elapsed = ",time() - t0,"s")

# File =  "fahy_H8.vtk"
# en = 5;
# vtkexportmesh(File, fes.conn, geom.values, FinEtools.MeshExportModule.H8;
#           scalars=[("Pressure_mode_$en", v[:,en])])
# println("Done")

  @test (fs[2]-343.08)/343.08 < 0.001

  true
end
end
using mmfahyH8example
mmfahyH8example.test()

module mmmmfahyH27example
using FinEtools
using Base.Test
function test()

# println("""
# Example from Sound and Structural Vibration, Second Edition: Radiation, Transmission and Response [Paperback]
# Frank J. Fahy, Paolo Gardonio, page 483.
#
# Hexahedral mesh.
# """)

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
fens,fes = H8toH27(fens,fes)

geom = NodalField(fens.xyz)
P = NodalField(zeros(size(fens.xyz,1),1))

numberdofs!(P)

femm = FEMMAcoust(GeoD(fes, GaussRule(3, 3)), MatAcoustFluid(bulk, rho))


S = acousticstiffness(femm, geom, P);
C = acousticmass(femm, geom, P);

d,v,nev,nconv = eigs(C+OmegaShift*S, S; nev=neigvs, which=:SM)
d = d - OmegaShift;
fs = real(sqrt.(complex(d)))/(2*pi)
# println("Eigenvalues: $fs [Hz]")
#
#
# println("Total time elapsed = ",time() - t0,"s")

# File =  "fahy_H8.vtk"
# en = 5;
# vtkexportmesh(File, fes.conn, geom.values, FinEtools.MeshExportModule.H8;
#           scalars=[("Pressure_mode_$en", v[:,en])])
# println("Done")

  @test (fs[2]-343.08)/343.08 < 0.001

  true
end
end
using mmmmfahyH27example
mmmmfahyH27example.test()

module mmmmmstraight_duct_H8_1
using FinEtools
using Base.Test
function test()
    t0  =  time()

    rho = 1.21*phun("kg/m^3");# mass density
    c  = 343.0*phun("m/s");# sound speed
    bulk =  c^2*rho;
    omega =  54.5901*phun("rev/s")
    vn0 =  -1.0*phun("m/s")
    Lx = 10.0*phun("m");# length of the box, millimeters
    Ly = 1.0*phun("m"); # length of the box, millimeters
    n = 20;#number of elements along the length

    # println("""
    #
    # Straight duct with anechoic termination.
    # Example from Boundary element acoustics: Fundamentals and computer codes, TW Wu, page 44.
    # Both real and imaginary components of the pressure should have amplitude of
    # rho*c = $(rho*c).
    #
    # Hexahedral mesh.
    # """)

    fens,fes  =  H8block(Lx,Ly,Ly,n,2,2); # Mesh
    bfes  =  meshboundary(fes)
    L0 = selectelem(fens,bfes,facing = true, direction = [-1.0 0.0 0.0])
    L10 = selectelem(fens,bfes,facing = true, direction = [+1.0 0.0 0.0])
    nLx = selectnode(fens,box = [0.0 Lx  0.0 0.0 0.0 0.0], inflate = Lx/1.0e5)

    geom  =  NodalField(fens.xyz)
    P  =  NodalField(zeros(Complex128,size(fens.xyz,1),1))

    numberdofs!(P)


    material = MatAcoustFluid(bulk,rho)
    femm  =  FEMMAcoust(GeoD(fes, GaussRule(3, 2)), material)

    S  =  acousticstiffness(femm, geom, P);
    C  =  acousticmass(femm, geom, P);


    E10femm  =  FEMMAcoustSurf(GeoD(subset(bfes,L10),GaussRule(2, 2)), material)
    D  =  acousticABC(E10femm, geom, P);

    E0femm  =  FEMMBase(GeoD(subset(bfes,L0), GaussRule(2,  2)))
    fi  =  ForceIntensity(-1.0im*omega*rho*vn0);
    F  =  distribloads(E0femm, geom, P, fi, 2);

    p = (-omega^2*S +omega*1.0im*D + C)\F
    scattersysvec!(P, p[:])

    # println("Pressure amplitude bounds: ")
    # println("  real $(minimum(real(P.values)))/$(maximum(real(P.values)))")
    # println("  imag $(minimum(imag(P.values)))/$(maximum(imag(P.values)))")
    #
    # println("Total time elapsed  =  ",time() - t0,"s")

    File  =   "straight_duct.vtk"
    scalars = real(P.values);
    vtkexportmesh(File, fes.conn, geom.values, FinEtools.MeshExportModule.H8;
    scalars = [("Pressure", scalars)])
    rm(File)
    # @async run(`"paraview.exe" $File`)

    ref=rho*c
    @test abs(minimum(real(P.values))-(-ref))/ref < 0.02
    @test abs(minimum(imag(P.values))-(-ref))/ref < 0.02
    @test abs(maximum(real(P.values))-(+ref))/ref < 0.02
    @test abs(maximum(imag(P.values))-(+ref))/ref < 0.02

    # plotly()
    # ix = sortperm(geom.values[nLx,1])
    # plot(geom.values[nLx,1][ix], real(P.values)[nLx][ix], color = :blue, label = "real")
    # plot!(geom.values[nLx,1][ix], imag(P.values)[nLx][ix], color = :red, label  =  "imag")
    # plot!(title = "Straight duct with anechoic termination",
    # xlabel = "x", ylabel = "Pressure")
    # gui()

    true

end
end
using mmmmmstraight_duct_H8_1
mmmmmstraight_duct_H8_1.test()

module mmmmmmmmmmsphere_dipole_1
using FinEtools
using Base.Test
function test()

# println("The interior sphere accelerates in the positive x-direction, generating
# positive pressure ahead of it, negative pressure behind
# ")
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

# println("""
#         Spherical domain filled with acoustic medium, no scatterer.
# Hexahedral mesh.
#         """)

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

# println("Pre-processing time elapsed  =  ",time() - t0,"s")

t1  =  time()


numberdofs!(P)

material = MatAcoustFluid(bulk,rho)
femm  =  FEMMAcoust(GeoD(fes, GaussRule(3, 2)), material)

S  =  acousticstiffness(femm, geom, P);
C  =  acousticmass(femm, geom, P);

abcfemm  =  FEMMAcoustSurf(GeoD(subset(bfes, louter), GaussRule(2, 2)), material)
D  =  acousticABC(abcfemm, geom, P);

# Inner sphere pressure loading
function dipole(dpdn, xyz, J, label)
    n = cross(J[:,1],J[:,2]);
    n = vec(n/norm(n));
    #println("$( (-1.0im)*dot(vec(k),n)*exp(-1.0im*dot(vec(k),vec(xyz))) )")
    dpdn[1] = -rho*a_amplitude*n[1]
end

fi  =  ForceIntensity(FCplxFlt, 1, dipole);
#F  =  F - distribloads(pfemm, nothing, geom, P, fi, 2);
dipfemm  =  FEMMAcoustSurf(GeoD(subset(bfes, linner), GaussRule(2, 2)), material)
F  = distribloads(dipfemm, geom, P, fi, 2);

K = lufact((1.0+0.0im)*(-omega^2*S + omega*1.0im*D + C)) # We fake a complex matrix here
p = K\F  #

scattersysvec!(P, p[:])

# println(" Minimum/maximum pressure= $(minimum(real(p)))/$(maximum(real(p)))")
@test abs(-minimum(real(p))-maximum(real(p))) < 1.0e-9
@test abs(3.110989923540772e-6-maximum(real(p))) < 1.0e-9
# println("Computing time elapsed  =  ",time() - t1,"s")
# println("Total time elapsed  =  ",time() - t0,"s")

File  =   "sphere_dipole_1.vtk"
vtkexportmesh(File, fes.conn, geom.values, FinEtools.MeshExportModule.H8;
 scalars = [( "realP", real(P.values))])
rm(File)
# @async run(`"paraview.exe" $File`)

true
end
end
using mmmmmmmmmmsphere_dipole_1
mmmmmmmmmmsphere_dipole_1.test()

module mmmmmmmstraight_duct_T10_examplemmm
using FinEtools
using Base.Test
function test()

t0  =  time()

rho = 1.21*phun("kg/m^3");# mass density
c  = 343.0*phun("m/s");# sound speed
bulk =  c^2*rho;
omega =  54.5901*phun("rev/s")
vn0 =  -1.0*phun("m/s")
Lx = 10.0*phun("m");# length of the box, millimeters
Ly = 1.0*phun("m"); # length of the box, millimeters
n = 20;#number of elements along the length

# println("""
#
# Straight duct with anechoic termination.
# Example from Boundary element acoustics: Fundamentals and computer codes, TW Wu, page 44.
# Both real and imaginary components of the pressure should have amplitude of
# rho*c = $(rho*c).
#
# Tetrahedral (quadratic) mesh.
# """)

fens,fes  =  T10block(Lx,Ly,Ly,n,2,2); # Mesh
bfes  =  meshboundary(fes)
L0 = selectelem(fens,bfes,facing = true, direction = [-1.0 0.0 0.0])
L10 = selectelem(fens,bfes,facing = true, direction = [+1.0 0.0 0.0])
nLx = selectnode(fens,box = [0.0 Lx  0.0 0.0 0.0 0.0], inflate = Lx/1.0e5)

geom  =  NodalField(fens.xyz)
P  =  NodalField(zeros(Complex128,size(fens.xyz,1),1))

numberdofs!(P)


material = MatAcoustFluid(bulk,rho)
femm  =  FEMMAcoust(GeoD(fes, TetRule(4)), material)

S  =  acousticstiffness(femm, geom, P);
C  =  acousticmass(femm, geom, P);


E10femm  =  FEMMAcoustSurf(GeoD(subset(bfes,L10), TriRule(3)), material)
D  =  acousticABC(E10femm, geom, P);

E0femm  =  FEMMBase(GeoD(subset(bfes,L0), TriRule(3)))
fi  =  ForceIntensity(-1.0im*omega*rho*vn0);
F  =  distribloads(E0femm, geom, P, fi, 2);

p = (-omega^2*S +omega*1.0im*D + C)\F
scattersysvec!(P, p[:])

# println("Pressure amplitude bounds: ")
# println("  real $(minimum(real(P.values)))/$(maximum(real(P.values)))")
# println("  imag $(minimum(imag(P.values)))/$(maximum(imag(P.values)))")
#
# println("Total time elapsed  =  ",time() - t0,"s")
@test abs(-414.0986190028914-minimum(imag(P.values))) < 1.e-6

File  =   "straight_duct.vtk"
scalars = real(P.values);
vtkexportmesh(File, fes.conn, geom.values, FinEtools.MeshExportModule.T10;
scalars = [("Pressure", scalars)])
# @async run(`"paraview.exe" $File`)
rm(File)

# plotly()
# ix = sortperm(geom.values[nLx,1])
# plot(geom.values[nLx,1][ix], real(P.values)[nLx][ix], color = :blue, label = "real")
# plot!(geom.values[nLx,1][ix], imag(P.values)[nLx][ix], color = :red, label  =  "imag")
# plot!(title = "Straight duct with anechoic termination",
# xlabel = "x", ylabel = "Pressure")
# gui()

true
end
end
using mmmmmmmstraight_duct_T10_examplemmm
mmmmmmmstraight_duct_T10_examplemmm.test()

module mmmmiintegrationmmmm
using FinEtools
using Base.Test
function test()
rho = 1000.*phun("kg/m^3");# mass density of water
c  = 1.4491e+3*phun("m/s");# sound speed in water
bulk =  c^2*rho;
rhos = 2500.*phun("kg/m^3");# mass density of the solid sphere
a_amplitude=1.*phun("mm/s^2");# amplitude of the  acceleration of the sphere
R = 5.0*phun("mm"); # radius of the interior sphere
Ro = 3*R # radius of the external sphere
# Mesh parameters, geometrical tolerance
nperR = 6;
nlayers=10;
tolerance = Ro/(nlayers)/100
# Hexahedral mesh of the solid sphere
fens,fes  =  H8spheren(R, nperR);
r(c) = c[[1, 4, 3, 2, 5, 8, 7, 6]]
fens1,fes1  =  mirrormesh(fens, fes, [-1.0, 0.0, 0.0], [0.0, 0.0, 0.0],
          renumb =  r);
fens,newfes1,fes2 =  mergemeshes(fens1, fes1, fens, fes, tolerance)
fes = cat(newfes1,fes2)
fens1,fes1  =  mirrormesh(fens, fes, [0.0, -1.0, 0.0], [0.0, 0.0, 0.0],
          renumb =  r);
fens,newfes1,fes2 =  mergemeshes(fens1, fes1, fens, fes, tolerance)
fes = cat(newfes1,fes2)
fens1,fes1  =  mirrormesh(fens, fes, [0.0, 0.0, -1.0], [0.0, 0.0, 0.0],
          renumb =  r);
fens,newfes1,fes2 =  mergemeshes(fens1, fes1, fens, fes, tolerance)
fes = cat(newfes1,fes2)

# Find the inertial properties of the solid sphere
geom  =  NodalField(fens.xyz)
femm  =  FEMMBase(GeoD(fes, GaussRule(3, 2)))
V = integratefunction(femm, geom, (x) ->  1.0)
  # println("V=$(V/phun("mm^3"))")
Sx = integratefunction(femm, geom, (x) ->  x[1])
  # println("Sx=$(Sx/phun("mm^4"))")
Sy = integratefunction(femm, geom, (x) ->  x[2])
  # println("Sy=$(Sy/phun("mm^4"))")
Sz = integratefunction(femm, geom, (x) ->  x[3])
  # println("Sz=$(Sz/phun("mm^4"))")
CG = vec([Sx Sy Sz]/V)
  # println("CG=$(CG/phun("mm"))")
function Iinteg(x)
  (norm(x-CG)^2*eye(3)-(x-CG)*(x-CG)')
end
Ixx = integratefunction(femm, geom, (x) ->  Iinteg(x)[1, 1])
Ixy = integratefunction(femm, geom, (x) ->  Iinteg(x)[1, 2])
Ixz = integratefunction(femm, geom, (x) ->  Iinteg(x)[1, 3])
Iyy = integratefunction(femm, geom, (x) ->  Iinteg(x)[2, 2])
Iyz = integratefunction(femm, geom, (x) ->  Iinteg(x)[2, 3])
Izz = integratefunction(femm, geom, (x) ->  Iinteg(x)[3, 3])
I = [Ixx Ixy Ixz; Ixy Iyy Iyz; Ixz Iyz Izz]
  # println("I=$(I/phun("mm^4"))")
@test abs(Ixx/phun("mm^4")-4.9809) < 0.001
mass =  V*rhos;
Inertia = I*rhos;
end
end
using mmmmiintegrationmmmm
mmmmiintegrationmmmm.test()
