using FinEtools
if VERSION >= v"0.7-"
    using SparseArrays
end

println("Rigid movable hemisphere in  water. Time-dependent simulation.
")
rho = 1000.*phun("kg/m^3");# mass density of water
c  = 1.4491e+3*phun("m/s");# sound speed in water
bulk =  c^2*rho;
rhos = 2500.*phun("kg/m^3");# mass density of the solid sphere
a_amplitude=1.*phun("mm/s^2");# amplitude of the  acceleration of the sphere
R = 5.0*phun("mm"); # radius of the interior sphere
Ro = 3*R # radius of the external sphere
P_amplitude = 1000*phun("Pa"); # pressure amplitude
frequency=200.; # Hz
omega=2*pi*frequency;
wavenumber=omega/c;
angl = 0.0
front_normal=vec([cos(pi/180*angl),sin(pi/180*angl),0]);
wavevector=wavenumber*front_normal;
tshift = 0.0
uincA=P_amplitude/rho/c/omega; # amplitude of the displacement
# This is the analytical solution of Hickling and Wang
uamplHW=P_amplitude/(rho*c^2)/(omega/c)*3/(2*rhos/rho+1);
# Mesh parameters, geometrical tolerance
nperR = 8;
nlayers=nperR;
tolerance = Ro/(nlayers)/100
# Time step
dt = 1.0/frequency/20;
tfinal = 9/frequency;
nsteps = round(tfinal/dt)+1;

t0  =  time()

# Hexahedral mesh of the solid sphere
fens,fes  =  H8spheren(R, nperR);
r(c) = c[[1, 4, 3, 2, 5, 8, 7, 6]]
fens1,fes1  =  mirrormesh(fens, fes, [-1.0, 0.0, 0.0], [0.0, 0.0, 0.0],
          renumb = r);
fens,newfes1,fes2 =  mergemeshes(fens1, fes1, fens, fes, tolerance)
fes = cat(newfes1,fes2)
fens1,fes1  =  mirrormesh(fens, fes, [0.0, -1.0, 0.0], [0.0, 0.0, 0.0],
          renumb = r);
fens,newfes1,fes2 =  mergemeshes(fens1, fes1, fens, fes, tolerance)
fes = cat(newfes1,fes2)
fens1,fes1  =  mirrormesh(fens, fes, [0.0, 0.0, -1.0], [0.0, 0.0, 0.0],
          renumb = r);
fens,newfes1,fes2 =  mergemeshes(fens1, fes1, fens, fes, tolerance)
fes = cat(newfes1,fes2)

# Retain only half for the hemisphere
Solidl = selectelem(fens, fes, box=[-Inf,Inf,-Inf,Inf,0,Inf], inflate=(Ro-R)/nlayers/100);
# This is the half sphere which is filled with water
fens1w, fes1w = deepcopy(fens), subset(fes, setdiff(collect(1:count(fes)), Solidl))

# Find the inertial properties of the solid sphere
geom  =  NodalField(fens.xyz)
femm  =  FEMMBase(IntegData(subset(fes, Solidl), GaussRule(3, 2)))
V = integratefunction(femm, geom, (x) ->  1.0)
  println("V=$(V/phun("mm^3"))")
Sx = integratefunction(femm, geom, (x) ->  x[1])
  println("Sx=$(Sx/phun("mm^4"))")
Sy = integratefunction(femm, geom, (x) ->  x[2])
  println("Sy=$(Sy/phun("mm^4"))")
Sz = integratefunction(femm, geom, (x) ->  x[3])
  println("Sz=$(Sz/phun("mm^4"))")
CG = vec([Sx Sy Sz]/V)
  println("CG=$(CG/phun("mm"))")
  I3 = [i==j ? one(FFlt) : zero(FFlt) for i=1:3, j=1:3]
function Iinteg(x)
  (norm(x-CG)^2*I3-(x-CG)*(x-CG)')
end
Ixx = integratefunction(femm, geom, (x) ->  Iinteg(x)[1, 1])
Ixy = integratefunction(femm, geom, (x) ->  Iinteg(x)[1, 2])
Ixz = integratefunction(femm, geom, (x) ->  Iinteg(x)[1, 3])
Iyy = integratefunction(femm, geom, (x) ->  Iinteg(x)[2, 2])
Iyz = integratefunction(femm, geom, (x) ->  Iinteg(x)[2, 3])
Izz = integratefunction(femm, geom, (x) ->  Iinteg(x)[3, 3])
It = [Ixx Ixy Ixz; Ixy Iyy Iyz; Ixz Iyz Izz]
  println("It=$(It/phun("mm^4"))")
mass =  V*rhos;
Inertia = It*rhos;

# The boundary surface of the solid sphere will now be extruded to form  the
# layers of fluid around it
# Outer spherical boundary
function dout(xyz)
    return xyz/norm(xyz)
end
bfes  =  meshboundary(fes)
louter = selectelem(fens, bfes, facing = true, direction = dout)
outer_fes = subset(bfes, louter);
fens, fes  = H8extrudeQ4(fens, outer_fes,
          nlayers, (xyz, layer)->(R+layer/nlayers*(Ro-R))*xyz/norm(xyz));
fens, newfes1, fes2 =  mergemeshes(fens, fes, fens1w, fes1w, tolerance)
fes = cat(newfes1, fes2)
connected = findunconnnodes(fens, fes);
fens, new_numbering = compactnodes(fens, connected);
fes = renumberconn!(fes, new_numbering);


# File  =   "Sphere1.vtk"
# vtkexportmesh(File, fes.conn, fens.xyz, FinEtools.MeshExportModule.H8)
# @async run(`"paraview.exe" $File`)

geom  =  NodalField(fens.xyz)
P  =  NodalField(zeros(FFlt,size(fens.xyz,1),1))

bfes  =  meshboundary(fes)

linner = selectelem(fens, bfes, distance = R, from = [0.0 0.0 0.0],
  inflate = tolerance)
louter = selectelem(fens, bfes, facing = true, direction = dout)

println("Pre-processing time elapsed  =  ",time() - t0,"s")

t1  =  time()

numberdofs!(P)

material = MatAcoustFluid(bulk,rho)
femm  =  FEMMAcoust(IntegData(fes, GaussRule(3, 2)), material)

S  =  acousticstiffness(femm, geom, P);
C  =  acousticmass(femm, geom, P);

abcfemm  =  FEMMAcoustSurf(IntegData(subset(bfes, louter), GaussRule(2, 2)), material)
D  =  acousticABC(abcfemm, geom, P);

# ABC surface pressure loading
function abcp(dpdn, xyz, J, label, t)
    n = cross(J[:,1],J[:,2]);
    n = vec(n/norm(n));
    arg = (-dot(xyz,wavevector)+omega*(t+tshift));
    dpdn[1] = P_amplitude*cos(arg)*(-dot(n,wavevector));
end


targetfemm  =  FEMMAcoustSurf(IntegData(subset(bfes, linner), GaussRule(2, 2)), material)

ForceF = GeneralField(zeros(3,1))
numberdofs!(ForceF)
TorqueF = GeneralField(zeros(3,1))
numberdofs!(TorqueF)

GF = pressure2resultantforce(targetfemm, geom, P, ForceF)
GT = pressure2resultanttorque(targetfemm, geom, P, TorqueF, CG)

H = transpose(GF)*(rho/mass)*GF + transpose(GT)*(sparse(rho*inv(Inertia)))*GT;

Ctild = C + H;

# Solve
P0 = deepcopy(P)
P0.values[:] = 0.0; # initially all pressure is zero
vP0 = gathersysvec(P0);
vP1 = zeros(vP0);
vQ0 = zeros(vP0);
vQ1 = zeros(vP0);
t = 0.0;
P1 = deepcopy(P0);
tTa1 = zeros(3)
tAa1 = zeros(3)
tTa0 = zeros(3)
tAa0 = zeros(3)
tTv1 = zeros(3)
tAv1 = zeros(3)
tTv0 = zeros(3)
tAv0 = zeros(3)
tTu1 = zeros(3)
tAu1 = zeros(3)
tTu0 = zeros(3)
tAu0 = zeros(3)
tTa_store = zeros(3, nsteps+1)
tAa_store = zeros(3, nsteps+1)
tTv_store = zeros(3, nsteps+1)
tAv_store = zeros(3, nsteps+1)
tTu_store = zeros(3, nsteps+1)
tAu_store = zeros(3, nsteps+1)
t_store = zeros(1, nsteps+1)

pinc = deepcopy(P)
pincdd = deepcopy(P)
pincv = zeros(P.nfreedofs)
pincddv = zeros(P.nfreedofs)
p1v = zeros(P.nfreedofs)
p1 = deepcopy(P)

function recalculate_incident!(t, pinc, pincdd)
  for j = 1:count(fens)
    arg=(-dot(vec(fens.xyz[j,:]),wavevector)+omega*(t+tshift));# wave along the wavevector
    pinc.values[j] = P_amplitude*sin(arg)
    pincdd.values[j] = -omega^2*P_amplitude*sin(arg)
  end
  return pinc, pincdd
end

pinc, pincdd = recalculate_incident!(t, pinc, pincdd)
L0 = -S*gathersysvec!(pincdd, pincddv) - Ctild*gathersysvec!(pinc, pincv);
fabcp(dpdn, xyz, J, label) = abcp(dpdn, xyz, J, label, t)
fi  =  ForceIntensity(FFlt, 1, fabcp);
La0 = distribloads(abcfemm, geom, P1, fi, 2);

A = lu((2.0/dt)*S + D + (dt/2.)*Ctild);

step = 1;
while step <= nsteps
  println("Time $t ($(step)/$(Int(round(tfinal/dt))+1))")
  # Store for output
  tTa_store[:, step] = tTa1
  tAa_store[:, step] = tAa1
  tTv_store[:, step] = tTv1
  tAv_store[:, step] = tAv1
  tTu_store[:, step] = tTu1
  tAu_store[:, step] = tAu1
  t_store[1, step] = t
  # New loads and update
  step = step  + 1;
  t = t+dt;
  pinc, pincdd = recalculate_incident!(t, pinc, pincdd)
  L1 = -S*gathersysvec!(pincdd, pincddv) - Ctild*gathersysvec!(pinc, pincv);
  fabcp(dpdn, xyz, J, label) = abcp(dpdn, xyz, J, label, t)
  fi  =  ForceIntensity(FFlt, 1, fabcp);
  La1 = distribloads(abcfemm, geom, P1, fi, 2);
  vQ1 = A\((2/dt)*(S*vQ0)-D*vQ0-Ctild*(2*vP0+(dt/2)*vQ0)+L0+L1+La0+La1);
  vP1 = vP0 + (dt/2)*(vQ0+vQ1);
  P1 = scattersysvec!(P1, vec(vP1));
  p1.values[:] = P1.values[:] + pinc.values[:]
  p1v = gathersysvec!(p1, p1v);
  Fresultant = GF*p1v;
  Tresultant = GT*p1v;
  # println("Fresultant = $Fresultant")
  tTa1 = Fresultant/mass; tAa1 = Inertia\Tresultant;
  # Update  the velocity and position of the rigid body
  tTv1 = tTv0 + dt/2*(tTa0+tTa1);
  tTu1 = tTu0 + dt/2*(tTv0+tTv1);
  tAv1 = tAv1 + dt/2*(tAa0+tAa1);
  tAu1 = tAu0 + dt/2*(tAv0+tAv1);
  # show(tTa1)
  # Reset for the next time step
  copy!(vP0, vP1);
  copy!(vQ0, vQ1);
  P0 = deepcopy(P1);
  copy!(L0, L1);
  copy!(La0, La1);
  copy!(tTa0, tTa1); copy!(tAa0, tAa1)
  copy!(tTv0, tTv1); copy!(tAv0, tAv1)
  copy!(tTu0, tTu1); copy!(tAu0, tAu1)
end
tTa_store[:, step] = tTa1
tAa_store[:, step] = tAa1
tTv_store[:, step] = tTv1
tAv_store[:, step] = tAv1
tTu_store[:, step] = tTu1
tAu_store[:, step] = tAu1
t_store[1, step] = t

function remove_bias(it, iy)
  local t = deepcopy(reshape(it,length(it),1))
  local y = deepcopy(reshape(iy,length(iy),1))
  t = t/mean(t);
  t = reshape(t,length(t),1);
  y = reshape(y,length(y),1);
  local B = [t ones(t)];
  local p = (B'*B)\(B'*y);
  return y-p[1]*t-p[2];
end

ux = remove_bias(t_store, tTu_store[1,:])

function find_peaks(t,y)
    p=Float64[];
    dy=diff(y);
    for j3=1:length(dy)-1
        if (dy[j3]*dy[j3+1]<0)
            push!(p, y[j3+1]);
        end
    end
    return p
end
function amplitude(t,y)
    p= find_peaks(t,y) ;
    A=mean(abs.(p)) ;
    return A
end
println("Normalized amplitude  = $(amplitude(t_store,ux)/uamplHW)")

using Plots
plotly()
plot(transpose(t_store), transpose(tTu_store)/phun("mm"),
 xlabel = "t", ylabel = "Displacement", label = "Displacements")
plot(transpose(t_store), transpose(tAu_store),
 xlabel = "t", ylabel = "Displacement", label = "Rotations")
# plot(transpose(t_store), (ux)/uamplHW,
#  xlabel = "t", ylabel = "Normalized displacement", label = "Rectified")
gui()

#
File  =   "hemisphere_final_td1.vtk"
vtkexportmesh(File, fes.conn, geom.values, FinEtools.MeshExportModule.H8;
 scalars = [( "realP", real(P1.values))])
@async run(`"paraview.exe" $File`)
#
# println("Computing time elapsed  =  ",time() - t1,"s")
# println("Total time elapsed  =  ",time() - t0,"s")
#
#
# true
