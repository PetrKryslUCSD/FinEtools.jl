using FinEtools
using FinEtools.MeshExportModule

# Thick elliptical plate with an elliptical hole is clamped on its exterior
# boundary and is loaded with transverse  pressure.
# This is a NAFEMS Benchmark, Test No. LE10.
# The plate is discretized with solid elements.
# Because of the symmetries of the geometry and load, only quarter of the plate is modeled.
# The $\sigma_y=\sigma_2$ at the point $P$ is to be determined. Since the
# target point is on the boundary of the domain it will not be an
# integration node as we use Gauss quadrature. The reference value is -5.38 MPa.

println("LE10NAFEMS: Transverse deflection of elliptical plate with elliptical hole."        )
t0 = time()

E = 210e3*phun("MEGA*PA");# 210e3 MPa
nu = 0.3;
qmagn = 1.0*phun("MEGA*PA");# transverse pressure
sigma_yP = -5.38*phun("MEGA*PA");# tensile stress at [2.0, 0.0] meters
Ae =3.25*phun("m"); # Major radius of the exterior ellipse
Be =2.75*phun("m"); # Minor radius of the exterior ellipse
Ai =2.0*phun("m"); # Major radius of the interior ellipse
Bi =1.0*phun("m"); # Minor radius of the interior ellipse
Thickness = 0.6*phun("m")
nc = 6; # number of elements per side
nr = 5; # number of elements per side
nt = 2; # number of elements through the thickness
# nc = 26; # number of elements per side
# nr = 25; # number of elements per side
# nt = 18; # number of elements through the thickness
tolerance = Thickness/nt/1000.; # Geometrical tolerance

fens,fes = Q4block(1.0, pi/2, nr, nc)
#
@assert nt % 2 == 0 "Number of elements through the thickness must be even"
fens,fes  = H8extrudeQ4(fens, fes,
  nt, (xyz, layer)->[xyz[1], xyz[2], (layer)/nt*Thickness]);

# Select the  boundary faces, on the boundary that is clamped,  and on the part
# of the boundary that is loaded with the transverse pressure
bdryfes = meshboundary(fes);
exteriorbfl = selectelem(fens, bdryfes, box=[1.0, 1.0, 0.0, pi/2, 0.0, Thickness], inflate=tolerance);
topbfl = selectelem(fens, bdryfes, box=[0.0, 1.0, 0.0, pi/2, Thickness, Thickness], inflate=tolerance);

# Reshape the generated block into the elliptical plate
for i=1:count(fens)
    r=fens.xyz[i,1]; a=fens.xyz[i,2]; z=fens.xyz[i,3]
    fens.xyz[i,:]=[(r*Ae+(1-r)*Ai)*cos(a) (r*Be+(1-r)*Bi)*sin(a) z];
end


geom = NodalField(fens.xyz)
u = NodalField(zeros(size(fens.xyz,1),3)) # displacement field

l12 =connectednodes(subset(bdryfes, exteriorbfl)) # external boundary
setebc!(u, l12, true, 1, 0.0)
setebc!(u, l12, true, 2, 0.0)
ll = selectnode(fens; box=[0.0, Inf, 0.0, Inf, Thickness/2.0, Thickness/2.0], inflate = tolerance)
l3 = intersect(ll, connectednodes(subset(bdryfes, exteriorbfl)))
setebc!(u, l3, true, 3, 0.0)
l1 =selectnode(fens; box=[0.0, 0.0, 0.0, Inf, 0.0, Thickness], inflate = tolerance)
setebc!(u,l1,true, 1, 0.0) # symmetry plane X = 0
l2 =selectnode(fens; box=[0.0, Inf, 0.0, 0.0, 0.0, Thickness], inflate = tolerance)
setebc!(u,l2,true, 2, 0.0) # symmetry plane Y = 0

applyebc!(u)
numberdofs!(u)

el1femm =  FEMMBase(GeoD(subset(bdryfes,topbfl), GaussRule(2, 2)))
function pfun(forceout::FVec{T}, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt) where {T}
    forceout .=  [0.0, 0.0, -qmagn]
    return forceout
end
fi = ForceIntensity(FFlt, 3, pfun);
F2 = distribloads(el1femm, geom, u, fi, 2);

# Note that the material object needs to be created with the proper
# model-dimension reduction in mind.  In this case that is the fully three-dimensional solid.
MR = DeforModelRed3D

material = MatDeforElastIso(MR, E, nu)

femm = FEMMDeforLinearMSH8(MR, GeoD(fes, GaussRule(3, 2)), material)

# The geometry field now needs to be associated with the FEMM
femm = associategeometry!(femm, geom)

K = stiffness(femm, geom, u)
K = cholfact(K)
U = K\(F2)
scattersysvec!(u, U[:])

nl = selectnode(fens, box=[Ai,Ai,0,0,Thickness,Thickness],inflate=tolerance);
thecorneru = zeros(FFlt,1,3)
gathervalues_asmat!(u, thecorneru, nl);
thecorneru = thecorneru/phun("mm")
println("displacement =$(thecorneru) [MM] as compared to reference [-0.030939, 0, -0.10488] [MM]")


fld= fieldfromintegpoints(femm, geom, u, :Cauchy, 2; reportat = :extraptrend)#
println("Sigma_y =$(fld.values[nl,1][1]/phun("MPa")) as compared to reference sigma_yP = $(sigma_yP/phun("MPa")) [MPa]")

println("$((nc, nr, nt)), $(fld.values[nl,1][1]/phun("MPa"))")

# File =  "LE10NAFEMS_sigmay.vtk"
# vtkexportmesh(File, fes.conn, geom.values,
#                FinEtools.MeshExportModule.H8; vectors=[("u", u.values)],
#                scalars=[("sigmay", fld.values)])
# @async run(`"paraview.exe" $File`)
# true


AE = AbaqusExporter("LE10NAFEMS_H8");
HEADING(AE, "LE10NAFEMS: Transverse deflection of elliptical plate with elliptical hole.");
PART(AE, "part1");
END_PART(AE);
ASSEMBLY(AE, "ASSEM1");
INSTANCE(AE, "INSTNC1", "PART1");
NODE(AE, fens.xyz);
ELEMENT(AE, "c3d8rh", "AllElements", 1, femm.geod.fes.conn)
ELEMENT(AE, "SFM3D4", "TractionElements",
1+count(femm.geod.fes), el1femm.geod.fes.conn)
NSET_NSET(AE, "l1", l1)
NSET_NSET(AE, "l2", l2)
NSET_NSET(AE, "l3", l3)
NSET_NSET(AE, "l12", l12)
ORIENTATION(AE, "GlobalOrientation", vec([1. 0 0]), vec([0 1. 0]));
SOLID_SECTION(AE, "elasticity", "GlobalOrientation", "AllElements", "Hourglassctl");
SURFACE_SECTION(AE, "TractionElements")
END_INSTANCE(AE);
END_ASSEMBLY(AE);
MATERIAL(AE, "elasticity")
ELASTIC(AE, E, nu)
SECTION_CONTROLS(AE, "Hourglassctl", "HOURGLASS=ENHANCED")
STEP_PERTURBATION_STATIC(AE)
BOUNDARY(AE, "ASSEM1.INSTNC1.l1", 1)
BOUNDARY(AE, "ASSEM1.INSTNC1.l2", 2)
BOUNDARY(AE, "ASSEM1.INSTNC1.l3", 3)
BOUNDARY(AE, "ASSEM1.INSTNC1.l12", 1)
BOUNDARY(AE, "ASSEM1.INSTNC1.l12", 2)
DLOAD(AE, "ASSEM1.INSTNC1.TractionElements", vec([0.0, 0.0, -qmagn]))
END_STEP(AE)
close(AE)
