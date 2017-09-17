using FinEtools
using FinEtools.AlgoDeforLinearModule

println("""
Laminated Strip Under Three-Point Bending
""")

# Determine the central transverse displacement in a simply-supported seven
# layer symmetric strip with a central line load. A 0/90/0/90/0/90/0
# material lay-up is specified with the center ply being four times as
# thick as the others.
# Reference: NAFEMS Report R0031, Test No.1, 17-Dec-1998.

# Because of the symmetries of the geometry and load, only the
# first-quadrant   (in XY) quarter of the plate is modeled.

# The coordinate system is centered at point E (at the difference with
# respect to the original benchmark definition).  The  load is applied
# along a curve passing through point C. The simple support is applied
# along the curve passing through point B.

# We realize the simple supports along the lines  A, B and the line load at
# point C  are illegal from the point of view of convergence.  No
# convergence can be hoped for as the stress underneath the load and above
# the simple supports  is infinite in the limit (these locations are stress
# singularities).   However, for relatively coarse meshes the results away
# from the singularities are still meaningful.

# The target quantities are displacement at the bottom surface at point E,
# the tensile axial stress at the same point,  and of the transverse shear
# stress at point D  in between the bottommost two layers (See figure 1).

t0 = time()
# Orthotropic material parameters of the material of the layers
E1s = 100.0*phun("GPa")
E2s = E3s = 5.0*phun("GPa")
nu12s = 0.4
nu13s = 0.3
nu23s = 0.3
G12s = 3.0*phun("GPa")
G13s = G23s = 2.0*phun("GPa")
CTE1 = 3.0e-6
CTE2 = 2.0e-5
CTE3 = 2.0e-5

AB = 30.0*phun("mm") # span between simple supports
OH = 10.0*phun("mm") # overhang
W = 10.0*phun("mm") # width of the strip

# Here we define the layout and the thicknesses of the layers.
angles = vec([0 90 0 90 0 90 0]);
ts = vec([0.1  0.1  0.1  0.4  0.1  0.1  0.1])*phun("mm"); # layer thicknesses
TH = sum(ts); # total thickness of the plate

tolerance = 0.0001*TH

# The line load is in the negative Z direction.
q0 = 10*phun("N/mm"); #    line load

# Reference deflection under the load is
wEref = -1.06*phun("mm");

# The reference tensile stress at the bottom of the lowest layer is
sigma11Eref = 684*phun("MPa");

# Because we model the first-quadrant quarter of the plate using
# coordinate axes centered  at the point E  the shear at the point D is
# positive instead of negative as in the benchmark where the coordinate
# system is located at the outer corner of the strip.
sigma13Dref=4.1*phun("MPa");

Refinement = 9
# We select 8 elements spanwise and 2 elements widthwise.  The overhang
# of the plate is given one element.
nL = Refinement * 4; nO = Refinement * 1; nW = Refinement * 1;

# Each layer is modeled with a single element.
nts= Refinement * ones(Int, length(angles));# number of elements per layer

xs = unique(vcat(collect(linspace(0,AB/2,nL+1)),
    collect(linspace(AB/2,AB/2+OH,nO+1))))
ys = collect(linspace(0,W/2,nW+1));

fens,fes = H8compositeplatex(xs, ys, ts, nts)


# This is the material  model
MR = DeforModelRed3D
material = MatDeforElastOrtho(MR,
    0.0, E1s, E2s, E3s,
    nu12s, nu13s, nu23s,
    G12s, G13s, G23s,
    CTE1, CTE2, CTE3)

# The material coordinate system function is defined as:
function updatecs!(csmatout::FFltMat, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)
    rotmat3!(csmatout, angles[fe_label]/180.0*pi* [0.0; 0.0; 1.0]);
end

# The vvolume integrals are evaluated using this rule
gr = GaussRule(3, 2)

# We will create two regions, one for the layers with 0°  orientation,
# and one for the layers with 90° orientation.
rl1 = vcat(selectelem(fens, fes, label = 1), selectelem(fens, fes, label = 3),
    selectelem(fens, fes, label = 5), selectelem(fens, fes, label = 7))
rl2 = vcat(selectelem(fens, fes, label = 2), selectelem(fens, fes, label = 4),
        selectelem(fens, fes, label = 6))
region1 = FDataDict("femm"=>FEMMDeforLinearMSH8(MR,
    GeoD(subset(fes, rl1), gr, CSys(3, 3, updatecs!)), material))
region2 = FDataDict("femm"=>FEMMDeforLinearMSH8(MR,
    GeoD(subset(fes, rl2), gr, CSys(3, 3, updatecs!)), material))

# File =  "NAFEMS-R0031-1-plate-r1.vtk"
# vtkexportmesh(File, region1["femm"].geod.fes.conn, fens.xyz, FinEtools.MeshExportModule.H8)
# # @async run(`"paraview.exe" $File`)
# File =  "NAFEMS-R0031-1-plate-r2.vtk"
# vtkexportmesh(File, region2["femm"].geod.fes.conn, fens.xyz, FinEtools.MeshExportModule.H8)
# @async run(`"paraview.exe" $File`)

# The essential boundary conditions are applied on the symmetry planes.
# First the plane X=0;...
lx0 = selectnode(fens, box=[0.0 0.0 -Inf Inf -Inf Inf], inflate=tolerance)
ex0 = FDataDict( "displacement"=>  0.0, "component"=> 1, "node_list"=>lx0 )
# ... and then the plane Y=0.
ly0 = selectnode(fens, box=[-Inf Inf 0.0 0.0 -Inf Inf], inflate=tolerance)
ey0 = FDataDict( "displacement"=>  0.0, "component"=> 2, "node_list"=>ly0 )
# The transverse displacement is fixed along the line  passing through
# point B. The nodes are fixed in the box along this line in the Z
# direction.
lz0 = selectnode(fens, box=[AB/2 AB/2 -Inf Inf -Inf Inf], inflate=tolerance)
ez0 = FDataDict( "displacement"=>  0.0, "component"=> 3, "node_list"=>lz0 )

# The traction boundary condition is applied  along  the edge of the
# mesh passing through point C at the top surface of the strip.   First
# we extract the boundary of the hexahedral mesh.
bfes = meshboundary(fes)
# From  the entire boundary we select those quadrilaterals that lie on the plane
# X = 0
xl = selectelem(fens, bfes, box = [0.0 0.0 -Inf Inf -Inf Inf], inflate=tolerance)
# Now we extract the boundary  of these selected quadrilaterals
bbfes = meshboundary(subset(bfes, xl))
# …  And from these  we extract the ones at the top
zl = selectelem(fens, bbfes, box = [0.0 0.0 -Inf Inf TH TH], inflate=tolerance)
# Note that  we have to apply only half of the line load given that
# were modeling  just one quarter of the geometry and we are splitting
# the line load  with the symmetry plane X=0. Also note that the
# quadrature rule is one-dimensional  since we are integrating along
# a curve.
Trac = FDataDict("traction_vector"=>vec([0.0; 0.0; -q0/2]),
    "femm"=>FEMMBase(GeoD(subset(bbfes, zl), GaussRule(1, 3))))

modeldata = FDataDict("fens"=>fens,
 "regions"=>[region1, region2],
 "essential_bcs"=>[ex0, ey0, ez0],
 "traction_bcs"=> [Trac]
 )
modeldata = AlgoDeforLinearModule.linearstatics(modeldata)

modeldata["postprocessing"] = FDataDict("file"=>"NAFEMS-R0031-1-plate")
modeldata = AlgoDeforLinearModule.exportdeformation(modeldata)

u = modeldata["u"]
geom = modeldata["geom"]

# The results of the displacement and stresses will be reported at
# nodes located at the appropriate points.
nE = selectnode(fens, box=[0.0 0.0 0.0 0.0 0.0 0.0], inflate=tolerance)
nC = selectnode(fens, box=[0.0 0.0 0.0 0.0 TH TH], inflate=tolerance)
nD = selectnode(fens, box=[0.0 0.0 0.0 0.0 ts[1] ts[1]], inflate=tolerance)
n0z = selectnode(fens, box=[0.0 0.0 0.0 0.0 0.0 TH], inflate=tolerance)
ix = sortperm(geom.values[n0z, 3])
# println("ix = $(ix)")

cdis = mean(u.values[nE, 3])
println("")
println("Normalized Center deflection: $(cdis/wEref)")

# # extrap = :extraptrend
# # extrap = :extrapmean
# inspectormeth = :averaging
extrap = :default
inspectormeth = :invdistance

modeldata["postprocessing"] = FDataDict("file"=>"NAFEMS-R0031-1-plate-sx",
    "quantity"=>:Cauchy, "component"=>1, "outputcsys"=>CSys(3),
     "nodevalmethod"=>inspectormeth, "reportat"=>extrap)
modeldata = AlgoDeforLinearModule.exportstress(modeldata)
s = modeldata["postprocessing"]["exported_fields"][1]
println("sx@E = $(s.values[nE]/phun("MPa")) [MPa]")

modeldata["postprocessing"] = FDataDict("file"=>"NAFEMS-R0031-1-plate-sxz",
"quantity"=>:Cauchy, "component"=>5, "outputcsys"=>CSys(3),
 "nodevalmethod"=>inspectormeth, "reportat"=>extrap)
modeldata = AlgoDeforLinearModule.exportstress(modeldata)
s = modeldata["postprocessing"]["exported_fields"][1]
println("sxz@D_1 = $(s.values[nD]/phun("MPa")) [MPa]")
s = modeldata["postprocessing"]["exported_fields"][2]
println("sxz@D_2 = $(s.values[nD]/phun("MPa")) [MPa]")


#
# s = fieldfromintegpoints(region1["femm"], geom, u, :Cauchy, 1;
#     outputcsys = CSys(3), nodevalmethod = inspectormeth, reportat = extrap)
# println("sx@E = $(s.values[nE]/phun("MPa")) [MPa]")
# sx_z = s.values[n0z]/phun("MPa")
# println("sx(z)_1 = $(sx_z)")
#
# s = fieldfromintegpoints(region1["femm"], geom, u, :Cauchy, 5;
#     outputcsys = CSys(3), nodevalmethod = inspectormeth, reportat = extrap)
# println("sxz@D_1 = $(s.values[nD]/phun("MPa")) [MPa]")
# sxz_z_1 = s.values[n0z]/phun("MPa")
# println("sxz(z)_1 = $(sxz_z_1)")
# s = fieldfromintegpoints(region2["femm"], geom, u, :Cauchy, 5;
#     outputcsys = CSys(3), nodevalmethod = inspectormeth, reportat = extrap)
# println("sxz@D_2 = $(s.values[nD]/phun("MPa")) [MPa]")
# sxz_z_2 = s.values[n0z]/phun("MPa")
# println("sxz(z)_2 = $(sxz_z_2)")

# function _inspector(idat, elnum, conn, xe,  out,  xq)
#     # xe = coordinates of the nodes of the element
#     # xq = coordinate of the quadrature point
#     println("@$(xq): $(out/1.0e6)")
#     return idat
# end
#
# felist = selectelem(fens, region1["femm"].geod.fes,
#     box=[0.0 0.0 0.0 0.0 0.0 0.0], inflate=tolerance, allin = false)
#
# inspectintegpoints(region1["femm"], geom, u, felist,
#     _inspector, 0, quantity=:Cauchy, outputcsys = CSys(3))
#
# femm = deepcopy(region1["femm"])
# femm.geod.fes = subset(femm.geod.fes, felist)
# associategeometry!(femm, geom)
# s = fieldfromintegpoints(femm, geom, u, :Cauchy, 5;
#     outputcsys = CSys(3), nodevalmethod = inspectormeth, reportat = extrap)
# println("sxz@D_1 = $(s.values[nD]/phun("MPa")) [MPa]")

# felist = selectelem(fens, region2["femm"].geod.fes,
#     box=[0.0 0.0 0.0 0.0 0.0 TH], inflate=tolerance, allin = false)
#
# inspectintegpoints(region2["femm"], geom, u, felist,
#     _inspector, 0, quantity=:Cauchy, outputcsys = CSys(3))


println("Done")
true
