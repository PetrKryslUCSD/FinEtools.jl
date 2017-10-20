using FinEtools
using FinEtools.AlgoDeforLinearModule
using FinEtools.MeshUtilModule

println("""
Pagano_3layer_cylindrical_bending_MSH8: mean-strain hexahedron
""")

# This example provides three-dimensional finite element model for  the
# transverse shear stress calculations. The problem consists of a one-, two- or
# three-layer plate subjected to a sinusoidal distributed load, as
# described by Pagano (1969). The resulting transverse shear and axial
# stresses through the thickness of the plate are compared to two existing
# analytical solutions by Pagano (1969). The first solution is derived from
# classical laminated plate theory (CPT), while the second is an exact
# solution from linear elasticity theory.


filebase = "Pagano_3layer_cylindrical_bending_MSH8"

t0 = time()
# Orthotropic material for the 3 layers
E1 = 25e6*phun("PSI"); E2 = 1e6*phun("PSI"); E3 = E2; 
G12 = 0.5e6*phun("PSI");  G13 = G12; G23 = 0.2e6*phun("PSI")
nu12 =  0.25; nu13 =  0.25; nu23 =  0.25;
Span_to_thickness = 4.0;
T = 2.5*phun("in"); # total thickness of the plate
L = Span_to_thickness*T; 
h = 1*phun("in");  # depth of the plate
q0 = 1*phun("PSI")
CTE1 =  CTE2 =  CTE3 = 0.0

# Here we define the layout and the thicknesses of the layers.
angles = vec([0.0 90.0 0.0]);
nLayers = length(angles)
ts = T/nLayers * ones(nLayers); # layer thicknesses

tolerance = 0.0001*T

# Reference deflection under the load is
wtopref = -3.789*phun("mm"); # From [1]
wbottomref = -2.16*phun("mm"); # Not given in [1]; guessed from the figure

# Select how find the mesh should be
Refinement = 1
nL, nh = Refinement*40, Refinement*1;
nts= Refinement * 8 * ones(Int, nLayers);# number of elements per layer

xs = collect(linspace(0.0, L, nL+1))
ys = collect(linspace(0.0, h, nh+1))

fens,fes = H8compositeplatex(xs, ys, ts, nts)

# This is the material  model
MR = DeforModelRed3D
skinmaterial = MatDeforElastOrtho(MR,
0.0, E1, E2, E3,
nu12, nu13, nu23,
G12, G13, G23,
CTE1, CTE2, CTE3)

# The material coordinate system function is defined as:
function updatecs!(csmatout::FFltMat, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)
    rotmat3!(csmatout, angles[fe_label]/180.0*pi* [0.0; 0.0; 1.0]);
end

# The vvolume integrals are evaluated using this rule
gr = GaussRule(3, 2)

# We will create 3 regions, one for each of the layers
regions = FDataDict[]
for layer = 1:nLayers
    rls = selectelem(fens, fes, label =  layer)
    push!(regions, FDataDict("femm"=>FEMMDeforLinearMSH8(MR,
    GeoD(subset(fes, rls), gr, CSys(3, 3, updatecs!)), skinmaterial)))
end

# File =  "Meyer_Piening_sandwich-r1.vtk"
# vtkexportmesh(File, skinregion["femm"].geod.fes.conn, fens.xyz, FinEtools.MeshExportModule.H8)
# # @async run(`"paraview.exe" $File`)


# The essential boundary conditions are applied to enforce the plane strain constraint.
ly0 = selectnode(fens, box=[-Inf Inf 0.0 0.0 -Inf Inf], inflate=tolerance)
lyh = selectnode(fens, box=[-Inf Inf h h -Inf Inf], inflate=tolerance)
ey = FDataDict("displacement"=>  0.0, "component"=> 2, "node_list"=>vcat(ly0, lyh))
# The transverse displacement is fixed at the two ends.
lz0 = selectnode(fens, box=[0.0 0.0 -Inf Inf -Inf Inf], inflate=tolerance)
lzL = selectnode(fens, box=[L L -Inf Inf -Inf Inf], inflate=tolerance)
ez = FDataDict("displacement"=>  0.0, "component"=> 3, "node_list"=>vcat(lz0, lzL))
ex = FDataDict("displacement"=>  0.0, "component"=> 1, "node_list"=>[1])

# The traction boundary condition is applied at the top of the plate.
bfes = meshboundary(fes)
function pfun(forceout::FVec{T}, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt) where {T}
    forceout[1] = 0.0
    forceout[2] = 0.0
    forceout[3] = -q0*sin(pi*XYZ[1]/L)
    return forceout
end
# From  the entire boundary we select those quadrilaterals that lie on the plane
# Z = thickness
tl = selectelem(fens, bfes, box = [-Inf Inf -Inf Inf T T], inflate=tolerance)
Trac = FDataDict("traction_vector"=>pfun,
"femm"=>FEMMBase(GeoD(subset(bfes, tl), GaussRule(2, 2))))

modeldata = FDataDict("fens"=>fens,
"regions"=>regions,
"essential_bcs"=>[ex, ey, ez],
"traction_bcs"=> [Trac]
)
modeldata = AlgoDeforLinearModule.linearstatics(modeldata)

modeldata["postprocessing"] = FDataDict("file"=>filebase * "-u")
modeldata = AlgoDeforLinearModule.exportdeformation(modeldata)

u = modeldata["u"]
geom = modeldata["geom"]

# The results of the displacement and stresses will be reported at
# nodes located at the appropriate points.
ntopcenter = selectnode(fens, box=[L/2 L/2 0.0 h T T], inflate=tolerance)
# ntopcenter = selectnode(fens, box=[0.0 0.0 0.0 0.0 TH TH], inflate=tolerance)
# ncenterline = selectnode(fens, box=[0.0 0.0 0.0 0.0 0.0 TH], inflate=tolerance)
# nintertop = selectnode(fens, box=[-Inf Inf 0.0 0.0 sum(ts[1:2]) sum(ts[1:2])], inflate=tolerance)
# ninterbot = selectnode(fens, box=[-Inf Inf 0.0 0.0 sum(ts[1:1]) sum(ts[1:1])], inflate=tolerance)

# zclo = sortperm(vec(geom.values[ncenterline, 3]))
# ncenterline = ncenterline[zclo]
# centerz = geom.values[ncenterline, 3]
# zclo = nothing

# xclotop = sortperm(vec(geom.values[nintertop, 1]))
# nintertop = nintertop[xclotop]
# topx = geom.values[nintertop, 1]
# xclobot = sortperm(vec(geom.values[ninterbot, 1]))
# ninterbot = ninterbot[xclobot]
# botx = geom.values[ninterbot, 1]
# xclotop = xclobot = nothing

# conninbotskin = intersect(connectednodes(botskinregion["femm"].geod.fes), ncenterline)
# connincore = intersect(connectednodes(coreregion["femm"].geod.fes), ncenterline)
# connintopskin = intersect(connectednodes(topskinregion["femm"].geod.fes), ncenterline)
# inbotskin = [n in conninbotskin for n in ncenterline]
# incore = [n in connincore for n in ncenterline]
# intopskin = [n in connintopskin for n in ncenterline]

println("")
println("Top Center deflection: $(mean(u.values[ntopcenter, 3], 1)/phun("in")) [in]")

# # # extrap = :extrapmean
# extrap = :extraptrend
# nodevalmeth = :averaging
# # extrap = :default
# # nodevalmeth = :invdistance

# # Normal stress in the X direction
# modeldata["postprocessing"] = FDataDict("file"=>filebase * "-sx",
# "quantity"=>:Cauchy, "component"=>1, "outputcsys"=>CSys(3),
# "nodevalmethod"=>nodevalmeth, "reportat"=>extrap)
# modeldata = AlgoDeforLinearModule.exportstress(modeldata)
# s = modeldata["postprocessing"]["exported"][1]["field"]
# sxbot = s.values[ncenterline, 1]
# s = modeldata["postprocessing"]["exported"][2]["field"]
# sxcore = s.values[ncenterline, 1]
# s = modeldata["postprocessing"]["exported"][3]["field"]
# sxtop = s.values[ncenterline, 1]

# # The graph data needs to be collected by going through each layer separately.
# # Some quantities may be discontinuous between layers.
# zs = vcat(  [z for (j,z) in enumerate(centerz) if inbotskin[j]],
# [z for (j,z) in enumerate(centerz) if incore[j]],
# [z for (j,z) in enumerate(centerz) if intopskin[j]]
# )
# sxs = vcat( [sxbot[j] for (j,z) in enumerate(centerz) if inbotskin[j]],
# [sxcore[j] for (j,z) in enumerate(centerz) if incore[j]],
# [sxtop[j] for (j,z) in enumerate(centerz) if intopskin[j]]
# )

# File = filebase * "-sx-$(extrap).CSV"
# savecsv(File, zs=vec(zs)/phun("mm"), sx=vec(sxs)/phun("MPa"))

# # @async run(`"paraview.exe" $File`)

# # Inter laminar stress between the skin and the core
# modeldata["postprocessing"] = FDataDict("file"=>filebase * "-sxz",
# "quantity"=>:Cauchy, "component"=>5, "outputcsys"=>CSys(3),
# "nodevalmethod"=>nodevalmeth, "reportat"=>extrap)
# modeldata = AlgoDeforLinearModule.exportstress(modeldata)
# s = modeldata["postprocessing"]["exported"][1]["field"]
# sxzskinbot = s.values[ninterbot, 1]
# s = modeldata["postprocessing"]["exported"][2]["field"]
# sxzcoretop = s.values[nintertop, 1]
# sxzcorebot = s.values[ninterbot, 1]
# s = modeldata["postprocessing"]["exported"][3]["field"]
# sxzskintop = s.values[nintertop, 1]

# File = filebase * "-sxz-$(extrap).CSV"
# savecsv(File, xstop=vec(topx)/phun("mm"),
# sxzskintop=vec(sxzskintop)/phun("MPa"),
# sxzcoretop=vec(sxzcoretop)/phun("MPa"),
# xsbot=vec(botx)/phun("mm"),
# sxzskinbot=vec(sxzskinbot)/phun("MPa"),
# sxzcorebot=vec(sxzcorebot)/phun("MPa"))

# @async run(`"paraview.exe" $File`)

println("Done")
true
