using FinEtools
using FinEtools.AlgoDeforLinearModule

println("""
Cantilever example.  Isotropic material.
@article{
author = {Krysl, P.},
title = {Mean-strain 8-node hexahedron with optimized energy-sampling stabilization},
journal = {Finite Elements in Analysis and Design},
volume = {108}, pages = {41-53}, DOI = {10.1016/j.finel.2015.09.008}, year = {2016}
}
""")

t0 = time()
# # Orthotropic material parameters of the external cylinder
# E1s = 130.0*phun("GPa")
# E2s = 5.0*phun("GPa")
# E3s = E2s
# nu12s = nu13s = 0.25
# nu23s = 0.0
# G12s = 10.0*phun("GPa")
# G13s = G12s
# G23s = 5.0*phun("GPa")
# CTE1 = 3.0e-6
# CTE2 = 2.0e-5
# CTE3 = 2.0e-5
# Isotropic material
E = 1.0e9*phun("Pa")
nu = 0.25
CTE = 0.0

# Reference value for  the vertical deflection of the tip
uz_ref =  -7.516310912734678e-06;

a = 90.0*phun("mm") # length of the cantilever
b = 10.0*phun("mm") # width of the cross-section
t = 20.0*phun("mm") # height of the cross-section
q0 = -1000.0*phun("Pa") # shear traction
dT = 0*phun("K") # temperature rise

tolerance = 0.00001*t

# Generate mesh
n = 24
na = n # number of elements lengthwise
nb = n # number of elements through the wwith
nt = n # number of elements through the thickness
xs = collect(linspace(0.0, a, na+1))
ys = collect(linspace(0.0, b, nb+1))
ts = collect(linspace(0.0, t, nt+1))
fens,fes = H8blockx(xs, ys, ts)
fens,fes = H8toH20(fens,fes)
bfes = meshboundary(fes)
# end cross-section surface  for the shear loading
sshearl = selectelem(fens, bfes; facing=true, direction = [+1.0 0.0 0.0])

MR = DeforModelRed3D
# externalmaterial = MatDeforElastOrtho(MR,
#   0.0, E1s, E2s, E3s,
#   nu12s, nu13s, nu23s,
#   G12s, G13s, G23s,
#   CTE1, CTE2, CTE3)
material = MatDeforElastIso(MR,
  0.0, E, nu, CTE)

# Material orientation matrix
csmat = [i==j ? one(FFlt) : zero(FFlt) for i=1:3, j=1:3]

function updatecs!(csmatout::FFltMat, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)
  copy!(csmatout, csmat)
end

gr = GaussRule(3, 2)

region = FDataDict("femm"=>FEMMDeforLinear(MR,
    FEMMBase(fes, gr), material))

lx0 = selectnode(fens, box=[0.0 0.0 -Inf Inf -Inf Inf], inflate=tolerance)

ex01 = FDataDict( "displacement"=>  0.0, "component"=> 1, "node_list"=>lx0 )
ex02 = FDataDict( "displacement"=>  0.0, "component"=> 2, "node_list"=>lx0 )
ex03 = FDataDict( "displacement"=>  0.0, "component"=> 3, "node_list"=>lx0 )

function getshr!(forceout::FFltVec, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)
  copy!(forceout, q0*[0.0; 0.0; 1.0])
end

Trac = FDataDict("traction_vector"=>getshr!,
    "femm"=>FEMMBase(subset(bfes, sshearl), GaussRule(2, 3)))

modeldata = FDataDict("fens"=>fens,
 "regions"=>[region],
 "essential_bcs"=>[ex01, ex02, ex03],
 "traction_bcs"=>[Trac],
 "temperature_change"=>FDataDict("temperature"=>dT)
 )
modeldata = AlgoDeforLinearModule.linearstatics(modeldata)

u = modeldata["u"]
geom = modeldata["geom"]

Tipl = selectnode(fens, box=[a a b b 0. 0.], inflate=tolerance)
utip = mean(u.values[Tipl, 3])
println(" Normalized deflection: $(utip/uz_ref)")
println("Solution: $(  time()-t0 )")

# File =  "NAFEMS-R0031-2-plate.vtk"
# vtkexportmesh(File, fes.conn, geom.values, FinEtools.MeshExportModule.H20;
#     scalars = [("Layer", fes.label)], vectors = [("displacement", u.values)])
# @async run(`"paraview.exe" $File`)

modeldata["postprocessing"] = FDataDict("file"=>"fiber_reinf_cant_iso",
  "outputcsys"=>CSys(updatecs!), "quantity"=>:Cauchy, "component"=>5)
# modeldata = AlgoDeforLinearModule.exportstress(modeldata)
modeldata = AlgoDeforLinearModule.exportdeformation(modeldata)

println("Done: $(  time()-t0 )")
true
