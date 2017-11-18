using FinEtools
using FinEtools.MeshExportModule

# using SymPy
#
# @vars uaa ubb ucc uab uac ubc ua ub uc ud
# @vars vaa vbb vcc vab vac vbc va vb vc vd
# @vars waa wbb wcc wab wac wbc wa wb wc wd
# @vars x y z
# u(x, y, z) = uaa*x^2 + ubb*y^2 + ucc*z^2 + uab*x*y + uac*x*z + ubc*y*z + ua*x + ub*y + uc*z + ud
# v(x, y, z) = vaa*x^2 + vbb*y^2 + vcc*z^2 + vab*x*y + vac*x*z + vbc*y*z + va*x + vb*y + vc*z + vd
# w(x, y, z) = waa*x^2 + wbb*y^2 + wcc*z^2 + wab*x*y + wac*x*z + wbc*y*z + wa*x + wb*y + wc*z + wd
#
# ex = diff(u(x, y, z), x)
# ey = diff(v(x, y, z), y)
# ez = diff(w(x, y, z), z)
# gxy = diff(v(x, y, z), x) + diff(u(x, y, z), y)
# gxz = diff(w(x, y, z), x) + diff(u(x, y, z), z)
# gyz = diff(w(x, y, z), y) + diff(v(x, y, z), z)
#
# display(ex)
# display(ey)
# display(ez)
# display(gxy)
# display(gxz)
# display(gyz)


E = 200.0
nu = 0.3;

uaa, ubb, ucc, uab, uac, ubc,  ua,  ub,  uc,  ud =
0.0, 0.0, 0.0, 0.0, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0
vaa, vbb, vcc, vab, vac, vbc,  va,  vb,  vc,  vd =
0.0, 0.0, 0.0, 0.3, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
waa, wbb, wcc, wab, wac, wbc,  wa,  wb,  wc,  wd =
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
u(x, y, z) = uaa*x^2 + ubb*y^2 + ucc*z^2 + uab*x*y + uac*x*z + ubc*y*z + ua*x + ub*y + uc*z + ud
v(x, y, z) = vaa*x^2 + vbb*y^2 + vcc*z^2 + vab*x*y + vac*x*z + vbc*y*z + va*x + vb*y + vc*z + vd
w(x, y, z) = waa*x^2 + wbb*y^2 + wcc*z^2 + wab*x*y + wac*x*z + wbc*y*z + wa*x + wb*y + wc*z + wd

ex(x, y, z) = ua + 2*uaa*x + uab*y + uac*z
ey(x, y, z) = vab*x + vb + 2*vbb*y + vbc*z
ez(x, y, z) = wac*x + wbc*y + wc + 2*wcc*z
gxy(x, y, z) = uab*x + ub + 2*ubb*y + ubc*z + va + 2*vaa*x + vab*y + vac*z
gxz(x, y, z) = uac*x + ubc*y + uc + 2*ucc*z + wa + 2*waa*x + wab*y + wac*z
gyz(x, y, z) = vac*x + vbc*y + vc + 2*vcc*z + wab*x + wb + 2*wbb*y + wbc*z

fens,fes = H8block(1.0, 1.0, 1.0, 1, 1, 1)
# Now we shift  the nodes so that the center of the element is of the origin
fens.xyz[:, 1] -= 0.5
fens.xyz[:, 2] -= 0.5
fens.xyz[:, 3] -= 0.5

geom = NodalField(fens.xyz)
Uf = NodalField(zeros(size(fens.xyz,1),3)) # displacement field

for ix = 1:count(fens)
    x, y, z = fens.xyz[ix, :]
    Uf.values[ix, :] = [u(x, y, z), v(x, y, z), w(x, y, z)]
end


numberdofs!(Uf)
MR = DeforModelRed3D

material = MatDeforElastIso(MR, E, nu)
println("D = $(material.D)")
println("")

epsln(x, y, z) = vec([ex(x, y, z), ey(x, y, z), ez(x, y, z), gxy(x, y, z), gxz(x, y, z), gyz(x, y, z)])
sig(x, y, z) = material.D * epsln(x, y, z)
sigs = zeros(6, count(fens))
for ix = 1:count(fens)
    x, y, z = fens.xyz[ix, :]
    sigs[:, ix] = sig(x, y, z)
end
println("sigs = $(sigs)")
println("")

femm = FEMMDeforLinearMSH8(MR, IntegData(fes, GaussRule(3, 2)), material)
femm = associategeometry!(femm, geom)

fld = fieldfromintegpoints(femm, geom, Uf, :Cauchy, 1; reportat = :extraptrend)
println("fld.values = $(fld.values)")

File =  "a.vtk"
vtkexportmesh(File, fes.conn, geom.values,
    FinEtools.MeshExportModule.H8; vectors=[("u", Uf.values)],
    scalars=[("sigmay", fld.values)])
@async run(`"paraview.exe" $File`)
