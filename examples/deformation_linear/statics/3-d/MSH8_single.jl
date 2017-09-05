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
0.1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
vaa, vbb, vcc, vab, vac, vbc,  va,  vb,  vc,  vd = 
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
waa, wbb, wcc, wab, wac, wbc,  wa,  wb,  wc,  wd =
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
u(x, y, z) = uaa*x^2 + ubb*y^2 + ucc*z^2 + uab*x*y + uac*x*z + ubc*y*z + ua*x + ub*y + uc*z + ud
v(x, y, z) = vaa*x^2 + vbb*y^2 + vcc*z^2 + vab*x*y + vac*x*z + vbc*y*z + va*x + vb*y + vc*z + vd
w(x, y, z) = waa*x^2 + wbb*y^2 + wcc*z^2 + wab*x*y + wac*x*z + wbc*y*z + wa*x + wb*y + wc*z + wd

fens,fes = H8block(1.0, 1.0, 1.0, 1, 1, 1)
geom = NodalField(fens.xyz)
Uf = NodalField(zeros(size(fens.xyz,1),3)) # displacement field

for ix = 1:count(fens)
    x, y, z = fens.xyz[ix, :]
    Uf.values[ix, :] = [u(x, y, z), v(x, y, z), w(x, y, z)]
end


numberdofs!(Uf)
MR = DeforModelRed3D

material = MatDeforElastIso(MR, E, nu)

femm = FEMMDeforLinearMSH8(MR, GeoD(fes, GaussRule(3, 2)), material)

fld = fieldfromintegpoints(femm, geom, Uf, :Cauchy, 1; tonode = :estimtrendpaper)
println("$(fld.values)")
