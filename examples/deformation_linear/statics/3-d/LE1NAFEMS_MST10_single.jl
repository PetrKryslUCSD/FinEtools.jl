using FinEtools
using FinEtools.MeshExportModule

using SymPy

@vars ua ub uc ud
@vars va vb vc vd
@vars wa wb wc wd
@vars x y z
u(x, y, z) = ua*x +  ub*y +  uc*z +  ud
v(x, y, z) = va*x +  vb*y +  vc*z +  vd
w(x, y, z) = wa*x +  wb*y +  wc*z +  wd

ex = diff(u(x, y, z), x)
ey = diff(v(x, y, z), y)
ez = diff(w(x, y, z), z)
gxy = diff(v(x, y, z), x) + diff(u(x, y, z), y)
gxz = diff(w(x, y, z), x) + diff(u(x, y, z), z)
gyz = diff(w(x, y, z), y) + diff(v(x, y, z), z)

display(ex)
display(ey)
display(ez)
display(gxy)
display(gxz)
display(gyz)
