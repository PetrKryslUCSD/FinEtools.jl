using FinEtools
using FinEtools.AlgoDeforLinearModule
using IterativeSolvers

println("""
Cantilever example.  Strongly orthotropic material. Orientation "y".
@article{
author = {Krysl, P.},
title = {Mean-strain 8-node hexahedron with optimized energy-sampling stabilization},
journal = {Finite Elements in Analysis and Design},
volume = {108}, pages = {41-53}, DOI = {10.1016/j.finel.2015.09.008}, year = {2016}
}
""")

t0 = time()
pu = ustring -> phun(ustring; system_of_units = :SIMM)

# # Orthotropic material
E1s = 100000.0*pu("GPa")
E2s = 1.0*pu("GPa")
E3s = E2s
nu23s = nu12s = nu13s = 0.25
G12s = 0.2*pu("GPa")
G23s = G13s = G12s
CTE1 = 0.0
CTE2 = 0.0
CTE3 = 0.0
# # Isotropic material
# E = 1.0e9*pu("Pa")
# nu = 0.25
# CTE = 0.0

# Reference value for  the vertical deflection of the tip
uz_ref = -1.027498445054843e-05*pu("m");

a = 90.0*pu("mm") # length of the cantilever
b = 10.0*pu("mm") # width of the cross-section
t = 20.0*pu("mm") # height of the cross-section
q0 = -1000.0*pu("Pa") # shear traction
dT = 0*pu("K") # temperature rise

tolerance = 0.00001*t

# Generate mesh
n = 20
na = n # number of elements lengthwise
nb = n # number of elements through the wwith
nt = n # number of elements through the thickness
xs = collect(linspace(0.0, a, na+1))
ys = collect(linspace(0.0, b, nb+1))
ts = collect(linspace(0.0, t, nt+1))
println("fens,fes = H8blockx(xs, ys, ts)")
@time fens,fes = H8blockx(xs, ys, ts)
println("fens,fes = H8toH20(fens,fes)")
@time fens,fes = H8toH20(fens,fes)
println("bfes = meshboundary(fes)")
@time bfes = meshboundary(fes)
# end cross-section surface  for the shear loading
sshearl = selectelem(fens, bfes; facing=true, direction = [+1.0 0.0 0.0])

MR = DeforModelRed3D
material = MatDeforElastOrtho(MR,
  0.0, E1s, E2s, E3s,
  nu12s, nu13s, nu23s,
  G12s, G13s, G23s,
  CTE1, CTE2, CTE3)
# material = MatDeforElastIso(MR,
#   0.0, E, nu, CTE)

# Material orientation matrix
csmat = zeros(3, 3)
rotmat3!(csmat, -45.0/180.0*pi*[0,1,0])

function updatecs!(csmatout::FFltMat, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)
  copy!(csmatout, csmat)
end

gr = GaussRule(3, 2)

femm = FEMMDeforLinear(MR,
    FEMMBase(fes, gr, CSys(updatecs!)), material)

lx0 = selectnode(fens, box=[0.0 0.0 -Inf Inf -Inf Inf], inflate=tolerance)

geom = NodalField(fens.xyz)
u = NodalField(zeros(size(fens.xyz,1),3)) # displacement field
nnodes(geom)

setebc!(u, lx0, true, 1, zeros(size(lx0)))
setebc!(u, lx0, true, 2, zeros(size(lx0)))
setebc!(u, lx0, true, 3, zeros(size(lx0)))
applyebc!(u)

S = connectionmatrix(femm.femmbase, nnodes(geom))
numberdofs!(u)

function getshr!(forceout::FFltVec, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)
  copy!(forceout, q0*[0.0; 0.0; 1.0])
end

Tracfemm = FEMMBase(subset(bfes, sshearl), GaussRule(2, 3))

println("K = stiffness(femm, geom, u)")
@time K = stiffness(femm, geom, u)
fi = ForceIntensity(FFlt, getshr!);
println("F =  distribloads(Tracfemm, geom, u, fi, 2);")
@time F =  distribloads(Tracfemm, geom, u, fi, 2);

println("K = cholfact(K)")
K = (K + K')/2;
@time K = cholfact(Symmetric(K))
println("U = K\\F")
@time U = K\F
# println("U = cg(K, F; tol=1e-3, maxiter=2000)")
# @time U = cg(K, F; tol=1e-3, maxiter=2000)
scattersysvec!(u, U[:])

Tipl = selectnode(fens, box=[a a b b 0. 0.], inflate=tolerance)
utip = mean(u.values[Tipl, 3])
println("Deflection $utip, normalized: $(utip/uz_ref)")
println("Solution: $(  time()-t0 )")

println("Done: $(  time()-t0 )")
true
