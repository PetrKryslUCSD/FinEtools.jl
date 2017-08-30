using FinEtools

println("""
Meshing  for an elliptical hole in a finite strip.
""")

t0 = time()

H = 100. # strip width
ell = 50. # distance between holes
a = 10. # radius of the hole
L = 200. # length of the strip
nL = 15
nH = 10
nR = 50
fens,fes = Q4elliphole(a, a, L/2, H/2, nL, nH, nR)
# Note that we have to number the quadrilaterals  in the mirrored mesh
# differently in order to generate positive areas.
fens2, fes2 = mirrormesh(fens, fes, [0.0; -1.0],
  [0.0; 0.0]; renumb = c -> c[end:-1:1])
fens, fes1, fes2 = mergemeshes(fens, fes,  fens2, fes2, a/1000.)
fes = cat(fes1, fes2)

geom = NodalField(fens.xyz)

File =  "elliptical.vtk"
vtkexportmesh(File, fes.conn, geom.values, FinEtools.MeshExportModule.Q4)
@async run(`"paraview.exe" $File`)

println("Done")
true
