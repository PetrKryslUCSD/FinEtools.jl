using FinEtools

println("""
Meshing  for a composite plate.
""")

t0 = time()

H = 100. # strip width
L = 200. # length of the strip
nL = 15
nH = 10
xs = collect(linspace(0.0, L, nL+1))
ys = collect(linspace(0.0, H, nH+1))
ts = [0.5; 2.0; 0.5; 2.0]
nts = [1; 2; 2; 1; ]
fens,fes = H8compositeplatex(xs, ys, ts, nts)

geom = NodalField(fens.xyz)

File =  "composite.vtk"
vtkexportmesh(File, fes.conn, geom.values, FinEtools.MeshExportModule.H8;
scalars = [( "Layer", fes.label)])
@async run(`"paraview.exe" $File`)

println("Done")
true
