using FinEtools

println("""
Meshing an L-shaped membrane: merging multiple meshes.
""")

t0 = time()

W= 100. # strip width
L = 200. # length of the strip
nL = 15
nW = 10
tolerance = W/nW/1.0e5
Meshes = Array{Tuple{FENodeSet, FESet},1}()
push!(Meshes, Q4quadrilateral([0.0 0.0; W W], nW, nW))
push!(Meshes, Q4quadrilateral([-L 0.0; 0.0 W], nL, nW))
push!(Meshes, Q4quadrilateral([0.0 -L; W 0.0], nW, nL))
fens, outputfes = mergenmeshes(Meshes, tolerance);
fes = cat(outputfes[1], cat(outputfes[2], outputfes[3]))

geom = NodalField(fens.xyz)

File =  "L_shape.vtk"
vtkexportmesh(File, fes.conn, geom.values, FinEtools.MeshExportModule.Q4);
@async run(`"paraview.exe" $File`)

println("Done")
true
