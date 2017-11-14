using FinEtools

println("""
Meshing, deforming  and smoothing
""")

A = 100. # strip width
N = 16
tolerance = A / N / 1.0e5
fens, fes = T3block(A, A, N, N)

bnl = connectednodes(meshboundary(fes))
for ixxxx = 1:length(bnl)
    x, y = fens.xyz[bnl[ixxxx], :]
    fens.xyz[bnl[ixxxx], 1] += A / N * sin(2 * pi * y / A)
    fens.xyz[bnl[ixxxx], 2] += -A / N * sin(2 * pi * x / A)
end

File =  "mesh_smoothing_before.vtk"
vtkexportmesh(File, fens, fes);
@async run(`"paraview.exe" $File`)
println("$(fens.xyz[Int(N^2 / 2), :] )")

fixedv = falses(count(fens))
fixedv[bnl] = true
fens = meshsmoothing(fens, fes; fixedv=fixedv, method=:taubin, npass=100)

println("$(fens.xyz[Int(N^2 / 2), :] )")

geom = NodalField(fens.xyz)

File =  "mesh_smoothing_after.vtk"
vtkexportmesh(File, fes.conn, geom.values, FinEtools.MeshExportModule.T3);
@async run(`"paraview.exe" $File`)

println("Done")
true
