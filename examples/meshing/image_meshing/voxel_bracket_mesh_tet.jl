using FinEtools

V = VoxelBoxVolume(Int, 6*[5,6,7], [4.0, 4.0, 5.0])

b1 = solidbox((0.0, 0.0, 0.0), (1.0, 4.0, 5.0))
b2 = solidbox((0.0, 0.0, 0.0), (4.0, 1.0, 5.0))
h1 = solidcylinder((2.0, 2.5, 2.5), (1.0, 0.0, 0.0), 00.75)
fillsolid!(V, differenceop(unionop(b1, b2), h1), 1)

fens, fes = T4voximg(V.data, vec([voxeldims(V)...]), [1])
fens = meshsmoothing(fens, fes; method = :laplace, npass = 5)
println("count(fes) = $(count(fes))")

File = "voxel_bracket_mesh_tet.vtk"
vtkexportmesh(File, fens, fes)
@async run(`"paraview.exe" $File`)
