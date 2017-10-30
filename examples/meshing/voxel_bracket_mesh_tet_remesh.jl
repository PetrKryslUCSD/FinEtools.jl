using FinEtools
using FinEtools.VoxelTetMeshingModule

V = VoxelBoxVolume(Int, 8*[5,6,7], [4.0, 4.0, 5.0])

b1 = solidbox((0.0, 0.0, 0.0), (1.0, 4.0, 5.0))
b2 = solidbox((0.0, 0.0, 0.0), (4.0, 1.0, 5.0))
h1 = solidcylinder((2.0, 2.5, 2.5), (1.0, 0.0, 0.0), 0.75)
fillsolid!(V, differenceop(unionop(b1, b2), h1), 1)

im = ImageMesher(V, zero(eltype(V.data)), eltype(V.data)[1])
mesh!(im)
println("Mesh size: initial = $(size(im.t,1))")

im.elementsizeweightfunctions = [ElementSizeWeightFunction(20.0, vec([0.0, 2.5, 2.5]), 1.0), ElementSizeWeightFunction(1.0, vec([0.0, 2.5, 2.5]), 3.5)]
for i = 1:12
    mesh!(im, 1.1)
    println("Mesh size: final = $(size(im.t,1))")
end

fens = FENodeSet(im.v)
fes = FESetT4(im.t)
setlabel!(fes, im.tmid)

File = "voxel_bracket_mesh_tet.vtk"
vtkexportmesh(File, fens, fes)
@async run(`"paraview.exe" $File`)
