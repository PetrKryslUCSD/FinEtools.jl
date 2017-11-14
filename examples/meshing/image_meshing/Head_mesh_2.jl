using FinEtools
using FinEtools.TetRemeshingModule
using NIfTI
using FinEtools.VoxelTetMeshingModule

File = "Head_256x256x126-256x256x252"
file = niread(File * ".nii")
bs = voxel_size(file.header) .* [size(file,1), size(file,2), size(file,3)-100]
vb = VoxelBoxVolume(file.raw[:, :, 101:end, 1], bs)
# vtkexport(File * ".vtk", vb)

im = ImageMesher(vb, zero(eltype(vb.data)), eltype(vb.data)[255])
mesh!(im)
println("Mesh size: initial = $(size(im.t,1))")
im.elementsizeweightfunctions = [ElementSizeWeightFunction(20.0, vec([112.0, 192.5, 26.5]), 45.0), ElementSizeWeightFunction(20.0, vec([157.0, 160.5, 0.5]), 45.0), ElementSizeWeightFunction(20.0, vec([56.0, 160.5, 0.5]), 45.0)]
for i = 1:12
    mesh!(im, 1.2)
    println("Mesh size: intermediate = $(size(im.t,1))")
end

fens = FENodeSet(im.v)
fes = FESetT4(im.t)
setlabel!(fes, im.tmid)

File = "Head_mesh.vtk"
vtkexportmesh(File, fens, fes)
@async run(`"paraview.exe" $File`)


