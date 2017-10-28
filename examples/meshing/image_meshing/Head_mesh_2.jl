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
mesh!(im)
println("Mesh size: final = $(size(im.t,1))")

fens = FENodeSet(im.v)
fes = FESetT4(im.t)
setlabel!(fes, im.tmid)

File = "Head_mesh.vtk"
vtkexportmesh(File, fens, fes)
@async run(`"paraview.exe" $File`)
