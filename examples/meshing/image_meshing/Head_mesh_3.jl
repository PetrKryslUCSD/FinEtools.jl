using FinEtools
using FinEtools.TetRemeshingModule
using FinEtools
using FinEtools.VoxelBoxModule
using FinEtools.VoxelTetMeshingModule
using FinEtools.MeshExportModule
using NIfTI

File = "Head_256x256x126-256x256x252"
file = niread(File * ".nii")
bs = voxel_size(file.header) .* [size(file,1)-100, size(file,2)-100, size(file,3)-100]
vb = VoxelBoxVolume(file.raw[101:end, 101:end, 101:end, 1], bs)
# vtkexport(File * ".vtk", vb)

im = ImageMesher(vb, zero(eltype(vb.data)), eltype(vb.data)[255])
mesh!(im)
println("Mesh size: initial = $(size(im.t,1))")
im.elementsizeweightfunctions = [ElementSizeWeightFunction(20.0, vec([112.0, 192.5, 26.5]), 45.0), ElementSizeWeightFunction(20.0, vec([157.0, 160.5, 0.5]), 45.0), ElementSizeWeightFunction(20.0, vec([56.0, 160.5, 0.5]), 45.0)]
for i = 1:5
    mesh!(im, 1.2)
    println("Mesh size: intermediate = $(size(im.t,1))")
end

fens = FENodeSet(im.v)
fes = FESetT4(im.t)
setlabel!(fes, im.tmid)

File = "Head_mesh_3.vtk"
vtkexportmesh(File, fens, fes)
@async run(`"paraview.exe" $File`)

ne = NASTRANExporter("Head_mesh_3.nas")
BEGIN_BULK(ne)
for i = 1:count(fens)
    GRID(ne, i, fens.xyz[i, :])
end
for i = 1:count(fes)
    CTETRA(ne, i, 1, fes.conn[i, :])
end
PSOLID(ne, 1, 1)
MAT1(ne, 1, 20.0e9, 0.3, 1800.0)
ENDDATA(ne)
close(ne)
