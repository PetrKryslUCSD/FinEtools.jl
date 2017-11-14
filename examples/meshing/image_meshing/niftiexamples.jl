using NIfTI
using FinEtools

# file = niread("example4d.nii", mmap=false)
# display(file.header)
# voxel_size(file.header)
# size(file)
# vb = VoxelBoxVolume(file.raw[:, :, :, 1], voxel_size(file.header))
# vtkexport("example4d.vtk", vb)

# file = niread("bigger_subset_MW_Ear_R_rotY-thinner-5mat-fin.nii")
# display(file.header)
# voxel_size(file.header)
# size(file)
# vb = VoxelBoxVolume(file.raw[:, :, :, 1], voxel_size(file.header))
# vtkexport("bigger_subset_MW_Ear_R_rotY-thinner-5mat-fin.vtk", vb)

File = "Head_256x256x126-256x256x252"
file = niread(File * ".nii")
display(file.header)
voxel_size(file.header)
size(file)
vb = VoxelBoxVolume(file.raw[:, :, :, 1], voxel_size(file.header))
vtkexport(File * ".vtk", vb)