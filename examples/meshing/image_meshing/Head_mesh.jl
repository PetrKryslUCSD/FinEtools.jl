using FinEtools
using FinEtools.TetRemeshingModule
using NIfTI

File = "Head_256x256x126-256x256x252"
file = niread(File * ".nii")
bs = voxel_size(file.header) .* [size(file,1), size(file,2), size(file,3)-119]
vb = VoxelBoxVolume(file.raw[:, :, 120:end, 1], bs)
vtkexport(File * ".vtk", vb)

fens, fes = T4voximg(vb.data, vec([voxeldims(vb)...]), eltype(vb.data)[255])
t = deepcopy(fes.conn);
v = deepcopy(fens.xyz);
tmid = ones(Int, size(t,1));

desired_ts = 4 * voxel_size(file.header)[1];
bfes = meshboundary(fes);
f = connectednodes(bfes);
bv = zeros(Bool, size(v,1));
# bv[f] = true;

println("Mesh size: initial = $(size(t,1))")
t0 = time()

t, v, tmid = TetRemeshingModule.coarsen(t, v, tmid; bv = bv, desired_ts = desired_ts);

println("Mesh size: final = $(size(t,1)) [$(time() - t0) sec]")

fens.xyz = deepcopy(v)
fes.conn = deepcopy(t)
setlabel!(fes, tmid)

fens = meshsmoothing(fens, fes; npass = 5)
println("count(fes) = $(count(fes))")

File = "Head_mesh.vtk"
vtkexportmesh(File, fens, fes)
@async run(`"paraview.exe" $File`)
