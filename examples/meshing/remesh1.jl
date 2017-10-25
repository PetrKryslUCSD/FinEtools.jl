
using FinEtools
using FinEtools.MeshExportModule
using Main.TetRemeshingModule
L= 0.3; 
W = 0.3;
a = 0.15;
nL=6; nW=6; na=6;

fens,fes = T4block(a,L,W,nL,nW,na,:a);
t = deepcopy(fes.conn);
v = deepcopy(fens.xyz);
tmid = ones(Int, size(t,1));

desired_ts =a;
bfes = meshboundary(fes);
f = connectednodes(bfes);
bv = zeros(Bool, size(v,1));
bv[f] = true;

println("Mesh size: initial = $(size(t,1))")

t, v, tmid = TetRemeshingModule.coarsen(t, v, tmid; bv = bv, desired_ts = desired_ts);

println("Mesh size: final = $(size(t,1))")

# File = "test1.vtk"
# MeshExportModule.vtkexportmesh(File, t, v, MeshExportModule.T4)
# @async run(`"paraview.exe" $File`)