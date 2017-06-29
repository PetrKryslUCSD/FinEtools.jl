using JFinEALE
using VoxelBox


V=VoxelBoxVolume(Int,23*[5,6,7],[4.0, 4.0, 5.0])
#show(V)

s1(x,y,z)=sphereobj(x,y,z,1.5,2.0,2.0,1.3)
s2(x,y,z)=sphereobj(x,y,z,1.5,1.0,2.0,1.3)
b1(x,y,z)=boxobj(x,y,z,1.5,1.0,2.0,1.3,1.3,1.3)
b2(x,y,z)=boxobj(x,y,z,1.5-1.3,1.0,2.0-1.3,1.3,1.3,1.3)
f(x,y,z)=max(max(s1(x,y,z),-s2(x,y,z)),-(b1(x,y,z)*b2(x,y,z)))
fillvolume!(V,f,1,0)
M,N,P=size(V.data)
fens,fes = H8voximg(V.data,V.dimensions./[M,N,P],[1,1])

bfes = meshboundary(fes)
fens,bfes = meshsmoothing(fens,bfes; npass=13, method=:taubin)

File =  "B4.vtk"
vtkexportmesh (File, bfes.conn, fens.xyz, JFinEALE.MeshExportModule.Q4)
 @async run(`"C:/Program Files (x86)/ParaView 4.3.1/bin/paraview.exe" $File`)

File =  "B8.vtk"
vtkexportmesh (File, fes.conn, fens.xyz, JFinEALE.MeshExportModule.H8)
 @async run(`"C:/Program Files (x86)/ParaView 4.3.1/bin/paraview.exe" $File`)

