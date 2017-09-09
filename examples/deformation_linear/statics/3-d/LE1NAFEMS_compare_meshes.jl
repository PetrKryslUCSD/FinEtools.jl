using FinEtools
using FinEtools.MeshExportModule
using FinEtools.MeshImportModule: import_ABAQUS

Thick0 = 0.1*phun("m")/2.0 # to account for the symmetry reduction

ref  = 0
Thickness = Thick0
tolerance = Thickness/2^ref/300.; # Geometrical tolerance


fens,fes = T10block(1.0, pi/2, Thickness, 2^ref*5, 2^ref*6, 1; orientation = :b)
for i=1:count(fens)
    t=fens.xyz[i,1]; a=fens.xyz[i,2]; z=fens.xyz[i,3]
    fens.xyz[i,:]=[(t*3.25+(1-t)*2)*cos(a), (t*2.75+(1-t)*1)*sin(a), z];
end
println("$((count(fens), count(fes)))")

output = import_ABAQUS("LE1AbaqusExport-C3D10HS-5-6-1.inp")
fens1, fes1 = output["fens"], output["fesets"][1]
println("$((count(fens1), count(fes1[1])))")

 fens, newfes1, fes2 = mergemeshes(fens,fes, fens1,fes1[1], tolerance)
 # fes = cat(fes2, newfes1)
 # println("$((count(fens), count(fes)))")

File =  "a1.vtk"
vtkexportmesh(File, newfes1.conn, fens.xyz,
               FinEtools.MeshExportModule.T10)
@async run(`"paraview.exe" $File`)
File =  "a2.vtk"
vtkexportmesh(File, fes2.conn, fens.xyz,
               FinEtools.MeshExportModule.T10)
@async run(`"paraview.exe" $File`)


#
# fens,fes = H8block(1.0, pi/2, Thickness, 2^ref*5, 2^ref*6, 1)
# for i=1:count(fens)
#     t=fens.xyz[i,1]; a=fens.xyz[i,2]; z=fens.xyz[i,3]
#     fens.xyz[i,:]=[(t*3.25+(1-t)*2)*cos(a), (t*2.75+(1-t)*1)*sin(a), z];
# end
# println("$((count(fens), count(fes)))")
#
# fens1,fes1 = import_ABAQUS("LE1AbaqusExport-C3D8S-80-96-1.inp")
# println("$((count(fens1), count(fes1[1])))")
#
#  fens, newfes1, fes2 = mergemeshes(fens,fes, fens1,fes1[1], tolerance)
#  fes = cat(fes2, newfes1)
#  println("$((count(fens), count(fes)))")
#
# File =  "a.vtk"
# vtkexportmesh(File, fes.conn, fens.xyz,
#                FinEtools.MeshExportModule.H8)
# @async run(`"paraview.exe" $File`)


true
