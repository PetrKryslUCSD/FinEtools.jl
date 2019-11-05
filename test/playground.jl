module mT4quartercylnm
using FinEtools
using FinEtools.MeshExportModule: VTK
using Test
function test()
	Radius, Length, nperradius, nL = 3.0, 5.0, 3, 5
	fens,fes = T4quartercyln(Radius, Length, nperradius, nL)
	@test  count(fes) == 270
	File = "mesh.vtk"
	VTK.vtkexportmesh(File, fens, fes)
	rm(File)
	# @async run(`"paraview.exe" $File`)
	bfes = meshboundary(fes)
	@test  count(bfes) == 156
	File = "bmesh.vtk"
	VTK.vtkexportmesh(File, fens, fes)
	rm(File)
	# @async run(`"paraview.exe" $File`)
end
end
using .mT4quartercylnm
mT4quartercylnm.test()