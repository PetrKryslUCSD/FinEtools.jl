module mmmmmiscellaneous2mmmmmm
using FinEtools
using Base.Test
function test()
  rho=1.21*1e-9;# mass density
  c =345.0*1000;# millimeters per second
  bulk= c^2*rho;
  Lx=1900.0;# length of the box, millimeters
  Ly=800.0; # length of the box, millimeters

  fens,fes = Q4block(Lx,Ly,3,2); # Mesh
  # show(fes.conn)

  bfes = meshboundary(fes)
  @test bfes.conn == [1 2; 5 1; 2 3; 3 4; 4 8; 9 5; 8 12; 10 9; 11 10; 12 11]
end
end
using mmmmmiscellaneous2mmmmmm
mmmmmiscellaneous2mmmmmm.test()

module mmmQ4blockneous2mmmmmm
using FinEtools
using FinEtools.MeshExportModule
using Base.Test
function test()
  rho=1.21*1e-9;# mass density
  c =345.0*1000;# millimeters per second
  bulk= c^2*rho;
  Lx=1900.0;# length of the box, millimeters
  Ly=800.0; # length of the box, millimeters

  fens,fes = Q8block(Lx,Ly,3,2); # Mesh
  # show(fes.conn)
  fens.xyz = xyz3(fens)
  fens.xyz[:, 3] .= (fens.xyz[:, 1].^2 + fens.xyz[:, 2].^2)/1.e3
  File = "mesh.vtk"
  MeshExportModule.vtkexportmesh(File, fens, fes)
  rm(File)
  # @async run(`"paraview.exe" $File`)
end
end
using mmmQ4blockneous2mmmmmm
mmmQ4blockneous2mmmmmm.test()

module mmmmmiimportexportmmmmm
using FinEtools
using FinEtools.MeshImportModule
using FinEtools.MeshExportModule
using Base.Test
function test()
  fens, fes = MeshImportModule.import_NASTRAN(dirname(@__FILE__) * "/" * "Slot-coarser.nas")
  # show(fes.conn[count(fes), :])
  File = "Slot-coarser.vtk"
  MeshExportModule.vtkexportmesh(File, fens, fes)
  rm(File)
  @test fes.conn[count(fes), :] == [143, 140, 144, 138, 361, 363, 176, 519, 781, 520]
  # @async run(`"paraview.exe" $File`)
end
end
using mmmmmiimportexportmmmmm
mmmmmiimportexportmmmmm.test()

module mmmLshapemmmmmmm
using FinEtools
using Base.Test
function test()
  # println("""
  # Meshing  for an L-shaped membrane.
  # """)

  t0 = time()

  W= 100. # strip width
  L = 200. # length of the strip
  nL = 5
  nW = 5
  tolerance = W/nW/1.0e5
  Meshes = Array{Tuple{FENodeSet, FESet},1}()
  push!(Meshes, Q4quadrilateral([0.0 0.0; W W], nW, nW))
  push!(Meshes, Q4quadrilateral([-L 0.0; 0.0 W], nL, nW))
  push!(Meshes, Q4quadrilateral([0.0 -L; W 0.0], nW, nL))
  fens, outputfes = mergenmeshes(Meshes, tolerance);
  fes = cat(outputfes[1], cat(outputfes[2], outputfes[3]))
  @test count(fens) == 96
  @test norm(vec(fes.conn[1,:]) - vec([1   2   8   7])) < 0.01
  geom = NodalField(fens.xyz)

  File =  "L_shape.vtk"
  vtkexportmesh(File, fes.conn, geom.values, FinEtools.MeshExportModule.Q4);
  # @async run(`"paraview.exe" $File`)
  rm(File)

  # println("Done")
  true
end
end
using mmmLshapemmmmmmm
mmmLshapemmmmmmm.test()
