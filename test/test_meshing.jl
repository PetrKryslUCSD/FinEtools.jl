module miscellaneous2mm
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
using miscellaneous2mm
miscellaneous2mm.test()

module mmmQ4blockneous2mm
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
using mmmQ4blockneous2mm
mmmQ4blockneous2mm.test()

module miimportexportm
using FinEtools
using FinEtools.MeshImportModule
using FinEtools.MeshExportModule
using Base.Test
function test()
  output = MeshImportModule.import_NASTRAN(dirname(@__FILE__) * "/" * "Slot-coarser.nas";
    allocationchunk = 13)
  # show(fes.conn[count(fes), :])
  File = "Slot-coarser.vtk"
  MeshExportModule.vtkexportmesh(File, output["fens"], output["fesets"][1])
  rm(File)
  @test output["fesets"][1].conn[count(output["fesets"][1]), :] ==
    [143, 140, 144, 138, 361, 363, 176, 519, 781, 520]
  # @async run(`"paraview.exe" $File`)
end
end
using miimportexportm
miimportexportm.test()

module miimportexportm2
using FinEtools
using FinEtools.MeshImportModule
using FinEtools.MeshExportModule
using Base.Test
function test()
  output = MeshImportModule.import_NASTRAN(dirname(@__FILE__) * "/" * "Slot-coarser-2.nas";
    allocationchunk = 13)
  # show(fes.conn[count(fes), :])
  File = "Slot-coarser.vtk"
  MeshExportModule.vtkexportmesh(File, output["fens"], output["fesets"][1])
  rm(File)
  @test output["fesets"][1].conn[count(output["fesets"][1]), :] ==
    [143, 140, 144, 138, 361, 363, 176, 519, 781, 520]
  # @async run(`"paraview.exe" $File`)
end
end
using miimportexportm2
miimportexportm2.test()

module mmmLshapemmm
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
using mmmLshapemmm
mmmLshapemmm.test()

module mAbaqusmmiimportmmm
using FinEtools
using FinEtools.MeshExportModule
using FinEtools.MeshImportModule
using Base.Test
function test()
  ## Solid cylinder/taper/sphere—-temperature loading; quadratic brick mesh

  # The mesh  will be created in a very coarse representation from the
  # key points in the drawing. The first coordinate is radial, the second coordinate is axial.
  rz=[1.     0.;#A
  1.4    0.;#B
  0.995184726672197   0.098017140329561;
  1.393258617341076 0.137223996461385;
  0.980785  0.195090;#
  1.37309939 0.27312645;
  0.956940335732209   0.290284677254462
  1.339716470025092 0.406398548156247
  0.9238795  0.38268;#C
  1.2124  0.7;#D
  0.7071  0.7071;#E
  1.1062  1.045;#F
  0.7071  (0.7071+1.79)/2;#(E+H)/2
  1.      1.39;#G
  0.7071  1.79;#H
  1.      1.79;#I
  ]*phun("M")
  tolerance =1.e-6*phun("M")

  # This is the quadrilateral mesh of the cross-section.   It will be modified and
  # refined as  we go.
  fens = FENodeSet(rz);
  fes = FESetQ4([1 2 4 3; 3 4 6 5; 5 6 8 7; 7 8 10 9; 9 10 12 11; 11 12 14 13; 13 14 16 15]);

  nref = 1;
  for ref = 1:nref
    fens,fes = Q4refine(fens,fes);
    list = selectnode(fens, distance=1.0+0.1/2^nref, from=[0. 0.], inflate=tolerance);
    fens.xyz[list,:] = FinEtools.MeshUtilModule.ontosphere(fens.xyz[list,:],1.0);
  end

  ##
  # The mesh is extruded by sweeping around the axis of symmetry.
  # Only a single layer of elements is generated of internal angle
  # |angslice|.
  nLayers = 7;
  angslice  = 5*pi/16;

  ##
  # First the mesh is extruded to a block whose third dimension
  # represents the angular coordinate.
  fens,fes = H8extrudeQ4(fens, fes, nLayers,
  (rz,k)->[rz[1],rz[2],0.0]-(k)/nLayers*[0.,0.,angslice]);
  # The block is now converted  to the axially symmetric geometry by using the
  # third (angular) coordinate  to sweep out  an axially symmetric domain. The
  # ccoordinates of the nodes at this point are |rza|,  radial distance,
  # Z-coordinate, angle.
  sweep(rza) = [-rza[1]*sin(rza[3]+angslice/2.0), rza[1]*cos(rza[3]+angslice/2.0), rza[2]]
  for j=1:size(fens.xyz,1)
    fens.xyz[j,:] = sweep(fens.xyz[j,:])
  end

  AE = AbaqusExporter("LE11NAFEMS_H8");
  HEADING(AE, "LE11NAFEMS: Linear bricks.");
  PART(AE, "part1");
  END_PART(AE);
  ASSEMBLY(AE, "ASSEM1");
  INSTANCE(AE, "INSTNC1", "PART1");
  NODE(AE, fens.xyz);
  ELEMENT(AE, "c3d8rh", "AllElements", 1, fes.conn)
  END_INSTANCE(AE);
  END_ASSEMBLY(AE);
  close(AE)


  output = MeshImportModule.import_ABAQUS("./LE11NAFEMS_H8.inp";
    allocationchunk = 11)
    fens, fes = output["fens"], output["fesets"][1]

  File = "LE11NAFEMS_H8.vtk"
  MeshExportModule.vtkexportmesh(File, fens, fes)
  # @async run(`"paraview.exe" $File`)
  try rm(File) catch end

  try rm(AE.filename) catch end
end
end
using mAbaqusmmiimportmmm
mAbaqusmmiimportmmm.test()

module mmsmoothingm1
using FinEtools
using Base.Test
function test()

    # println("""
    # Meshing, deforming  and smoothing
    # """)

    A= 100. # strip width
    N = 16
    tolerance = A/N/1.0e5
    fens,fes = T3block(A, A, N, N)

    bnl = connectednodes(meshboundary(fes))
    for ixxxx = 1:length(bnl)
        x, y = fens.xyz[bnl[ixxxx], :]
        fens.xyz[bnl[ixxxx], 1] += A/N*sin(2*pi*y/A)
        fens.xyz[bnl[ixxxx], 2] += -A/N*sin(2*pi*x/A)
    end

    File =  "mesh_smoothing_before.vtk"
    vtkexportmesh(File, fens, fes);
    # @async run(`"paraview.exe" $File`)
    try rm(File) catch end
    before = [50.0, 43.75]
    # println("$(fens.xyz[Int(N^2/2), :] )")
    @test norm(before-fens.xyz[Int(N^2/2), :]) < 1.e-4

    fixedv = falses(count(fens))
    fixedv[bnl] = true
    fens = meshsmoothing(fens, fes; fixedv = fixedv, method = :laplace, npass = 100)

    after = [50.0438, 44.0315]
    # println("$(fens.xyz[Int(N^2/2), :] )")
    @test norm(fens.xyz[Int(N^2/2), :]-after) < 1.e-4

    geom = NodalField(fens.xyz)

    File =  "mesh_smoothing_after.vtk"
    vtkexportmesh(File, fes.conn, geom.values, FinEtools.MeshExportModule.T3);
    # @async run(`"paraview.exe" $File`)
    try rm(File) catch end

    # println("Done")
    true
end
end
using mmsmoothingm1
mmsmoothingm1.test()


module mmsmoothingm2
using FinEtools
using Base.Test
function test()

    # println("""
    # Meshing, deforming  and smoothing
    # """)

    A= 100. # strip width
    N = 16
    tolerance = A/N/1.0e5
    fens,fes = T3block(A, A, N, N)

    bnl = connectednodes(meshboundary(fes))
    for ixxxx = 1:length(bnl)
        x, y = fens.xyz[bnl[ixxxx], :]
        fens.xyz[bnl[ixxxx], 1] += A/N*sin(2*pi*y/A)
        fens.xyz[bnl[ixxxx], 2] += -A/N*sin(2*pi*x/A)
    end

    File =  "mesh_smoothing_before.vtk"
    vtkexportmesh(File, fens, fes);
    # @async run(`"paraview.exe" $File`)
    try rm(File) catch end
    before = [50.0, 43.75]
    # println("$(fens.xyz[Int(N^2/2), :] )")
    @test norm(before-fens.xyz[Int(N^2/2), :]) < 1.e-4

    fixedv = falses(count(fens))
    fixedv[bnl] = true
    fens = meshsmoothing(fens, fes; fixedv = fixedv, method = :taubin, npass = 100)

    after = [50.0059, 43.6281]
    # println("$(fens.xyz[Int(N^2/2), :] )")
    @test norm(fens.xyz[Int(N^2/2), :]-after) < 1.e-4

    geom = NodalField(fens.xyz)

    File =  "mesh_smoothing_after.vtk"
    vtkexportmesh(File, fes.conn, geom.values, FinEtools.MeshExportModule.T3);
    # @async run(`"paraview.exe" $File`)
    try rm(File) catch end

    # println("Done")
    true
end
end
using mmsmoothingm2
mmsmoothingm2.test()

module mmbracketmm
using FinEtools
using Base.Test
function test()

    V = VoxelBoxVolume(Int, 6*[5,6,7], [4.0, 4.0, 5.0])

    b1 = solidbox((0.0, 0.0, 0.0), (1.0, 4.0, 5.0))
    b2 = solidbox((0.0, 0.0, 0.0), (4.0, 1.0, 5.0))
    h1 = solidcylinder((2.0, 2.5, 2.5), (1.0, 0.0, 0.0), 00.75)
    fillsolid!(V, differenceop(unionop(b1, b2), h1), 1)

    fens, fes = H8voximg(V.data, vec([voxeldims(V)...]), [1])
    fens = meshsmoothing(fens, fes; npass = 5)
    # println("count(fes) = $(count(fes))")
    @test count(fes) == 16337

    File = "voxel_bracket_mesh.vtk"
    vtkexportmesh(File, fens, fes)
    # @async run(`"paraview.exe" $File`)
    try rm(File) catch end
end
end
using mmbracketmm
mmbracketmm.test()

module mmmMergem1
using FinEtools
using Base.Test
function test()
    h = 0.05*phun("M");
    l = 10*h;
    Rmed = h/0.2;
    psi   = 0;    # Cylinder
    nh = 12; nl  = 40; nc = 120;
    nh = 6; nl  = 20; nc = 60;
    # nh = 3; nl  = 8; nc = 30;
    tolerance = h/nh/100;

    t0 = time()
    MR = DeforModelRed3D
    fens,fes  = H8block(h,l,2.0*pi,nh,nl,nc)
    # Shape into a cylinder
    R = zeros(3, 3)
    for i = 1:count(fens)
        x, y, z = fens.xyz[i,:];
        rotmat3!(R, [0, z, 0])
        Q = [cos(psi*pi/180) sin(psi*pi/180) 0;
            -sin(psi*pi/180) cos(psi*pi/180) 0;
            0 0 1]
        fens.xyz[i,:] = reshape([x+Rmed-h/2, y-l/2, 0], 1, 3)*Q*R;
    end
    # println("  before merging  = $(count(fens))")
    @test count(fens) == 8967
    # File =  "unit_cube_modes.vtk"
    # vtkexportmesh(File, fens, fes)
    # @async run(`"paraview.exe" $File`)

    candidates = selectnode(fens, box = boundingbox([Rmed-h -Inf 0.0; Rmed+h +Inf 0.0]), inflate = tolerance)
    fens,fes = mergenodes(fens, fes,  tolerance, candidates);
    @test count(fens) == 8820

    # fens,fes = mergenodes(fens, fes,  tolerance);
    # println("  after merging  = $(count(fens))")
    #
    # println("Mesh generation ($(time() - t0) sec)")
    #
    # File =  "unit_cube_modes.vtk"
    # vtkexportmesh(File, fens, fes)
    # @async run(`"paraview.exe" $File`)
end
end
using mmmMergem1
mmmMergem1.test()


module mmmMergem2
using FinEtools
using Base.Test
function test()
    h = 0.05*phun("M");
    l = 10*h;
    Rmed = h/0.2;
    psi   = 0;    # Cylinder
    nh = 12; nl  = 40; nc = 120;
    nh = 6; nl  = 20; nc = 60;
    # nh = 3; nl  = 8; nc = 30;
    tolerance = h/nh/100;

    t0 = time()
    MR = DeforModelRed3D
    fens,fes  = H8block(h,l,2.0*pi,nh,nl,nc)
    # Shape into a cylinder
    R = zeros(3, 3)
    for i = 1:count(fens)
        x, y, z = fens.xyz[i,:];
        rotmat3!(R, [0, z, 0])
        Q = [cos(psi*pi/180) sin(psi*pi/180) 0;
            -sin(psi*pi/180) cos(psi*pi/180) 0;
            0 0 1]
        fens.xyz[i,:] = reshape([x+Rmed-h/2, y-l/2, 0], 1, 3)*Q*R;
    end
    # println("  before merging  = $(count(fens))")
    @test count(fens) == 8967
    # File =  "unit_cube_modes.vtk"
    # vtkexportmesh(File, fens, fes)
    # @async run(`"paraview.exe" $File`)

    candidates = selectnode(fens, plane = [0.0 0.0 1.0 0.0], inflate = tolerance)
    fens,fes = mergenodes(fens, fes,  tolerance, candidates);
    @test count(fens) == 8820

    # fens,fes = mergenodes(fens, fes,  tolerance);
    # println("  after merging  = $(count(fens))")
    #
    # println("Mesh generation ($(time() - t0) sec)")
    #
    # File =  "unit_cube_modes.vtk"
    # vtkexportmesh(File, fens, fes)
    # @async run(`"paraview.exe" $File`)
end
end
using mmmMergem2
mmmMergem2.test()


module mmmMergem3
using FinEtools
using Base.Test
function test()
    h = 0.05*phun("M");
    l = 10*h;
    Rmed = h/0.2;
    psi   = 0;    # Cylinder
    nh = 12; nl  = 40; nc = 120;
    nh = 6; nl  = 20; nc = 60;
    # nh = 3; nl  = 8; nc = 30;
    tolerance = h/nh/100;

    t0 = time()
    MR = DeforModelRed3D
    fens,fes  = H8block(h,l,2.0*pi,nh,nl,nc)
    # Shape into a cylinder
    R = zeros(3, 3)
    for i = 1:count(fens)
        x, y, z = fens.xyz[i,:];
        rotmat3!(R, [0, z, 0])
        Q = [cos(psi*pi/180) sin(psi*pi/180) 0;
            -sin(psi*pi/180) cos(psi*pi/180) 0;
            0 0 1]
        fens.xyz[i,:] = reshape([x+Rmed-h/2, y-l/2, 0], 1, 3)*Q*R;
    end
    # println("  before merging  = $(count(fens))")
    @test count(fens) == 8967
    # File =  "unit_cube_modes.vtk"
    # vtkexportmesh(File, fens, fes)
    # @async run(`"paraview.exe" $File`)

    candidates = selectnode(fens, plane = [0.0 0.0 1.0 0.0], thickness = h/1000)
    fens,fes = mergenodes(fens, fes,  tolerance, candidates);
    @test count(fens) == 8820

    # fens,fes = mergenodes(fens, fes,  tolerance);
    # println("  after merging  = $(count(fens))")
    #
    # println("Mesh generation ($(time() - t0) sec)")
    #
    # File =  "unit_cube_modes.vtk"
    # vtkexportmesh(File, fens, fes)
    # @async run(`"paraview.exe" $File`)
end
end
using mmmMergem3
mmmMergem3.test()


module mAbaqusmmiimport_1m
using FinEtools
using FinEtools.MeshExportModule
using FinEtools.MeshImportModule
using Base.Test
function test()



  ## Solid cylinder/taper/sphere—-temperature loading; quadratic brick mesh

  # The mesh  will be created in a very coarse representation from the
  # key points in the drawing. The first coordinate is radial, the second coordinate is axial.
  rz=[1.     0.;#A
  1.4    0.;#B
  0.995184726672197   0.098017140329561;
  1.393258617341076 0.137223996461385;
  0.980785  0.195090;#
  1.37309939 0.27312645;
  0.956940335732209   0.290284677254462
  1.339716470025092 0.406398548156247
  0.9238795  0.38268;#C
  1.2124  0.7;#D
  0.7071  0.7071;#E
  1.1062  1.045;#F
  0.7071  (0.7071+1.79)/2;#(E+H)/2
  1.      1.39;#G
  0.7071  1.79;#H
  1.      1.79;#I
  ]*phun("M")
  tolerance =1.e-6*phun("M")

  # This is the quadrilateral mesh of the cross-section.   It will be modified and
  # refined as  we go.
  fens = FENodeSet(rz);
  fes = FESetQ4([1 2 4 3; 3 4 6 5; 5 6 8 7; 7 8 10 9; 9 10 12 11; 11 12 14 13; 13 14 16 15]);

  nref = 1;
  for ref = 1:nref
    fens,fes = Q4refine(fens,fes);
    list = selectnode(fens, distance=1.0+0.1/2^nref, from=[0. 0.], inflate=tolerance);
    fens.xyz[list,:] = FinEtools.MeshUtilModule.ontosphere(fens.xyz[list,:],1.0);
  end

  ##
  # The mesh is extruded by sweeping around the axis of symmetry.
  # Only a single layer of elements is generated of internal angle
  # |angslice|.
  nLayers = 7;
  angslice  = 5*pi/16;

  ##
  # First the mesh is extruded to a block whose third dimension
  # represents the angular coordinate.
  fens,fes = H8extrudeQ4(fens, fes, nLayers,
  (rz,k)->[rz[1],rz[2],0.0]-(k)/nLayers*[0.,0.,angslice]);
  # The block is now converted  to the axially symmetric geometry by using the
  # third (angular) coordinate  to sweep out  an axially symmetric domain. The
  # ccoordinates of the nodes at this point are |rza|,  radial distance,
  # Z-coordinate, angle.
  sweep(rza) = [-rza[1]*sin(rza[3]+angslice/2.0), rza[1]*cos(rza[3]+angslice/2.0), rza[2]]
  for j=1:size(fens.xyz,1)
    fens.xyz[j,:] = sweep(fens.xyz[j,:])
  end

  AE = AbaqusExporter("LE11NAFEMS_H8");
  HEADING(AE, "LE11NAFEMS: Linear bricks.");
  PART(AE, "part1");
  END_PART(AE);
  ASSEMBLY(AE, "ASSEM1");
  INSTANCE(AE, "INSTNC1", "PART1");
  NODE(AE, fens.xyz);
  ELEMENT(AE, "c3d8rh", "AllElements", 1, fes.conn)
  ELSET_ELSET(AE, "AllSolids", collect(1:count(fes)))
  END_INSTANCE(AE);
  END_ASSEMBLY(AE);
  close(AE)


  output = MeshImportModule.import_ABAQUS("./LE11NAFEMS_H8.inp";
    allocationchunk = 12)
    fens, fes = output["fens"], output["fesets"][1]

  File = "LE11NAFEMS_H8.vtk"
  MeshExportModule.vtkexportmesh(File, fens, fes)
  # @async run(`"paraview.exe" $File`)
  try rm(File) catch end
  try rm(AE.filename) catch end
end
end
using mAbaqusmmiimport_1m
mAbaqusmmiimport_1m.test()


module mmfflood1
using FinEtools
using Base.Test
function test()
    h = 0.05*phun("M");
    l = 10*h;
    Rmed = h/0.2;
    psi   = 0;    # Cylinder
    nh = 12; nl  = 40; nc = 120;
    nh = 6; nl  = 20; nc = 60;
    # nh = 3; nl  = 8; nc = 30;
    tolerance = h/nh/100;

    fens,fes  = H8block(h,l,2.0*pi,nh,nl,nc)
    # Shape into a cylinder
    R = zeros(3, 3)
    for i = 1:count(fens)
        x, y, z = fens.xyz[i,:];
        rotmat3!(R, [0, z, 0])
        Q = [cos(psi*pi/180) sin(psi*pi/180) 0;
            -sin(psi*pi/180) cos(psi*pi/180) 0;
            0 0 1]
        fens.xyz[i,:] = reshape([x+Rmed-h/2, y-l/2, 0], 1, 3)*Q*R;
    end
    # println("  before merging  = $(count(fens))")
    @test count(fens) == 8967
    # File =  "unit_cube_modes.vtk"
    # vtkexportmesh(File, fens, fes)
    # @async run(`"paraview.exe" $File`)

    candidates = selectnode(fens, box = boundingbox([Rmed-h -Inf 0.0; Rmed+h +Inf 0.0]), inflate = tolerance)
    fens,fes = mergenodes(fens, fes,  tolerance, candidates);
    @test count(fens) == 8820

    bfes = meshboundary(fes)
    startnode = bfes.conn[1,1]
    lb = selectelem(fens, bfes, flood=true, startnode=startnode)
    # println("$(length(lb))")
    @test length(lb) == 3120

    bbfes = meshboundary(subset(bfes, lb))
    @test count(bbfes) == 0
    # File =  "unit_cube_modes.vtk"
    # vtkexportmesh(File, fens, subset(bfes, lb))
    # @async run(`"paraview.exe" $File`)
    # try rm(File) catch end
end
end
using mmfflood1
mmfflood1.test()


module mmfflood2
using FinEtools
using Base.Test
function test()
    h = 0.05*phun("M");
    l = 10*h;
    Rmed = h/0.2;
    psi   = 0;    # Cylinder
    # nh = 12; nl  = 40; nc = 120;
    # nh = 6; nl  = 20; nc = 60;
    nh = 3; nl  = 8; nc = 30;
    nh = 2; nl  = 2; nc = 10;
    tolerance = h/nh/100;

    fens,fes  = H20block(h,l,2.0*pi,nh,nl,nc)
    # Shape into a cylinder
    R = zeros(3, 3)
    for i = 1:count(fens)
        x, y, z = fens.xyz[i,:];
        rotmat3!(R, [0, z, 0])
        Q = [cos(psi*pi/180) sin(psi*pi/180) 0;
            -sin(psi*pi/180) cos(psi*pi/180) 0;
            0 0 1]
        fens.xyz[i,:] = reshape([x+Rmed-h/2, y-l/2, 0], 1, 3)*Q*R;
    end
    # println("  before merging  = $(count(fens))")
    # @test count(fens) == 4025
    # File =  "unit_cube_modes.vtk"
    # vtkexportmesh(File, fens, fes)
    # @async run(`"paraview.exe" $File`)

    candidates = selectnode(fens, box = boundingbox([Rmed-h -Inf 0.0; Rmed+h +Inf 0.0]), inflate = tolerance)
    fens,fes = mergenodes(fens, fes,  tolerance, candidates);
    # @test count(fens) == 3930
    # @test count(fens) == 8820

    bfes = meshboundary(fes)
    startnode = bfes.conn[1,1]
    lb = selectelem(fens, bfes, flood=true, startnode=startnode)
    # println("$(lb)")
    @test length(lb) == 80

    bbfes = meshboundary(subset(bfes, lb))
    @test count(bbfes) == 0
    # File =  "unit_cube_modes.vtk"
    # vtkexportmesh(File, fens, subset(bfes, lb))
    # @async run(`"paraview.exe" $File`)
end
end
using mmfflood2
mmfflood2.test()

module mt4orientation2
using FinEtools
using FinEtools.MeshImportModule
using FinEtools.MeshExportModule
using Base.Test
function test()
    xs = collect(linspace(0.0, 2.0, 5))
    ys = collect(linspace(0.0, 1.0, 6).^2)
    zs = collect(linspace(0.0, 1.0, 3))
    fens, fes = T4blockx(xs, ys, zs, :a)
    @test count(fes) == 240
    fens, fes = T4blockx(xs, ys, zs, :b)
    @test count(fes) == 240
    fens, fes = T4blockx(xs, ys, zs, :ca)
    @test count(fes) == 200
    fens, fes = T4blockx(xs, ys, zs, :cb)
    @test count(fes) == 200
    # show(fes.conn[count(fes), :])
    # File = "Slot-coarser.vtk"
    # MeshExportModule.vtkexportmesh(File, fens, fes)
    # rm(File)
    # @test fes.conn[count(fes), :] == [143, 140, 144, 138, 361, 363, 176, 519, 781, 520]
    # @async run(`"paraview.exe" $File`)
end
end
using mt4orientation2
mt4orientation2.test()

module mmAbaqusexport3
using FinEtools
using FinEtools.MeshExportModule
using FinEtools.MeshImportModule
using Base.Test
function test()
  ## Solid cylinder/taper/sphere—-temperature loading; quadratic brick mesh

  # The mesh  will be created in a very coarse representation from the
  # key points in the drawing. The first coordinate is radial, the second coordinate is axial.
  rz=[1.     0.;#A
  1.4    0.;#B
  0.995184726672197   0.098017140329561;
  1.393258617341076 0.137223996461385;
  0.980785  0.195090;#
  1.37309939 0.27312645;
  0.956940335732209   0.290284677254462
  1.339716470025092 0.406398548156247
  0.9238795  0.38268;#C
  1.2124  0.7;#D
  0.7071  0.7071;#E
  1.1062  1.045;#F
  0.7071  (0.7071+1.79)/2;#(E+H)/2
  1.      1.39;#G
  0.7071  1.79;#H
  1.      1.79;#I
  ]*phun("M")
  tolerance =1.e-6*phun("M")

  # This is the quadrilateral mesh of the cross-section.   It will be modified and
  # refined as  we go.
  fens = FENodeSet(rz);
  fes = FESetQ4([1 2 4 3; 3 4 6 5; 5 6 8 7; 7 8 10 9; 9 10 12 11; 11 12 14 13; 13 14 16 15]);

  nref = 1;
  for ref = 1:nref
    fens,fes = Q4refine(fens,fes);
    list = selectnode(fens, distance=1.0+0.1/2^nref, from=[0. 0.], inflate=tolerance);
    fens.xyz[list,:] = FinEtools.MeshUtilModule.ontosphere(fens.xyz[list,:],1.0);
  end

  ##
  # The mesh is extruded by sweeping around the axis of symmetry.
  # Only a single layer of elements is generated of internal angle
  # |angslice|.
  nLayers = 7;
  angslice  = 5*pi/16;

  ##
  # First the mesh is extruded to a block whose third dimension
  # represents the angular coordinate.
  fens,fes = H8extrudeQ4(fens, fes, nLayers,
  (rz,k)->[rz[1],rz[2],0.0]-(k)/nLayers*[0.,0.,angslice]);
  # The block is now converted  to the axially symmetric geometry by using the
  # third (angular) coordinate  to sweep out  an axially symmetric domain. The
  # ccoordinates of the nodes at this point are |rza|,  radial distance,
  # Z-coordinate, angle.
  sweep(rza) = [-rza[1]*sin(rza[3]+angslice/2.0), rza[1]*cos(rza[3]+angslice/2.0), rza[2]]
  for j=1:size(fens.xyz,1)
    fens.xyz[j,:] = sweep(fens.xyz[j,:])
  end

  AE = AbaqusExporter("LE11NAFEMS_H8_B.inp");
  HEADING(AE, "LE11NAFEMS: Linear bricks.");
  PART(AE, "part1");
  END_PART(AE);
  ASSEMBLY(AE, "ASSEM1");
  INSTANCE(AE, "INSTNC1", "PART1");
  NODE(AE, fens.xyz);
  ELEMENT(AE, "c3d8rh", "AllElements", 1, fes.conn)
  END_INSTANCE(AE);
  END_ASSEMBLY(AE);
  STEP_PERTURBATION_BUCKLE(AE, 1)
  CLOAD(AE, 1, 2, 10.0)
  EL_PRINT(AE, "AllElements", "S")
  ENERGY_PRINT(AE);
  END_STEP(AE);
  close(AE)


  output = MeshImportModule.import_ABAQUS("./LE11NAFEMS_H8_B.inp";
    allocationchunk = 11)
    fens, fes = output["fens"], output["fesets"][1]

  File = "LE11NAFEMS_H8.vtk"
  MeshExportModule.vtkexportmesh(File, fens, fes)
  # @async run(`"paraview.exe" $File`)
  # try rm(File) catch end

end
end
using mmAbaqusexport3
mmAbaqusexport3.test()


module mttriangles13
using FinEtools
using FinEtools.MeshExportModule
using Base.Test
function test()
    xs = collect(linspace( 1.0, 3.0, 4))
    ys = collect(linspace(-1.0, 3.0, 4))
    fens, fes = T3blockx(xs, ys, :a)
    @test count(fes) == 3*3*2
    fens, fes = T3blockx(xs, ys, :b)
    @test count(fes) == 3*3*2

    # File = "playground.vtk"
    # MeshExportModule.vtkexportmesh(File, fens, fes)
    # @async run(`"paraview.exe" $File`)
    # try rm(File) catch end

end
end
using mttriangles13
mttriangles13.test()

module mttriangles14
using FinEtools
using FinEtools.MeshExportModule
using Base.Test
function test()
    xs = collect(linspace( 1.0, 3.0, 4))
    ys = collect(linspace(-1.0, 3.0, 5))
    fens, fes = Q4blockx(xs, ys)
    fens, fes = Q4toT3(fens, fes)
    @test count(fes) == 3*4*2
    fens, fes = Q4blockx(xs, ys)
    fens, fes = Q4toT3(fens, fes, :alternate)
    @test count(fes) == 3*4*2

    # File = "playground.vtk"
    # MeshExportModule.vtkexportmesh(File, fens, fes)
    # @async run(`"paraview.exe" $File`)
    # try rm(File) catch end

end
end
using mttriangles14
mttriangles14.test()

module mselecte1
using FinEtools

using Base.Test
function test()

    # println("""
    # The initially twisted cantilever beam is one of the standard test
    # problems for verifying the finite-element accuracy [1]. The beam is
    #   clamped at one end and loaded either with unit in-plane or
    #   unit out-of-plane force at the other. The centroidal axis of the beam is
    #   straight at the undeformed  configuration, while its cross-sections are
    #   twisted about the centroidal axis from 0 at the clamped end to pi/2 at
    #   the free end.
    #
    #   Reference:
    #   Zupan D, Saje M (2004) On "A proposed standard set of problems to test
    #   finite element accuracy": the twisted beam. Finite Elements in Analysis
    #   and Design 40: 1445-1451.
    #   """)
    W = 1.1;
    L = 12.;
    t =  0.32;
    nl = 2; nt = 1; nw = 1; ref = 7;
    tolerance  = t/1000;

    fens,fes  = H8block(L,W,t, nl*ref,nw*ref,nt*ref)

    # Reshape into a twisted beam shape
    for i = 1:count(fens)
        a = fens.xyz[i,1]/L*(pi/2); y = fens.xyz[i,2]-(W/2); z = fens.xyz[i,3]-(t/2);
        fens.xyz[i,:] = [fens.xyz[i,1],y*cos(a)-z*sin(a),y*sin(a)+z*cos(a)];
    end

    # Clamped end of the beam
    l1  = selectnode(fens; box = [0 0 -100*W 100*W W W], inflate  =  tolerance)
    # display(l1)
    @test isempty(l1)

    # # Traction on the opposite edge
    boundaryfes  =   meshboundary(fes);
    Toplist   = selectelem(fens,boundaryfes, facing = true,
    direction = [-1.0 0.0 0.0], dotmin = 0.999);
    # display(Toplist)
    @test length(Toplist) == 49
    Toplist   = selectelem(fens,boundaryfes, box = [L L -Inf Inf -Inf Inf],
    allin = true, inflate  =  tolerance);
    # display(Toplist)
    @test length(Toplist) == 49
    Toplist   = selectelem(fens,boundaryfes, box = [L L -Inf Inf -Inf Inf],
    allin = false, inflate  =  tolerance);
    # display(Toplist)
    @test length(Toplist) == 77

    # File =  "b.vtk"
    # vtkexportmesh(File, subset(boundaryfes, Toplist).conn, fens.xyz,
    #     FinEtools.MeshExportModule.Q4)
    # File =  "a.vtk"
    # vtkexportmesh(File, fes.conn, fens.xyz,  FinEtools.MeshExportModule.H8)
    # @async run(`"paraview.exe" $File`)
end
end
using mselecte1
mselecte1.test()


module mselecte2
using FinEtools

using Base.Test
function test()

    W = 10.0;
    L = 12.;
    nl = 2; nw = 1; ref = 7;
    tolerance  = W/1000;

    fens,fes  = Q4block(L, W, nl*ref, nw*ref)

    boundaryfes  =   meshboundary(fes);
    alist   = selectelem(fens,boundaryfes, facing = true,
    direction = [0.0 -1.0], dotmin = 0.999);
    # display(alist)
    @test length(alist) == 14

    # File =  "b.vtk"
    # vtkexportmesh(File, subset(boundaryfes, alist).conn, fens.xyz,
    #     FinEtools.MeshExportModule.L2)
    # File =  "a.vtk"
    # vtkexportmesh(File, fes.conn, fens.xyz,  FinEtools.MeshExportModule.Q4)
    # @async run(`"paraview.exe" $File`)
end
end
using mselecte2
mselecte2.test()

module mMeyer_Piening_1
using FinEtools
using FinEtools.MeshUtilModule
using Base.Test
function test()

    # Reference results from:
    # [1] Application of the Elasticity Solution
    # to Linear Sandwich Beam, Plate
    # and Shell Analyses
    # H.-R. MEYER -PIENING
    # Journal of SANDWICH STRUCTURES AND MATERIALS , Vol. 6—July 2004

    Lx = 5.0*phun("mm") # length  of loaded rectangle
    Ly = 20.0*phun("mm") # length  of loaded rectangle
    Sx = 100.0*phun("mm") # span of the plate
    Sy = 200.0*phun("mm") # span of the plate

    # Here we define the layout and the thicknesses of the layers.
    angles = vec([0.0 0.0 0.0]);
    ts = vec([0.5  11.4  0.1])*phun("mm"); # layer thicknesses
    TH = sum(ts); # total thickness of the plate

    tolerance = 0.0001*TH

    # Select how find the mesh should be
    Refinement = 5
    nL = Refinement * 1;
    nSx = nL + Refinement * 4;
    nSy = 2 * nSx;

    # Each layer is modeled with a single element.
    nts= Refinement * [1, 2, 1];# number of elements per layer
    strength = 1.5
    xs = unique(vcat(reverse(collect(MeshUtilModule.gradedspace(Lx/2, 0.0, nL+1, strength))),
        collect(MeshUtilModule.gradedspace(Lx/2, Sx/2, nSx-nL+1, strength))))
    ys = unique(vcat(reverse(collect(MeshUtilModule.gradedspace(Ly/2, 0.0, nL+1, strength))),
        collect(MeshUtilModule.gradedspace(Ly/2, Sy/2, nSy-nL+1, strength))))

    fens,fes = H8compositeplatex(xs, ys, ts, nts)
    @test count(fens) == 27846
    @test count(fes) == 25000
    rls = selectelem(fens, fes, label = 1)
    @test length(rls) == 6250
end
end
using mMeyer_Piening_1
mMeyer_Piening_1.test()


module mMeyer_Piening_2
using FinEtools
using FinEtools.MeshUtilModule
using Base.Test
function test()

    # Reference results from:
    # [1] Application of the Elasticity Solution
    # to Linear Sandwich Beam, Plate
    # and Shell Analyses
    # H.-R. MEYER -PIENING
    # Journal of SANDWICH STRUCTURES AND MATERIALS , Vol. 6—July 2004

    Lx = 5.0*phun("mm") # length  of loaded rectangle
    Ly = 20.0*phun("mm") # length  of loaded rectangle
    Sx = 100.0*phun("mm") # span of the plate
    Sy = 200.0*phun("mm") # span of the plate

    # Here we define the layout and the thicknesses of the layers.
    angles = vec([0.0 0.0 0.0]);
    ts = vec([0.5  11.4  0.1])*phun("mm"); # layer thicknesses
    TH = sum(ts); # total thickness of the plate

    tolerance = 0.0001*TH

    # Select how find the mesh should be
    Refinement = 5
    nL = Refinement * 1;
    nSx = nL + Refinement * 4;
    nSy = 2 * nSx;

    # Each layer is modeled with a single element.
    nts= Refinement * [1, 2, 1];# number of elements per layer
    strength = 1.5
    xs = unique(vcat(reverse(collect(MeshUtilModule.gradedspace(Lx/2, 0.0, nL+1, strength))),
        collect(MeshUtilModule.gradedspace(Lx/2, Sx/2, nSx-nL+1, strength))))
    ys = unique(vcat(reverse(collect(MeshUtilModule.gradedspace(Ly/2, 0.0, nL+1, strength))),
        collect(MeshUtilModule.gradedspace(Ly/2, Sy/2, nSy-nL+1, strength))))

    fens,fes = T10compositeplatex(xs, ys, ts, nts)

    @test count(fens) == 211191
    @test count(fes) == 150000
    rls = selectelem(fens, fes, label = 3)
    @test length(rls) == 37500
end
end
using mMeyer_Piening_2
mMeyer_Piening_2.test()


module mMeyer_Piening_3
using FinEtools
using FinEtools.MeshUtilModule
using Base.Test
function test()

    # Reference results from:
    # [1] Application of the Elasticity Solution
    # to Linear Sandwich Beam, Plate
    # and Shell Analyses
    # H.-R. MEYER -PIENING
    # Journal of SANDWICH STRUCTURES AND MATERIALS , Vol. 6—July 2004

    Lx = 5.0*phun("mm") # length  of loaded rectangle
    Ly = 20.0*phun("mm") # length  of loaded rectangle
    Sx = 100.0*phun("mm") # span of the plate
    Sy = 200.0*phun("mm") # span of the plate

    # Here we define the layout and the thicknesses of the layers.
    angles = vec([0.0 0.0 0.0]);
    ts = vec([0.5  11.4  0.1])*phun("mm"); # layer thicknesses
    TH = sum(ts); # total thickness of the plate

    tolerance = 0.0001*TH

    # Select how find the mesh should be
    Refinement = 5
    nL = Refinement * 1;
    nSx = nL + Refinement * 4;
    nSy = 2 * nSx;

    # Each layer is modeled with a single element.
    nts= Refinement * [1, 2, 1];# number of elements per layer
    strength = 1.5
    xs = unique(vcat(reverse(collect(MeshUtilModule.gradedspace(Lx/2, 0.0, nL+1, strength))),
        collect(MeshUtilModule.gradedspace(Lx/2, Sx/2, nSx-nL+1, strength))))
    ys = unique(vcat(reverse(collect(MeshUtilModule.gradedspace(Ly/2, 0.0, nL+1, strength))),
        collect(MeshUtilModule.gradedspace(Ly/2, Sy/2, nSy-nL+1, strength))))

    fens,fes = T10compositeplatex(xs, ys, ts, nts, :b)

    @test count(fens) == 211191
    @test count(fes) == 150000
    rls = selectelem(fens, fes, label = 3)
    @test length(rls) == 37500
end
end
using mMeyer_Piening_3
mMeyer_Piening_3.test()

module mmtetblocksmm
using FinEtools
using Base.Test
function test()
    A = 5.0*phun("mm") # length  of loaded rectangle
    B = 20.0*phun("mm") # length  of loaded rectangle
    C = 100.0*phun("mm") # span of the plate

    # Select how find the mesh should be
    Refinement = 5
    nA, nB, nC = Refinement * 1, Refinement * 2, Refinement * 1;
    xs = reshape(collect(linspace(0.0, A, nA + 1)), nA + 1, 1)
    ys = reshape(collect(linspace(0.0, A, nB + 1)), nB + 1, 1)
    zs = reshape(collect(linspace(0.0, A, nC + 1)), nC + 1, 1)
    fens,fes = T10blockx(xs, ys, zs, :b)
    # println("$(count(fens))")
    # println("$(count(fes))")
    @test count(fens) == 2541
    @test count(fes) == 1500

    fens,fes = T10blockx(xs, ys, zs, :a)
    # println("$(count(fens))")
    # println("$(count(fes))")
    @test count(fens) == 2541
    @test count(fes) == 1500


    fens,fes = T10blockx(vec(xs), vec(ys), vec(zs), :a)
    # println("$(count(fens))")
    # println("$(count(fes))")
    @test count(fens) == 2541
    @test count(fes) == 1500


    fens,fes = T4blockx(xs, ys, zs, :ca)
    # println("$(count(fens))")
    # println("$(count(fes))")
    @test count(fens) == 396
    @test count(fes) == 1250

    fens,fes = T4blockx(vec(xs), vec(ys), vec(zs), :ca)
    @test count(fens) == 396
    @test count(fes) == 1250

    return true
end
end
using mmtetblocksmm
mmtetblocksmm.test()

module momap2para1
using FinEtools
using Base.Test
function test()
    X = [-1.0 -1.0; 2.0 0.4; 1.0 2.3; -2.0 1.0]
    fens = FENodeSet(X);
    fes = FESetQ4(reshape([1 2 3 4], 1, 4))
    pt = [0.1 0.2]
    pc, success = map2parametric(fes, fens.xyz[fes.conn[1, :], :], vec(pt))
    N = bfun(fes,  pc)
    pt1 = transpose(N) * fens.xyz[fes.conn[1, :], :]
    # println("pt = $(pt)")
    # println("pt1 = $(pt1)")
    @test norm(pt - pt1) < 1.0e-7
end
end
using momap2para1
momap2para1.test()

module momap2para2
using FinEtools
using Base.Test
function test()
    X = [-1.0 -1.0; 2.0 0.4; 1.0 2.3; -2.0 1.0]
    fens = FENodeSet(X);
    fes = FESetQ4(reshape([1 2 3 4], 1, 4))
    pt = [3.1 0.2]
    pc, success = map2parametric(fes, fens.xyz[fes.conn[1, :], :], vec(pt);
        Tolerance = 0.000001, maxiter =7)
    # println("pc = $(pc)")
    N = bfun(fes,  pc)
    pt1 = transpose(N) * fens.xyz[fes.conn[1, :], :]
    # println("pt = $(pt)")
    # println("pt1 = $(pt1)")
    @test norm(pt - pt1) < 1.0e-7
    isin = inparametric(fes, pc)
    # println("isin = $(isin)")
    @test !isin
end
end
using momap2para2
momap2para2.test()

module momap2para3
using FinEtools
using FinEtools.MeshSelectionModule: vselect
using FinEtools.MeshExportModule
using Base.Test



function test()
    A = 50.0*phun("m") # length  of loaded rectangle
    B = 200.0*phun("m") # length  of loaded rectangle
    C = 100.0*phun("m") # span of the plate

    # Select how find the mesh should be
    Refinement = 2
    nA, nB, nC = Refinement * 1, Refinement * 2, Refinement * 4;
    xs = reshape(collect(linspace(0.0, A, nA + 1)), nA + 1, 1)
    ys = reshape(collect(linspace(0.0, B, nB + 1)), nB + 1, 1)
    zs = reshape(collect(linspace(0.0, C, nC + 1)), nC + 1, 1)
    fensc,fesc = T10blockx(xs, ys, zs, :b)
    fc = NodalField(zeros(count(fensc), 1))
    for i = 1:count(fensc)
        x, y, z = fensc.xyz[i, :]
        fc.values[i, :] = sin(2*x/A) * cos(6.5*y/B) * sin(3*z/C-1.0)
    end
    File = "momap2para3-coarse.vtk"
    MeshExportModule.vtkexportmesh(File, fensc, fesc; scalars = [("fc", fc.values)])
    # @async run(`"paraview.exe" $File`)
    try rm(File) catch end

    Refinement = Refinement + 1
    nA, nB, nC = Refinement * 1, Refinement * 2, Refinement * 4;
    xs = reshape(collect(linspace(0.0, A, nA + 1)), nA + 1, 1)
    ys = reshape(collect(linspace(0.0, B, nB + 1)), nB + 1, 1)
    zs = reshape(collect(linspace(0.0, C, nC + 1)), nC + 1, 1)
    fensf,fesf = T10blockx(xs, ys, zs, :b)
    ff = NodalField(zeros(count(fensf), 1))
    tolerance = min(A/nA, B/nB, C/nC)/1000.0

    referenceff = NodalField(zeros(count(fensf), 1))
    for i = 1:count(fensf)
        x, y, z = fensf.xyz[i, :]
        referenceff.values[i, :] = sin(2*x/A) * cos(6.5*y/B) * sin(3*z/C-1.0)
    end
    File = "momap2para3-reference.vtk"
    MeshExportModule.vtkexportmesh(File, fensf, fesf; scalars = [("referenceff", referenceff.values)])
    # @async run(`"paraview.exe" $File`)
    try rm(File) catch end

    ff = transferfield!(ff, fensf, fesf, fc, fensc, fesc, tolerance)
    File = "momap2para3-fine.vtk"
    MeshExportModule.vtkexportmesh(File, fensf, fesf; scalars = [("ff", ff.values)])
    # @async run(`"paraview.exe" $File`)
    try rm(File) catch end

    diffff = NodalField(referenceff.values - ff.values)
    femm  = FEMMBase(GeoD(fesf, SimplexRule(3, 4)))
    geom = NodalField(fensf.xyz)
    error = integratefieldfunction(femm, geom, diffff, (x, v) -> norm(v), 0.0)
    ref = integratefieldfunction(femm, geom, referenceff, (x, v) -> norm(v), 0.0)
    # println("error/ref = $(error/ref)")
    @test abs(error/ref - 0.02541369940759616) < 1.0e-4
end
end
using momap2para3
momap2para3.test()

module momap2para4
using FinEtools
using FinEtools.MeshSelectionModule: vselect
using FinEtools.MeshExportModule
using Base.Test



function test()
    A = 50.0*phun("m") # length  of loaded rectangle
    B = 200.0*phun("m") # length  of loaded rectangle
    C = 100.0*phun("m") # span of the plate

    # Select how find the mesh should be
    Refinement = 2
    nA, nB, nC = Refinement * 1, Refinement * 6, Refinement * 4;
    xs = reshape(collect(linspace(0.0, A, nA + 1)), nA + 1, 1)
    ys = reshape(collect(linspace(0.0, B, nB + 1)), nB + 1, 1)
    zs = reshape(collect(linspace(0.0, C, nC + 1)), nC + 1, 1)
    fensc,fesc = T10blockx(xs, ys, zs, :b)

    centroidpc = centroidparametric(fesc)
    N = bfun(fesc, centroidpc)
    NT = transpose(N)

    fc = ElementalField(zeros(count(fesc), 1))
    for i = 1:count(fesc)
        c = view(fesc.conn, i, :)
        centroid = NT * fensc.xyz[c, :]
        x, y, z = centroid
        fc.values[i, :] = sin(2*x/A) * cos(6.5*y/B) * sin(3*z/C-1.0)
    end
    File = "momap2para3-coarse.vtk"
    MeshExportModule.vtkexportmesh(File, fensc, fesc; scalars = [("fc", fc.values)])
    # @async run(`"paraview.exe" $File`)
    try rm(File) catch end

    Refinement = Refinement + 1
    nA, nB, nC = Refinement * 1, Refinement * 6, Refinement * 4;
    xs = reshape(collect(linspace(0.0, A, nA + 1)), nA + 1, 1)
    ys = reshape(collect(linspace(0.0, B, nB + 1)), nB + 1, 1)
    zs = reshape(collect(linspace(0.0, C, nC + 1)), nC + 1, 1)
    fensf,fesf = T10blockx(xs, ys, zs, :b)
    tolerance = min(A/nA, B/nB, C/nC)/1000.0

    ff = ElementalField(zeros(count(fesf), 1))
    referenceff = ElementalField(zeros(count(fesf), 1))
    for i = 1:count(fesf)
        c = view(fesf.conn, i, :)
        centroid = NT * fensf.xyz[c, :]
        x, y, z = centroid
        referenceff.values[i, :] = sin(2*x/A) * cos(6.5*y/B) * sin(3*z/C-1.0)
    end
    File = "momap2para3-reference.vtk"
    MeshExportModule.vtkexportmesh(File, fensf, fesf; scalars = [("referenceff", referenceff.values)])
    # @async run(`"paraview.exe" $File`)
    try rm(File) catch end

    ff = transferfield!(ff, fensf, fesf, fc, fensc, fesc, tolerance)
    File = "momap2para3-fine.vtk"
    MeshExportModule.vtkexportmesh(File, fensf, fesf; scalars = [("ff", ff.values)])
    # @async run(`"paraview.exe" $File`)
    try rm(File) catch end

    diffff = ElementalField(referenceff.values - ff.values)
    femm  = FEMMBase(GeoD(fesf, SimplexRule(3, 4)))
    geom = NodalField(fensf.xyz)
    error = integratefieldfunction(femm, geom, diffff, (x, v) -> norm(v), 0.0)
    ref = integratefieldfunction(femm, geom, referenceff, (x, v) -> norm(v), 0.0)
    # println("error/ref = $(error/ref)")
    @test abs(error/ref - 0.19808425992688541) < 1.0e-4
end
end
using momap2para4
momap2para4.test()

module momap2para5
using FinEtools
using FinEtools.MeshSelectionModule: vselect
using FinEtools.MeshExportModule
using Base.Test



function test()
    A = 50.0*phun("m") # length  of loaded rectangle
    B = 200.0*phun("m") # length  of loaded rectangle
    C = 100.0*phun("m") # span of the plate

Meshing = T4blockx
    # Select how find the mesh should be
    Refinement = 2
    nA, nB, nC = Refinement * 1, Refinement * 6, Refinement * 4;
    xs = reshape(collect(linspace(0.0, A, nA + 1)), nA + 1, 1)
    ys = reshape(collect(linspace(0.0, B, nB + 1)), nB + 1, 1)
    zs = reshape(collect(linspace(0.0, C, nC + 1)), nC + 1, 1)
    fensc,fesc = Meshing(xs, ys, zs, :b)

    centroidpc = centroidparametric(fesc)
    N = bfun(fesc, centroidpc)
    NT = transpose(N)

    fc = ElementalField(zeros(count(fesc), 1))
    for i = 1:count(fesc)
        c = view(fesc.conn, i, :)
        centroid = NT * fensc.xyz[c, :]
        x, y, z = centroid
        fc.values[i, :] = sin(2*x/A) * cos(6.5*y/B) * sin(3*z/C-1.0)
    end
    File = "momap2para3-coarse.vtk"
    MeshExportModule.vtkexportmesh(File, fensc, fesc; scalars = [("fc", fc.values)])
    # @async run(`"paraview.exe" $File`)
    try rm(File) catch end

    Refinement = Refinement + 1
    nA, nB, nC = Refinement * 1, Refinement * 6, Refinement * 4;
    xs = reshape(collect(linspace(0.0, A, nA + 1)), nA + 1, 1)
    ys = reshape(collect(linspace(0.0, B, nB + 1)), nB + 1, 1)
    zs = reshape(collect(linspace(0.0, C, nC + 1)), nC + 1, 1)
    fensf,fesf = Meshing(xs, ys, zs, :b)
    tolerance = min(A/nA, B/nB, C/nC)/1000.0

    ff = ElementalField(zeros(count(fesf), 1))
    referenceff = ElementalField(zeros(count(fesf), 1))
    for i = 1:count(fesf)
        c = view(fesf.conn, i, :)
        centroid = NT * fensf.xyz[c, :]
        x, y, z = centroid
        referenceff.values[i, :] = sin(2*x/A) * cos(6.5*y/B) * sin(3*z/C-1.0)
    end
    File = "momap2para3-reference.vtk"
    MeshExportModule.vtkexportmesh(File, fensf, fesf; scalars = [("referenceff", referenceff.values)])
    # @async run(`"paraview.exe" $File`)
    try rm(File) catch end

    ff = transferfield!(ff, fensf, fesf, fc, fensc, fesc, tolerance)
    File = "momap2para3-fine.vtk"
    MeshExportModule.vtkexportmesh(File, fensf, fesf; scalars = [("ff", ff.values)])
    # @async run(`"paraview.exe" $File`)
    try rm(File) catch end

    diffff = ElementalField(referenceff.values - ff.values)
    femm  = FEMMBase(GeoD(fesf, SimplexRule(3, 4)))
    geom = NodalField(fensf.xyz)
    error = integratefieldfunction(femm, geom, diffff, (x, v) -> norm(v), 0.0)
    ref = integratefieldfunction(femm, geom, referenceff, (x, v) -> norm(v), 0.0)
    # println("error/ref = $(error/ref)")
    @test abs(error/ref - 0.19808425992688541) < 1.0e-4
end
end
using momap2para5
momap2para5.test()


module momap2para6
using FinEtools
using FinEtools.MeshSelectionModule: vselect
using FinEtools.MeshExportModule
using Base.Test



function test()
    A = 50.0*phun("m") # length  of loaded rectangle
    B = 200.0*phun("m") # length  of loaded rectangle
    C = 100.0*phun("m") # span of the plate

Meshing = H20blockx
    # Select how find the mesh should be
    Refinement = 2
    nA, nB, nC = Refinement * 1, Refinement * 6, Refinement * 4;
    xs = collect(linspace(0.0, A, nA + 1))
    ys = collect(linspace(0.0, B, nB + 1))
    zs = collect(linspace(0.0, C, nC + 1))
    fensc,fesc = Meshing(xs, ys, zs)

    centroidpc = centroidparametric(fesc)
    N = bfun(fesc, centroidpc)
    NT = transpose(N)

    fc = ElementalField(zeros(count(fesc), 1))
    for i = 1:count(fesc)
        c = view(fesc.conn, i, :)
        centroid = NT * fensc.xyz[c, :]
        x, y, z = centroid
        fc.values[i, :] = sin(2*x/A) * cos(6.5*y/B) * sin(3*z/C-1.0)
    end
    File = "momap2para3-coarse.vtk"
    MeshExportModule.vtkexportmesh(File, fensc, fesc; scalars = [("fc", fc.values)])
    # @async run(`"paraview.exe" $File`)
    try rm(File) catch end

    Refinement = Refinement + 1
    nA, nB, nC = Refinement * 1, Refinement * 6, Refinement * 4;
    xs = collect(linspace(0.0, A, nA + 1))
    ys = collect(linspace(0.0, B, nB + 1))
    zs = collect(linspace(0.0, C, nC + 1))
    fensf,fesf = Meshing(xs, ys, zs)
    tolerance = min(A/nA, B/nB, C/nC)/1000.0

    ff = ElementalField(zeros(count(fesf), 1))
    referenceff = ElementalField(zeros(count(fesf), 1))
    for i = 1:count(fesf)
        c = view(fesf.conn, i, :)
        centroid = NT * fensf.xyz[c, :]
        x, y, z = centroid
        referenceff.values[i, :] = sin(2*x/A) * cos(6.5*y/B) * sin(3*z/C-1.0)
    end
    File = "momap2para3-reference.vtk"
    MeshExportModule.vtkexportmesh(File, fensf, fesf; scalars = [("referenceff", referenceff.values)])
    # @async run(`"paraview.exe" $File`)
    try rm(File) catch end

    ff = transferfield!(ff, fensf, fesf, fc, fensc, fesc, tolerance)
    File = "momap2para3-fine.vtk"
    MeshExportModule.vtkexportmesh(File, fensf, fesf; scalars = [("ff", ff.values)])
    # @async run(`"paraview.exe" $File`)
    try rm(File) catch end

    diffff = ElementalField(referenceff.values - ff.values)
    femm  = FEMMBase(GeoD(fesf, GaussRule(3, 3)))
    geom = NodalField(fensf.xyz)
    error = integratefieldfunction(femm, geom, diffff, (x, v) -> norm(v), 0.0)
    ref = integratefieldfunction(femm, geom, referenceff, (x, v) -> norm(v), 0.0)
    # println("error/ref = $(error/ref)")
    @test abs(error/ref - 0.3602708839379691) < 1.0e-4
end
end
using momap2para6
momap2para6.test()

module momap2para7
using FinEtools
using FinEtools.MeshSelectionModule: vselect
using FinEtools.MeshExportModule
using Base.Test



function test()
    A = 50.0*phun("m") # length  of loaded rectangle
    B = 200.0*phun("m") # length  of loaded rectangle
    C = 100.0*phun("m") # span of the plate

Meshing = H8blockx
    # Select how find the mesh should be
    Refinement = 2
    nA, nB, nC = Refinement * 1, Refinement * 6, Refinement * 4;
    xs = collect(linspace(0.0, A, nA + 1))
    ys = collect(linspace(0.0, B, nB + 1))
    zs = collect(linspace(0.0, C, nC + 1))
    fensc,fesc = Meshing(xs, ys, zs)

    centroidpc = centroidparametric(fesc)
    N = bfun(fesc, centroidpc)
    NT = transpose(N)

    fc = ElementalField(zeros(count(fesc), 1))
    for i = 1:count(fesc)
        c = view(fesc.conn, i, :)
        centroid = NT * fensc.xyz[c, :]
        x, y, z = centroid
        fc.values[i, :] = sin(2*x/A) * cos(6.5*y/B) * sin(3*z/C-1.0)
    end
    File = "momap2para3-coarse.vtk"
    MeshExportModule.vtkexportmesh(File, fensc, fesc; scalars = [("fc", fc.values)])
    # @async run(`"paraview.exe" $File`)
    try rm(File) catch end

    Refinement = Refinement + 1
    nA, nB, nC = Refinement * 1, Refinement * 6, Refinement * 4;
    xs = collect(linspace(0.0, A, nA + 1))
    ys = collect(linspace(0.0, B, nB + 1))
    zs = collect(linspace(0.0, C, nC + 1))
    fensf,fesf = Meshing(xs, ys, zs)
    tolerance = min(A/nA, B/nB, C/nC)/1000.0

    ff = ElementalField(zeros(count(fesf), 1))
    referenceff = ElementalField(zeros(count(fesf), 1))
    for i = 1:count(fesf)
        c = view(fesf.conn, i, :)
        centroid = NT * fensf.xyz[c, :]
        x, y, z = centroid
        referenceff.values[i, :] = sin(2*x/A) * cos(6.5*y/B) * sin(3*z/C-1.0)
    end
    File = "momap2para3-reference.vtk"
    MeshExportModule.vtkexportmesh(File, fensf, fesf; scalars = [("referenceff", referenceff.values)])
    # @async run(`"paraview.exe" $File`)
    try rm(File) catch end

    ff = transferfield!(ff, fensf, fesf, fc, fensc, fesc, tolerance)
    File = "momap2para3-fine.vtk"
    MeshExportModule.vtkexportmesh(File, fensf, fesf; scalars = [("ff", ff.values)])
    # @async run(`"paraview.exe" $File`)
    try rm(File) catch end

    diffff = ElementalField(referenceff.values - ff.values)
    femm  = FEMMBase(GeoD(fesf, GaussRule(3, 3)))
    geom = NodalField(fensf.xyz)
    error = integratefieldfunction(femm, geom, diffff, (x, v) -> norm(v), 0.0)
    ref = integratefieldfunction(femm, geom, referenceff, (x, v) -> norm(v), 0.0)
    # println("error/ref = $(error/ref)")
    @test abs(error/ref - 0.36027088393796847) < 1.0e-4
end
end
using momap2para7
momap2para7.test()

module momap2para8
using FinEtools
using FinEtools.MeshSelectionModule: vselect
using FinEtools.MeshExportModule
using Base.Test



function test()
    A = 50.0*phun("m") # length  of loaded rectangle
    B = 200.0*phun("m") # length  of loaded rectangle
    C = 100.0*phun("m") # span of the plate

Meshing = H27blockx
    # Select how find the mesh should be
    Refinement = 2
    nA, nB, nC = Refinement * 1, Refinement * 6, Refinement * 4;
    xs = collect(linspace(0.0, A, nA + 1))
    ys = collect(linspace(0.0, B, nB + 1))
    zs = collect(linspace(0.0, C, nC + 1))
    fensc,fesc = Meshing(xs, ys, zs)

    centroidpc = centroidparametric(fesc)
    N = bfun(fesc, centroidpc)
    NT = transpose(N)

    fc = ElementalField(zeros(count(fesc), 1))
    for i = 1:count(fesc)
        c = view(fesc.conn, i, :)
        centroid = NT * fensc.xyz[c, :]
        x, y, z = centroid
        fc.values[i, :] = sin(2*x/A) * cos(6.5*y/B) * sin(3*z/C-1.0)
    end
    # File = "momap2para3-coarse.vtk"
    # MeshExportModule.vtkexportmesh(File, fensc, fesc; scalars = [("fc", fc.values)])
    # # @async run(`"paraview.exe" $File`)
    # try rm(File) catch end

    Refinement = Refinement + 1
    nA, nB, nC = Refinement * 1, Refinement * 6, Refinement * 4;
    xs = collect(linspace(0.0, A, nA + 1))
    ys = collect(linspace(0.0, B, nB + 1))
    zs = collect(linspace(0.0, C, nC + 1))
    fensf,fesf = Meshing(xs, ys, zs)
    tolerance = min(A/nA, B/nB, C/nC)/1000.0

    ff = ElementalField(zeros(count(fesf), 1))
    referenceff = ElementalField(zeros(count(fesf), 1))
    for i = 1:count(fesf)
        c = view(fesf.conn, i, :)
        centroid = NT * fensf.xyz[c, :]
        x, y, z = centroid
        referenceff.values[i, :] = sin(2*x/A) * cos(6.5*y/B) * sin(3*z/C-1.0)
    end
    # File = "momap2para3-reference.vtk"
    # MeshExportModule.vtkexportmesh(File, fensf, fesf; scalars = [("referenceff", referenceff.values)])
    # # @async run(`"paraview.exe" $File`)
    # try rm(File) catch end

    ff = transferfield!(ff, fensf, fesf, fc, fensc, fesc, tolerance)
    # File = "momap2para3-fine.vtk"
    # MeshExportModule.vtkexportmesh(File, fensf, fesf; scalars = [("ff", ff.values)])
    # # @async run(`"paraview.exe" $File`)
    # try rm(File) catch end

    diffff = ElementalField(referenceff.values - ff.values)
    femm  = FEMMBase(GeoD(fesf, GaussRule(3, 3)))
    geom = NodalField(fensf.xyz)
    error = integratefieldfunction(femm, geom, diffff, (x, v) -> norm(v), 0.0)
    ref = integratefieldfunction(femm, geom, referenceff, (x, v) -> norm(v), 0.0)
    # println("error/ref = $(error/ref)")
    @test abs(error/ref - 0.3602708839379695) < 1.0e-4
end
end
using momap2para8
momap2para8.test()


module momap2para9
using FinEtools
using FinEtools.MeshSelectionModule: vselect
using FinEtools.MeshExportModule
using Base.Test



function test()
    A = 50.0*phun("m") # length  of loaded rectangle
    B = 200.0*phun("m") # length  of loaded rectangle
    C = 100.0*phun("m") # span of the plate

Meshing = Q4blockx
    # Select how find the mesh should be
    Refinement = 2
    nA, nB, nC = Refinement * 1, Refinement * 6, Refinement * 4;
    xs = collect(linspace(0.0, A, nA + 1))
    ys = collect(linspace(0.0, B, nB + 1))
    fensc,fesc = Meshing(xs, ys)

    centroidpc = centroidparametric(fesc)
    N = bfun(fesc, centroidpc)
    NT = transpose(N)

    fc = ElementalField(zeros(count(fesc), 1))
    for i = 1:count(fesc)
        c = view(fesc.conn, i, :)
        centroid = NT * fensc.xyz[c, :]
        x, y = centroid
        fc.values[i, :] = sin(2*x/A) * cos(6.5*y/B)
    end
    # File = "momap2para3-coarse.vtk"
    # MeshExportModule.vtkexportmesh(File, fensc, fesc; scalars = [("fc", fc.values)])
    # # @async run(`"paraview.exe" $File`)
    # try rm(File) catch end

    Refinement = Refinement + 1
    nA, nB, nC = Refinement * 1, Refinement * 6, Refinement * 4;
    xs = collect(linspace(0.0, A, nA + 1))
    ys = collect(linspace(0.0, B, nB + 1))
    zs = collect(linspace(0.0, C, nC + 1))
    fensf,fesf = Meshing(xs, ys)
    tolerance = min(A/nA, B/nB, C/nC)/1000.0

    ff = ElementalField(zeros(count(fesf), 1))
    referenceff = ElementalField(zeros(count(fesf), 1))
    for i = 1:count(fesf)
        c = view(fesf.conn, i, :)
        centroid = NT * fensf.xyz[c, :]
        x, y = centroid
        referenceff.values[i, :] = sin(2*x/A) * cos(6.5*y/B)
    end
    # File = "momap2para3-reference.vtk"
    # MeshExportModule.vtkexportmesh(File, fensf, fesf; scalars = [("referenceff", referenceff.values)])
    # # @async run(`"paraview.exe" $File`)
    # try rm(File) catch end

    ff = transferfield!(ff, fensf, fesf, fc, fensc, fesc, tolerance)
    # File = "momap2para3-fine.vtk"
    # MeshExportModule.vtkexportmesh(File, fensf, fesf; scalars = [("ff", ff.values)])
    # # @async run(`"paraview.exe" $File`)
    # try rm(File) catch end

    diffff = ElementalField(referenceff.values - ff.values)
    femm  = FEMMBase(GeoD(fesf, GaussRule(2, 3)))
    geom = NodalField(fensf.xyz)
    error = integratefieldfunction(femm, geom, diffff, (x, v) -> norm(v), 0.0)
    ref = integratefieldfunction(femm, geom, referenceff, (x, v) -> norm(v), 0.0)
    # println("error/ref = $(error/ref)")
    @test abs(error/ref - 0.32348898713315494) < 1.0e-4
end
end
using momap2para9
momap2para9.test()


module momap2para10
using FinEtools
using FinEtools.MeshSelectionModule: vselect
using FinEtools.MeshExportModule
using Base.Test



function test()
    A = 50.0*phun("m") # length  of loaded rectangle
    B = 200.0*phun("m") # length  of loaded rectangle
    C = 100.0*phun("m") # span of the plate

Meshing = Q8blockx
    # Select how find the mesh should be
    Refinement = 2
    nA, nB, nC = Refinement * 1, Refinement * 6, Refinement * 4;
    xs = collect(linspace(0.0, A, nA + 1))
    ys = collect(linspace(0.0, B, nB + 1))
    fensc,fesc = Meshing(xs, ys)

    centroidpc = centroidparametric(fesc)
    N = bfun(fesc, centroidpc)
    NT = transpose(N)

    fc = ElementalField(zeros(count(fesc), 1))
    for i = 1:count(fesc)
        c = view(fesc.conn, i, :)
        centroid = NT * fensc.xyz[c, :]
        x, y = centroid
        fc.values[i, :] = sin(2*x/A) * cos(6.5*y/B)
    end
    # File = "momap2para3-coarse.vtk"
    # MeshExportModule.vtkexportmesh(File, fensc, fesc; scalars = [("fc", fc.values)])
    # # @async run(`"paraview.exe" $File`)
    # try rm(File) catch end

    Refinement = Refinement + 1
    nA, nB, nC = Refinement * 1, Refinement * 6, Refinement * 4;
    xs = collect(linspace(0.0, A, nA + 1))
    ys = collect(linspace(0.0, B, nB + 1))
    zs = collect(linspace(0.0, C, nC + 1))
    fensf,fesf = Meshing(xs, ys)
    tolerance = min(A/nA, B/nB, C/nC)/1000.0

    ff = ElementalField(zeros(count(fesf), 1))
    referenceff = ElementalField(zeros(count(fesf), 1))
    for i = 1:count(fesf)
        c = view(fesf.conn, i, :)
        centroid = NT * fensf.xyz[c, :]
        x, y = centroid
        referenceff.values[i, :] = sin(2*x/A) * cos(6.5*y/B)
    end
    # File = "momap2para3-reference.vtk"
    # MeshExportModule.vtkexportmesh(File, fensf, fesf; scalars = [("referenceff", referenceff.values)])
    # # @async run(`"paraview.exe" $File`)
    # try rm(File) catch end

    ff = transferfield!(ff, fensf, fesf, fc, fensc, fesc, tolerance)
    # File = "momap2para3-fine.vtk"
    # MeshExportModule.vtkexportmesh(File, fensf, fesf; scalars = [("ff", ff.values)])
    # # @async run(`"paraview.exe" $File`)
    # try rm(File) catch end

    diffff = ElementalField(referenceff.values - ff.values)
    femm  = FEMMBase(GeoD(fesf, GaussRule(2, 3)))
    geom = NodalField(fensf.xyz)
    error = integratefieldfunction(femm, geom, diffff, (x, v) -> norm(v), 0.0)
    ref = integratefieldfunction(femm, geom, referenceff, (x, v) -> norm(v), 0.0)
    # println("error/ref = $(error/ref)")
    @test abs(error/ref - 0.32348898713315494) < 1.0e-4
end
end
using momap2para10
momap2para10.test()


module momap2para11
using FinEtools
using FinEtools.MeshSelectionModule: vselect
using FinEtools.MeshExportModule
using Base.Test



function test()
    A = 50.0*phun("m") # length  of loaded rectangle
    B = 200.0*phun("m") # length  of loaded rectangle
    C = 100.0*phun("m") # span of the plate

Meshing = T3blockx
    # Select how find the mesh should be
    Refinement = 2
    nA, nB, nC = Refinement * 1, Refinement * 6, Refinement * 4;
    xs = collect(linspace(0.0, A, nA + 1))
    ys = collect(linspace(0.0, B, nB + 1))
    fensc,fesc = Meshing(xs, ys)

    centroidpc = centroidparametric(fesc)
    N = bfun(fesc, centroidpc)
    NT = transpose(N)

    fc = ElementalField(zeros(count(fesc), 1))
    for i = 1:count(fesc)
        c = view(fesc.conn, i, :)
        centroid = NT * fensc.xyz[c, :]
        x, y = centroid
        fc.values[i, :] = sin(2*x/A) * cos(6.5*y/B)
    end
    File = "momap2para3-coarse.vtk"
    MeshExportModule.vtkexportmesh(File, fensc, fesc; scalars = [("fc", fc.values)])
    # @async run(`"paraview.exe" $File`)
    try rm(File) catch end

    Refinement = Refinement + 1
    nA, nB, nC = Refinement * 1, Refinement * 6, Refinement * 4;
    xs = collect(linspace(0.0, A, nA + 1))
    ys = collect(linspace(0.0, B, nB + 1))
    zs = collect(linspace(0.0, C, nC + 1))
    fensf,fesf = Meshing(xs, ys)
    tolerance = min(A/nA, B/nB, C/nC)/1000.0

    ff = ElementalField(zeros(count(fesf), 1))
    referenceff = ElementalField(zeros(count(fesf), 1))
    for i = 1:count(fesf)
        c = view(fesf.conn, i, :)
        centroid = NT * fensf.xyz[c, :]
        x, y = centroid
        referenceff.values[i, :] = sin(2*x/A) * cos(6.5*y/B)
    end
    File = "momap2para3-reference.vtk"
    MeshExportModule.vtkexportmesh(File, fensf, fesf; scalars = [("referenceff", referenceff.values)])
    # @async run(`"paraview.exe" $File`)
    try rm(File) catch end

    ff = transferfield!(ff, fensf, fesf, fc, fensc, fesc, tolerance)
    File = "momap2para3-fine.vtk"
    MeshExportModule.vtkexportmesh(File, fensf, fesf; scalars = [("ff", ff.values)])
    # @async run(`"paraview.exe" $File`)
    try rm(File) catch end

    diffff = ElementalField(referenceff.values - ff.values)
    femm  = FEMMBase(GeoD(fesf, SimplexRule(2, 3)))
    geom = NodalField(fensf.xyz)
    error = integratefieldfunction(femm, geom, diffff, (x, v) -> norm(v), 0.0)
    ref = integratefieldfunction(femm, geom, referenceff, (x, v) -> norm(v), 0.0)
    println("error/ref = $(error/ref)")
    @test abs(error/ref -  0.19654730310399787) < 1.0e-4
end
end
using momap2para11
momap2para11.test()


module momap2para12
using FinEtools
using FinEtools.MeshSelectionModule: vselect
using FinEtools.MeshExportModule
using Base.Test



function test()
    A = 50.0*phun("m") # length  of loaded rectangle
    B = 200.0*phun("m") # length  of loaded rectangle
    C = 100.0*phun("m") # span of the plate

Meshing = T6blockx
    # Select how find the mesh should be
    Refinement = 2
    nA, nB, nC = Refinement * 1, Refinement * 6, Refinement * 4;
    xs = collect(linspace(0.0, A, nA + 1))
    ys = collect(linspace(0.0, B, nB + 1))
    fensc,fesc = Meshing(xs, ys)

    centroidpc = centroidparametric(fesc)
    N = bfun(fesc, centroidpc)
    NT = transpose(N)

    fc = ElementalField(zeros(count(fesc), 1))
    for i = 1:count(fesc)
        c = view(fesc.conn, i, :)
        centroid = NT * fensc.xyz[c, :]
        x, y = centroid
        fc.values[i, :] = sin(2*x/A) * cos(6.5*y/B)
    end
    File = "momap2para12-coarse.vtk"
    MeshExportModule.vtkexportmesh(File, fensc, fesc; scalars = [("fc", fc.values)])
    # @async run(`"paraview.exe" $File`)
    # try rm(File) catch end

    Refinement = Refinement + 1
    nA, nB, nC = Refinement * 1, Refinement * 6, Refinement * 4;
    xs = collect(linspace(0.0, A, nA + 1))
    ys = collect(linspace(0.0, B, nB + 1))
    zs = collect(linspace(0.0, C, nC + 1))
    fensf,fesf = Meshing(xs, ys)
    tolerance = min(A/nA, B/nB, C/nC)/1000.0

    ff = ElementalField(zeros(count(fesf), 1))
    referenceff = ElementalField(zeros(count(fesf), 1))
    for i = 1:count(fesf)
        c = view(fesf.conn, i, :)
        centroid = NT * fensf.xyz[c, :]
        x, y = centroid
        referenceff.values[i, :] = sin(2*x/A) * cos(6.5*y/B)
    end
    File = "momap2para12-reference.vtk"
    MeshExportModule.vtkexportmesh(File, fensf, fesf; scalars = [("referenceff", referenceff.values)])
    # @async run(`"paraview.exe" $File`)
    # try rm(File) catch end

    ff = transferfield!(ff, fensf, fesf, fc, fensc, fesc, tolerance)
    File = "momap2para12-fine.vtk"
    MeshExportModule.vtkexportmesh(File, fensf, fesf; scalars = [("ff", ff.values)])
    # @async run(`"paraview.exe" $File`)
    # try rm(File) catch end

    diffff = ElementalField(referenceff.values - ff.values)
    femm  = FEMMBase(GeoD(fesf, SimplexRule(2, 3)))
    geom = NodalField(fensf.xyz)
    error = integratefieldfunction(femm, geom, diffff, (x, v) -> norm(v), 0.0)
    ref = integratefieldfunction(femm, geom, referenceff, (x, v) -> norm(v), 0.0)
    # println("error/ref = $(error/ref)")
    @test abs(error/ref -  0.19654730310399787) < 1.0e-4
end
end
using momap2para12
momap2para12.test()


module momap2para13
using FinEtools
using FinEtools.MeshSelectionModule: vselect
using FinEtools.MeshExportModule
using Base.Test



function test()
    A = 50.0*phun("m") # length  of loaded rectangle
    B = 200.0*phun("m") # length  of loaded rectangle
    C = 100.0*phun("m") # span of the plate

Meshing = L2blockx
    # Select how find the mesh should be
    Refinement = 2
    nA, nB, nC = Refinement * 1, Refinement * 6, Refinement * 4;
    xs = collect(linspace(0.0, A, nA + 1))
    fensc,fesc = Meshing(xs)

    centroidpc = centroidparametric(fesc)
    N = bfun(fesc, centroidpc)
    NT = transpose(N)

    fc = ElementalField(zeros(count(fesc), 1))
    for i = 1:count(fesc)
        c = view(fesc.conn, i, :)
        centroid = NT * fensc.xyz[c, :]
        x = centroid
        fc.values[i, :] = sin.(2*x/A)
    end
    File = "momap2para12-coarse.vtk"
    MeshExportModule.vtkexportmesh(File, fensc, fesc; scalars = [("fc", fc.values)])
    # @async run(`"paraview.exe" $File`)
    # try rm(File) catch end

    Refinement = Refinement + 1
    nA, nB, nC = Refinement * 1, Refinement * 6, Refinement * 4;
    xs = collect(linspace(0.0, A, nA + 1))
    fensf,fesf = Meshing(xs)
    tolerance = min(A/nA, B/nB, C/nC)/1000.0

    ff = ElementalField(zeros(count(fesf), 1))
    referenceff = ElementalField(zeros(count(fesf), 1))
    for i = 1:count(fesf)
        c = view(fesf.conn, i, :)
        centroid = NT * fensf.xyz[c, :]
        x = centroid
        referenceff.values[i, :] = sin.(2*x/A)
    end
    File = "momap2para12-reference.vtk"
    MeshExportModule.vtkexportmesh(File, fensf, fesf; scalars = [("referenceff", referenceff.values)])
    # @async run(`"paraview.exe" $File`)
    # try rm(File) catch end

    ff = transferfield!(ff, fensf, fesf, fc, fensc, fesc, tolerance)
    File = "momap2para12-fine.vtk"
    MeshExportModule.vtkexportmesh(File, fensf, fesf; scalars = [("ff", ff.values)])
    # @async run(`"paraview.exe" $File`)
    # try rm(File) catch end

    diffff = ElementalField(referenceff.values - ff.values)
    femm  = FEMMBase(GeoD(fesf, GaussRule(1, 3)))
    geom = NodalField(fensf.xyz)
    error = integratefieldfunction(femm, geom, diffff, (x, v) -> norm(v), 0.0)
    ref = integratefieldfunction(femm, geom, referenceff, (x, v) -> norm(v), 0.0)
    # println("error/ref = $(error/ref)")
    @test abs(error/ref -  0.23860709149331033) < 1.0e-4
end
end
using momap2para13
momap2para13.test()

module momap2para14
using FinEtools
using FinEtools.MeshSelectionModule: vselect
using FinEtools.MeshExportModule
using Base.Test



function test()
    A = 50.0*phun("m") # length  of loaded rectangle
    B = 200.0*phun("m") # length  of loaded rectangle
    C = 100.0*phun("m") # span of the plate

Meshing = L3blockx
    # Select how find the mesh should be
    Refinement = 2
    nA, nB, nC = Refinement * 1, Refinement * 6, Refinement * 4;
    xs = collect(linspace(0.0, A, nA + 1))
    fensc,fesc = Meshing(xs)
    println("fensc = $(fensc)")
println("fesc = $(fesc)")

    centroidpc = centroidparametric(fesc)
    N = bfun(fesc, centroidpc)
    NT = transpose(N)

    fc = ElementalField(zeros(count(fesc), 1))
    for i = 1:count(fesc)
        c = view(fesc.conn, i, :)
        centroid = NT * fensc.xyz[c, :]
        x = centroid
        fc.values[i, :] = sin.(2*x/A)
    end
    File = "momap2para12-coarse.vtk"
    MeshExportModule.vtkexportmesh(File, fensc, fesc; scalars = [("fc", fc.values)])
    # @async run(`"paraview.exe" $File`)
    # try rm(File) catch end

    Refinement = Refinement + 1
    nA, nB, nC = Refinement * 1, Refinement * 6, Refinement * 4;
    xs = collect(linspace(0.0, A, nA + 1))
    fensf,fesf = Meshing(xs)
    tolerance = min(A/nA, B/nB, C/nC)/1000.0

    ff = ElementalField(zeros(count(fesf), 1))
    referenceff = ElementalField(zeros(count(fesf), 1))
    for i = 1:count(fesf)
        c = view(fesf.conn, i, :)
        centroid = NT * fensf.xyz[c, :]
        x = centroid
        referenceff.values[i, :] = sin.(2*x/A)
    end
    File = "momap2para12-reference.vtk"
    MeshExportModule.vtkexportmesh(File, fensf, fesf; scalars = [("referenceff", referenceff.values)])
    # @async run(`"paraview.exe" $File`)
    # try rm(File) catch end

    ff = transferfield!(ff, fensf, fesf, fc, fensc, fesc, tolerance)
    File = "momap2para12-fine.vtk"
    MeshExportModule.vtkexportmesh(File, fensf, fesf; scalars = [("ff", ff.values)])
    # @async run(`"paraview.exe" $File`)
    # try rm(File) catch end

    diffff = ElementalField(referenceff.values - ff.values)
    femm  = FEMMBase(GeoD(fesf, GaussRule(1, 3)))
    geom = NodalField(fensf.xyz)
    error = integratefieldfunction(femm, geom, diffff, (x, v) -> norm(v), 0.0)
    ref = integratefieldfunction(femm, geom, referenceff, (x, v) -> norm(v), 0.0)
    # println("error/ref = $(error/ref)")
    @test abs(error/ref -  0.23860709149331033) < 1.0e-4
end
end
using momap2para14
momap2para14.test()
