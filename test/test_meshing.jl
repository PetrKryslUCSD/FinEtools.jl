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
  fens, fes = MeshImportModule.import_NASTRAN(dirname(@__FILE__) * "/" * "Slot-coarser.nas";
    allocationchunk = 13)
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

module mmmmmAbaqusmmiimportmmm
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


  fens, fesarray = MeshImportModule.import_ABAQUS("./LE11NAFEMS_H8.inp";
    allocationchunk = 11)

  File = "LE11NAFEMS_H8.vtk"
  MeshExportModule.vtkexportmesh(File, fens, fesarray[1])
  # @async run(`"paraview.exe" $File`)
  try rm(File) catch end

end
end
using mmmmmAbaqusmmiimportmmm
mmmmmAbaqusmmiimportmmm.test()

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


  fens, fesarray = MeshImportModule.import_ABAQUS("./LE11NAFEMS_H8.inp";
    allocationchunk = 12)

  File = "LE11NAFEMS_H8.vtk"
  MeshExportModule.vtkexportmesh(File, fens, fesarray[1])
  # @async run(`"paraview.exe" $File`)
  try rm(File) catch end

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
