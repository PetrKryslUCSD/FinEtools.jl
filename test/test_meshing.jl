module mexpmeshx1ba3
using FinEtools
using FinEtools.MeshExportModule.MESH: write_MESH
using FinEtools.MeshImportModule
using LinearAlgebra: norm
using Test
function test()
    rho = 1.21 * 1e-9# mass density
    c = 345.0 * 1000# millimeters per second
    bulk = c^2 * rho
    Lx = 1900.0# length of the box, millimeters
    Ly = 800.0 # length of the box, millimeters

    fens, fes = Q4block(Lx, Ly, 4, 2) # Mesh
    setlabel!(fes, 3)

    write_MESH("q4-4-2", fens, fes)

    output = MeshImportModule.import_MESH("q4-4-2.mesh")
    @test count(output["fens"]) == count(fens)
    @test count(output["fesets"][1]) == count(fes)
    @test norm(output["fesets"][1].label - fes.label) == 0

    try
        rm("q4-4-2-conn.dat")
    catch
    end
    try
        rm("q4-4-2-label.dat")
    catch
    end
    try
        rm("q4-4-2-xyz.dat")
    catch
    end
    try
        rm("q4-4-2.mesh")
    catch
    end
    true
end
end
using .mexpmeshx1ba3
mexpmeshx1ba3.test()

module mexpmeshx1bb1
using FinEtools
using FinEtools.MeshExportModule.MESH: write_MESH
using FinEtools.MeshImportModule
using LinearAlgebra: norm
using Test
function test()
    rho = 1.21 * 1e-9# mass density
    c = 345.0 * 1000# millimeters per second
    bulk = c^2 * rho
    Lx = 1900.0# length of the box, millimeters
    Ly = 800.0 # length of the box, millimeters

    fens, fes = Q4block(Lx, Ly, 4, 2) # Mesh
    setlabel!(fes, 3)

    l = selectnode(fens, nearestto = [0.0, 0.0])
    @test l[1] == 1
    # @show fens.xyz[l, :]

    l = selectnode(fens, farthestfrom = [0.0, 0.0])
    @test l[1] == 15
    # @show fens.xyz[l, :]
    true
end
end
using .mexpmeshx1bb1
mexpmeshx1bb1.test()

# module mimptm2
# using FinEtools
# using FinEtools.MeshImportModule
# using FinEtools.MeshExportModule: VTK
# using Test
# function test()
#  output = MeshImportModule.import_GMSH(dirname(@__FILE__) * "/" * "cylinder.msh";
#    allocationchunk = 13)
#  @show length(output["fesets"])
#  @show typeof.(output["fesets"])
#  fens, fes = output["fens"], output["fesets"][1]
#    # show(fes.conn[count(fes), :])
#    File = "cylinder.vtk"
#    VTK.vtkexportmesh(File, fens, fes)
#    # rm(File)

#    # @test fes.conn[1,:] == NTuple{4,Int64}[(1, 2, 3, 4),]
#    # @test fes.conn[2,:] == NTuple{4,Int64}[(2, 5, 6, 3),]
#    # @async run(`"paraview.exe" $File`)
#  end
# end
# using .mimptm2
# mimptm2.test()

# module mimptm1
# using FinEtools
# using FinEtools.MeshImportModule
# using FinEtools.MeshExportModule: VTK
# using Test
# function test()
#   output = MeshImportModule.import_GMSH(dirname(@__FILE__) * "/" * "twoquads.msh";
#     allocationchunk = 13)
#   fens, fes = output["fens"], output["fesets"][1]
#   # show(fes.conn[count(fes), :])
#   File = "twoquads.vtk"
#   VTK.vtkexportmesh(File, fens, fes)
#   # rm(File)

#   @test fes.conn[1,:] == NTuple{4,Int64}[(1, 2, 3, 4),]
#   @test fes.conn[2,:] == NTuple{4,Int64}[(2, 5, 6, 3),]
#   # @async run(`"paraview.exe" $File`)
# end
# end
# using .mimptm1
# mimptm1.test()

module miscellaneous2mm
using FinEtools
using Test
function test()
    rho = 1.21 * 1e-9# mass density
    c = 345.0 * 1000# millimeters per second
    bulk = c^2 * rho
    Lx = 1900.0# length of the box, millimeters
    Ly = 800.0 # length of the box, millimeters

    fens, fes = Q4block(Lx, Ly, 3, 2) # Mesh
    # show(fes.conn)

    bfes = meshboundary(fes)
    @test bfes.conn == Tuple{Int64,Int64}[
        (1, 2),
        (5, 1),
        (2, 3),
        (3, 4),
        (4, 8),
        (9, 5),
        (8, 12),
        (10, 9),
        (11, 10),
        (12, 11),
    ]
end
end
using .miscellaneous2mm
miscellaneous2mm.test()

module mmmQ4blockneous2mm
using FinEtools
using FinEtools.MeshExportModule: VTK
using Test
function test()
    rho = 1.21 * 1e-9# mass density
    c = 345.0 * 1000# millimeters per second
    bulk = c^2 * rho
    Lx = 1900.0# length of the box, millimeters
    Ly = 800.0 # length of the box, millimeters

    fens, fes = Q8block(Lx, Ly, 3, 2) # Mesh
    # show(fes.conn)
    fens.xyz = xyz3(fens)
    fens.xyz[:, 3] .= (fens.xyz[:, 1] .^ 2 + fens.xyz[:, 2] .^ 2) / 1.e3
    File = "mesh.vtk"
    VTK.vtkexportmesh(File, fens, fes)
    rm(File)
    # @async run(`"paraview.exe" $File`)
end
end
using .mmmQ4blockneous2mm
mmmQ4blockneous2mm.test()

module miimportexportm
using FinEtools
using FinEtools.MeshImportModule
using FinEtools.MeshExportModule: VTK
using Test
function test()
    output = MeshImportModule.import_NASTRAN(
        dirname(@__FILE__) * "/" * "Slot-coarser.nas";
        allocationchunk = 13,
    )
    # show(fes.conn[count(fes), :])
    File = "Slot-coarser.vtk"
    VTK.vtkexportmesh(File, output["fens"], output["fesets"][1])
    rm(File)
    @test output["fesets"][1].conn[count(output["fesets"][1]), :] ==
          NTuple{10,Int64}[(143, 140, 144, 138, 361, 363, 176, 519, 781, 520)]
    # @async run(`"paraview.exe" $File`)
end
end
using .miimportexportm
miimportexportm.test()

module miimportexportm2
using FinEtools
using FinEtools.MeshImportModule
using FinEtools.MeshExportModule: VTK
using Test
function test()
    output = MeshImportModule.import_NASTRAN(
        dirname(@__FILE__) * "/" * "Slot-coarser-2.nas";
        allocationchunk = 13,
    )
    # show(fes.conn[count(fes), :])
    File = "Slot-coarser.vtk"
    VTK.vtkexportmesh(File, output["fens"], output["fesets"][1])
    rm(File)
    @test output["fesets"][1].conn[count(output["fesets"][1]), :] ==
          NTuple{10,Int64}[(143, 140, 144, 138, 361, 363, 176, 519, 781, 520)]

    # @async run(`"paraview.exe" $File`)
end
end
using .miimportexportm2
miimportexportm2.test()

module mmmLshapemmm
using FinEtools
using FinEtools.MeshExportModule: VTK
using Test
function test()
    # println("""
    # Meshing  for an L-shaped membrane.
    # """)

    t0 = time()

    W = 100.0 # strip width
    L = 200.0 # length of the strip
    nL = 5
    nW = 5
    tolerance = W / nW / 1.0e5
    Meshes = Array{Tuple{FENodeSet,AbstractFESet},1}()
    push!(Meshes, Q4quadrilateral([0.0 0.0; W W], nW, nW))
    push!(Meshes, Q4quadrilateral([-L 0.0; 0.0 W], nL, nW))
    push!(Meshes, Q4quadrilateral([0.0 -L; W 0.0], nW, nL))
    fens, outputfes = mergenmeshes(Meshes, tolerance)
    fes = cat(outputfes[1], cat(outputfes[2], outputfes[3]))
    @test count(fens) == 96
    @test fes.conn == NTuple{4,Int64}[
        (1, 2, 8, 7),
        (7, 8, 14, 13),
        (13, 14, 20, 19),
        (19, 20, 26, 25),
        (25, 26, 32, 31),
        (2, 3, 9, 8),
        (8, 9, 15, 14),
        (14, 15, 21, 20),
        (20, 21, 27, 26),
        (26, 27, 33, 32),
        (3, 4, 10, 9),
        (9, 10, 16, 15),
        (15, 16, 22, 21),
        (21, 22, 28, 27),
        (27, 28, 34, 33),
        (4, 5, 11, 10),
        (10, 11, 17, 16),
        (16, 17, 23, 22),
        (22, 23, 29, 28),
        (28, 29, 35, 34),
        (5, 6, 12, 11),
        (11, 12, 18, 17),
        (17, 18, 24, 23),
        (23, 24, 30, 29),
        (29, 30, 36, 35),
        (37, 38, 43, 42),
        (42, 43, 48, 47),
        (47, 48, 53, 52),
        (52, 53, 58, 57),
        (57, 58, 63, 62),
        (38, 39, 44, 43),
        (43, 44, 49, 48),
        (48, 49, 54, 53),
        (53, 54, 59, 58),
        (58, 59, 64, 63),
        (39, 40, 45, 44),
        (44, 45, 50, 49),
        (49, 50, 55, 54),
        (54, 55, 60, 59),
        (59, 60, 65, 64),
        (40, 41, 46, 45),
        (45, 46, 51, 50),
        (50, 51, 56, 55),
        (55, 56, 61, 60),
        (60, 61, 66, 65),
        (41, 1, 7, 46),
        (46, 7, 13, 51),
        (51, 13, 19, 56),
        (56, 19, 25, 61),
        (61, 25, 31, 66),
        (67, 68, 74, 73),
        (73, 74, 80, 79),
        (79, 80, 86, 85),
        (85, 86, 92, 91),
        (91, 92, 2, 1),
        (68, 69, 75, 74),
        (74, 75, 81, 80),
        (80, 81, 87, 86),
        (86, 87, 93, 92),
        (92, 93, 3, 2),
        (69, 70, 76, 75),
        (75, 76, 82, 81),
        (81, 82, 88, 87),
        (87, 88, 94, 93),
        (93, 94, 4, 3),
        (70, 71, 77, 76),
        (76, 77, 83, 82),
        (82, 83, 89, 88),
        (88, 89, 95, 94),
        (94, 95, 5, 4),
        (71, 72, 78, 77),
        (77, 78, 84, 83),
        (83, 84, 90, 89),
        (89, 90, 96, 95),
        (95, 96, 6, 5),
    ]
    geom = NodalField(fens.xyz)

    File = "L_shape.vtk"
    vtkexportmesh(File, connasarray(fes), geom.values, VTK.Q4)
    # @async run(`"paraview.exe" $File`)
    rm(File)
    vtkexportmesh(File, fes.conn, geom.values, VTK.Q4)
    # @async run(`"paraview.exe" $File`)
    rm(File)

    # println("Done")
    true
end
end
using .mmmLshapemmm
mmmLshapemmm.test()

module mAbaqusmmiimportmmm
using FinEtools
using FinEtools.MeshExportModule: VTK
using FinEtools.MeshImportModule
using Test
function test()
    ## Solid cylinder/taper/sphere—-temperature loading; quadratic brick mesh

    # The mesh  will be created in a very coarse representation from the
    # key points in the drawing. The first coordinate is radial, the second coordinate is axial.
    rz =
        [
            1.0 0.0#A
            1.4 0.0#B
            0.995184726672197 0.098017140329561
            1.393258617341076 0.137223996461385
            0.980785 0.195090#
            1.37309939 0.27312645
            0.956940335732209 0.290284677254462
            1.339716470025092 0.406398548156247
            0.9238795 0.38268#C
            1.2124 0.7#D
            0.7071 0.7071#E
            1.1062 1.045#F
            0.7071 (0.7071+1.79)/2#(E+H)/2
            1.0 1.39#G
            0.7071 1.79#H
            1.0 1.79#I
        ] * phun("M")
    tolerance = 1.e-6 * phun("M")

    # This is the quadrilateral mesh of the cross-section.   It will be modified and
    # refined as  we go.
    fens = FENodeSet(rz)
    fes =
        FESetQ4([1 2 4 3; 3 4 6 5; 5 6 8 7; 7 8 10 9; 9 10 12 11; 11 12 14 13; 13 14 16 15])

    nref = 1
    for ref in 1:nref
        fens, fes = Q4refine(fens, fes)
        list = selectnode(
            fens,
            distance = 1.0 + 0.1 / 2^nref,
            from = [0.0 0.0],
            inflate = tolerance,
        )
        fens.xyz[list, :] = FinEtools.MeshUtilModule.ontosphere(fens.xyz[list, :], 1.0)
    end

    ##
    # The mesh is extruded by sweeping around the axis of symmetry.
    # Only a single layer of elements is generated of internal angle
    # |angslice|.
    nLayers = 7
    angslice = 5 * pi / 16

    ##
    # First the mesh is extruded to a block whose third dimension
    # represents the angular coordinate.
    fens, fes = H8extrudeQ4(
        fens,
        fes,
        nLayers,
        (rz, k) -> [rz[1], rz[2], 0.0] - (k) / nLayers * [0.0, 0.0, angslice],
    )
    # The block is now converted  to the axially symmetric geometry by using the
    # third (angular) coordinate  to sweep out  an axially symmetric domain. The
    # ccoordinates of the nodes at this point are |rza|,  radial distance,
    # Z-coordinate, angle.
    sweep(rza) = [
        -rza[1] * sin(rza[3] + angslice / 2.0),
        rza[1] * cos(rza[3] + angslice / 2.0),
        rza[2],
    ]
    for j in axes(fens.xyz, 1)
        fens.xyz[j, :] = sweep(fens.xyz[j, :])
    end

    AE = AbaqusExporter("LE11NAFEMS_H8")
    HEADING(AE, "LE11NAFEMS: Linear bricks.")
    PART(AE, "part1")
    END_PART(AE)
    ASSEMBLY(AE, "ASSEM1")
    INSTANCE(AE, "INSTNC1", "PART1")
    NODE(AE, fens.xyz)
    ELEMENT(AE, "c3d8rh", "AllElements", 1, connasarray(fes))
    END_INSTANCE(AE)
    END_ASSEMBLY(AE)
    close(AE)


    output = MeshImportModule.import_ABAQUS("./LE11NAFEMS_H8.inp"; allocationchunk = 11)
    fens, fes = output["fens"], output["fesets"][1]

    File = "LE11NAFEMS_H8.vtk"
    VTK.vtkexportmesh(File, fens, fes)
    # @async run(`"paraview.exe" $File`)
    try
        rm(File)
    catch
    end

    try
        rm(AE.filename)
    catch
    end

    AE = AbaqusExporter("LE11NAFEMS_H8")
    HEADING(AE, "LE11NAFEMS: Linear bricks.")
    PART(AE, "part1")
    END_PART(AE)
    ASSEMBLY(AE, "ASSEM1")
    INSTANCE(AE, "INSTNC1", "PART1")
    NODE(AE, fens.xyz)
    ELEMENT(AE, "c3d8rh", "AllElements", 1, fes.conn)
    END_INSTANCE(AE)
    END_ASSEMBLY(AE)
    close(AE)


    output = MeshImportModule.import_ABAQUS("./LE11NAFEMS_H8.inp"; allocationchunk = 11)
    fens, fes = output["fens"], output["fesets"][1]

    File = "LE11NAFEMS_H8.vtk"
    VTK.vtkexportmesh(File, fens, fes)
    # @async run(`"paraview.exe" $File`)
    try
        rm(File)
    catch
    end

    try
        rm(AE.filename)
    catch
    end

    AE = AbaqusExporter("LE11NAFEMS_H8")
    HEADING(AE, "LE11NAFEMS: Linear bricks.")
    PART(AE, "part1")
    END_PART(AE)
    ASSEMBLY(AE, "ASSEM1")
    INSTANCE(AE, "INSTNC1", "PART1")
    NODE(AE, fens.xyz)
    ELEMENT(AE, "c3d8rh", "AllElements", fes.conn)
    END_INSTANCE(AE)
    END_ASSEMBLY(AE)
    close(AE)


    output = MeshImportModule.import_ABAQUS("./LE11NAFEMS_H8.inp"; allocationchunk = 11)
    fens, fes = output["fens"], output["fesets"][1]

    File = "LE11NAFEMS_H8.vtk"
    VTK.vtkexportmesh(File, fens, fes)
    # @async run(`"paraview.exe" $File`)
    try
        rm(File)
    catch
    end

    try
        rm(AE.filename)
    catch
    end
end
end
using .mAbaqusmmiimportmmm
mAbaqusmmiimportmmm.test()

module mmsmoothingm1
using FinEtools
using FinEtools.MeshExportModule: VTK
using Test
import LinearAlgebra: norm
function test()

    # println("""
    # Meshing, deforming  and smoothing
    # """)

    A = 100.0 # strip width
    N = 16
    tolerance = A / N / 1.0e5
    fens, fes = T3block(A, A, N, N)

    bnl = connectednodes(meshboundary(fes))
    for ixxxx in eachindex(bnl)
        x, y = fens.xyz[bnl[ixxxx], :]
        fens.xyz[bnl[ixxxx], 1] += A / N * sin(2 * pi * y / A)
        fens.xyz[bnl[ixxxx], 2] += -A / N * sin(2 * pi * x / A)
    end

    File = "mesh_smoothing_before.vtk"
    vtkexportmesh(File, fens, fes)
    # @async run(`"paraview.exe" $File`)
    try
        rm(File)
    catch
    end
    before = [50.0, 43.75]
    # println("$(fens.xyz[Int(N^2/2), :] )")
    @test norm(before - fens.xyz[Int(N^2 / 2), :]) < 1.e-4

    fixedv = falses(count(fens))
    fixedv[bnl] .= true
    fens = meshsmoothing(fens, fes; fixedv = fixedv, method = :laplace, npass = 100)

    after = [50.0438, 44.0315]
    # println("$(fens.xyz[Int(N^2/2), :] )")
    @test norm(fens.xyz[Int(N^2 / 2), :] - after) < 1.e-4

    geom = NodalField(fens.xyz)

    File = "mesh_smoothing_after.vtk"
    VTK.vtkexportmesh(File, connasarray(fes), geom.values, VTK.T3)
    # @async run(`"paraview.exe" $File`)
    try
        rm(File)
    catch
    end

    # println("Done")
    true
end
end
using .mmsmoothingm1
mmsmoothingm1.test()


module mmsmoothingm2
using FinEtools
using FinEtools.MeshExportModule: VTK
using Test
import LinearAlgebra: norm
function test()

    # println("""
    # Meshing, deforming  and smoothing
    # """)

    A = 100.0 # strip width
    N = 16
    tolerance = A / N / 1.0e5
    fens, fes = T3block(A, A, N, N)

    bnl = connectednodes(meshboundary(fes))
    for ixxxx in eachindex(bnl)
        x, y = fens.xyz[bnl[ixxxx], :]
        fens.xyz[bnl[ixxxx], 1] += A / N * sin(2 * pi * y / A)
        fens.xyz[bnl[ixxxx], 2] += -A / N * sin(2 * pi * x / A)
    end

    File = "mesh_smoothing_before.vtk"
    VTK.vtkexportmesh(File, fens, fes)
    # @async run(`"paraview.exe" $File`)
    try
        rm(File)
    catch
    end
    before = [50.0, 43.75]
    # println("$(fens.xyz[Int(N^2/2), :] )")
    @test norm(before - fens.xyz[Int(N^2 / 2), :]) < 1.e-4

    fixedv = falses(count(fens))
    fixedv[bnl] .= true
    fens = meshsmoothing(fens, fes; fixedv = fixedv, method = :taubin, npass = 100)

    after = [50.0059, 43.6281]
    # println("$(fens.xyz[Int(N^2/2), :] )")
    @test norm(fens.xyz[Int(N^2 / 2), :] - after) < 1.e-4

    geom = NodalField(fens.xyz)

    File = "mesh_smoothing_after.vtk"
    VTK.vtkexportmesh(File, connasarray(fes), geom.values, VTK.T3)
    # @async run(`"paraview.exe" $File`)
    try
        rm(File)
    catch
    end

    # println("Done")
    true
end
end
using .mmsmoothingm2
mmsmoothingm2.test()



module mAbaqusmmiimport_1m
using FinEtools
using FinEtools.MeshExportModule: VTK
using FinEtools.MeshImportModule
using Test
function test()



    ## Solid cylinder/taper/sphere—-temperature loading; quadratic brick mesh

    # The mesh  will be created in a very coarse representation from the
    # key points in the drawing. The first coordinate is radial, the second coordinate is axial.
    rz =
        [
            1.0 0.0#A
            1.4 0.0#B
            0.995184726672197 0.098017140329561
            1.393258617341076 0.137223996461385
            0.980785 0.195090#
            1.37309939 0.27312645
            0.956940335732209 0.290284677254462
            1.339716470025092 0.406398548156247
            0.9238795 0.38268#C
            1.2124 0.7#D
            0.7071 0.7071#E
            1.1062 1.045#F
            0.7071 (0.7071+1.79)/2#(E+H)/2
            1.0 1.39#G
            0.7071 1.79#H
            1.0 1.79#I
        ] * phun("M")
    tolerance = 1.e-6 * phun("M")

    # This is the quadrilateral mesh of the cross-section.   It will be modified and
    # refined as  we go.
    fens = FENodeSet(rz)
    fes =
        FESetQ4([1 2 4 3; 3 4 6 5; 5 6 8 7; 7 8 10 9; 9 10 12 11; 11 12 14 13; 13 14 16 15])

    nref = 1
    for ref in 1:nref
        fens, fes = Q4refine(fens, fes)
        list = selectnode(
            fens,
            distance = 1.0 + 0.1 / 2^nref,
            from = [0.0 0.0],
            inflate = tolerance,
        )
        fens.xyz[list, :] = FinEtools.MeshUtilModule.ontosphere(fens.xyz[list, :], 1.0)
    end

    ##
    # The mesh is extruded by sweeping around the axis of symmetry.
    # Only a single layer of elements is generated of internal angle
    # |angslice|.
    nLayers = 7
    angslice = 5 * pi / 16

    ##
    # First the mesh is extruded to a block whose third dimension
    # represents the angular coordinate.
    fens, fes = H8extrudeQ4(
        fens,
        fes,
        nLayers,
        (rz, k) -> [rz[1], rz[2], 0.0] - (k) / nLayers * [0.0, 0.0, angslice],
    )
    # The block is now converted  to the axially symmetric geometry by using the
    # third (angular) coordinate  to sweep out  an axially symmetric domain. The
    # ccoordinates of the nodes at this point are |rza|,  radial distance,
    # Z-coordinate, angle.
    sweep(rza) = [
        -rza[1] * sin(rza[3] + angslice / 2.0),
        rza[1] * cos(rza[3] + angslice / 2.0),
        rza[2],
    ]
    for j in axes(fens.xyz, 1)
        fens.xyz[j, :] = sweep(fens.xyz[j, :])
    end

    AE = AbaqusExporter("LE11NAFEMS_H8")
    HEADING(AE, "LE11NAFEMS: Linear bricks.")
    PART(AE, "part1")
    END_PART(AE)
    ASSEMBLY(AE, "ASSEM1")
    INSTANCE(AE, "INSTNC1", "PART1")
    NODE(AE, fens.xyz)
    ELEMENT(AE, "c3d8rh", "AllElements", 1, connasarray(fes))
    ELSET_ELSET(AE, "AllSolids", collect(eachindex(fes)))
    END_INSTANCE(AE)
    END_ASSEMBLY(AE)
    close(AE)


    output = MeshImportModule.import_ABAQUS("./LE11NAFEMS_H8.inp"; allocationchunk = 12)
    fens, fes = output["fens"], output["fesets"][1]

    File = "LE11NAFEMS_H8.vtk"
    VTK.vtkexportmesh(File, fens, fes)
    # @async run(`"paraview.exe" $File`)
    try
        rm(File)
    catch
    end
    try
        rm(AE.filename)
    catch
    end
end
end
using .mAbaqusmmiimport_1m
mAbaqusmmiimport_1m.test()


module mmfflood1
using FinEtools
using Test
function test()
    h = 0.05 * phun("M")
    l = 10 * h
    Rmed = h / 0.2
    psi = 0    # Cylinder
    nh = 12
    nl = 40
    nc = 120
    nh = 6
    nl = 20
    nc = 60
    # nh = 3; nl  = 8; nc = 30;
    tolerance = h / nh / 100

    fens, fes = H8block(h, l, 2.0 * pi, nh, nl, nc)
    # Shape into a cylinder
    R = zeros(3, 3)
    for i in eachindex(fens)
        x, y, z = fens.xyz[i, :]
        rotmat3!(R, [0, z, 0])
        Q = [
            cos(psi * pi / 180) sin(psi * pi / 180) 0
            -sin(psi * pi / 180) cos(psi * pi / 180) 0
            0 0 1
        ]
        fens.xyz[i, :] = reshape([x + Rmed - h / 2, y - l / 2, 0], 1, 3) * Q * R
    end
    # println("  before merging  = $(count(fens))")
    @test count(fens) == 8967
    # File =  "unit_cube_modes.vtk"
    # vtkexportmesh(File, fens, fes)
    # @async run(`"paraview.exe" $File`)

    candidates = selectnode(
        fens,
        box = boundingbox([Rmed-h -Inf 0.0; Rmed+h +Inf 0.0]),
        inflate = tolerance,
    )
    fens, fes = mergenodes(fens, fes, tolerance, candidates)
    @test count(fens) == 8820

    bfes = meshboundary(fes)
    startnode = bfes.conn[1][1]
    lb = selectelem(fens, bfes, flood = true, startnode = startnode)
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
using .mmfflood1
mmfflood1.test()


module mmfflood2
using FinEtools
using Test
function test()
    h = 0.05 * phun("M")
    l = 10 * h
    Rmed = h / 0.2
    psi = 0    # Cylinder
    # nh = 12; nl  = 40; nc = 120;
    # nh = 6; nl  = 20; nc = 60;
    nh = 3
    nl = 8
    nc = 30
    nh = 2
    nl = 2
    nc = 10
    tolerance = h / nh / 100

    fens, fes = H20block(h, l, 2.0 * pi, nh, nl, nc)
    # Shape into a cylinder
    R = zeros(3, 3)
    for i in eachindex(fens)
        x, y, z = fens.xyz[i, :]
        rotmat3!(R, [0, z, 0])
        Q = [
            cos(psi * pi / 180) sin(psi * pi / 180) 0
            -sin(psi * pi / 180) cos(psi * pi / 180) 0
            0 0 1
        ]
        fens.xyz[i, :] = reshape([x + Rmed - h / 2, y - l / 2, 0], 1, 3) * Q * R
    end
    # println("  before merging  = $(count(fens))")
    # @test count(fens) == 4025
    # File =  "unit_cube_modes.vtk"
    # vtkexportmesh(File, fens, fes)
    # @async run(`"paraview.exe" $File`)

    candidates = selectnode(
        fens,
        box = boundingbox([Rmed-h -Inf 0.0; Rmed+h +Inf 0.0]),
        inflate = tolerance,
    )
    fens, fes = mergenodes(fens, fes, tolerance, candidates)
    # @test count(fens) == 3930
    # @test count(fens) == 8820

    bfes = meshboundary(fes)
    startnode = bfes.conn[1][1]
    lb = selectelem(fens, bfes, flood = true, startnode = startnode)
    # println("$(lb)")
    @test length(lb) == 80

    bbfes = meshboundary(subset(bfes, lb))
    @test count(bbfes) == 0
    # File =  "unit_cube_modes.vtk"
    # vtkexportmesh(File, fens, subset(bfes, lb))
    # @async run(`"paraview.exe" $File`)
end
end
using .mmfflood2
mmfflood2.test()

module mt4orientation2
using FinEtools
using FinEtools.MeshImportModule
using FinEtools.MeshExportModule: VTK
using Test
function test()
    xs = collect(linearspace(0.0, 2.0, 5))
    ys = collect(linearspace(0.0, 1.0, 6) .^ 2)
    zs = collect(linearspace(0.0, 1.0, 3))
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
    # VTK.vtkexportmesh(File, fens, fes)
    # rm(File)
    # @test fes.conn[count(fes), :] == [143, 140, 144, 138, 361, 363, 176, 519, 781, 520]
    # @async run(`"paraview.exe" $File`)
end
end
using .mt4orientation2
mt4orientation2.test()

module mmAbaqusexport3
using FinEtools
using FinEtools.MeshExportModule
using FinEtools.MeshImportModule
using Test
function test()
    ## Solid cylinder/taper/sphere—-temperature loading; quadratic brick mesh

    # The mesh  will be created in a very coarse representation from the
    # key points in the drawing. The first coordinate is radial, the second coordinate is axial.
    rz =
        [
            1.0 0.0#A
            1.4 0.0#B
            0.995184726672197 0.098017140329561
            1.393258617341076 0.137223996461385
            0.980785 0.195090#
            1.37309939 0.27312645
            0.956940335732209 0.290284677254462
            1.339716470025092 0.406398548156247
            0.9238795 0.38268#C
            1.2124 0.7#D
            0.7071 0.7071#E
            1.1062 1.045#F
            0.7071 (0.7071+1.79)/2#(E+H)/2
            1.0 1.39#G
            0.7071 1.79#H
            1.0 1.79#I
        ] * phun("M")
    tolerance = 1.e-6 * phun("M")

    # This is the quadrilateral mesh of the cross-section.   It will be modified and
    # refined as  we go.
    fens = FENodeSet(rz)
    fes =
        FESetQ4([1 2 4 3; 3 4 6 5; 5 6 8 7; 7 8 10 9; 9 10 12 11; 11 12 14 13; 13 14 16 15])

    nref = 1
    for ref in 1:nref
        fens, fes = Q4refine(fens, fes)
        list = selectnode(
            fens,
            distance = 1.0 + 0.1 / 2^nref,
            from = [0.0 0.0],
            inflate = tolerance,
        )
        fens.xyz[list, :] = FinEtools.MeshUtilModule.ontosphere(fens.xyz[list, :], 1.0)
    end

    ##
    # The mesh is extruded by sweeping around the axis of symmetry.
    # Only a single layer of elements is generated of internal angle
    # |angslice|.
    nLayers = 7
    angslice = 5 * pi / 16

    ##
    # First the mesh is extruded to a block whose third dimension
    # represents the angular coordinate.
    fens, fes = H8extrudeQ4(
        fens,
        fes,
        nLayers,
        (rz, k) -> [rz[1], rz[2], 0.0] - (k) / nLayers * [0.0, 0.0, angslice],
    )
    # The block is now converted  to the axially symmetric geometry by using the
    # third (angular) coordinate  to sweep out  an axially symmetric domain. The
    # ccoordinates of the nodes at this point are |rza|,  radial distance,
    # Z-coordinate, angle.
    sweep(rza) = [
        -rza[1] * sin(rza[3] + angslice / 2.0),
        rza[1] * cos(rza[3] + angslice / 2.0),
        rza[2],
    ]
    for j in axes(fens.xyz, 1)
        fens.xyz[j, :] = sweep(fens.xyz[j, :])
    end

    AE = AbaqusExporter("LE11NAFEMS_H8_B.inp")
    HEADING(AE, "LE11NAFEMS: Linear bricks.")
    PART(AE, "part1")
    END_PART(AE)
    ASSEMBLY(AE, "ASSEM1")
    INSTANCE(AE, "INSTNC1", "PART1")
    NODE(AE, fens.xyz)
    ELEMENT(AE, "c3d8rh", "AllElements", 1, connasarray(fes))
    END_INSTANCE(AE)
    END_ASSEMBLY(AE)
    STEP_PERTURBATION_BUCKLE(AE, 1)
    CLOAD(AE, 1, 2, 10.0)
    EL_PRINT(AE, "AllElements", "S")
    ENERGY_PRINT(AE)
    END_STEP(AE)
    close(AE)


    output = MeshImportModule.import_ABAQUS("./LE11NAFEMS_H8_B.inp"; allocationchunk = 11)
    fens, fes = output["fens"], output["fesets"][1]

    File = "LE11NAFEMS_H8.vtk"
    vtkexportmesh(File, fens, fes)
    # @async run(`"paraview.exe" $File`)
    try
        rm(File)
    catch
    end

    try
        rm(AE.filename)
    catch
    end
end
end
using .mmAbaqusexport3
mmAbaqusexport3.test()


module mttriangles13
using FinEtools
using FinEtools.MeshExportModule
using Test
function test()
    xs = collect(linearspace(1.0, 3.0, 4))
    ys = collect(linearspace(-1.0, 3.0, 4))
    fens, fes = T3blockx(xs, ys, :a)
    @test count(fes) == 3 * 3 * 2
    fens, fes = T3blockx(xs, ys, :b)
    @test count(fes) == 3 * 3 * 2

    # File = "playground.vtk"
    # VTK.vtkexportmesh(File, fens, fes)
    # @async run(`"paraview.exe" $File`)
    # try rm(File) catch end

end
end
using .mttriangles13
mttriangles13.test()

module mttriangles14
using FinEtools
using FinEtools.MeshExportModule
using Test
function test()
    xs = collect(linearspace(1.0, 3.0, 4))
    ys = collect(linearspace(-1.0, 3.0, 5))
    fens, fes = Q4blockx(xs, ys)
    fens, fes = Q4toT3(fens, fes)
    @test count(fes) == 3 * 4 * 2
    fens, fes = Q4blockx(xs, ys)
    fens, fes = Q4toT3(fens, fes, :alternate)
    @test count(fes) == 3 * 4 * 2

    # File = "playground.vtk"
    # VTK.vtkexportmesh(File, fens, fes)
    # @async run(`"paraview.exe" $File`)
    # try rm(File) catch end

end
end
using .mttriangles14
mttriangles14.test()

module mselecte1
using FinEtools

using Test
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
    W = 1.1
    L = 12.0
    t = 0.32
    nl = 2
    nt = 1
    nw = 1
    ref = 7
    tolerance = t / 1000

    fens, fes = H8block(L, W, t, nl * ref, nw * ref, nt * ref)

    # Reshape into a twisted beam shape
    for i in eachindex(fens)
        a = fens.xyz[i, 1] / L * (pi / 2)
        y = fens.xyz[i, 2] - (W / 2)
        z = fens.xyz[i, 3] - (t / 2)
        fens.xyz[i, :] = [fens.xyz[i, 1], y * cos(a) - z * sin(a), y * sin(a) + z * cos(a)]
    end

    # Clamped end of the beam
    l1 = selectnode(fens; box = [0 0 -100 * W 100 * W W W], inflate = tolerance)
    # display(l1)
    @test isempty(l1)

    # # Traction on the opposite edge
    boundaryfes = meshboundary(fes)
    Toplist = selectelem(
        fens,
        boundaryfes,
        facing = true,
        direction = [-1.0 0.0 0.0],
        dotmin = 0.999,
    )
    # display(Toplist)
    @test length(Toplist) == 49
    Toplist = selectelem(
        fens,
        boundaryfes,
        box = [L L -Inf Inf -Inf Inf],
        allin = true,
        inflate = tolerance,
    )
    # display(Toplist)
    @test length(Toplist) == 49
    Toplist = selectelem(
        fens,
        boundaryfes,
        box = [L L -Inf Inf -Inf Inf],
        allin = false,
        inflate = tolerance,
    )
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
using .mselecte1
mselecte1.test()


module mselecte2
using FinEtools

using Test
function test()

    W = 10.0
    L = 12.0
    nl = 2
    nw = 1
    ref = 7
    tolerance = W / 1000

    fens, fes = Q4block(L, W, nl * ref, nw * ref)

    boundaryfes = meshboundary(fes)
    alist =
        selectelem(fens, boundaryfes, facing = true, direction = [0.0 -1.0], dotmin = 0.999)
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
using .mselecte2
mselecte2.test()

module mMeyer_Piening_1
using FinEtools
using FinEtools.MeshUtilModule
using Test
function test()

    # Reference results from:
    # [1] Application of the Elasticity Solution
    # to Linear Sandwich Beam, Plate
    # and Shell Analyses
    # H.-R. MEYER -PIENING
    # Journal of SANDWICH STRUCTURES AND MATERIALS , Vol. 6—July 2004

    Lx = 5.0 * phun("mm") # length  of loaded rectangle
    Ly = 20.0 * phun("mm") # length  of loaded rectangle
    Sx = 100.0 * phun("mm") # span of the plate
    Sy = 200.0 * phun("mm") # span of the plate

    # Here we define the layout and the thicknesses of the layers.
    angles = vec([0.0 0.0 0.0])
    ts = vec([0.5 11.4 0.1]) * phun("mm") # layer thicknesses
    TH = sum(ts) # total thickness of the plate

    tolerance = 0.0001 * TH

    # Select how find the mesh should be
    Refinement = 5
    nL = Refinement * 1
    nSx = nL + Refinement * 4
    nSy = 2 * nSx

    # Each layer is modeled with a single element.
    nts = Refinement * [1, 2, 1]# number of elements per layer
    strength = 1.5
    xs = unique(
        vcat(
            reverse(collect(MeshUtilModule.gradedspace(Lx / 2, 0.0, nL + 1, strength))),
            collect(MeshUtilModule.gradedspace(Lx / 2, Sx / 2, nSx - nL + 1, strength)),
        ),
    )
    ys = unique(
        vcat(
            reverse(collect(MeshUtilModule.gradedspace(Ly / 2, 0.0, nL + 1, strength))),
            collect(MeshUtilModule.gradedspace(Ly / 2, Sy / 2, nSy - nL + 1, strength)),
        ),
    )

    fens, fes = H8layeredplatex(xs, ys, ts, nts)
    @test count(fens) == 27846
    @test count(fes) == 25000
    rls = selectelem(fens, fes, label = 1)
    @test length(rls) == 6250
end
end
using .mMeyer_Piening_1
mMeyer_Piening_1.test()


module mMeyer_Piening_2
using FinEtools
using FinEtools.MeshUtilModule
using Test
function test()

    # Reference results from:
    # [1] Application of the Elasticity Solution
    # to Linear Sandwich Beam, Plate
    # and Shell Analyses
    # H.-R. MEYER -PIENING
    # Journal of SANDWICH STRUCTURES AND MATERIALS , Vol. 6—July 2004

    Lx = 5.0 * phun("mm") # length  of loaded rectangle
    Ly = 20.0 * phun("mm") # length  of loaded rectangle
    Sx = 100.0 * phun("mm") # span of the plate
    Sy = 200.0 * phun("mm") # span of the plate

    # Here we define the layout and the thicknesses of the layers.
    angles = vec([0.0 0.0 0.0])
    ts = vec([0.5 11.4 0.1]) * phun("mm") # layer thicknesses
    TH = sum(ts) # total thickness of the plate

    tolerance = 0.0001 * TH

    # Select how find the mesh should be
    Refinement = 5
    nL = Refinement * 1
    nSx = nL + Refinement * 4
    nSy = 2 * nSx

    # Each layer is modeled with a single element.
    nts = Refinement * [1, 2, 1]# number of elements per layer
    strength = 1.5
    xs = unique(
        vcat(
            reverse(collect(MeshUtilModule.gradedspace(Lx / 2, 0.0, nL + 1, strength))),
            collect(MeshUtilModule.gradedspace(Lx / 2, Sx / 2, nSx - nL + 1, strength)),
        ),
    )
    ys = unique(
        vcat(
            reverse(collect(MeshUtilModule.gradedspace(Ly / 2, 0.0, nL + 1, strength))),
            collect(MeshUtilModule.gradedspace(Ly / 2, Sy / 2, nSy - nL + 1, strength)),
        ),
    )

    fens, fes = T10layeredplatex(xs, ys, ts, nts)

    @test count(fens) == 211191
    @test count(fes) == 150000
    rls = selectelem(fens, fes, label = 3)
    @test length(rls) == 37500
end
end
using .mMeyer_Piening_2
mMeyer_Piening_2.test()


module mMeyer_Piening_3
using FinEtools
using FinEtools.MeshUtilModule
using Test
function test()

    # Reference results from:
    # [1] Application of the Elasticity Solution
    # to Linear Sandwich Beam, Plate
    # and Shell Analyses
    # H.-R. MEYER -PIENING
    # Journal of SANDWICH STRUCTURES AND MATERIALS , Vol. 6—July 2004

    Lx = 5.0 * phun("mm") # length  of loaded rectangle
    Ly = 20.0 * phun("mm") # length  of loaded rectangle
    Sx = 100.0 * phun("mm") # span of the plate
    Sy = 200.0 * phun("mm") # span of the plate

    # Here we define the layout and the thicknesses of the layers.
    angles = vec([0.0 0.0 0.0])
    ts = vec([0.5 11.4 0.1]) * phun("mm") # layer thicknesses
    TH = sum(ts) # total thickness of the plate

    tolerance = 0.0001 * TH

    # Select how find the mesh should be
    Refinement = 5
    nL = Refinement * 1
    nSx = nL + Refinement * 4
    nSy = 2 * nSx

    # Each layer is modeled with a single element.
    nts = Refinement * [1, 2, 1]# number of elements per layer
    strength = 1.5
    xs = unique(
        vcat(
            reverse(collect(MeshUtilModule.gradedspace(Lx / 2, 0.0, nL + 1, strength))),
            collect(MeshUtilModule.gradedspace(Lx / 2, Sx / 2, nSx - nL + 1, strength)),
        ),
    )
    ys = unique(
        vcat(
            reverse(collect(MeshUtilModule.gradedspace(Ly / 2, 0.0, nL + 1, strength))),
            collect(MeshUtilModule.gradedspace(Ly / 2, Sy / 2, nSy - nL + 1, strength)),
        ),
    )

    fens, fes = T10layeredplatex(xs, ys, ts, nts, :b)

    @test count(fens) == 211191
    @test count(fes) == 150000
    rls = selectelem(fens, fes, label = 3)
    @test length(rls) == 37500
end
end
using .mMeyer_Piening_3
mMeyer_Piening_3.test()

module mmtetblocksmm
using FinEtools
using Test
function test()
    A = 5.0 * phun("mm") # length  of loaded rectangle
    B = 20.0 * phun("mm") # length  of loaded rectangle
    C = 100.0 * phun("mm") # span of the plate

    # Select how find the mesh should be
    Refinement = 5
    nA, nB, nC = Refinement * 1, Refinement * 2, Refinement * 1
    xs = reshape(collect(linearspace(0.0, A, nA + 1)), nA + 1, 1)
    ys = reshape(collect(linearspace(0.0, A, nB + 1)), nB + 1, 1)
    zs = reshape(collect(linearspace(0.0, A, nC + 1)), nC + 1, 1)
    fens, fes = T10blockx(xs, ys, zs, :b)
    # println("$(count(fens))")
    # println("$(count(fes))")
    @test count(fens) == 2541
    @test count(fes) == 1500

    fens, fes = T10blockx(xs, ys, zs, :a)
    # println("$(count(fens))")
    # println("$(count(fes))")
    @test count(fens) == 2541
    @test count(fes) == 1500


    fens, fes = T10blockx(vec(xs), vec(ys), vec(zs), :a)
    # println("$(count(fens))")
    # println("$(count(fes))")
    @test count(fens) == 2541
    @test count(fes) == 1500


    fens, fes = T4blockx(xs, ys, zs, :ca)
    # println("$(count(fens))")
    # println("$(count(fes))")
    @test count(fens) == 396
    @test count(fes) == 1250

    fens, fes = T4blockx(vec(xs), vec(ys), vec(zs), :ca)
    @test count(fens) == 396
    @test count(fes) == 1250

    return true
end
end
using .mmtetblocksmm
mmtetblocksmm.test()

module momap2para1
using FinEtools
using Test
import LinearAlgebra: norm
function test()
    X = [-1.0 -1.0; 2.0 0.4; 1.0 2.3; -2.0 1.0]
    fens = FENodeSet(X)
    fes = FESetQ4(reshape([1 2 3 4], 1, 4))
    pt = [0.1 0.2]
    pc, success = map2parametric(fes, fens.xyz[[i for i in fes.conn[1]], :], vec(pt))
    N = bfun(fes, pc)
    pt1 = transpose(N) * fens.xyz[[i for i in fes.conn[1]], :]
    # println("pt = $(pt)")
    # println("pt1 = $(pt1)")
    @test norm(pt - pt1) < 1.0e-7
end
end
using .momap2para1
momap2para1.test()

module momap2para2
using FinEtools
using Test
import LinearAlgebra: norm
function test()
    X = [-1.0 -1.0; 2.0 0.4; 1.0 2.3; -2.0 1.0]
    fens = FENodeSet(X)
    fes = FESetQ4(reshape([1 2 3 4], 1, 4))
    pt = [3.1 0.2]
    pc, success = map2parametric(
        fes,
        fens.xyz[[i for i in fes.conn[1]], :],
        vec(pt);
        tolerance = 0.000001,
        maxiter = 7,
    )
    # println("pc = $(pc)")
    N = bfun(fes, pc)
    pt1 = transpose(N) * fens.xyz[[i for i in fes.conn[1]], :]
    # println("pt = $(pt)")
    # println("pt1 = $(pt1)")
    @test norm(pt - pt1) < 1.0e-7
    isin = inparametric(fes, pc)
    # println("isin = $(isin)")
    @test !isin
end
end
using .momap2para2
momap2para2.test()

module momap2para3
using FinEtools
using FinEtools.MeshSelectionModule: vselect
using FinEtools.MeshExportModule
using Test
import LinearAlgebra: norm
function test()
    A = 50.0 * phun("m") # length  of loaded rectangle
    B = 200.0 * phun("m") # length  of loaded rectangle
    C = 100.0 * phun("m") # span of the plate

    # Select how find the mesh should be
    Refinement = 2
    nA, nB, nC = Refinement * 1, Refinement * 2, Refinement * 4
    xs = reshape(collect(linearspace(0.0, A, nA + 1)), nA + 1, 1)
    ys = reshape(collect(linearspace(0.0, B, nB + 1)), nB + 1, 1)
    zs = reshape(collect(linearspace(0.0, C, nC + 1)), nC + 1, 1)
    fensc, fesc = T10blockx(xs, ys, zs, :b)
    fc = NodalField(zeros(count(fensc), 1))
    for i in eachindex(fensc)
        x, y, z = fensc.xyz[i, :]
        fc.values[i, :] .= sin(2 * x / A) * cos(6.5 * y / B) * sin(3 * z / C - 1.0)
    end
    File = "momap2para3-coarse.vtk"
    vtkexportmesh(File, fensc, fesc; scalars = [("fc", fc.values)])
    # @async run(`"paraview.exe" $File`)
    try
        rm(File)
    catch
    end

    Refinement = Refinement + 1
    nA, nB, nC = Refinement * 1, Refinement * 2, Refinement * 4
    xs = reshape(collect(linearspace(0.0, A, nA + 1)), nA + 1, 1)
    ys = reshape(collect(linearspace(0.0, B, nB + 1)), nB + 1, 1)
    zs = reshape(collect(linearspace(0.0, C, nC + 1)), nC + 1, 1)
    fensf, fesf = T10blockx(xs, ys, zs, :b)
    ff = NodalField(zeros(count(fensf), 1))
    tolerance = min(A / nA, B / nB, C / nC) / 1000.0

    referenceff = NodalField(zeros(count(fensf), 1))
    for i in eachindex(fensf)
        x, y, z = fensf.xyz[i, :]
        referenceff.values[i, :] .= sin(2 * x / A) * cos(6.5 * y / B) * sin(3 * z / C - 1.0)
    end
    File = "momap2para3-reference.vtk"
    vtkexportmesh(File, fensf, fesf; scalars = [("referenceff", referenceff.values)])
    # @async run(`"paraview.exe" $File`)
    try
        rm(File)
    catch
    end

    ff = transferfield!(ff, fensf, fesf, fc, fensc, fesc, tolerance)
    File = "momap2para3-fine.vtk"
    vtkexportmesh(File, fensf, fesf; scalars = [("ff", ff.values)])
    # @async run(`"paraview.exe" $File`)
    try
        rm(File)
    catch
    end

    diffff = NodalField(referenceff.values - ff.values)
    femm = FEMMBase(IntegDomain(fesf, SimplexRule(3, 4)))
    geom = NodalField(fensf.xyz)
    error = integratefieldfunction(femm, geom, diffff, (x, v) -> norm(v); initial = 0.0)
    ref = integratefieldfunction(femm, geom, referenceff, (x, v) -> norm(v); initial = 0.0)
    # println("error/ref = $(error/ref)")
    @test abs(error / ref - 0.02541369940759616) < 1.0e-4
end
end
using .momap2para3
momap2para3.test()

module momap2para4
using FinEtools
using FinEtools.MeshSelectionModule: vselect
using FinEtools.MeshExportModule
using Test
import LinearAlgebra: norm
function test()
    A = 50.0 * phun("m") # length  of loaded rectangle
    B = 200.0 * phun("m") # length  of loaded rectangle
    C = 100.0 * phun("m") # span of the plate

    # Select how find the mesh should be
    Refinement = 2
    nA, nB, nC = Refinement * 1, Refinement * 6, Refinement * 4
    xs = reshape(collect(linearspace(0.0, A, nA + 1)), nA + 1, 1)
    ys = reshape(collect(linearspace(0.0, B, nB + 1)), nB + 1, 1)
    zs = reshape(collect(linearspace(0.0, C, nC + 1)), nC + 1, 1)
    fensc, fesc = T10blockx(xs, ys, zs, :b)

    centroidpc = centroidparametric(fesc)
    N = bfun(fesc, centroidpc)
    NT = transpose(N)

    fc = ElementalField(zeros(count(fesc), 1))
    for i in eachindex(fesc)
        c = [k for k in fesc.conn[i]]
        centroid = NT * fensc.xyz[c, :]
        x, y, z = centroid
        fc.values[i, :] .= sin(2 * x / A) * cos(6.5 * y / B) * sin(3 * z / C - 1.0)
    end
    File = "momap2para3-coarse.vtk"
    vtkexportmesh(File, fensc, fesc; scalars = [("fc", fc.values)])
    # @async run(`"paraview.exe" $File`)
    try
        rm(File)
    catch
    end

    Refinement = Refinement + 1
    nA, nB, nC = Refinement * 1, Refinement * 6, Refinement * 4
    xs = reshape(collect(linearspace(0.0, A, nA + 1)), nA + 1, 1)
    ys = reshape(collect(linearspace(0.0, B, nB + 1)), nB + 1, 1)
    zs = reshape(collect(linearspace(0.0, C, nC + 1)), nC + 1, 1)
    fensf, fesf = T10blockx(xs, ys, zs, :b)
    tolerance = min(A / nA, B / nB, C / nC) / 1000.0

    ff = ElementalField(zeros(count(fesf), 1))
    referenceff = ElementalField(zeros(count(fesf), 1))
    for i in eachindex(fesf)
        c = [k for k in fesf.conn[i]]
        centroid = NT * fensf.xyz[c, :]
        x, y, z = centroid
        referenceff.values[i, :] .= sin(2 * x / A) * cos(6.5 * y / B) * sin(3 * z / C - 1.0)
    end
    File = "momap2para3-reference.vtk"
    vtkexportmesh(File, fensf, fesf; scalars = [("referenceff", referenceff.values)])
    # @async run(`"paraview.exe" $File`)
    try
        rm(File)
    catch
    end

    ff = transferfield!(ff, fensf, fesf, fc, fensc, fesc, tolerance)
    File = "momap2para3-fine.vtk"
    vtkexportmesh(File, fensf, fesf; scalars = [("ff", ff.values)])
    # @async run(`"paraview.exe" $File`)
    try
        rm(File)
    catch
    end

    diffff = ElementalField(referenceff.values - ff.values)
    femm = FEMMBase(IntegDomain(fesf, SimplexRule(3, 4)))
    geom = NodalField(fensf.xyz)
    error = integratefieldfunction(femm, geom, diffff, (x, v) -> norm(v); initial = 0.0)
    ref = integratefieldfunction(femm, geom, referenceff, (x, v) -> norm(v); initial = 0.0)
    # println("error/ref = $(error/ref)")
    @test abs(error / ref - 0.19808425992688541) < 1.0e-4
end
end
using .momap2para4
momap2para4.test()

module momap2para5
using FinEtools
using FinEtools.MeshSelectionModule: vselect
using FinEtools.MeshExportModule
using Test
import LinearAlgebra: norm
function test()
    A = 50.0 * phun("m") # length  of loaded rectangle
    B = 200.0 * phun("m") # length  of loaded rectangle
    C = 100.0 * phun("m") # span of the plate

    Meshing = T4blockx
    # Select how find the mesh should be
    Refinement = 2
    nA, nB, nC = Refinement * 1, Refinement * 6, Refinement * 4
    xs = reshape(collect(linearspace(0.0, A, nA + 1)), nA + 1, 1)
    ys = reshape(collect(linearspace(0.0, B, nB + 1)), nB + 1, 1)
    zs = reshape(collect(linearspace(0.0, C, nC + 1)), nC + 1, 1)
    fensc, fesc = Meshing(xs, ys, zs, :b)

    centroidpc = centroidparametric(fesc)
    N = bfun(fesc, centroidpc)
    NT = transpose(N)

    fc = ElementalField(zeros(count(fesc), 1))
    for i in eachindex(fesc)
        c = [k for k in fesc.conn[i]]
        centroid = NT * fensc.xyz[c, :]
        x, y, z = centroid
        fc.values[i, :] .= sin(2 * x / A) * cos(6.5 * y / B) * sin(3 * z / C - 1.0)
    end
    File = "momap2para3-coarse.vtk"
    vtkexportmesh(File, fensc, fesc; scalars = [("fc", fc.values)])
    # @async run(`"paraview.exe" $File`)
    try
        rm(File)
    catch
    end

    Refinement = Refinement + 1
    nA, nB, nC = Refinement * 1, Refinement * 6, Refinement * 4
    xs = reshape(collect(linearspace(0.0, A, nA + 1)), nA + 1, 1)
    ys = reshape(collect(linearspace(0.0, B, nB + 1)), nB + 1, 1)
    zs = reshape(collect(linearspace(0.0, C, nC + 1)), nC + 1, 1)
    fensf, fesf = Meshing(xs, ys, zs, :b)
    tolerance = min(A / nA, B / nB, C / nC) / 1000.0

    ff = ElementalField(zeros(count(fesf), 1))
    referenceff = ElementalField(zeros(count(fesf), 1))
    for i in eachindex(fesf)
        c = [k for k in fesf.conn[i]]
        centroid = NT * fensf.xyz[c, :]
        x, y, z = centroid
        referenceff.values[i, :] .= sin(2 * x / A) * cos(6.5 * y / B) * sin(3 * z / C - 1.0)
    end
    File = "momap2para3-reference.vtk"
    vtkexportmesh(File, fensf, fesf; scalars = [("referenceff", referenceff.values)])
    # @async run(`"paraview.exe" $File`)
    try
        rm(File)
    catch
    end

    ff = transferfield!(ff, fensf, fesf, fc, fensc, fesc, tolerance)
    File = "momap2para3-fine.vtk"
    vtkexportmesh(File, fensf, fesf; scalars = [("ff", ff.values)])
    # @async run(`"paraview.exe" $File`)
    try
        rm(File)
    catch
    end

    diffff = ElementalField(referenceff.values - ff.values)
    femm = FEMMBase(IntegDomain(fesf, SimplexRule(3, 4)))
    geom = NodalField(fensf.xyz)
    error = integratefieldfunction(femm, geom, diffff, (x, v) -> norm(v); initial = 0.0)
    ref = integratefieldfunction(femm, geom, referenceff, (x, v) -> norm(v); initial = 0.0)
    # println("error/ref = $(error/ref)")
    @test abs(error / ref - 0.19808425992688541) < 1.0e-4
end
end
using .momap2para5
momap2para5.test()


module momap2para6
using FinEtools
using FinEtools.MeshSelectionModule: vselect
using FinEtools.MeshExportModule
using Test
import LinearAlgebra: norm
function test()
    A = 50.0 * phun("m") # length  of loaded rectangle
    B = 200.0 * phun("m") # length  of loaded rectangle
    C = 100.0 * phun("m") # span of the plate

    Meshing = H20blockx
    # Select how find the mesh should be
    Refinement = 2
    nA, nB, nC = Refinement * 1, Refinement * 6, Refinement * 4
    xs = collect(linearspace(0.0, A, nA + 1))
    ys = collect(linearspace(0.0, B, nB + 1))
    zs = collect(linearspace(0.0, C, nC + 1))
    fensc, fesc = Meshing(xs, ys, zs)

    centroidpc = centroidparametric(fesc)
    N = bfun(fesc, centroidpc)
    NT = transpose(N)

    fc = ElementalField(zeros(count(fesc), 1))
    for i in eachindex(fesc)
        c = [k for k in fesc.conn[i]]
        centroid = NT * fensc.xyz[c, :]
        x, y, z = centroid
        fc.values[i, :] .= sin(2 * x / A) * cos(6.5 * y / B) * sin(3 * z / C - 1.0)
    end
    File = "momap2para3-coarse.vtk"
    vtkexportmesh(File, fensc, fesc; scalars = [("fc", fc.values)])
    # @async run(`"paraview.exe" $File`)
    try
        rm(File)
    catch
    end

    Refinement = Refinement + 1
    nA, nB, nC = Refinement * 1, Refinement * 6, Refinement * 4
    xs = collect(linearspace(0.0, A, nA + 1))
    ys = collect(linearspace(0.0, B, nB + 1))
    zs = collect(linearspace(0.0, C, nC + 1))
    fensf, fesf = Meshing(xs, ys, zs)
    tolerance = min(A / nA, B / nB, C / nC) / 1000.0

    ff = ElementalField(zeros(count(fesf), 1))
    referenceff = ElementalField(zeros(count(fesf), 1))
    for i in eachindex(fesf)
        c = [k for k in fesf.conn[i]]
        centroid = NT * fensf.xyz[c, :]
        x, y, z = centroid
        referenceff.values[i, :] .= sin(2 * x / A) * cos(6.5 * y / B) * sin(3 * z / C - 1.0)
    end
    File = "momap2para3-reference.vtk"
    vtkexportmesh(File, fensf, fesf; scalars = [("referenceff", referenceff.values)])
    # @async run(`"paraview.exe" $File`)
    try
        rm(File)
    catch
    end

    ff = transferfield!(ff, fensf, fesf, fc, fensc, fesc, tolerance)
    File = "momap2para3-fine.vtk"
    vtkexportmesh(File, fensf, fesf; scalars = [("ff", ff.values)])
    # @async run(`"paraview.exe" $File`)
    try
        rm(File)
    catch
    end

    diffff = ElementalField(referenceff.values - ff.values)
    femm = FEMMBase(IntegDomain(fesf, GaussRule(3, 3)))
    geom = NodalField(fensf.xyz)
    error = integratefieldfunction(femm, geom, diffff, (x, v) -> norm(v); initial = 0.0)
    ref = integratefieldfunction(femm, geom, referenceff, (x, v) -> norm(v); initial = 0.0)
    # println("error/ref = $(error/ref)")
    @test abs(error / ref - 0.3602708839379691) < 1.0e-4
end
end
using .momap2para6
momap2para6.test()

module momap2para7
using FinEtools
using FinEtools.MeshSelectionModule: vselect
using FinEtools.MeshExportModule
using Test
import LinearAlgebra: norm
function test()
    A = 50.0 * phun("m") # length  of loaded rectangle
    B = 200.0 * phun("m") # length  of loaded rectangle
    C = 100.0 * phun("m") # span of the plate

    Meshing = H8blockx
    # Select how find the mesh should be
    Refinement = 2
    nA, nB, nC = Refinement * 1, Refinement * 6, Refinement * 4
    xs = collect(linearspace(0.0, A, nA + 1))
    ys = collect(linearspace(0.0, B, nB + 1))
    zs = collect(linearspace(0.0, C, nC + 1))
    fensc, fesc = Meshing(xs, ys, zs)

    centroidpc = centroidparametric(fesc)
    N = bfun(fesc, centroidpc)
    NT = transpose(N)

    fc = ElementalField(zeros(count(fesc), 1))
    for i in eachindex(fesc)
        c = [k for k in fesc.conn[i]]
        centroid = NT * fensc.xyz[c, :]
        x, y, z = centroid
        fc.values[i, :] .= sin(2 * x / A) * cos(6.5 * y / B) * sin(3 * z / C - 1.0)
    end
    File = "momap2para3-coarse.vtk"
    vtkexportmesh(File, fensc, fesc; scalars = [("fc", fc.values)])
    # @async run(`"paraview.exe" $File`)
    try
        rm(File)
    catch
    end

    Refinement = Refinement + 1
    nA, nB, nC = Refinement * 1, Refinement * 6, Refinement * 4
    xs = collect(linearspace(0.0, A, nA + 1))
    ys = collect(linearspace(0.0, B, nB + 1))
    zs = collect(linearspace(0.0, C, nC + 1))
    fensf, fesf = Meshing(xs, ys, zs)
    tolerance = min(A / nA, B / nB, C / nC) / 1000.0

    ff = ElementalField(zeros(count(fesf), 1))
    referenceff = ElementalField(zeros(count(fesf), 1))
    for i in eachindex(fesf)
        c = [k for k in fesf.conn[i]]
        centroid = NT * fensf.xyz[c, :]
        x, y, z = centroid
        referenceff.values[i, :] .= sin(2 * x / A) * cos(6.5 * y / B) * sin(3 * z / C - 1.0)
    end
    File = "momap2para3-reference.vtk"
    vtkexportmesh(File, fensf, fesf; scalars = [("referenceff", referenceff.values)])
    # @async run(`"paraview.exe" $File`)
    try
        rm(File)
    catch
    end

    ff = transferfield!(ff, fensf, fesf, fc, fensc, fesc, tolerance)
    File = "momap2para3-fine.vtk"
    vtkexportmesh(File, fensf, fesf; scalars = [("ff", ff.values)])
    # @async run(`"paraview.exe" $File`)
    try
        rm(File)
    catch
    end

    diffff = ElementalField(referenceff.values - ff.values)
    femm = FEMMBase(IntegDomain(fesf, GaussRule(3, 3)))
    geom = NodalField(fensf.xyz)
    error = integratefieldfunction(femm, geom, diffff, (x, v) -> norm(v); initial = 0.0)
    ref = integratefieldfunction(femm, geom, referenceff, (x, v) -> norm(v); initial = 0.0)
    # println("error/ref = $(error/ref)")
    @test abs(error / ref - 0.36027088393796847) < 1.0e-4
end
end
using .momap2para7
momap2para7.test()

module momap2para8
using FinEtools
using FinEtools.MeshSelectionModule: vselect
using FinEtools.MeshExportModule
using Test
import LinearAlgebra: norm
function test()
    A = 50.0 * phun("m") # length  of loaded rectangle
    B = 200.0 * phun("m") # length  of loaded rectangle
    C = 100.0 * phun("m") # span of the plate

    Meshing = H27blockx
    # Select how find the mesh should be
    Refinement = 2
    nA, nB, nC = Refinement * 1, Refinement * 6, Refinement * 4
    xs = collect(linearspace(0.0, A, nA + 1))
    ys = collect(linearspace(0.0, B, nB + 1))
    zs = collect(linearspace(0.0, C, nC + 1))
    fensc, fesc = Meshing(xs, ys, zs)

    centroidpc = centroidparametric(fesc)
    N = bfun(fesc, centroidpc)
    NT = transpose(N)

    fc = ElementalField(zeros(count(fesc), 1))
    for i in eachindex(fesc)
        c = [k for k in fesc.conn[i]]
        centroid = NT * fensc.xyz[c, :]
        x, y, z = centroid
        fc.values[i, :] .= sin(2 * x / A) * cos(6.5 * y / B) * sin(3 * z / C - 1.0)
    end
    # File = "momap2para3-coarse.vtk"
    # VTK.vtkexportmesh(File, fensc, fesc; scalars = [("fc", fc.values)])
    # # @async run(`"paraview.exe" $File`)
    # try rm(File) catch end

    Refinement = Refinement + 1
    nA, nB, nC = Refinement * 1, Refinement * 6, Refinement * 4
    xs = collect(linearspace(0.0, A, nA + 1))
    ys = collect(linearspace(0.0, B, nB + 1))
    zs = collect(linearspace(0.0, C, nC + 1))
    fensf, fesf = Meshing(xs, ys, zs)
    tolerance = min(A / nA, B / nB, C / nC) / 1000.0

    ff = ElementalField(zeros(count(fesf), 1))
    referenceff = ElementalField(zeros(count(fesf), 1))
    for i in eachindex(fesf)
        c = [k for k in fesf.conn[i]]
        centroid = NT * fensf.xyz[c, :]
        x, y, z = centroid
        referenceff.values[i, :] .= sin(2 * x / A) * cos(6.5 * y / B) * sin(3 * z / C - 1.0)
    end
    # File = "momap2para3-reference.vtk"
    # VTK.vtkexportmesh(File, fensf, fesf; scalars = [("referenceff", referenceff.values)])
    # # @async run(`"paraview.exe" $File`)
    # try rm(File) catch end

    ff = transferfield!(ff, fensf, fesf, fc, fensc, fesc, tolerance)
    # File = "momap2para3-fine.vtk"
    # VTK.vtkexportmesh(File, fensf, fesf; scalars = [("ff", ff.values)])
    # # @async run(`"paraview.exe" $File`)
    # try rm(File) catch end

    diffff = ElementalField(referenceff.values - ff.values)
    femm = FEMMBase(IntegDomain(fesf, GaussRule(3, 3)))
    geom = NodalField(fensf.xyz)
    error = integratefieldfunction(femm, geom, diffff, (x, v) -> norm(v); initial = 0.0)
    ref = integratefieldfunction(femm, geom, referenceff, (x, v) -> norm(v); initial = 0.0)
    # println("error/ref = $(error/ref)")
    @test abs(error / ref - 0.3602708839379695) < 1.0e-4
end
end
using .momap2para8
momap2para8.test()


module momap2para9
using FinEtools
using FinEtools.MeshSelectionModule: vselect
using FinEtools.MeshExportModule
using Test
import LinearAlgebra: norm
function test()
    A = 50.0 * phun("m") # length  of loaded rectangle
    B = 200.0 * phun("m") # length  of loaded rectangle
    C = 100.0 * phun("m") # span of the plate

    Meshing = Q4blockx
    # Select how find the mesh should be
    Refinement = 2
    nA, nB, nC = Refinement * 1, Refinement * 6, Refinement * 4
    xs = collect(linearspace(0.0, A, nA + 1))
    ys = collect(linearspace(0.0, B, nB + 1))
    fensc, fesc = Meshing(xs, ys)

    centroidpc = centroidparametric(fesc)
    N = bfun(fesc, centroidpc)
    NT = transpose(N)

    fc = ElementalField(zeros(count(fesc), 1))
    for i in eachindex(fesc)
        c = [k for k in fesc.conn[i]]
        centroid = NT * fensc.xyz[c, :]
        x, y = centroid
        fc.values[i, :] .= sin(2 * x / A) * cos(6.5 * y / B)
    end
    # File = "momap2para3-coarse.vtk"
    # VTK.vtkexportmesh(File, fensc, fesc; scalars = [("fc", fc.values)])
    # # @async run(`"paraview.exe" $File`)
    # try rm(File) catch end

    Refinement = Refinement + 1
    nA, nB, nC = Refinement * 1, Refinement * 6, Refinement * 4
    xs = collect(linearspace(0.0, A, nA + 1))
    ys = collect(linearspace(0.0, B, nB + 1))
    zs = collect(linearspace(0.0, C, nC + 1))
    fensf, fesf = Meshing(xs, ys)
    tolerance = min(A / nA, B / nB, C / nC) / 1000.0

    ff = ElementalField(zeros(count(fesf), 1))
    referenceff = ElementalField(zeros(count(fesf), 1))
    for i in eachindex(fesf)
        c = [k for k in fesf.conn[i]]
        centroid = NT * fensf.xyz[c, :]
        x, y = centroid
        referenceff.values[i, :] .= sin(2 * x / A) * cos(6.5 * y / B)
    end
    # File = "momap2para3-reference.vtk"
    # VTK.vtkexportmesh(File, fensf, fesf; scalars = [("referenceff", referenceff.values)])
    # # @async run(`"paraview.exe" $File`)
    # try rm(File) catch end

    ff = transferfield!(ff, fensf, fesf, fc, fensc, fesc, tolerance)
    # File = "momap2para3-fine.vtk"
    # VTK.vtkexportmesh(File, fensf, fesf; scalars = [("ff", ff.values)])
    # # @async run(`"paraview.exe" $File`)
    # try rm(File) catch end

    diffff = ElementalField(referenceff.values - ff.values)
    femm = FEMMBase(IntegDomain(fesf, GaussRule(2, 3)))
    geom = NodalField(fensf.xyz)
    error = integratefieldfunction(femm, geom, diffff, (x, v) -> norm(v); initial = 0.0)
    ref = integratefieldfunction(femm, geom, referenceff, (x, v) -> norm(v); initial = 0.0)
    # println("error/ref = $(error/ref)")
    @test abs(error / ref - 0.32348898713315494) < 1.0e-4
end
end
using .momap2para9
momap2para9.test()


module momap2para10
using FinEtools
using FinEtools.MeshSelectionModule: vselect
using FinEtools.MeshExportModule
using Test
import LinearAlgebra: norm
function test()
    A = 50.0 * phun("m") # length  of loaded rectangle
    B = 200.0 * phun("m") # length  of loaded rectangle
    C = 100.0 * phun("m") # span of the plate

    Meshing = Q8blockx
    # Select how find the mesh should be
    Refinement = 2
    nA, nB, nC = Refinement * 1, Refinement * 6, Refinement * 4
    xs = collect(linearspace(0.0, A, nA + 1))
    ys = collect(linearspace(0.0, B, nB + 1))
    fensc, fesc = Meshing(xs, ys)

    centroidpc = centroidparametric(fesc)
    N = bfun(fesc, centroidpc)
    NT = transpose(N)

    fc = ElementalField(zeros(count(fesc), 1))
    for i in eachindex(fesc)
        c = [k for k in fesc.conn[i]]
        centroid = NT * fensc.xyz[c, :]
        x, y = centroid
        fc.values[i, :] .= sin(2 * x / A) * cos(6.5 * y / B)
    end
    # File = "momap2para3-coarse.vtk"
    # VTK.vtkexportmesh(File, fensc, fesc; scalars = [("fc", fc.values)])
    # # @async run(`"paraview.exe" $File`)
    # try rm(File) catch end

    Refinement = Refinement + 1
    nA, nB, nC = Refinement * 1, Refinement * 6, Refinement * 4
    xs = collect(linearspace(0.0, A, nA + 1))
    ys = collect(linearspace(0.0, B, nB + 1))
    zs = collect(linearspace(0.0, C, nC + 1))
    fensf, fesf = Meshing(xs, ys)
    tolerance = min(A / nA, B / nB, C / nC) / 1000.0

    ff = ElementalField(zeros(count(fesf), 1))
    referenceff = ElementalField(zeros(count(fesf), 1))
    for i in eachindex(fesf)
        c = [k for k in fesf.conn[i]]
        centroid = NT * fensf.xyz[c, :]
        x, y = centroid
        referenceff.values[i, :] .= sin(2 * x / A) * cos(6.5 * y / B)
    end
    # File = "momap2para3-reference.vtk"
    # VTK.vtkexportmesh(File, fensf, fesf; scalars = [("referenceff", referenceff.values)])
    # # @async run(`"paraview.exe" $File`)
    # try rm(File) catch end

    ff = transferfield!(ff, fensf, fesf, fc, fensc, fesc, tolerance)
    # File = "momap2para3-fine.vtk"
    # VTK.vtkexportmesh(File, fensf, fesf; scalars = [("ff", ff.values)])
    # # @async run(`"paraview.exe" $File`)
    # try rm(File) catch end

    diffff = ElementalField(referenceff.values - ff.values)
    femm = FEMMBase(IntegDomain(fesf, GaussRule(2, 3)))
    geom = NodalField(fensf.xyz)
    error = integratefieldfunction(femm, geom, diffff, (x, v) -> norm(v); initial = 0.0)
    ref = integratefieldfunction(femm, geom, referenceff, (x, v) -> norm(v); initial = 0.0)
    # println("error/ref = $(error/ref)")
    @test abs(error / ref - 0.32348898713315494) < 1.0e-4
end
end
using .momap2para10
momap2para10.test()


module momap2para11
using FinEtools
using FinEtools.MeshSelectionModule: vselect
using FinEtools.MeshExportModule
using Test
import LinearAlgebra: norm
function test()
    A = 50.0 * phun("m") # length  of loaded rectangle
    B = 200.0 * phun("m") # length  of loaded rectangle
    C = 100.0 * phun("m") # span of the plate

    Meshing = T3blockx
    # Select how find the mesh should be
    Refinement = 2
    nA, nB, nC = Refinement * 1, Refinement * 6, Refinement * 4
    xs = collect(linearspace(0.0, A, nA + 1))
    ys = collect(linearspace(0.0, B, nB + 1))
    fensc, fesc = Meshing(xs, ys)

    centroidpc = centroidparametric(fesc)
    N = bfun(fesc, centroidpc)
    NT = transpose(N)

    fc = ElementalField(zeros(count(fesc), 1))
    for i in eachindex(fesc)
        c = [k for k in fesc.conn[i]]
        centroid = NT * fensc.xyz[c, :]
        x, y = centroid
        fc.values[i, :] .= sin(2 * x / A) * cos(6.5 * y / B)
    end
    File = "momap2para3-coarse.vtk"
    vtkexportmesh(File, fensc, fesc; scalars = [("fc", fc.values)])
    # @async run(`"paraview.exe" $File`)
    try
        rm(File)
    catch
    end

    Refinement = Refinement + 1
    nA, nB, nC = Refinement * 1, Refinement * 6, Refinement * 4
    xs = collect(linearspace(0.0, A, nA + 1))
    ys = collect(linearspace(0.0, B, nB + 1))
    zs = collect(linearspace(0.0, C, nC + 1))
    fensf, fesf = Meshing(xs, ys)
    tolerance = min(A / nA, B / nB, C / nC) / 1000.0

    ff = ElementalField(zeros(count(fesf), 1))
    referenceff = ElementalField(zeros(count(fesf), 1))
    for i in eachindex(fesf)
        c = [k for k in fesf.conn[i]]
        centroid = NT * fensf.xyz[c, :]
        x, y = centroid
        referenceff.values[i, :] .= sin(2 * x / A) * cos(6.5 * y / B)
    end
    File = "momap2para3-reference.vtk"
    vtkexportmesh(File, fensf, fesf; scalars = [("referenceff", referenceff.values)])
    # @async run(`"paraview.exe" $File`)
    try
        rm(File)
    catch
    end

    ff = transferfield!(ff, fensf, fesf, fc, fensc, fesc, tolerance)
    File = "momap2para3-fine.vtk"
    vtkexportmesh(File, fensf, fesf; scalars = [("ff", ff.values)])
    # @async run(`"paraview.exe" $File`)
    try
        rm(File)
    catch
    end

    diffff = ElementalField(referenceff.values - ff.values)
    femm = FEMMBase(IntegDomain(fesf, SimplexRule(2, 3)))
    geom = NodalField(fensf.xyz)
    error = integratefieldfunction(femm, geom, diffff, (x, v) -> norm(v); initial = 0.0)
    ref = integratefieldfunction(femm, geom, referenceff, (x, v) -> norm(v); initial = 0.0)
    # println("error/ref = $(error/ref)")
    @test abs(error / ref - 0.19654730310399787) < 1.0e-4
end
end
using .momap2para11
momap2para11.test()


module momap2para12
using FinEtools
using FinEtools.MeshSelectionModule: vselect
using FinEtools.MeshExportModule
using Test
import LinearAlgebra: norm
function test()
    A = 50.0 * phun("m") # length  of loaded rectangle
    B = 200.0 * phun("m") # length  of loaded rectangle
    C = 100.0 * phun("m") # span of the plate

    Meshing = T6blockx
    # Select how find the mesh should be
    Refinement = 2
    nA, nB, nC = Refinement * 1, Refinement * 6, Refinement * 4
    xs = collect(linearspace(0.0, A, nA + 1))
    ys = collect(linearspace(0.0, B, nB + 1))
    fensc, fesc = Meshing(xs, ys)

    centroidpc = centroidparametric(fesc)
    N = bfun(fesc, centroidpc)
    NT = transpose(N)

    fc = ElementalField(zeros(count(fesc), 1))
    for i in eachindex(fesc)
        c = [k for k in fesc.conn[i]]
        centroid = NT * fensc.xyz[c, :]
        x, y = centroid
        fc.values[i, :] .= sin(2 * x / A) * cos(6.5 * y / B)
    end
    File = "momap2para12-coarse.vtk"
    vtkexportmesh(File, fensc, fesc; scalars = [("fc", fc.values)])
    # @async run(`"paraview.exe" $File`)
    try
        rm(File)
    catch
    end

    Refinement = Refinement + 1
    nA, nB, nC = Refinement * 1, Refinement * 6, Refinement * 4
    xs = collect(linearspace(0.0, A, nA + 1))
    ys = collect(linearspace(0.0, B, nB + 1))
    zs = collect(linearspace(0.0, C, nC + 1))
    fensf, fesf = Meshing(xs, ys)
    tolerance = min(A / nA, B / nB, C / nC) / 1000.0

    ff = ElementalField(zeros(count(fesf), 1))
    referenceff = ElementalField(zeros(count(fesf), 1))
    for i in eachindex(fesf)
        c = [k for k in fesf.conn[i]]
        centroid = NT * fensf.xyz[c, :]
        x, y = centroid
        referenceff.values[i, :] .= sin(2 * x / A) * cos(6.5 * y / B)
    end
    File = "momap2para12-reference.vtk"
    vtkexportmesh(File, fensf, fesf; scalars = [("referenceff", referenceff.values)])
    # @async run(`"paraview.exe" $File`)
    try
        rm(File)
    catch
    end

    ff = transferfield!(ff, fensf, fesf, fc, fensc, fesc, tolerance)
    File = "momap2para12-fine.vtk"
    vtkexportmesh(File, fensf, fesf; scalars = [("ff", ff.values)])
    # @async run(`"paraview.exe" $File`)
    try
        rm(File)
    catch
    end

    diffff = ElementalField(referenceff.values - ff.values)
    femm = FEMMBase(IntegDomain(fesf, SimplexRule(2, 3)))
    geom = NodalField(fensf.xyz)
    error = integratefieldfunction(femm, geom, diffff, (x, v) -> norm(v); initial = 0.0)
    ref = integratefieldfunction(femm, geom, referenceff, (x, v) -> norm(v); initial = 0.0)
    # println("error/ref = $(error/ref)")
    @test abs(error / ref - 0.19654730310399787) < 1.0e-4
end
end
using .momap2para12
momap2para12.test()


module momap2para13
using FinEtools
using FinEtools.MeshSelectionModule: vselect
using FinEtools.MeshExportModule
using Test
import LinearAlgebra: norm, cross
function test()
    A = 50.0 * phun("m") # length  of loaded rectangle
    B = 200.0 * phun("m") # length  of loaded rectangle
    C = 100.0 * phun("m") # span of the plate

    Meshing = L2blockx
    # Select how find the mesh should be
    Refinement = 2
    nA, nB, nC = Refinement * 1, Refinement * 6, Refinement * 4
    xs = collect(linearspace(0.0, A, nA + 1))
    fensc, fesc = Meshing(xs)

    centroidpc = centroidparametric(fesc)
    N = bfun(fesc, centroidpc)
    NT = transpose(N)

    fc = ElementalField(zeros(count(fesc), 1))
    for i in eachindex(fesc)
        c = [k for k in fesc.conn[i]]
        centroid = NT * fensc.xyz[c, :]
        x = centroid
        fc.values[i, :] = sin.(2 * x / A)
    end
    File = "momap2para12-coarse.vtk"
    vtkexportmesh(File, fensc, fesc; scalars = [("fc", fc.values)])
    # @async run(`"paraview.exe" $File`)
    try
        rm(File)
    catch
    end

    Refinement = Refinement + 1
    nA, nB, nC = Refinement * 1, Refinement * 6, Refinement * 4
    xs = collect(linearspace(0.0, A, nA + 1))
    fensf, fesf = Meshing(xs)
    tolerance = min(A / nA, B / nB, C / nC) / 1000.0

    ff = ElementalField(zeros(count(fesf), 1))
    referenceff = ElementalField(zeros(count(fesf), 1))
    for i in eachindex(fesf)
        c = [k for k in fesf.conn[i]]
        centroid = NT * fensf.xyz[c, :]
        x = centroid
        referenceff.values[i, :] = sin.(2 * x / A)
    end
    File = "momap2para12-reference.vtk"
    vtkexportmesh(File, fensf, fesf; scalars = [("referenceff", referenceff.values)])
    # @async run(`"paraview.exe" $File`)
    try
        rm(File)
    catch
    end

    ff = transferfield!(ff, fensf, fesf, fc, fensc, fesc, tolerance)
    File = "momap2para12-fine.vtk"
    vtkexportmesh(File, fensf, fesf; scalars = [("ff", ff.values)])
    # @async run(`"paraview.exe" $File`)
    try
        rm(File)
    catch
    end

    diffff = ElementalField(referenceff.values - ff.values)
    femm = FEMMBase(IntegDomain(fesf, GaussRule(1, 3)))
    geom = NodalField(fensf.xyz)
    error = integratefieldfunction(femm, geom, diffff, (x, v) -> norm(v); initial = 0.0)
    ref = integratefieldfunction(femm, geom, referenceff, (x, v) -> norm(v); initial = 0.0)
    # println("error/ref = $(error/ref)")
    @test abs(error / ref - 0.23860709149331033) < 1.0e-4
end
end
using .momap2para13
momap2para13.test()

module momap2para14
using FinEtools
using FinEtools.MeshSelectionModule: vselect
using FinEtools.MeshExportModule
using Test
import LinearAlgebra: norm, cross
function test()
    A = 50.0 * phun("m") # length  of loaded rectangle
    B = 200.0 * phun("m") # length  of loaded rectangle
    C = 100.0 * phun("m") # span of the plate

    Meshing = L3blockx
    # Select how find the mesh should be
    Refinement = 2
    nA, nB, nC = Refinement * 1, Refinement * 6, Refinement * 4
    xs = collect(linearspace(0.0, A, nA + 1))
    fensc, fesc = Meshing(xs)
    #     println("fensc = $(fensc)")
    # println("fesc = $(fesc)")

    centroidpc = centroidparametric(fesc)
    N = bfun(fesc, centroidpc)
    NT = transpose(N)

    fc = ElementalField(zeros(count(fesc), 1))
    for i in eachindex(fesc)
        c = [k for k in fesc.conn[i]]
        centroid = NT * fensc.xyz[c, :]
        x = centroid
        fc.values[i, :] = sin.(2 * x / A)
    end
    File = "momap2para12-coarse.vtk"
    vtkexportmesh(File, fensc, fesc; scalars = [("fc", fc.values)])
    # @async run(`"paraview.exe" $File`)
    try
        rm(File)
    catch
    end

    Refinement = Refinement + 1
    nA, nB, nC = Refinement * 1, Refinement * 6, Refinement * 4
    xs = collect(linearspace(0.0, A, nA + 1))
    fensf, fesf = Meshing(xs)
    tolerance = min(A / nA, B / nB, C / nC) / 1000.0

    ff = ElementalField(zeros(count(fesf), 1))
    referenceff = ElementalField(zeros(count(fesf), 1))
    for i in eachindex(fesf)
        c = [k for k in fesf.conn[i]]
        centroid = NT * fensf.xyz[c, :]
        x = centroid
        referenceff.values[i, :] = sin.(2 * x / A)
    end
    File = "momap2para12-reference.vtk"
    vtkexportmesh(File, fensf, fesf; scalars = [("referenceff", referenceff.values)])
    # @async run(`"paraview.exe" $File`)
    try
        rm(File)
    catch
    end

    ff = transferfield!(ff, fensf, fesf, fc, fensc, fesc, tolerance)
    File = "momap2para12-fine.vtk"
    vtkexportmesh(File, fensf, fesf; scalars = [("ff", ff.values)])
    # @async run(`"paraview.exe" $File`)
    try
        rm(File)
    catch
    end

    diffff = ElementalField(referenceff.values - ff.values)
    femm = FEMMBase(IntegDomain(fesf, GaussRule(1, 3)))
    geom = NodalField(fensf.xyz)
    error = integratefieldfunction(femm, geom, diffff, (x, v) -> norm(v); initial = 0.0)
    ref = integratefieldfunction(femm, geom, referenceff, (x, v) -> norm(v); initial = 0.0)
    # println("error/ref = $(error/ref)")
    @test abs(error / ref - 0.23860709149331033) < 1.0e-4
end
end
using .momap2para14
momap2para14.test()

module mpointmm1
using FinEtools
using Test
function test()
    x = [i == j ? one(Float64) : zero(Float64) for i in 1:4, j in 1:3]
    fes = FESetP1(reshape([1 2 4 3], 4, 1))
    pt = x[4, :]
    pc, success = map2parametric(fes, reshape(x[4, :], 1, 3), pt)
    @test success
    @test pc[1] == 00.0
    pc, success = map2parametric(fes, reshape(x[2, :], 1, 3), pt)
    @test !success

    c = centroidparametric(fes)
    @test c[1] == 0.0

    b = inparametric(fes, [0.0], tolerance = 0.001)
    @test b
end
end
using .mpointmm1
mpointmm1.test()

module ml314
using FinEtools
using FinEtools.MeshSelectionModule: vselect
using FinEtools.MeshExportModule
using Test
function test()
    A = 50.0 * phun("m") # length  of loaded rectangle
    B = 200.0 * phun("m") # length  of loaded rectangle
    C = 100.0 * phun("m") # span of the plate

    Meshing = L3blockx
    # Select how find the mesh should be
    Refinement = 2
    nA, nB, nC = Refinement * 1, Refinement * 6, Refinement * 4
    xs = collect(linearspace(0.0, A, nA + 1))
    fensc, fesc = Meshing(xs)
    bfes = meshboundary(fesc)
    @test count(bfes) == 2


end
end
using .ml314
ml314.test()

module mtetmeshedges1
using FinEtools
using FinEtools.MeshTetrahedronModule
using Test
import LinearAlgebra: norm, cross
t = [
    1 17 13 14
    6 5 2 18
    18 2 14 17
    5 18 17 2
    2 1 14 17
    5 17 1 2
    13 29 25 26
    18 17 14 30
    30 14 26 29
    17 30 29 14
    14 13 26 29
    17 29 13 14
    25 41 37 38
    30 29 26 42
    42 26 38 41
    29 42 41 26
    26 25 38 41
    29 41 25 26
    37 53 49 50
    42 41 38 54
    54 38 50 53
    41 54 53 38
    38 37 50 53
    41 53 37 38
    5 21 17 18
    10 9 6 22
    22 6 18 21
    9 22 21 6
    6 5 18 21
    9 21 5 6
    17 33 29 30
    22 21 18 34
    34 18 30 33
    21 34 33 18
    18 17 30 33
    21 33 17 18
    29 45 41 42
    34 33 30 46
    46 30 42 45
    33 46 45 30
    30 29 42 45
    33 45 29 30
    41 57 53 54
    46 45 42 58
    58 42 54 57
    45 58 57 42
    42 41 54 57
    45 57 41 42
    2 18 14 15
    7 6 3 19
    19 3 15 18
    6 19 18 3
    3 2 15 18
    6 18 2 3
    14 30 26 27
    19 18 15 31
    31 15 27 30
    18 31 30 15
    15 14 27 30
    18 30 14 15
    26 42 38 39
    31 30 27 43
    43 27 39 42
    30 43 42 27
    27 26 39 42
    30 42 26 27
    38 54 50 51
    43 42 39 55
    55 39 51 54
    42 55 54 39
    39 38 51 54
    42 54 38 39
    6 22 18 19
    11 10 7 23
    23 7 19 22
    10 23 22 7
    7 6 19 22
    10 22 6 7
    18 34 30 31
    23 22 19 35
    35 19 31 34
    22 35 34 19
    19 18 31 34
    22 34 18 19
    30 46 42 43
    35 34 31 47
    47 31 43 46
    34 47 46 31
    31 30 43 46
    34 46 30 31
    42 58 54 55
    47 46 43 59
    59 43 55 58
    46 59 58 43
    43 42 55 58
    46 58 42 43
    3 19 15 16
    8 7 4 20
    20 4 16 19
    7 20 19 4
    4 3 16 19
    7 19 3 4
    15 31 27 28
    20 19 16 32
    32 16 28 31
    19 32 31 16
    16 15 28 31
    19 31 15 16
    27 43 39 40
    32 31 28 44
    44 28 40 43
    31 44 43 28
    28 27 40 43
    31 43 27 28
    39 55 51 52
    44 43 40 56
    56 40 52 55
    43 56 55 40
    40 39 52 55
    43 55 39 40
    7 23 19 20
    12 11 8 24
    24 8 20 23
    11 24 23 8
    8 7 20 23
    11 23 7 8
    19 35 31 32
    24 23 20 36
    36 20 32 35
    23 36 35 20
    20 19 32 35
    23 35 19 20
    31 47 43 44
    36 35 32 48
    48 32 44 47
    35 48 47 32
    32 31 44 47
    35 47 31 32
    43 59 55 56
    48 47 44 60
    60 44 56 59
    47 60 59 44
    44 43 56 59
    47 59 43 44
]

e = [
    1 2
    1 5
    1 13
    1 14
    1 17
    2 3
    2 5
    2 6
    2 14
    2 15
    2 17
    2 18
    3 4
    3 6
    3 7
    3 15
    3 16
    3 18
    3 19
    4 7
    4 8
    4 16
    4 19
    4 20
    5 6
    5 9
    5 17
    5 18
    5 21
    6 7
    6 9
    6 10
    6 18
    6 19
    6 21
    6 22
    7 8
    7 10
    7 11
    7 19
    7 20
    7 22
    7 23
    8 11
    8 12
    8 20
    8 23
    8 24
    9 10
    9 21
    9 22
    10 11
    10 22
    10 23
    11 12
    11 23
    11 24
    12 24
    13 14
    13 17
    13 25
    13 26
    13 29
    14 15
    14 17
    14 18
    14 26
    14 27
    14 29
    14 30
    15 16
    15 18
    15 19
    15 27
    15 28
    15 30
    15 31
    16 19
    16 20
    16 28
    16 31
    16 32
    17 18
    17 21
    17 29
    17 30
    17 33
    18 19
    18 21
    18 22
    18 30
    18 31
    18 33
    18 34
    19 20
    19 22
    19 23
    19 31
    19 32
    19 34
    19 35
    20 23
    20 24
    20 32
    20 35
    20 36
    21 22
    21 33
    21 34
    22 23
    22 34
    22 35
    23 24
    23 35
    23 36
    24 36
    25 26
    25 29
    25 37
    25 38
    25 41
    26 27
    26 29
    26 30
    26 38
    26 39
    26 41
    26 42
    27 28
    27 30
    27 31
    27 39
    27 40
    27 42
    27 43
    28 31
    28 32
    28 40
    28 43
    28 44
    29 30
    29 33
    29 41
    29 42
    29 45
    30 31
    30 33
    30 34
    30 42
    30 43
    30 45
    30 46
    31 32
    31 34
    31 35
    31 43
    31 44
    31 46
    31 47
    32 35
    32 36
    32 44
    32 47
    32 48
    33 34
    33 45
    33 46
    34 35
    34 46
    34 47
    35 36
    35 47
    35 48
    36 48
    37 38
    37 41
    37 49
    37 50
    37 53
    38 39
    38 41
    38 42
    38 50
    38 51
    38 53
    38 54
    39 40
    39 42
    39 43
    39 51
    39 52
    39 54
    39 55
    40 43
    40 44
    40 52
    40 55
    40 56
    41 42
    41 45
    41 53
    41 54
    41 57
    42 43
    42 45
    42 46
    42 54
    42 55
    42 57
    42 58
    43 44
    43 46
    43 47
    43 55
    43 56
    43 58
    43 59
    44 47
    44 48
    44 56
    44 59
    44 60
    45 46
    45 57
    45 58
    46 47
    46 58
    46 59
    47 48
    47 59
    47 60
    48 60
    49 50
    49 53
    50 51
    50 53
    50 54
    51 52
    51 54
    51 55
    52 55
    52 56
    53 54
    53 57
    54 55
    54 57
    54 58
    55 56
    55 58
    55 59
    56 59
    56 60
    57 58
    58 59
    59 60
]

function tetmeshedges(t::Array{Int,2})
    ec = [
        1 2
        2 3
        3 1
        4 1
        4 2
        4 3
    ]
    e = vcat(
        t[:, ec[1, :]],
        t[:, ec[2, :]],
        t[:, ec[3, :]],
        t[:, ec[4, :]],
        t[:, ec[5, :]],
        t[:, ec[6, :]],
    )
    e = sort(e, dims = 2)
    ix = sortperm(e[:, 1])
    e = e[ix, :]
    ue = deepcopy(e)
    i = 1
    n = 1
    while n <= size(e, 1)
        c = ue[n, 1]
        m = n + 1
        while m <= size(e, 1)
            if (ue[m, 1] != c)
                break
            end
            m = m + 1
        end
        us = unique(ue[n:m-1, 2], dims = 1)
        ls = length(us)
        e[i:i+ls-1, 1] .= c
        e[i:i+ls-1, 2] = sort(us)
        i = i + ls
        n = m
    end
    e = e[1:i-1, :]
end

function test()
    mye = tetmeshedges(t)
    # println("mye = $(mye)")
    @test norm(e - mye) == 0
    mye = MeshTetrahedronModule.T4meshedges(t)
    # println("mye = $(mye)")
    @test norm(e - mye) == 0
end
end
using .mtetmeshedges1
mtetmeshedges1.test()

module minterior2boundary1
using FinEtools
using Test
import LinearAlgebra: norm, cross
function test()
    t = [
        1 17 13 14
        6 5 2 18
        18 2 14 17
        5 18 17 2
        2 1 14 17
        5 17 1 2
        13 29 25 26
        18 17 14 30
        30 14 26 29
        17 30 29 14
        14 13 26 29
        17 29 13 14
        25 41 37 38
        30 29 26 42
        42 26 38 41
        29 42 41 26
        26 25 38 41
        29 41 25 26
        37 53 49 50
        42 41 38 54
        54 38 50 53
        41 54 53 38
        38 37 50 53
        41 53 37 38
        5 21 17 18
        10 9 6 22
        22 6 18 21
        9 22 21 6
        6 5 18 21
        9 21 5 6
        17 33 29 30
        22 21 18 34
        34 18 30 33
        21 34 33 18
        18 17 30 33
        21 33 17 18
        29 45 41 42
        34 33 30 46
        46 30 42 45
        33 46 45 30
        30 29 42 45
        33 45 29 30
        41 57 53 54
        46 45 42 58
        58 42 54 57
        45 58 57 42
        42 41 54 57
        45 57 41 42
        2 18 14 15
        7 6 3 19
        19 3 15 18
        6 19 18 3
        3 2 15 18
        6 18 2 3
        14 30 26 27
        19 18 15 31
        31 15 27 30
        18 31 30 15
        15 14 27 30
        18 30 14 15
        26 42 38 39
        31 30 27 43
        43 27 39 42
        30 43 42 27
        27 26 39 42
        30 42 26 27
        38 54 50 51
        43 42 39 55
        55 39 51 54
        42 55 54 39
        39 38 51 54
        42 54 38 39
        6 22 18 19
        11 10 7 23
        23 7 19 22
        10 23 22 7
        7 6 19 22
        10 22 6 7
        18 34 30 31
        23 22 19 35
        35 19 31 34
        22 35 34 19
        19 18 31 34
        22 34 18 19
        30 46 42 43
        35 34 31 47
        47 31 43 46
        34 47 46 31
        31 30 43 46
        34 46 30 31
        42 58 54 55
        47 46 43 59
        59 43 55 58
        46 59 58 43
        43 42 55 58
        46 58 42 43
        3 19 15 16
        8 7 4 20
        20 4 16 19
        7 20 19 4
        4 3 16 19
        7 19 3 4
        15 31 27 28
        20 19 16 32
        32 16 28 31
        19 32 31 16
        16 15 28 31
        19 31 15 16
        27 43 39 40
        32 31 28 44
        44 28 40 43
        31 44 43 28
        28 27 40 43
        31 43 27 28
        39 55 51 52
        44 43 40 56
        56 40 52 55
        43 56 55 40
        40 39 52 55
        43 55 39 40
        7 23 19 20
        12 11 8 24
        24 8 20 23
        11 24 23 8
        8 7 20 23
        11 23 7 8
        19 35 31 32
        24 23 20 36
        36 20 32 35
        23 36 35 20
        20 19 32 35
        23 35 19 20
        31 47 43 44
        36 35 32 48
        48 32 44 47
        35 48 47 32
        32 31 44 47
        35 47 31 32
        43 59 55 56
        48 47 44 60
        60 44 56 59
        47 60 59 44
        44 43 56 59
        47 59 43 44
    ]
    bt = [
        5 2 1
        2 14 1
        5 1 17
        1 14 13
        1 13 17
        6 3 2
        3 15 2
        6 2 5
        2 15 14
        7 4 3
        4 16 3
        7 3 6
        3 16 15
        8 4 7
        8 20 4
        20 16 4
        9 6 5
        9 5 21
        5 17 21
        10 7 6
        10 6 9
        11 8 7
        11 7 10
        12 8 11
        12 24 8
        24 20 8
        10 9 22
        9 21 22
        11 10 23
        10 22 23
        12 11 24
        11 23 24
        14 26 13
        17 13 29
        13 26 25
        13 25 29
        15 27 14
        14 27 26
        16 28 15
        15 28 27
        20 32 16
        32 28 16
        21 17 33
        17 29 33
        24 36 20
        36 32 20
        22 21 34
        21 33 34
        23 22 35
        22 34 35
        24 23 36
        23 35 36
        26 38 25
        29 25 41
        25 38 37
        25 37 41
        27 39 26
        26 39 38
        28 40 27
        27 40 39
        32 44 28
        44 40 28
        33 29 45
        29 41 45
        36 48 32
        48 44 32
        34 33 46
        33 45 46
        35 34 47
        34 46 47
        36 35 48
        35 47 48
        38 50 37
        41 37 53
        37 50 49
        37 49 53
        39 51 38
        38 51 50
        40 52 39
        39 52 51
        44 56 40
        56 52 40
        45 41 57
        41 53 57
        48 60 44
        60 56 44
        46 45 58
        45 57 58
        47 46 59
        46 58 59
        48 47 60
        47 59 60
        53 49 50
        54 50 51
        54 53 50
        55 51 52
        55 54 51
        56 55 52
        57 53 54
        58 54 55
        58 57 54
        59 55 56
        59 58 55
        60 59 56
    ]
    nbt = interior2boundary(t, [1 3 2; 1 2 4; 2 3 4; 1 4 3])
    # println("nbt = $(nbt)")
    @test norm(bt - nbt) == 0
end
end
using .minterior2boundary1
minterior2boundary1.test()


module miscellan3m
using FinEtools
using Test
function test()
    Lx = 1900.0# length of the box, millimeters
    Ly = 800.0 # length of the box, millimeters

    fens, fes = Q4block(Lx, Ly, 3, 2) # Mesh
    # show(fes.conn)

    fen2fe1 = FENodeToFEMap(connasarray(fes), count(fens))
    # display(connasarray(fes))
    # display(fen2fe1)
    fen2fe2 = FENodeToFEMap(fes.conn, count(fens))
    # display(fes.conn)
    # display(fen2fe2)
    identical = true
    for i in eachindex(fen2fe1.map)
        identical = identical && (fen2fe1.map[i] == fen2fe2.map[i])
        # if !(fen2fe1.map[i] == fen2fe2.map[i])
        #     display(fen2fe1.map[i])
        #     display(fen2fe2.map[i])
        # end
    end
    @test identical
end
end
using .miscellan3m
miscellan3m.test()

module mt4refine1
using FinEtools
using FinEtools.MeshExportModule
using Test
function test()
    xs = collect(linearspace(0.0, pi / 2, 5))
    ys = collect(linearspace(0.0, 1.0, 6) .^ 2)
    zs = collect(linearspace(0.0, 1.0, 3))
    fens, fes = T4blockx(xs, ys, zs, :a)
    # for i = eachindex(fens)
    #     a, y, z = fens.xyz[i,:]
    #     fens.xyz[i,1] = sin(a) * (y + 0.5)
    #     fens.xyz[i,2] = cos(a) * (y + 0.5)
    #     fens.xyz[i,3] = z
    # end
    @test count(fes) == 240
    bfes = meshboundary(fes)
    @test count(bfes) == 2 * 2 * (4 * 5 + 5 * 2 + 4 * 2)
    fens, fes = T4refine(fens, fes)
    for i in eachindex(fens)
        a, y, z = fens.xyz[i, :]
        fens.xyz[i, 1] = sin(a) * (y + 0.5)
        fens.xyz[i, 2] = cos(a) * (y + 0.5)
        fens.xyz[i, 3] = z
    end
    @test count(fes) == 240 * 8
    bfes = meshboundary(fes)
    @test count(bfes) == 4 * 2 * 2 * (4 * 5 + 5 * 2 + 4 * 2)

    geom = NodalField(fens.xyz)
    femm = FEMMBase(IntegDomain(fes, SimplexRule(3, 4)))
    V = integratefunction(femm, geom, (x) -> 1.0)

    File = "Refine-T4-a.vtk"
    vtkexportmesh(File, fens, bfes)
    rm(File)
    # @async run(`"paraview.exe" $File`)
end
end
using .mt4refine1
mt4refine1.test()


module mt4refine2
using FinEtools
using FinEtools.MeshExportModule
using Test
function test()
    xs = collect(linearspace(0.0, pi / 2, 5))
    ys = collect(linearspace(0.0, 1.0, 6) .^ 2)
    zs = collect(linearspace(0.0, 1.0, 3))
    fens, fes = T4blockx(xs, ys, zs, :a)
    fens, fes = T4toT10(fens, fes)
    # for i = eachindex(fens)
    #     a, y, z = fens.xyz[i,:]
    #     fens.xyz[i,1] = sin(a) * (y + 0.5)
    #     fens.xyz[i,2] = cos(a) * (y + 0.5)
    #     fens.xyz[i,3] = z
    # end
    @test count(fes) == 240
    bfes = meshboundary(fes)
    @test count(bfes) == 2 * 2 * (4 * 5 + 5 * 2 + 4 * 2)
    fens, fes = T10refine(fens, fes)
    for i in eachindex(fens)
        a, y, z = fens.xyz[i, :]
        fens.xyz[i, 1] = sin(a) * (y + 0.5)
        fens.xyz[i, 2] = cos(a) * (y + 0.5)
        fens.xyz[i, 3] = z
    end
    @test count(fes) == 240 * 8
    bfes = meshboundary(fes)
    @test count(bfes) == 4 * 2 * 2 * (4 * 5 + 5 * 2 + 4 * 2)

    geom = NodalField(fens.xyz)
    femm = FEMMBase(IntegDomain(fes, SimplexRule(3, 4)))
    V = integratefunction(femm, geom, (x) -> 1.0)

    File = "Refine-T10-a.vtk"
    vtkexportmesh(File, fens, bfes)
    rm(File)
    # @async run(`"paraview.exe" $File`)
end
end
using .mt4refine2
mt4refine2.test()


module mimportexportm1
using FinEtools
using FinEtools.MeshImportModule
using FinEtools.MeshExportModule
using Test
function test()
    output = MeshImportModule.import_NASTRAN(
        dirname(@__FILE__) * "/" * "cylinder.nas";
        allocationchunk = 13,
        expectfixedformat = true,
    )
    @test count(output["fens"]) == 1406
    @test count(output["fesets"][1]) == 829
    # show(fes.conn[count(fes), :])
    File = "cylinder.vtk"
    vtkexportmesh(File, output["fens"], output["fesets"][1])
    rm(File)
    #   @test output["fesets"][1].conn[count(output["fesets"][1]), :] == NTuple{10,Int64}[(143, 140, 144, 138, 361, 363, 176, 519, 781, 520)]
    # @async run(`"paraview.exe" $File`)
end
end
using .mimportexportm1
mimportexportm1.test()




module momap2para61
using FinEtools
using FinEtools.MeshSelectionModule: vselect
using FinEtools.MeshExportModule
using Test
import LinearAlgebra: norm
function test()
    A = 50.0 * phun("m") # length  of loaded rectangle
    B = 200.0 * phun("m") # length  of loaded rectangle
    C = 100.0 * phun("m") # span of the plate

    # Select how find the mesh should be
    Refinement = 2
    nA, nB, nC = Refinement * 1, Refinement * 2, Refinement * 4
    xs = reshape(collect(linearspace(0.0, A, nA + 1)), nA + 1, 1)
    ys = reshape(collect(linearspace(0.0, B, nB + 1)), nB + 1, 1)
    zs = reshape(collect(linearspace(0.0, C, nC + 1)), nC + 1, 1)
    fensc, fesc = T10blockx(xs, ys, zs, :b)
    fc = NodalField(zeros(count(fensc), 1))
    for i in eachindex(fensc)
        x, y, z = fensc.xyz[i, :]
        fc.values[i, :] .= sin(2 * x / A) * cos(6.5 * y / B) * sin(3 * z / C - 1.0)
    end
    File = "momap2para61-coarse.vtk"
    vtkexportmesh(File, fensc, fesc; scalars = [("fc", fc.values)])
    # @async run(`"paraview.exe" $File`)
    try
        rm(File)
    catch
    end

    Refinement = Refinement + 1
    nA, nB, nC = Refinement * 1, Refinement * 2, Refinement * 4
    xs = reshape(collect(linearspace(0.0, A, nA + 1)), nA + 1, 1)
    ys = reshape(collect(linearspace(0.0, B, nB + 1)), nB + 1, 1)
    zs = reshape(collect(linearspace(0.0, C, nC + 1)), nC + 1, 1)
    fensf, fesf = T10blockx(xs, ys, zs, :b)
    ff = NodalField(zeros(count(fensf), 1))
    tolerance = min(A / nA, B / nB, C / nC) / 1000.0

    referenceff = NodalField(zeros(count(fensf), 1))
    for i in eachindex(fensf)
        x, y, z = fensf.xyz[i, :]
        referenceff.values[i, :] .= sin(2 * x / A) * cos(6.5 * y / B) * sin(3 * z / C - 1.0)
    end
    File = "momap2para61-reference.vtk"
    vtkexportmesh(File, fensf, fesf; scalars = [("referenceff", referenceff.values)])
    # @async run(`"paraview.exe" $File`)
    try
        rm(File)
    catch
    end

    ff = transferfield!(ff, fensf, fesf, fc, fensc, fesc, tolerance)
    File = "momap2para61-fine.vtk"
    vtkexportmesh(File, fensf, fesf; scalars = [("ff", ff.values), ("ffcopy", ff.values)])
    # @async run(`"paraview.exe" $File`)
    try
        rm(File)
    catch
    end

    diffff = NodalField(referenceff.values - ff.values)
    femm = FEMMBase(IntegDomain(fesf, SimplexRule(3, 4)))
    geom = NodalField(fensf.xyz)
    error = integratefieldfunction(femm, geom, diffff, (x, v) -> norm(v); initial = 0.0)
    ref = integratefieldfunction(femm, geom, referenceff, (x, v) -> norm(v); initial = 0.0)
    # println("error/ref = $(error/ref)")
    @test abs(error / ref - 0.02541369940759616) < 1.0e-4
end
end
using .momap2para61
momap2para61.test()

module momap2para6378
using FinEtools
using FinEtools.MeshSelectionModule: vselect
using FinEtools.MeshExportModule
using Test
import LinearAlgebra: norm
function test()
    A = 50.0 * phun("m") # length  of loaded rectangle
    B = 200.0 * phun("m") # length  of loaded rectangle
    C = 100.0 * phun("m") # span of the plate

    # Select how find the mesh should be
    Refinement = 2
    nA, nB, nC = Refinement * 1, Refinement * 6, Refinement * 4
    xs = reshape(collect(linearspace(0.0, A, nA + 1)), nA + 1, 1)
    ys = reshape(collect(linearspace(0.0, B, nB + 1)), nB + 1, 1)
    zs = reshape(collect(linearspace(0.0, C, nC + 1)), nC + 1, 1)
    fensc, fesc = T10blockx(xs, ys, zs, :b)

    centroidpc = centroidparametric(fesc)
    N = bfun(fesc, centroidpc)
    NT = transpose(N)

    fc = ElementalField(zeros(count(fesc), 1))
    for i in eachindex(fesc)
        c = [k for k in fesc.conn[i]]
        centroid = NT * fensc.xyz[c, :]
        x, y, z = centroid
        fc.values[i, :] .= sin(2 * x / A) * cos(6.5 * y / B) * sin(3 * z / C - 1.0)
    end
    File = "momap2para3-coarse.vtk"
    vtkexportmesh(File, fensc, fesc; scalars = [("fc", fc.values)])
    # @async run(`"paraview.exe" $File`)
    try
        rm(File)
    catch
    end

    Refinement = Refinement + 1
    nA, nB, nC = Refinement * 1, Refinement * 6, Refinement * 4
    xs = reshape(collect(linearspace(0.0, A, nA + 1)), nA + 1, 1)
    ys = reshape(collect(linearspace(0.0, B, nB + 1)), nB + 1, 1)
    zs = reshape(collect(linearspace(0.0, C, nC + 1)), nC + 1, 1)
    fensf, fesf = T10blockx(xs, ys, zs, :b)
    tolerance = min(A / nA, B / nB, C / nC) / 1000.0

    ff = ElementalField(zeros(count(fesf), 1))
    referenceff = ElementalField(zeros(count(fesf), 1))
    for i in eachindex(fesf)
        c = [k for k in fesf.conn[i]]
        centroid = NT * fensf.xyz[c, :]
        x, y, z = centroid
        referenceff.values[i, :] .= sin(2 * x / A) * cos(6.5 * y / B) * sin(3 * z / C - 1.0)
    end
    File = "momap2para3-reference.vtk"
    vtkexportmesh(File, fensf, fesf; scalars = [("referenceff", referenceff.values)])
    # @async run(`"paraview.exe" $File`)
    try
        rm(File)
    catch
    end

    ff = transferfield!(ff, fensf, fesf, fc, fensc, fesc, tolerance)
    File = "momap2para3-fine.vtk"
    vtkexportmesh(File, fensf, fesf; scalars = [("ff", ff.values), ("ffcopy", ff.values)])
    # @async run(`"paraview.exe" $File`)
    try
        rm(File)
    catch
    end

    diffff = ElementalField(referenceff.values - ff.values)
    femm = FEMMBase(IntegDomain(fesf, SimplexRule(3, 4)))
    geom = NodalField(fensf.xyz)
    error = integratefieldfunction(femm, geom, diffff, (x, v) -> norm(v); initial = 0.0)
    ref = integratefieldfunction(femm, geom, referenceff, (x, v) -> norm(v); initial = 0.0)
    # println("error/ref = $(error/ref)")
    @test abs(error / ref - 0.19808425992688541) < 1.0e-4
end
end
using .momap2para6378
momap2para6378.test()

module mesh_Q4spheren
using FinEtools
using FinEtools.MeshExportModule
using Test
function test()
    rex = 2.0 #external radius
    nr = 13

    fens, fes = Q4spheren(rex, nr)

    # # Postprocessing
    # vtkexportmesh("sphere.vtk", fes.conn, fens.xyz,  FinEtools.MeshExportModule.Q4)

    # @show count(fens), count(fes)
    @test count(fens) == 169
    @test count(fes) == 147
end

end
using .mesh_Q4spheren
mesh_Q4spheren.test()

module mesh_Q4circlen
using FinEtools
using FinEtools.MeshExportModule
using Test
function test()
    rex = 2.0 #external radius
    nr = 14

    fens, fes = Q4circlen(rex, nr)

    # # Postprocessing
    # vtkexportmesh("circle.vtk", fes.conn, fens.xyz,  FinEtools.MeshExportModule.Q4)

    # @show count(fens), count(fes)
    @test count(fens) == 169
    @test count(fes) == 147
end

end
using .mesh_Q4circlen
mesh_Q4circlen.test()

module mesh_H8cylindern
using FinEtools
using FinEtools.MeshExportModule
using Test
function test()
    Radius, Length, nperradius, nL = 1.0, 2.0, 7, 9

    fens, fes = H8cylindern(Radius, Length, nperradius, nL)

    # # Postprocessing
    # vtkexportmesh("cylinder.vtk", fes.conn, fens.xyz,  FinEtools.MeshExportModule.H8)

    # @show count(fens), count(fes)all
    @test (count(fens), count(fes)) == (2090, 1728)
end

end
using .mesh_H8cylindern
mesh_H8cylindern.test()

module mesh_H8cylindern_1
using FinEtools
using FinEtools.MeshExportModule
using Test
function test()
    Radius, Length, nperradius, nL = 1.0, 2.0, 17, 3

    fens, fes = H8cylindern(Radius, Length, nperradius, nL)

    # # Postprocessing
    # vtkexportmesh("cylinder.vtk", fes.conn, fens.xyz,  FinEtools.MeshExportModule.H8)

    geom = NodalField(fens.xyz)

    femm = FEMMBase(IntegDomain(fes, GaussRule(3, 2)))
    V = integratefunction(femm, geom, (x) -> 1.0)

    @test abs(V - pi * Radius^2 * Length) / (pi * Radius^2 * Length) < 1.30e-3
end

end
using .mesh_H8cylindern_1
mesh_H8cylindern_1.test()

module mesh_rcm_1
using FinEtools
using FinEtools.MeshExportModule
using SymRCM: adjgraph, nodedegrees, symrcm
using Test
using SparseArrays

function test()
    conn = [
        9 1 8 4
        1 3 2 8
        8 2 7 5
        2 6 7 7
    ]
    nfens = 9
    ag = adjgraph(conn, nfens)
    nd = nodedegrees(ag)
    @test ag == Array{Int64,1}[
        [9, 8, 4, 3, 2],
        [1, 3, 8, 7, 5, 6],
        [1, 2, 8],
        [9, 1, 8],
        [8, 2, 7],
        [2, 7],
        [8, 2, 5, 6],
        [9, 1, 4, 3, 2, 7, 5],
        [1, 8, 4],
    ]
    @test nd == [5, 6, 3, 3, 3, 2, 4, 7, 3]
    numbering = symrcm(ag, nd)
    @test numbering == [4, 9, 5, 8, 3, 1, 7, 2, 6]


    A = sprand(nfens, nfens, 1 / nfens)
    A = A + A'
    # display(spy(A))
    # display(spy(A[numbering, numbering]))
end

end
using .mesh_rcm_1
mesh_rcm_1.test()


module mesh_rcm_2
using FinEtools
using FinEtools.MeshExportModule
using SymRCM: adjgraph, nodedegrees, symrcm
using Test
using SparseArrays

function test()
    nfens = 29

    # A = sprand(nfens, nfens, 1/nfens)
    # A = A+A'
    A = spzeros(nfens, nfens)
    A[5, 1] = 0.559559
    A[24, 1] = 0.079212
    A[19, 2] = 0.459102
    A[8, 3] = 0.844709
    A[16, 3] = 0.206808
    A[24, 4] = 0.82036
    A[1, 5] = 0.559559
    A[7, 5] = 0.595562
    A[28, 6] = 0.151036
    A[5, 7] = 0.595562
    A[3, 8] = 0.844709
    A[24, 8] = 0.47093
    A[21, 9] = 0.673707
    A[22, 9] = 0.492159
    A[24, 9] = 0.10736
    A[19, 10] = 0.99992
    A[18, 11] = 0.573561
    A[22, 11] = 0.803174
    A[17, 12] = 0.127183
    A[28, 12] = 0.644722
    A[29, 14] = 0.839357
    A[3, 16] = 0.206808
    A[20, 16] = 0.0470789
    A[12, 17] = 0.127183
    A[11, 18] = 0.573561
    A[2, 19] = 0.459102
    A[10, 19] = 0.99992
    A[16, 20] = 0.0470789
    A[26, 20] = 0.661536
    A[9, 21] = 0.673707
    A[9, 22] = 0.492159
    A[11, 22] = 0.803174
    A[24, 23] = 0.373656
    A[1, 24] = 0.079212
    A[4, 24] = 0.82036
    A[8, 24] = 0.47093
    A[9, 24] = 0.10736
    A[23, 24] = 0.373656
    A[20, 26] = 0.661536
    A[6, 28] = 0.151036
    A[12, 28] = 0.644722
    A[14, 29] = 0.839357

    ag = adjgraph(A; sortbydeg = true)
    nd = nodedegrees(ag)
    numbering = symrcm(ag, nd)
    # display(spy(A))
    # display(spy(A[numbering, numbering]))
end

end
using .mesh_rcm_2
mesh_rcm_2.test()


module mesh_rcm_3
using FinEtools
using FinEtools.MeshExportModule
using SymRCM: adjgraph, nodedegrees, symrcm
using Test
using SparseArrays

function test()
    nfens = 19

    # A1 = sprand(nfens, nfens, 1/nfens)
    # @show A1 = A1+A1'
    A1 = spzeros(nfens, nfens)
    A1[7, 1] = 0.783374
    A1[18, 1] = 0.6411
    A1[8, 2] = 0.66032
    A1[13, 2] = 0.552169
    A1[5, 3] = 0.522678
    A1[11, 4] = 0.244274
    A1[3, 5] = 0.522678
    A1[19, 5] = 0.870687
    A1[15, 6] = 0.254443
    A1[17, 6] = 0.138423
    A1[1, 7] = 0.783374
    A1[8, 7] = 0.274651
    A1[11, 7] = 0.255421
    A1[15, 7] = 0.961861
    A1[2, 8] = 0.66032
    A1[7, 8] = 0.274651
    A1[11, 8] = 0.0421145
    A1[4, 11] = 0.244274
    A1[7, 11] = 0.255421
    A1[8, 11] = 0.0421145
    A1[12, 11] = 0.610131
    A1[16, 11] = 0.678996
    A1[11, 12] = 0.610131
    A1[19, 12] = 0.510702
    A1[2, 13] = 0.552169
    A1[18, 13] = 0.0696182
    A1[14, 14] = 0.213021
    A1[15, 14] = 0.516788
    A1[6, 15] = 0.254443
    A1[7, 15] = 0.961861
    A1[14, 15] = 0.516788
    A1[17, 15] = 0.34131
    A1[11, 16] = 0.678996
    A1[6, 17] = 0.138423
    A1[15, 17] = 0.34131
    A1[1, 18] = 0.6411
    A1[13, 18] = 0.0696182
    A1[5, 19] = 0.870687
    A1[12, 19] = 0.510702

    A = vcat(
        hcat(A1, 0 * A1, 0 * A1),
        hcat(0 * A1, 0 * A1, 0 * A1),
        hcat(0 * A1, 0 * A1, A1),
    )
    ag = adjgraph(A; sortbydeg = true)
    nd = nodedegrees(ag)
    numbering = symrcm(ag, nd)
    # display(spy(A))
    # display(spy(A[numbering, numbering]))
    @test numbering[end-4:end] == [22, 21, 20, 10, 9]
end

end
using .mesh_rcm_3
mesh_rcm_3.test()

module mesh_rcm_4
using FinEtools
using FinEtools.MeshExportModule
using SymRCM: adjgraph, nodedegrees, symrcm
using Test
using SparseArrays
using LinearAlgebra: norm

function test()
    nfens = 15000

    A = sprand(nfens, nfens, 1 / nfens)
    A = A + A'
    ag = adjgraph(A; sortbydeg = true)
    nd = nodedegrees(ag)
    numbering = symrcm(ag, nd)
    # display(spy(A))
    # display(spy(A[numbering, numbering]))
    b = rand(nfens)
    inumbering = deepcopy(numbering)
    inumbering[numbering] = 1:nfens
    x = A * b
    xp = A[numbering, numbering] * b[numbering]
    xp[numbering]
    @test norm(x - xp[inumbering]) < 1.0e-5 * nfens
end

end
using .mesh_rcm_4
mesh_rcm_4.test()

module mesh_triangle_conversion_4
using FinEtools
using FinEtools.MeshExportModule
using Test
using SparseArrays
using LinearAlgebra: norm, I

function test()
    rho = 1.21 * 1e-9# mass density
    c = 345.0 * 1000# millimeters per second
    bulk = c^2 * rho
    Lx = 1900.0# length of the box, millimeters
    Ly = 800.0 # length of the box, millimeters

    fens, fes = T3block(Lx, Ly, 3, 2) # Mesh
    @test (count(fens), count(fes)) == (12, 12)
    fens, fes = T3toT6(fens, fes)
    @test (count(fens), count(fes)) == (35, 12)
end

end
using .mesh_triangle_conversion_4
mesh_triangle_conversion_4.test()

module mmT4ttH80
using FinEtools
using Test
function test()
    a, b, h, na, nb, nh = 1.0, 2.0, 3.0, 2, 3, 4
    fens, fes = T4block(a, b, h, na, nb, nh)
    fens, fes = T4toH8(fens, fes)

    geom = NodalField(fens.xyz)
    u = NodalField(zeros(size(fens.xyz, 1), 3)) # displacement field

    # File =  "T4toH8-mesh.vtk"
    # vtkexportmesh(File, fes.conn, fens.xyz, FinEtools.MeshExportModule.H8)
    # @async run(`"paraview.exe" $File`)

    # @show count(fes)
    @test count(fes) == 576
end

end
using .mmT4ttH80
mmT4ttH80.test()

module mmT4ttH81
using FinEtools
using Test
function test()
    a, b, h, na, nb, nh = 1.0, 2.0, 3.0, 2, 3, 4
    fens, fes = T4block(a, b, h, na, nb, nh)
    fens, fes = T4toH8(fens, fes)

    geom = NodalField(fens.xyz)
    u = NodalField(zeros(size(fens.xyz, 1), 3)) # displacement field

    bfes = meshboundary(fes)
    # @show count(bfes)
    # File =  "T4toH8-mesh.vtk"
    # vtkexportmesh(File, bfes.conn, fens.xyz, FinEtools.MeshExportModule.Q4)
    # @async run(`"paraview.exe" $File`)

    @test count(bfes) == 312
end

end
using .mmT4ttH81
mmT4ttH81.test()

module mmT4refine20a
using FinEtools
using Test
function test()
    a, b, h, na, nb, nh = 1.0, 2.0, 3.0, 2, 3, 4
    fens, fes = T4block(a, b, h, na, nb, nh)
    fens, fes = T4refine20(fens, fes)

    geom = NodalField(fens.xyz)
    u = NodalField(zeros(size(fens.xyz, 1), 3)) # displacement field

    # File =  "T4refine20-mesh.vtk"
    # vtkexportmesh(File, fes.conn, fens.xyz, FinEtools.MeshExportModule.T4)
    # @async run(`"paraview.exe" $File`)

    # @show count(fes)
    @test count(fes) == 2880
end

end
using .mmT4refine20a
mmT4refine20a.test()

module mmT4refine20b
using FinEtools
using Test
function test()
    a, b, h, na, nb, nh = 1.0, 2.0, 3.0, 2, 3, 4
    fens, fes = T4block(a, b, h, na, nb, nh)
    fens, fes = T4refine20(fens, fes)

    geom = NodalField(fens.xyz)
    u = NodalField(zeros(size(fens.xyz, 1), 3)) # displacement field

    bfes = meshboundary(fes)
    # @show count(bfes)
    # File =  "T4refine20b-mesh.vtk"
    # vtkexportmesh(File, bfes.conn, fens.xyz, FinEtools.MeshExportModule.T3)
    # @async run(`"paraview.exe" $File`)

    # @show count(fes)
    @test count(bfes) == 624
end

end
using .mmT4refine20b
mmT4refine20b.test()

module mq4patchtest
using FinEtools
using Test
function test()
    # println("Q4. Plane stress.")

    E = 1.0
    nu = 1.0 / 3
    alpha, beta, gamma, delta, eta, phi =
        1.0 / 30, 1.0 / 34, -1.0 / 21, -1.0 / 51, -1.0 / 26, -1.0 / 35
    ux(x, y) = alpha + beta * x + gamma * y
    uy(x, y) = delta + eta * x + phi * y

    fens = FENodeSet(
        [1.0 -0.3; 2.3 -0.3; 2.3 0.95; 1.0 0.95; 1.4 0.05; 1.9 -0.03; 1.7 0.5; 1.3 0.6],
    )
    fes = FESetQ4([1 2 6 5; 6 2 3 7; 7 3 4 8; 8 4 1 5; 5 6 7 8])

    geom = NodalField(fens.xyz)
    u = NodalField(zeros(size(fens.xyz, 1), 2)) # displacement field

    # Apply prescribed displacements to exterior nodes
    for i in 1:4
        setebc!(u, [i], 1, ux(fens.xyz[i, :]...))
        setebc!(u, [i], 2, uy(fens.xyz[i, :]...))
    end

    applyebc!(u)
    numberdofs!(u)

    # for i in 5:8
    # 	uexact = [ux(fens.xyz[i, :]...), uy(fens.xyz[i, :]...)]
    # 	println("u.values[$i, :] = $(u.values[i, :]), uexact = [$(uexact)]")
    # end

    AE = AbaqusExporter("q4_stress_export")
    HEADING(AE, "q4_stress_export")
    COMMENT(AE, "")
    PART(AE, "part1")
    END_PART(AE)
    ASSEMBLY(AE, "ASSEM1")
    INSTANCE(AE, "INSTNC1", "PART1")
    NODE(AE, fens.xyz)
    ELEMENT(AE, "CPS4", "AllElements", connasarray(fes))
    NSET_NSET(AE, "clamped", 1:4)
    ORIENTATION(AE, "GlobalOrientation", vec([1.0 0 0]), vec([0 1.0 0]))
    SOLID_SECTION(AE, "elasticity", "GlobalOrientation", "AllElements", 1.0)
    END_INSTANCE(AE)
    END_ASSEMBLY(AE)
    MATERIAL(AE, "elasticity")
    ELASTIC(AE, E, nu)
    STEP_PERTURBATION_STATIC(AE)
    BOUNDARY(
        AE,
        "ASSEM1.INSTNC1",
        1:4,
        fill(true, 4, 2),
        [[ux(fens.xyz[i, :]...) for i in 1:4] [uy(fens.xyz[i, :]...) for i in 1:4]],
    )
    END_STEP(AE)
    close(AE)
    s = readlines("q4_stress_export.inp")
    @test length(s) == 49
    try
        rm("q4_stress_export.inp")
    catch
    end

    true

end
end
using .mq4patchtest
mq4patchtest.test()


module mq4patchtest1
using FinEtools
using Test
function test()
    # println("Q4. Plane stress.")

    E = 1.0
    nu = 1.0 / 3
    alpha, beta, gamma, delta, eta, phi =
        1.0 / 30, 1.0 / 34, -1.0 / 21, -1.0 / 51, -1.0 / 26, -1.0 / 35
    ux(x, y) = alpha + beta * x + gamma * y
    uy(x, y) = delta + eta * x + phi * y

    fens = FENodeSet(
        [1.0 -0.3; 2.3 -0.3; 2.3 0.95; 1.0 0.95; 1.4 0.05; 1.9 -0.03; 1.7 0.5; 1.3 0.6],
    )
    fes = FESetQ4([1 2 6 5; 6 2 3 7; 7 3 4 8; 8 4 1 5; 5 6 7 8])

    geom = NodalField(fens.xyz)
    u = NodalField(zeros(size(fens.xyz, 1), 2)) # displacement field

    # Apply prescribed displacements to exterior nodes
    for i in 1:4
        setebc!(u, [i], 1, ux(fens.xyz[i, :]...))
        setebc!(u, [i], 2, uy(fens.xyz[i, :]...))
    end

    applyebc!(u)
    numberdofs!(u)

    # for i in 5:8
    # 	uexact = [ux(fens.xyz[i, :]...), uy(fens.xyz[i, :]...)]
    # 	println("u.values[$i, :] = $(u.values[i, :]), uexact = [$(uexact)]")
    # end

    AE = AbaqusExporter("q4_stress_export1")
    HEADING(AE, "q4_stress_export")
    COMMENT(AE, "")
    PART(AE, "part1")
    END_PART(AE)
    ASSEMBLY(AE, "ASSEM1")
    INSTANCE(AE, "INSTNC1", "PART1")
    NODE(AE, fens.xyz)
    ELEMENT(AE, "CPS4", "AllElements", connasarray(fes))
    NSET_NSET(AE, "clamped", 1:4)
    ORIENTATION(AE, "GlobalOrientation", vec([1.0 0 0]), vec([0 1.0 0]))
    SOLID_SECTION(AE, "elasticity", "GlobalOrientation", "AllElements", 1.0)
    END_INSTANCE(AE)
    END_ASSEMBLY(AE)
    MATERIAL(AE, "elasticity")
    ELASTIC(AE, E, nu)
    STEP_PERTURBATION_STATIC(AE)
    BOUNDARY(AE, "clamped", 1, -1.0)
    END_STEP(AE)
    close(AE)
    s = readlines("q4_stress_export1.inp")
    @test length(s) == 42
    try
        rm("q4_stress_export1.inp")
    catch
    end

    true

end
end
using .mq4patchtest1
mq4patchtest1.test()

module mq4patchtest2
using FinEtools
using Test
function test()
    # println("Q4. Plane stress.")

    E = 1.0
    nu = 1.0 / 3
    alpha, beta, gamma, delta, eta, phi =
        1.0 / 30, 1.0 / 34, -1.0 / 21, -1.0 / 51, -1.0 / 26, -1.0 / 35
    ux(x, y) = alpha + beta * x + gamma * y
    uy(x, y) = delta + eta * x + phi * y

    fens = FENodeSet(
        [1.0 -0.3; 2.3 -0.3; 2.3 0.95; 1.0 0.95; 1.4 0.05; 1.9 -0.03; 1.7 0.5; 1.3 0.6],
    )
    fes = FESetQ4([1 2 6 5; 6 2 3 7; 7 3 4 8; 8 4 1 5; 5 6 7 8])

    geom = NodalField(fens.xyz)
    u = NodalField(zeros(size(fens.xyz, 1), 2)) # displacement field

    # Apply prescribed displacements to exterior nodes
    for i in 1:4
        setebc!(u, [i], 1, ux(fens.xyz[i, :]...))
        setebc!(u, [i], 2, uy(fens.xyz[i, :]...))
    end

    applyebc!(u)
    numberdofs!(u)

    # for i in 5:8
    # 	uexact = [ux(fens.xyz[i, :]...), uy(fens.xyz[i, :]...)]
    # 	println("u.values[$i, :] = $(u.values[i, :]), uexact = [$(uexact)]")
    # end

    AE = AbaqusExporter("q4_stress_export2")
    HEADING(AE, "q4_stress_export2")
    COMMENT(AE, "")
    PART(AE, "part1")
    END_PART(AE)
    ASSEMBLY(AE, "ASSEM1")
    INSTANCE(AE, "INSTNC1", "PART1")
    NODE(AE, fens.xyz)
    ELEMENT(AE, "CPS4", "AllElements", connasarray(fes))
    NSET_NSET(AE, "clamped", 1:4)
    ORIENTATION(AE, "GlobalOrientation", vec([1.0 0 0]), vec([0 1.0 0]))
    SOLID_SECTION(AE, "elasticity", "GlobalOrientation", "AllElements", 1.0)
    END_INSTANCE(AE)
    END_ASSEMBLY(AE)
    MATERIAL(AE, "elasticity")
    ELASTIC(AE, E, nu)
    STEP_PERTURBATION_STATIC(AE)
    BOUNDARY(AE, "ASSEM1.INSTNC1", u.is_fixed, u.values)
    END_STEP(AE)
    close(AE)
    s = readlines("q4_stress_export2.inp")
    @test length(s) == 49
    try
        rm("q4_stress_export2.inp")
    catch
    end

    true

end
end
using .mq4patchtest2
mq4patchtest2.test()

module mq4patchtest3a
using FinEtools
using Test
function test()
    # println("Q4. Plane stress.")

    E = 1.0
    nu = 1.0 / 3
    alpha, beta, gamma, delta, eta, phi =
        1.0 / 30, 1.0 / 34, -1.0 / 21, -1.0 / 51, -1.0 / 26, -1.0 / 35
    ux(x, y) = alpha + beta * x + gamma * y
    uy(x, y) = delta + eta * x + phi * y

    fens = FENodeSet(
        [1.0 -0.3; 2.3 -0.3; 2.3 0.95; 1.0 0.95; 1.4 0.05; 1.9 -0.03; 1.7 0.5; 1.3 0.6],
    )
    fes = FESetQ4([1 2 6 5; 6 2 3 7; 7 3 4 8; 8 4 1 5; 5 6 7 8])

    geom = NodalField(fens.xyz)
    u = NodalField(zeros(size(fens.xyz, 1), 2)) # displacement field

    # Apply prescribed displacements to exterior nodes
    for i in 1:4
        setebc!(u, [i], 1, ux(fens.xyz[i, :]...))
        setebc!(u, [i], 2, uy(fens.xyz[i, :]...))
    end

    applyebc!(u)
    numberdofs!(u)

    # for i in 5:8
    # 	uexact = [ux(fens.xyz[i, :]...), uy(fens.xyz[i, :]...)]
    # 	println("u.values[$i, :] = $(u.values[i, :]), uexact = [$(uexact)]")
    # end

    AE = AbaqusExporter("q4_stress_export3a")
    HEADING(AE, "q4_stress_export3a")
    COMMENT(AE, "")
    PART(AE, "part1")
    END_PART(AE)
    ASSEMBLY(AE, "ASSEM1")
    INSTANCE(AE, "INSTNC1", "PART1")
    NODE(AE, fens.xyz)
    ELEMENT(AE, "CPS4", "AllElements", 1, connasarray(fes))
    NSET_NSET(AE, "clamped", 1:4)
    ORIENTATION(AE, "GlobalOrientation", vec([1.0 0 0]), vec([0 1.0 0]))
    SOLID_SECTION(AE, "elasticity", "GlobalOrientation", "AllElements", 1.0)
    END_INSTANCE(AE)
    END_ASSEMBLY(AE)
    MATERIAL(AE, "elasticity")
    ELASTIC(AE, E, E, E, nu, nu, nu, E / 2 / (1 + nu), E / 2 / (1 + nu), E / 2 / (1 + nu))
    EXPANSION(AE, 0.0)
    DENSITY(AE, 1.0)
    STEP_PERTURBATION_STATIC(AE)
    BOUNDARY(AE, "ASSEM1.INSTNC1", u.is_fixed, u.values)
    END_STEP(AE)
    close(AE)
    s = readlines("q4_stress_export3a.inp")
    @test length(s) == 54
    try
        rm("q4_stress_export3a.inp")
    catch
    end

    true

end
end
using .mq4patchtest3a
mq4patchtest3a.test()

module mT4quartercylnm
using FinEtools
using FinEtools.MeshExportModule: VTK
using Test
function test()
    Radius, Length, nperradius, nL = 3.0, 5.0, 3, 5
    fens, fes = T4quartercyln(Radius, Length, nperradius, nL)
    @test count(fes) == 270
    File = "t4cyl-mesh.vtk"
    VTK.vtkexportmesh(File, fens, fes)
    # rm(File)
    # @async run(`"paraview.exe" $File`)
    bfes = meshboundary(fes)
    @test count(bfes) == 156
    File = "t4cyl-bmesh.vtk"
    VTK.vtkexportmesh(File, fens, fes)
    # rm(File)
    # @async run(`"paraview.exe" $File`)
end
end
using .mT4quartercylnm
mT4quartercylnm.test()

module mt10convt41
using FinEtools
using FinEtools.MeshImportModule
using FinEtools.MeshExportModule: VTK
using Test
function test()
    xs = collect(linearspace(0.0, 2.0, 5))
    ys = collect(linearspace(0.0, 1.0, 6) .^ 2)
    zs = collect(linearspace(0.0, 1.0, 3))
    fens, fes = T4blockx(xs, ys, zs, :a)
    @test count(fes) == 240
    fens, fes = T4toT10(fens, fes)
    @test count(fes) == 240
    fens, fes = T10toT4(fens, fes)
    @test count(fes) == 240
    true
end
end
using .mt10convt41
mt10convt41.test()

module msels2
using FinEtools
using FinEtools.MeshSelectionModule: vselect
using LinearAlgebra: norm, I
# using BenchmarkTools
using Test
function test()
    box = fill(0.0, 6)
    initbox!(box, [0.0, 0.0, 0.0])
    xyz = [
        240.0 0.0 0.0
        120.00000000000001 120.0 0.0
        0.0 1.469576158976824e-14 0.0
        119.99999999999997 -120.0 0.0
    ]
    tolerance = 0.03
    clampedn = vselect(xyz; box = box, inflate = tolerance)
    @test clampedn[1] == 3
    true
end
end
using .msels2
msels2.test()

module mt4vox1
using FinEtools
using Test
function test()
    img = fill(1, 30, 20, 10)
    fens, fes = T4voximg(img, Float64[1.0, 1.1, 1.3], [1])
    @test count(fes) == 30 * 20 * 10 * 5
    true
end
end
using .mt4vox1
mt4vox1.test()

module mt6block1
using FinEtools
using FinEtools.MeshExportModule: VTK
using Test
function test()
    rho = 1.21 * 1e-9# mass density
    c = 345.0 * 1000# millimeters per second
    bulk = c^2 * rho
    Lx = 1900.0# length of the box, millimeters
    Ly = 800.0 # length of the box, millimeters

    fens, fes = T6block(Lx, Ly, 3, 2) # Mesh
    # @show count(fens), count(fes)
    @test (count(fens), count(fes)) == (35, 12)
    fens.xyz = xyz3(fens)
    fens.xyz[:, 3] .= (fens.xyz[:, 1] .^ 2 + fens.xyz[:, 2] .^ 2) / 1.e3
    File = "mesh.vtk"
    VTK.vtkexportmesh(File, fens, fes)
    rm(File)
    # @async run(`"paraview.exe" $File`)
end
end
using .mt6block1
mt6block1.test()

module mt6annul1
using FinEtools
using FinEtools.MeshExportModule: VTK
using Test
function test()
    rin, rex, nr, nc, Angl, orientation::Symbol =
        100.0, 200.0, 3, 6, pi / 3, :a
    fens, fes =
        T6annulus(rin, rex, nr, nc, Angl, orientation::Symbol)
    # @show count(fens), count(fes)
    @test (count(fens), count(fes)) == (91, 36)
    fens.xyz = xyz3(fens)
    fens.xyz[:, 3] .= (fens.xyz[:, 1] .^ 2 + fens.xyz[:, 2] .^ 2) / 1.e3
    File = "mesh.vtk"
    VTK.vtkexportmesh(File, fens, fes)
    rm(File)
    # @async run(`"paraview.exe" $File`)
end
end
using .mt6annul1
mt6annul1.test()

module mt3annul1
using FinEtools
using FinEtools.MeshExportModule: VTK
using Test
function test()
    rin, rex, nr, nc, Angl, orientation::Symbol =
        100.0, 200.0, 3, 6, pi / 3, :a
    fens, fes =
        T3annulus(rin, rex, nr, nc, Angl, orientation::Symbol)
    # @show count(fens), count(fes)
    @test (count(fens), count(fes)) == (28, 36)
    fens.xyz = xyz3(fens)
    fens.xyz[:, 3] .= (fens.xyz[:, 1] .^ 2 + fens.xyz[:, 2] .^ 2) / 1.e3
    File = "mesh.vtk"
    VTK.vtkexportmesh(File, fens, fes)
    rm(File)
    # @async run(`"paraview.exe" $File`)
end
end
using .mt3annul1
mt3annul1.test()

module mt3refi1
using FinEtools
using FinEtools.MeshExportModule: VTK
using Test
function test()
    rin, rex, nr, nc, Angl, orientation::Symbol =
        100.0, 200.0, 3, 6, pi / 3, :a
    fens, fes =
        T3annulus(rin, rex, nr, nc, Angl, orientation::Symbol)
    fens, fes = T3refine(fens, fes)
    # @show count(fens), count(fes)
    @test (count(fens), count(fes)) == (91, 144)
    fens.xyz = xyz3(fens)
    fens.xyz[:, 3] .= (fens.xyz[:, 1] .^ 2 + fens.xyz[:, 2] .^ 2) / 1.e3
    File = "mesh.vtk"
    VTK.vtkexportmesh(File, fens, fes)
    # rm(File)
    # @async run(`"paraview.exe" $File`)
end
end
using .mt3refi1
mt3refi1.test()

module mh8refix
using FinEtools
using FinEtools.MeshExportModule: VTK
using Test
function test()
    fens, fes = H8block(3.1, 10.1, 7.7, 3, 9, 6)
    fens, fes = H8refine(fens, fes)
    @test (count(fens), count(fes)) == (1729, 1296)
    # File = "mesh.vtk"
    # VTK.vtkexportmesh(File, fens, fes)
    # rm(File)
    # @async run(`"/c/Program\ Files/ParaView\ 5.7.0-Windows-Python3.7-msvc2015-64bit/bin/paraview" $File`)
end
end
using .mh8refix
mh8refix.test()

module mh8sphere1
using FinEtools
using FinEtools.MeshExportModule: VTK
using Test
function test()
    fens, fes = H8sphere(3.1, 1)
    @test (count(fens), count(fes)) == (65, 32)
    # File = "mesh.vtk"
    # VTK.vtkexportmesh(File, fens, fes)
    # rm(File)
    # @async run(`"/c/Program\ Files/ParaView\ 5.7.0-Windows-Python3.7-msvc2015-64bit/bin/paraview" $File`)
end
end
using .mh8sphere1
mh8sphere1.test()

module mh8spheren1
using FinEtools
using FinEtools.MeshExportModule: VTK
using Test
function test()
    fens, fes = H8spheren(3.1, 5)
    @test (count(fens), count(fes)) == (175, 108)
    File = "mesh.vtk"
    VTK.vtkexportmesh(File, fens, fes)
    # rm(File)
    # @async run(`"/c/Program\ Files/ParaView\ 5.7.0-Windows-Python3.7-msvc2015-64bit/bin/paraview" $File`)
end
end
using .mh8spheren1
mh8spheren1.test()

module mh8voximg1
using FinEtools
using FinEtools.MeshExportModule: VTK
using Test
function test()
    img = fill(1, 30, 20, 10)
    fens, fes = H8voximg(img, Float64[1.0, 1.1, 1.3], [1])
    @test (count(fens), count(fes)) == (7161, 6000)
    # File = "mesh.vtk"
    # VTK.vtkexportmesh(File, fens, fes)
    # rm(File)
    # @async run(`"/c/Program\ Files/ParaView\ 5.7.0-Windows-Python3.7-msvc2015-64bit/bin/paraview" $File`)
    geom = NodalField(fens.xyz)
    femm = FEMMBase(IntegDomain(fes, GaussRule(3, 2)))
    # Test the volume
    V = integratefunction(femm, geom, (x) -> 1.0)
    @test abs(V - (*).(30, 20, 10) * (*).(1.0, 1.1, 1.3)) / V < 1.0e-5
    true
end
end
using .mh8voximg1
mh8voximg1.test()

module mh8hex2pt
using FinEtools
using FinEtools.MeshExportModule: VTK
using Test
function test()
    fens, fes = H8hexahedron([-1.0 -1.0 -1.0; 1.0 1.0 1.0], 2, 3, 4)
    @test (count(fens), count(fes)) == (60, 24)
    # File = "mesh.vtk"
    # VTK.vtkexportmesh(File, fens, fes)
    # rm(File)
    # @async run(`"/c/Program\ Files/ParaView\ 5.7.0-Windows-Python3.7-msvc2015-64bit/bin/paraview" $File`)
    geom = NodalField(fens.xyz)
    femm = FEMMBase(IntegDomain(fes, GaussRule(3, 2)))
    # Test the volume
    V = integratefunction(femm, geom, (x) -> 1.0)
    @test abs(V - 8.0) / V < 1.0e-5
    true
end
end
using .mh8hex2pt
mh8hex2pt.test()

module mexpmulvf1
using FinEtools
using FinEtools.MeshExportModule: VTK
using Test
function test()
    rin, rex, nr, nc, Angl, orientation::Symbol =
        100.0, 200.0, 13, 16, pi / 3, :a
    fens, fes =
        T3annulus(rin, rex, nr, nc, Angl, orientation::Symbol)
    d = (fens.xyz[:, 1] .^ 2 + fens.xyz[:, 2] .^ 2)
    File = "mesh.vtk"
    # Export of multiple scalar fields
    result = VTK.vtkexportmesh(File, fens, fes; scalars = [("d", d), ("invd", 1 ./ d)])
    @test result == true
    rm(File)
    # @async run(`"paraview.exe" $File`)
end
end
using .mexpmulvf1
mexpmulvf1.test()

module mexpstl1
using FinEtools
using FinEtools.MeshExportModule: STL
using Test
function test()
    rin, rex, nr, nc, Angl, orientation::Symbol =
        1.0, 2.0, 7, 16, pi / 3, :a
    fens, fes =
        T3annulus(rin, rex, nr, nc, Angl, orientation::Symbol)
    fens.xyz = xyz3(fens)
    filename = "mesh.stl"
    e = STLExporter(filename::AbstractString)
    solid(e, "thesolid")
    for i in eachindex(fes)
        c = fes.conn[i]
        facet(e, fens.xyz[c[1], :], fens.xyz[c[2], :], fens.xyz[c[3], :])
    end
    endsolid(e, "thesolid")
    close(e)
    @test filesize(filename) > 0
    rm(filename)
    # @async run(`"paraview.exe" $File`)
end
end
using .mexpstl1
mexpstl1.test()


module mexpvecv1
using FinEtools
using FinEtools.MeshExportModule: VTK
using Test
function test()
    rin, rex, nr, nc, Angl, orientation::Symbol =
        100.0, 200.0, 13, 16, pi / 3, :a
    fens, fes =
        T3annulus(rin, rex, nr, nc, Angl, orientation::Symbol)
    d = (fens.xyz[:, 1] .^ 2 + fens.xyz[:, 2] .^ 2)
    File = "mesh.vtk"
    # Export of multiple scalar fields
    result = VTK.vtkexportmesh(File, fens, fes; scalars = [("d", d), ("invd", 1 ./ d)])
    @test result == true
    rm(File)

    vectors = [("grad", deepcopy(fens.xyz))]
    filename = "vectors.vtk"
    result = VTK.vtkexportmesh(filename, fens, fes; vectors = vectors)
    @test result == true
    rm(filename)
end
end
using .mexpvecv1
mexpvecv1.test()

module mexpvecv2
using FinEtools
using FinEtools.MeshExportModule: VTK
using Test
function test()
    rin, rex, nr, nc, Angl, orientation::Symbol =
        100.0, 200.0, 3, 6, pi / 3, :a
    fens, fes =
        T3annulus(rin, rex, nr, nc, Angl, orientation::Symbol)
    fens.xyz = xyz3(fens)
    d = (fens.xyz[:, 1] .^ 2 + fens.xyz[:, 2] .^ 2)
    File = "mesh.vtk"
    # Export of multiple scalar fields
    result = VTK.vtkexportmesh(File, fens, fes; scalars = [("d", d), ("invd", 1 ./ d)])
    @test result == true
    rm(File)

    points = [vec(fens.xyz[idx, :]) for idx in eachindex(fens)]
    vectors = [("grad", [vec(fens.xyz[idx, [2, 1, 3]]) for idx in eachindex(fens)])]
    filename = "vectors.vtk"
    result = VTK.vtkexportvectors(filename, points, vectors)
    @test result == true
    rm(filename)
end
end
using .mexpvecv2
mexpvecv2.test()

module mt4nastranexp1
using FinEtools
using FinEtools.MeshExportModule
using FinEtools.MeshImportModule
using Test
function test()
    xs = collect(linearspace(0.0, pi / 2, 5))
    ys = collect(linearspace(0.0, 1.0, 6) .^ 2)
    zs = collect(linearspace(0.0, 1.0, 3))
    fens, fes = T4blockx(xs, ys, zs, :a)
    fens, fes = T4toT10(fens, fes)
    @test count(fes) == 240
    bfes = meshboundary(fes)
    @test count(bfes) == 2 * 2 * (4 * 5 + 5 * 2 + 4 * 2)
    fens, fes = T10refine(fens, fes)
    for i in eachindex(fens)
        a, y, z = fens.xyz[i, :]
        fens.xyz[i, 1] = sin(a) * (y + 0.5)
        fens.xyz[i, 2] = cos(a) * (y + 0.5)
        fens.xyz[i, 3] = z
    end
    @test count(fes) == 240 * 8
    bfes = meshboundary(fes)
    @test count(bfes) == 4 * 2 * 2 * (4 * 5 + 5 * 2 + 4 * 2)

    geom = NodalField(fens.xyz)
    femm = FEMMBase(IntegDomain(fes, SimplexRule(3, 4)))

    File = "Refine-T10-a.vtk"
    vtkexportmesh(File, fens, bfes)
    try
        rm(File)
    catch
    end
    # @async run(`"paraview.exe" $File`)
    e = NASTRANExporter("Refine-T10-a.nas")
    BEGIN_BULK(e)
    for i in eachindex(fens)
        GRID(e, i, fens.xyz[i, :])
    end
    conn = connasarray(fes)
    for i in eachindex(fes)
        CTETRA(e, i, 1, conn[i, :])
    end
    ENDDATA(e)
    close(e)

    output = MeshImportModule.import_NASTRAN(
        "Refine-T10-a.nas";
        allocationchunk = 13,
        expectfixedformat = false,
    )
    @test count(output["fens"]) == 3213
    @test count(output["fesets"][1]) == 1920
    rm("Refine-T10-a.nas")
end
end
using .mt4nastranexp1
mt4nastranexp1.test()

module mefexp22
using FinEtools
using LinearAlgebra
using FinEtools.MeshExportModule: VTK
using FinEtools.AlgoBaseModule: evalconvergencestudy
using Test
function test()
    Lx = 1900.0# length of the box, millimeters
    Ly = 800.0 # length of the box, millimeters
    nx, ny = 10, 9
    th = 7.0
    f(x) = cos(0.93 * pi * x[1] / Lx) + sin(1.7 * pi * x[2] / Ly)
    g(x) = -cos(0.93 * pi * x[1] / Lx) + 2 * sin(1.7 * pi * x[2] / Ly)

    fens, fes = Q4block(Lx, Ly, nx, ny) # Mesh
    femm = FEMMBase(IntegDomain(fes, GaussRule(2, 2), th))
    geom = NodalField(fens.xyz)
    psi = NodalField([f(fens.xyz[idx, :]) for idx in eachindex(fens)])
    numberdofs!(psi)
    psie = ElementalField([
        f(
            (
                fens.xyz[fes.conn[idx][1], :] +
                fens.xyz[fes.conn[idx][2], :] +
                fens.xyz[fes.conn[idx][3], :] +
                fens.xyz[fes.conn[idx][4], :]
            ) ./ 4,
        ) for idx in eachindex(fes)
    ])
    phie = ElementalField([
        g(
            (
                fens.xyz[fes.conn[idx][1], :] +
                fens.xyz[fes.conn[idx][2], :] +
                fens.xyz[fes.conn[idx][3], :] +
                fens.xyz[fes.conn[idx][4], :]
            ) ./ 4,
        ) for idx in eachindex(fes)
    ])
    File = "mesh.vtk"
    # Export of multiple scalar fields
    result = VTK.vtkexportmesh(
        File,
        fens,
        fes;
        scalars = [("psie", psie.values), ("phie", phie.values)],
    )
    @test result == true
    rm(File)
    true
end
end
using .mefexp22
mefexp22.test()

module mefexp23
using FinEtools
using LinearAlgebra
using FinEtools.MeshExportModule: VTK
using FinEtools.AlgoBaseModule: evalconvergencestudy
using Test
function test()
    Lx = 1900.0# length of the box, millimeters
    Ly = 800.0 # length of the box, millimeters
    nx, ny = 10, 9
    th = 7.0
    f(x) = cos(0.93 * pi * x[1] / Lx) + sin(1.7 * pi * x[2] / Ly)
    g(x) = -cos(0.93 * pi * x[1] / Lx) + 2 * sin(1.7 * pi * x[2] / Ly)

    fens, fes = Q4block(Lx, Ly, nx, ny) # Mesh
    femm = FEMMBase(IntegDomain(fes, GaussRule(2, 2), th))
    geom = NodalField(fens.xyz)
    psi = NodalField(fill(0.0, count(fens), 2))
    for i in eachindex(fens)
        psi.values[i, 1] = f(fens.xyz[i, :])
        psi.values[i, 2] = g(fens.xyz[i, :])
    end
    numberdofs!(psi)
    # psie  = ElementalField(
    #     [f((fens.xyz[fes.conn[idx][1], :]+fens.xyz[fes.conn[idx][2], :]+fens.xyz[fes.conn[idx][3], :]+fens.xyz[fes.conn[idx][4], :])./4) 
    #     for idx in eachindex(fes)])
    # phie  = ElementalField(
    #     [g((fens.xyz[fes.conn[idx][1], :]+fens.xyz[fes.conn[idx][2], :]+fens.xyz[fes.conn[idx][3], :]+fens.xyz[fes.conn[idx][4], :])./4) 
    #     for idx in eachindex(fes)])
    File = "mesh.vtk"
    # Export of multiple scalar fields
    result = VTK.vtkexportmesh(File, fens, fes; scalars = [("psi", psi.values)])
    @test result == true
    rm(File)
    true
end
end
using .mefexp23
mefexp23.test()

module mefexp24
using FinEtools
using LinearAlgebra
using FinEtools.MeshExportModule: VTK
using FinEtools.AlgoBaseModule: evalconvergencestudy
using Test
function test()
    Lx = 1900.0# length of the box, millimeters
    Ly = 800.0 # length of the box, millimeters
    nx, ny = 10, 9
    th = 7.0
    f(x) = cos(0.93 * pi * x[1] / Lx) + sin(1.7 * pi * x[2] / Ly)
    g(x) = -cos(0.93 * pi * x[1] / Lx) + 2 * sin(1.7 * pi * x[2] / Ly)

    fens, fes = Q4block(Lx, Ly, nx, ny) # Mesh
    femm = FEMMBase(IntegDomain(fes, GaussRule(2, 2), th))
    geom = NodalField(fens.xyz)
    psie = ElementalField(fill(0.0, count(fes), 2))
    for i in eachindex(fes)
        psie.values[i, 1] = f(
            (
                fens.xyz[fes.conn[i][1], :] +
                fens.xyz[fes.conn[i][2], :] +
                fens.xyz[fes.conn[i][3], :] +
                fens.xyz[fes.conn[i][4], :]
            ) ./ 4,
        )
        psie.values[i, 2] = g(
            (
                fens.xyz[fes.conn[i][1], :] +
                fens.xyz[fes.conn[i][2], :] +
                fens.xyz[fes.conn[i][3], :] +
                fens.xyz[fes.conn[i][4], :]
            ) ./ 4,
        )
    end
    numberdofs!(psie)
    File = "mesh.vtk"
    # Export of multiple scalar fields
    result = VTK.vtkexportmesh(File, fens, fes; scalars = [("psie", psie.values)])
    @test result == true
    rm(File)
    true
end
end
using .mefexp24
mefexp24.test()

module mtetvt1
using FinEtools
using FinEtools.MeshTetrahedronModule: tetv
using Test
function reftetv(X)
    local one6th = 1.0 / 6
    # @assert size(X, 1) == 4
    # @assert size(X, 2) == 3
    @inbounds let
        A1 = X[2, 1] - X[1, 1]
        A2 = X[2, 2] - X[1, 2]
        A3 = X[2, 3] - X[1, 3]
        B1 = X[3, 1] - X[1, 1]
        B2 = X[3, 2] - X[1, 2]
        B3 = X[3, 3] - X[1, 3]
        C1 = X[4, 1] - X[1, 1]
        C2 = X[4, 2] - X[1, 2]
        C3 = X[4, 3] - X[1, 3]
        return one6th * (
            (-A3 * B2 + A2 * B3) * C1 +
            (A3 * B1 - A1 * B3) * C2 +
            (-A2 * B1 + A1 * B2) * C3
        )
    end
end
function test()
    X = fill(0.0, 4, 3)
    xs = collect(linearspace(0.0, pi / 2, 5))
    ys = collect(linearspace(0.0, 1.0, 6) .^ 2)
    zs = collect(linearspace(0.0, 1.0, 7) .^ 2)
    fens, fes = T4blockx(xs, ys, zs, :a)
    match = true
    for i in eachindex(fes)
        X[1, :] = fens.xyz[fes.conn[i][1], :]
        X[2, :] = fens.xyz[fes.conn[i][2], :]
        X[3, :] = fens.xyz[fes.conn[i][3], :]
        X[4, :] = fens.xyz[fes.conn[i][4], :]
        match = match && (reftetv(X) == tetv(X))
    end
    @test match

    X = Float64[
        0 4 3
        9 2 4
        6 1 7
        0 1 5
    ] # for these points the volume is 10.0
    @test tetv((X')...) == tetv(X)
    @test tetv(X, 1, 2, 3, 4) == tetv(X)
end
end
using .mtetvt1
mtetvt1.test()

module mh8ellip1
using FinEtools
using FinEtools.MeshExportModule: VTK
using Test
function test()
    xradius,
    yradius,
    L,
    H,
    T,
    nL,
    nH,
    nW,
    nT = 0.3, 0.45, 0.9, 0.8, 0.1, 5, 4, 3, 2
    fens, fes = H8elliphole(
        xradius,
        yradius,
        L,
        H,
        T,
        nL,
        nH,
        nW,
        nT,
    )
    # @show count(fens), count(fes)
    @test (count(fens), count(fes)) == (120, 54)
    File = "mesh.vtk"
    VTK.vtkexportmesh(File, fens, fes)
    rm(File)
end
end
using .mh8ellip1
mh8ellip1.test()

module mh8sphere1x1
using FinEtools
using LinearAlgebra
using FinEtools.MeshExportModule: VTK
using Test
function test()
    fens, fes = H8sphere(3.1, 1)
    @test (count(fens), count(fes)) == (65, 32)
    # File = "mesh.vtk"
    # VTK.vtkexportmesh(File, fens, fes)
    # rm(File)

    bfes = meshboundary(fes)

    list = selectelem(fens, bfes, plane = [0.0, -1.0, 0.0, 0.0], thickness = 0.001)
    @test length(list) == 12
    # File = "bmesh.vtk"
    # VTK.vtkexportmesh(File, fens, subset(bfes, list))
    # rm(File)

    list = selectelem(
        fens,
        bfes,
        facing = true,
        direction = function (x)x / norm(x) end,
        dotmin = 0.99,
    )
    @test length(list) == 12
    # File = "bmesh.vtk"
    # VTK.vtkexportmesh(File, fens, subset(bfes, list))
    # rm(File)
end
end
using .mh8sphere1x1
mh8sphere1x1.test()

module mh8sphere1x2
using FinEtools
using LinearAlgebra
using FinEtools.MeshExportModule: VTK
using Test
function test()
    fens, fes = H8sphere(3.1, 1)
    @test (count(fens), count(fes)) == (65, 32)
    meshes = Array{Tuple{FENodeSet,AbstractFESet},1}()
    push!(meshes, (fens, fes))
    fens2, fesa = mergenmeshes(meshes, 0.0001)
    @test (count(fens2), count(fesa[1])) == (65, 32)


    fens1, fes1 = mirrormesh(fens, fes, [0, -1.0, 0], [0.0, 0, 0])
    meshes = Array{Tuple{FENodeSet,AbstractFESet},1}()
    push!(meshes, (fens, fes))
    push!(meshes, (fens1, fes1))
    fens, fesa = mergenmeshes(meshes, 0.0001)
    fes = cat(fesa[1], fesa[2])
    @test (count(fens), count(fes)) == (111, 64)

    # File = "mesh.vtk"
    # VTK.vtkexportmesh(File, fens, fes)
    # rm(File)
    true
end
end
using .mh8sphere1x2
mh8sphere1x2.test()

module mq8ann1
using FinEtools
using FinEtools.MeshExportModule: VTK
using Test
function test()
    rin, rex, nr, nc, Angl, orientation::Symbol =
        100.0, 200.0, 13, 16, pi / 3, :a
    fens, fes = Q8annulus(rin, rex, nr, nc, Angl)
    d = (fens.xyz[:, 1] .^ 2 + fens.xyz[:, 2] .^ 2)
    File = "mesh.vtk"
    # Export of multiple scalar fields
    result = VTK.vtkexportmesh(File, fens, fes; scalars = [("d", d), ("invd", 1 ./ d)])
    @test result == true
    rm(File)

    vectors = [("grad", deepcopy(fens.xyz))]
    filename = "vectors.vtk"
    result = VTK.vtkexportmesh(filename, fens, fes; vectors = vectors)
    @test result == true
    rm(filename)
end
end
using .mq8ann1
mq8ann1.test()

module mexpmeshx1a2
using FinEtools
using FinEtools.MeshExportModule.MESH: write_MESH
using Test
function test()
    rho = 1.21 * 1e-9# mass density
    c = 345.0 * 1000# millimeters per second
    bulk = c^2 * rho
    Lx = 1900.0# length of the box, millimeters
    Ly = 800.0 # length of the box, millimeters

    fens, fes = Q4block(Lx, Ly, 4, 2) # Mesh

    write_MESH("q4-4-2", fens, fes)
    try
        rm("q4-4-2-conn.dat")
    catch
    end
    try
        rm("q4-4-2-xyz.dat")
    catch
    end
    try
        rm("q4-4-2.mesh")
    catch
    end
    true
end
end
using .mexpmeshx1a2
mexpmeshx1a2.test()

module mmmfieldma2
using FinEtools
using Test
import LinearAlgebra: norm
using SymRCM: adjgraph, nodedegrees, symrcm
using SparseArrays
function test()
    N = 10000
    nz = 0.01
    S = sprand(N, N, nz)
    ag = adjgraph(S)
    nd = nodedegrees(ag)
    p = symrcm(ag, nd)
end
end
using .mmmfieldma2
mmmfieldma2.test()


module mmmfieldma3
using FinEtools
using Test
import LinearAlgebra: norm
using SymRCM: adjgraph, nodedegrees, symrcm
using SparseArrays
# using UnicodePlots
function test()
    # The "Bucky" matrix from Matlab documentation
    B = [
        0 1 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
        1 0 1 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
        0 1 0 1 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
        0 0 1 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
        1 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
        1 0 0 0 0 0 1 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
        0 0 0 0 0 1 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
        0 0 0 0 0 0 1 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
        0 0 0 0 0 0 0 1 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
        0 0 0 0 0 1 0 0 1 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
        0 1 0 0 0 0 0 0 0 0 0 1 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
        0 0 0 0 0 0 0 0 0 1 1 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
        0 0 0 0 0 0 0 0 0 0 0 1 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
        0 0 0 0 0 0 0 0 0 0 0 0 1 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
        0 0 0 0 0 0 0 0 0 0 1 0 0 1 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
        0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 1 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 1 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
        0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 1 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
        0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 1 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
        0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 1 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
        0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 1 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
        0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
        0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 1 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
        0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 1 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 1 0 0 1 0 0 0 0 0 0 0 0 0 0
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 1 0 0 0 0 0 0 0 0 0 0 0 0
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 1 0 0 0 0 0 0 0 0 0 0 0
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 1 1 0 0 0 0 0 0 0 0 0
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 1 0 0 0 0 0 0 0 0 0 1 0
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 1 0 0 1 0 0 0 0 0
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 1 0 0 0 0 0 0 0
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 1 0 0 0 0 0 0
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 1 0 0 0 0 0
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 1 0 0 0 0 0 1
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 1
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 1 0 0
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 1 0 1 0
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 1 0 1
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 1 0
    ]
    B = sparse(B)
    # pl = spy(B, canvas = DotCanvas)
    # display(pl)
    ag = adjgraph(B)
    nd = nodedegrees(ag)
    p = symrcm(ag, nd)
    # pl = spy(B[p, p], canvas = DotCanvas)
    # display(pl)
end
end
using .mmmfieldma3
mmmfieldma3.test()

module mAbaqusmnset1
using FinEtools
using FinEtools.MeshExportModule: VTK
using FinEtools.MeshImportModule
using Test
function test()
    output = MeshImportModule.import_ABAQUS("./testsets.inp"; allocationchunk = 11)
    fens, fes = output["fens"], output["fesets"][1]

    nsets = output["nsets"]
    @test nsets["THE-REGION"] == collect(1:358)
    @test nsets["TOP-FACE"] == [
        2,
        33,
        34,
        35,
        36,
        37,
        38,
        39,
        40,
        41,
        42,
        43,
        44,
        45,
        46,
        47,
        48,
        49,
        50,
        51,
        52,
        53,
        54,
        55,
        56,
        57,
        58,
        59,
        60,
        61,
        62,
        168,
        169,
        170,
        171,
        172,
        173,
        174,
        175,
        176,
        177,
        178,
        179,
        180,
        181,
        182,
        183,
        184,
        185,
        186,
        187,
        188,
        189,
        190,
        191,
        192,
        193,
        194,
        195,
        196,
        197,
        198,
        199,
        200,
        201,
        202,
        203,
        204,
        205,
        206,
        207,
        208,
        209,
        210,
        211,
        212,
        213,
        214,
        215,
        216,
        217,
        218,
        219,
        220,
        221,
        222,
        223,
        224,
        225,
        226,
        227,
        228,
        229,
        230,
        231,
        232,
        233,
        234,
        235,
        236,
        237,
        238,
        239,
        240,
        241,
    ]
    @test nsets["BOTTOM-FACE"] == [
        1,
        3,
        4,
        5,
        6,
        7,
        8,
        9,
        10,
        11,
        12,
        13,
        14,
        15,
        16,
        17,
        18,
        19,
        20,
        21,
        22,
        23,
        24,
        25,
        26,
        27,
        28,
        29,
        30,
        31,
        32,
        63,
        64,
        65,
        66,
        67,
        68,
        69,
        70,
        71,
        72,
        73,
        74,
        75,
        76,
        77,
        78,
        79,
        80,
        81,
        82,
        83,
        84,
        85,
        86,
        87,
        88,
        89,
        90,
        91,
        92,
        93,
        94,
        95,
        96,
        97,
        98,
        99,
        100,
        101,
        102,
        103,
        104,
        105,
        106,
        107,
        108,
        109,
        110,
        111,
        112,
        113,
        114,
        115,
        116,
        117,
        118,
        119,
        120,
        121,
        122,
        123,
        124,
        125,
        126,
        127,
        128,
        129,
        130,
        131,
        132,
        133,
        134,
        135,
        136,
    ]
    @test nsets["THE-SIDE"] == [
        1,
        2,
        3,
        4,
        5,
        6,
        7,
        8,
        9,
        10,
        11,
        12,
        13,
        14,
        15,
        16,
        17,
        18,
        19,
        20,
        21,
        22,
        23,
        24,
        25,
        26,
        27,
        28,
        29,
        30,
        31,
        32,
        33,
        34,
        35,
        36,
        37,
        38,
        39,
        40,
        41,
        42,
        43,
        44,
        45,
        46,
        47,
        48,
        49,
        50,
        51,
        52,
        53,
        54,
        55,
        56,
        57,
        58,
        59,
        60,
        61,
        62,
        137,
        138,
        139,
        140,
        141,
        142,
        143,
        144,
        145,
        146,
        147,
        148,
        149,
        150,
        151,
        152,
        153,
        154,
        155,
        156,
        157,
        158,
        159,
        160,
        161,
        162,
        163,
        164,
        165,
        166,
        167,
    ]
    # File = "LE11NAFEMS_H8.vtk"
    # VTK.vtkexportmesh(File, fens, fes)
    # # @async run(`"paraview.exe" $File`)
    # try rm(File) catch end

    # try rm(AE.filename) catch end
end
end
using .mAbaqusmnset1
mAbaqusmnset1.test()

module mef2nf2para9
using FinEtools
using FinEtools.MeshSelectionModule: vselect
using FinEtools.MeshExportModule
using FinEtools.MeshExportModule.VTK
using FinEtools.FEMMBaseModule: field_elem_to_nodal!
using Test
import LinearAlgebra: norm
function test()
    A = 1500.0 * phun("m") # length  of loaded rectangle
    B = 2000.0 * phun("m") # length  of loaded rectangle
    C = 1000.0 * phun("m") # span of the plate

    Meshing = Q4blockx
    # Select how find the mesh should be
    Refinement = 4
    nA, nB, nC = Refinement * 4, Refinement * 6, Refinement * 4
    xs = collect(linearspace(0.0, A, nA + 1))
    ys = collect(linearspace(0.0, B, nB + 1))
    fens, fes = Meshing(xs, ys)

    centroidpc = centroidparametric(fes)
    N = bfun(fes, centroidpc)
    NT = transpose(N)

    ef = ElementalField(zeros(count(fes), 1))
    for i in eachindex(fes)
        c = [k for k in fes.conn[i]]
        centroid = NT * fens.xyz[c, :]
        x, y = centroid
        ef.values[i, :] .= sin(2 * x / A) * cos(6.5 * y / B)
    end
    File = "mef2nf2para9-coarse-ef.vtk"
    VTK.vtkexportmesh(File, fens, fes; scalars = [("ef", ef.values)])
    # @async run(`"paraview.exe" $File`)
    try
        rm(File)
    catch
    end

    femm = FEMMBase(IntegDomain(fes, GaussRule(2, 4)))
    geom = NodalField(fens.xyz)
    nf = NodalField(zeros(count(fens), 1))
    field_elem_to_nodal!(femm, geom, ef, nf)
    File = "mef2nf2para9-coarse-nf.vtk"
    VTK.vtkexportmesh(File, fens, fes; scalars = [("nf", nf.values)])
    try
        rm(File)
    catch
    end
    @test true
end
end
using .mef2nf2para9
mef2nf2para9.test()

module mnf2ef2para9
using FinEtools
using FinEtools.MeshSelectionModule: vselect
using FinEtools.MeshExportModule
using FinEtools.MeshExportModule.VTK
using FinEtools.FEMMBaseModule: field_elem_to_nodal!
using FinEtools.FEMMBaseModule: field_nodal_to_elem!
using Test
import LinearAlgebra: norm
function test()
    A = 1500.0 * phun("m") # length  of loaded rectangle
    B = 2000.0 * phun("m") # length  of loaded rectangle
    C = 1000.0 * phun("m") # span of the plate

    Meshing = Q4blockx
    # Select how find the mesh should be
    Refinement = 4
    nA, nB, nC = Refinement * 4, Refinement * 6, Refinement * 4
    xs = collect(linearspace(0.0, A, nA + 1))
    ys = collect(linearspace(0.0, B, nB + 1))
    fens, fes = Meshing(xs, ys)

    centroidpc = centroidparametric(fes)
    N = bfun(fes, centroidpc)
    NT = transpose(N)

    nf = NodalField(zeros(count(fens), 1))
    for i in eachindex(fens)
        x, y = fens.xyz[i, :]
        nf.values[i, :] .= sin(2 * x / A) * cos(6.5 * y / B)
    end
    File = "mnf2ef2para9-coarse-nf.vtk"
    VTK.vtkexportmesh(File, fens, fes; scalars = [("nf", nf.values)])
    try
        rm(File)
    catch
    end
    # @async run(`"paraview.exe" $File`)
    # try rm(File) catch end

    femm = FEMMBase(IntegDomain(fes, GaussRule(2, 4)))
    geom = NodalField(fens.xyz)
    ef = ElementalField(zeros(count(fes), 1))
    field_nodal_to_elem!(femm, geom, nf, ef)
    File = "mnf2ef2para9-coarse-ef.vtk"
    VTK.vtkexportmesh(File, fens, fes; scalars = [("ef", ef.values)])
    try
        rm(File)
    catch
    end
    @test true
end
end
using .mnf2ef2para9
mnf2ef2para9.test()

module mef2nf2para10
using FinEtools
using FinEtools.MeshSelectionModule: vselect
using FinEtools.MeshExportModule
using FinEtools.MeshExportModule.VTK
using FinEtools.FEMMBaseModule: field_elem_to_nodal!
using Test
import LinearAlgebra: norm
function test()
    A = 1500.0 * phun("m") # length  of loaded rectangle
    B = 2000.0 * phun("m") # length  of loaded rectangle
    C = 1000.0 * phun("m") # span of the plate

    Meshing = T3blockx
    # Select how find the mesh should be
    Refinement = 4
    nA, nB, nC = Refinement * 4, Refinement * 4, Refinement * 1
    xs = collect(linearspace(0.0, A, nA + 1))
    ys = collect(linearspace(0.0, B, nB + 1))
    fens, fes = Meshing(xs, ys)
    # fens,fes = T3toT6(fens,fes)	

    centroidpc = centroidparametric(fes)
    N = bfun(fes, centroidpc)
    NT = transpose(N)

    ef = ElementalField(zeros(count(fes), 1))
    for i in eachindex(fes)
        c = [k for k in fes.conn[i]]
        centroid = NT * fens.xyz[c, :]
        x, y = centroid
        ef.values[i, :] .= sin(2 * x / A) * cos(6.5 * y / B)
    end
    File = "mef2nf2para10-coarse-ef.vtk"
    VTK.vtkexportmesh(File, fens, fes; scalars = [("ef", ef.values)])
    # @async run(`"paraview.exe" $File`)
    try
        rm(File)
    catch
    end

    femm = FEMMBase(IntegDomain(fes, TriRule(3)))
    geom = NodalField(fens.xyz)
    nf = NodalField(zeros(count(fens), 1))
    field_elem_to_nodal!(femm, geom, ef, nf)
    File = "mef2nf2para10-coarse-nf.vtk"
    VTK.vtkexportmesh(File, fens, fes; scalars = [("nf", nf.values)])
    try
        rm(File)
    catch
    end
    @test true
end
end
using .mef2nf2para10
mef2nf2para10.test()

module mnf2ef2para10
using FinEtools
using FinEtools.MeshSelectionModule: vselect
using FinEtools.MeshExportModule
using FinEtools.MeshExportModule.VTK
using FinEtools.FEMMBaseModule: field_elem_to_nodal!
using FinEtools.FEMMBaseModule: field_nodal_to_elem!
using Test
import LinearAlgebra: norm
function test()
    A = 1500.0 * phun("m") # length  of loaded rectangle
    B = 2000.0 * phun("m") # length  of loaded rectangle
    C = 1000.0 * phun("m") # span of the plate

    Meshing = T3blockx
    # Select how find the mesh should be
    Refinement = 4
    nA, nB, nC = Refinement * 4, Refinement * 6, Refinement * 4
    xs = collect(linearspace(0.0, A, nA + 1))
    ys = collect(linearspace(0.0, B, nB + 1))
    fens, fes = Meshing(xs, ys)
    fens, fes = T3toT6(fens, fes)

    centroidpc = centroidparametric(fes)
    N = bfun(fes, centroidpc)
    NT = transpose(N)

    nf = NodalField(zeros(count(fens), 1))
    for i in eachindex(fens)
        x, y = fens.xyz[i, :]
        nf.values[i, :] .= sin(2 * x / A) * cos(6.5 * y / B)
    end
    File = "mnf2ef2para10-coarse-nf.vtk"
    VTK.vtkexportmesh(File, fens, fes; scalars = [("nf", nf.values)])
    try
        rm(File)
    catch
    end
    # @async run(`"paraview.exe" $File`)
    # try rm(File) catch end

    femm = FEMMBase(IntegDomain(fes, TriRule(3)))
    geom = NodalField(fens.xyz)
    ef = ElementalField(zeros(count(fes), 1))
    field_nodal_to_elem!(femm, geom, nf, ef)
    File = "mnf2ef2para10-coarse-ef.vtk"
    VTK.vtkexportmesh(File, fens, fes; scalars = [("ef", ef.values)])
    try
        rm(File)
    catch
    end
    @test true
end
end
using .mnf2ef2para10
mnf2ef2para10.test()


module mef2nf2para11
using FinEtools
using FinEtools.MeshSelectionModule: vselect
using FinEtools.MeshExportModule
using FinEtools.MeshExportModule.VTK
using FinEtools.FEMMBaseModule: field_elem_to_nodal!
using Test
import LinearAlgebra: norm
function test()
    A = 1500.0 * phun("m") # length  of loaded rectangle
    B = 2000.0 * phun("m") # length  of loaded rectangle
    C = 1000.0 * phun("m") # span of the plate

    Meshing = T3blockx
    # Select how find the mesh should be
    Refinement = 4
    nA, nB, nC = Refinement * 2, Refinement * 2, Refinement * 1
    xs = collect(linearspace(0.0, A, nA + 1))
    ys = collect(linearspace(0.0, B, nB + 1))
    fens, fes = Meshing(xs, ys)
    fens, fes = T3toT6(fens, fes)

    centroidpc = centroidparametric(fes)
    N = bfun(fes, centroidpc)
    NT = transpose(N)

    ef = ElementalField(zeros(count(fes), 1))
    for i in eachindex(fes)
        c = [k for k in fes.conn[i]]
        centroid = NT * fens.xyz[c, :]
        x, y = centroid
        ef.values[i, :] .= sin(2 * x / A) * cos(6.5 * y / B)
    end
    File = "mef2nf2para11-coarse-ef.vtk"
    VTK.vtkexportmesh(File, fens, fes; scalars = [("ef", ef.values)])
    try
        rm(File)
    catch
    end
    # @async run(`"paraview.exe" $File`)
    # try rm(File) catch end

    femm = FEMMBase(IntegDomain(fes, TriRule(3)))
    geom = NodalField(fens.xyz)
    nf = NodalField(zeros(count(fens), 1))
    field_elem_to_nodal!(femm, geom, ef, nf)
    File = "mef2nf2para11-coarse-nf.vtk"
    VTK.vtkexportmesh(File, fens, fes; scalars = [("nf", nf.values)])
    try
        rm(File)
    catch
    end
    @test true
end
end
using .mef2nf2para11
mef2nf2para11.test()

module mnf2ef2para11
using FinEtools
using FinEtools.MeshSelectionModule: vselect
using FinEtools.MeshExportModule
using FinEtools.MeshExportModule.VTK
using FinEtools.FEMMBaseModule: field_elem_to_nodal!
using FinEtools.FEMMBaseModule: field_nodal_to_elem!
using Test
import LinearAlgebra: norm
function test()
    A = 1500.0 * phun("m") # length  of loaded rectangle
    B = 2000.0 * phun("m") # length  of loaded rectangle
    C = 1000.0 * phun("m") # span of the plate

    Meshing = T3blockx
    # Select how find the mesh should be
    Refinement = 4
    nA, nB, nC = Refinement * 3, Refinement * 2, Refinement * 3
    xs = collect(linearspace(0.0, A, nA + 1))
    ys = collect(linearspace(0.0, B, nB + 1))
    fens, fes = Meshing(xs, ys)
    fens, fes = T3toT6(fens, fes)

    centroidpc = centroidparametric(fes)
    N = bfun(fes, centroidpc)
    NT = transpose(N)

    nf = NodalField(zeros(count(fens), 1))
    for i in eachindex(fens)
        x, y = fens.xyz[i, :]
        nf.values[i, :] .= sin(2 * x / A) * cos(6.5 * y / B)
    end
    File = "mnf2ef2para11-coarse-nf.vtk"
    VTK.vtkexportmesh(File, fens, fes; scalars = [("nf", nf.values)])
    try
        rm(File)
    catch
    end
    # @async run(`"paraview.exe" $File`)
    # try rm(File) catch end

    femm = FEMMBase(IntegDomain(fes, TriRule(3)))
    geom = NodalField(fens.xyz)
    ef = ElementalField(zeros(count(fes), 1))
    field_nodal_to_elem!(femm, geom, nf, ef)
    File = "mnf2ef2para11-coarse-ef.vtk"
    VTK.vtkexportmesh(File, fens, fes; scalars = [("ef", ef.values)])
    try
        rm(File)
    catch
    end
    @test true
end
end
using .mnf2ef2para11
mnf2ef2para11.test()

module mef2nf2para12
using FinEtools
using FinEtools.MeshSelectionModule: vselect
using FinEtools.MeshExportModule
using FinEtools.MeshExportModule.VTK
using FinEtools.FEMMBaseModule: field_elem_to_nodal!
using FinEtools.FEMMBaseModule: field_nodal_to_elem!
using Test
import LinearAlgebra: norm
function test()
    A = 1500.0 * phun("m") # length  of loaded rectangle
    B = 2000.0 * phun("m") # length  of loaded rectangle
    C = 1000.0 * phun("m") # span of the plate

    Meshing = T3blockx
    # Select how find the mesh should be
    Refinement = 4
    nA, nB, nC = Refinement * 5, Refinement * 5, Refinement * 1
    xs = collect(linearspace(0.0, A, nA + 1))
    ys = collect(linearspace(0.0, B, nB + 1))
    fens, fes = Meshing(xs, ys)
    # fens,fes = T3toT6(fens,fes)	

    centroidpc = centroidparametric(fes)
    N = bfun(fes, centroidpc)
    NT = transpose(N)

    ef = ElementalField(zeros(count(fes), 1))
    for i in eachindex(fes)
        c = [k for k in fes.conn[i]]
        centroid = NT * fens.xyz[c, :]
        x, y = centroid
        ef.values[i, :] .= sin(2 * x / A) * cos(6.5 * y / B)
    end
    File = "mef2nf2para12-coarse-ef.vtk"
    VTK.vtkexportmesh(File, fens, fes; scalars = [("ef", ef.values)])
    try
        rm(File)
    catch
    end
    # @async run(`"paraview.exe" $File`)
    # try rm(File) catch end

    femm = FEMMBase(IntegDomain(fes, TriRule(3)))
    geom = NodalField(fens.xyz)
    nf = NodalField(zeros(count(fens), 1))
    field_elem_to_nodal!(femm, geom, ef, nf)

    for sm in 1:3
        field_nodal_to_elem!(femm, geom, nf, ef)
        field_elem_to_nodal!(femm, geom, ef, nf)
    end

    File = "mef2nf2para12-coarse-nf-smoothed.vtk"
    VTK.vtkexportmesh(File, fens, fes; scalars = [("nf", nf.values)])
    try
        rm(File)
    catch
    end
    File = "mef2nf2para12-coarse-ef-smoothed.vtk"
    VTK.vtkexportmesh(File, fens, fes; scalars = [("ef", ef.values)])
    try
        rm(File)
    catch
    end
    @test true
end
end
using .mef2nf2para12
mef2nf2para12.test()

module mef2nf2maxa9
using FinEtools
using FinEtools.MeshSelectionModule: vselect
using FinEtools.MeshExportModule
using FinEtools.MeshExportModule.VTK
using FinEtools.FEMMBaseModule: field_elem_to_nodal!
using Test
import LinearAlgebra: norm
function test()
    A = 1500.0 * phun("m") # length  of loaded rectangle
    B = 2000.0 * phun("m") # length  of loaded rectangle
    C = 1000.0 * phun("m") # span of the plate

    Meshing = Q4blockx
    # Select how find the mesh should be
    Refinement = 4
    nA, nB, nC = Refinement * 4, Refinement * 6, Refinement * 4
    xs = collect(linearspace(0.0, A, nA + 1))
    ys = collect(linearspace(0.0, B, nB + 1))
    fens, fes = Meshing(xs, ys)

    centroidpc = centroidparametric(fes)
    N = bfun(fes, centroidpc)
    NT = transpose(N)

    ef = ElementalField(zeros(count(fes), 1))
    for i in eachindex(fes)
        c = [k for k in fes.conn[i]]
        centroid = NT * fens.xyz[c, :]
        x, y = centroid
        ef.values[i, :] .= sin(2 * x / A) * cos(6.5 * y / B)
    end
    File = "mef2nf2maxa9-coarse-ef.vtk"
    VTK.vtkexportmesh(File, fens, fes; scalars = [("ef", ef.values)])
    # @async run(`"paraview.exe" $File`)
    try
        rm(File)
    catch
    end

    femm = FEMMBase(IntegDomain(fes, GaussRule(2, 4)))
    geom = NodalField(fens.xyz)
    nf = NodalField(zeros(count(fens), 1))
    field_elem_to_nodal!(femm, geom, ef, nf; kind = :max)
    File = "mef2nf2maxa9-coarse-nf.vtk"
    VTK.vtkexportmesh(File, fens, fes; scalars = [("nf", nf.values)])
    try
        rm(File)
    catch
    end
    @test true
end
end
using .mef2nf2maxa9
mef2nf2maxa9.test()

module mnf2ef2maxa9
using FinEtools
using FinEtools.MeshSelectionModule: vselect
using FinEtools.MeshExportModule
using FinEtools.MeshExportModule.VTK
using FinEtools.FEMMBaseModule: field_elem_to_nodal!
using FinEtools.FEMMBaseModule: field_nodal_to_elem!
using Test
import LinearAlgebra: norm
function test()
    A = 1500.0 * phun("m") # length  of loaded rectangle
    B = 2000.0 * phun("m") # length  of loaded rectangle
    C = 1000.0 * phun("m") # span of the plate

    Meshing = Q4blockx
    # Select how find the mesh should be
    Refinement = 4
    nA, nB, nC = Refinement * 4, Refinement * 6, Refinement * 4
    xs = collect(linearspace(0.0, A, nA + 1))
    ys = collect(linearspace(0.0, B, nB + 1))
    fens, fes = Meshing(xs, ys)

    centroidpc = centroidparametric(fes)
    N = bfun(fes, centroidpc)
    NT = transpose(N)

    nf = NodalField(zeros(count(fens), 1))
    for i in eachindex(fens)
        x, y = fens.xyz[i, :]
        nf.values[i, :] .= sin(2 * x / A) * cos(6.5 * y / B)
    end
    File = "mnf2ef2maxa9-coarse-nf.vtk"
    VTK.vtkexportmesh(File, fens, fes; scalars = [("nf", nf.values)])
    try
        rm(File)
    catch
    end
    # @async run(`"paraview.exe" $File`)
    # try rm(File) catch end

    femm = FEMMBase(IntegDomain(fes, GaussRule(2, 4)))
    geom = NodalField(fens.xyz)
    ef = ElementalField(zeros(count(fes), 1))
    field_nodal_to_elem!(femm, geom, nf, ef; kind = :max)
    File = "mnf2ef2maxa9-coarse-ef.vtk"
    VTK.vtkexportmesh(File, fens, fes; scalars = [("ef", ef.values)])
    try
        rm(File)
    catch
    end
    @test true
end
end
using .mnf2ef2maxa9
mnf2ef2maxa9.test()

module mef2nf2maxa10
using FinEtools
using FinEtools.MeshSelectionModule: vselect
using FinEtools.MeshExportModule
using FinEtools.MeshExportModule.VTK
using FinEtools.FEMMBaseModule: field_elem_to_nodal!
using Test
import LinearAlgebra: norm
function test()
    A = 1500.0 * phun("m") # length  of loaded rectangle
    B = 2000.0 * phun("m") # length  of loaded rectangle
    C = 1000.0 * phun("m") # span of the plate

    Meshing = T3blockx
    # Select how find the mesh should be
    Refinement = 4
    nA, nB, nC = Refinement * 4, Refinement * 4, Refinement * 1
    xs = collect(linearspace(0.0, A, nA + 1))
    ys = collect(linearspace(0.0, B, nB + 1))
    fens, fes = Meshing(xs, ys)
    # fens,fes = T3toT6(fens,fes)   

    centroidpc = centroidparametric(fes)
    N = bfun(fes, centroidpc)
    NT = transpose(N)

    ef = ElementalField(zeros(count(fes), 1))
    for i in eachindex(fes)
        c = [k for k in fes.conn[i]]
        centroid = NT * fens.xyz[c, :]
        x, y = centroid
        ef.values[i, :] .= sin(2 * x / A) * cos(6.5 * y / B)
    end
    File = "mef2nf2maxa10-coarse-ef.vtk"
    VTK.vtkexportmesh(File, fens, fes; scalars = [("ef", ef.values)])
    # @async run(`"paraview.exe" $File`)
    try
        rm(File)
    catch
    end

    femm = FEMMBase(IntegDomain(fes, TriRule(3)))
    geom = NodalField(fens.xyz)
    nf = NodalField(zeros(count(fens), 1))
    field_elem_to_nodal!(femm, geom, ef, nf; kind = :max)
    File = "mef2nf2maxa10-coarse-nf.vtk"
    VTK.vtkexportmesh(File, fens, fes; scalars = [("nf", nf.values)])
    try
        rm(File)
    catch
    end
    @test true
end
end
using .mef2nf2maxa10
mef2nf2maxa10.test()

module mnf2ef2maxa10
using FinEtools
using FinEtools.MeshSelectionModule: vselect
using FinEtools.MeshExportModule
using FinEtools.MeshExportModule.VTK
using FinEtools.FEMMBaseModule: field_elem_to_nodal!
using FinEtools.FEMMBaseModule: field_nodal_to_elem!
using Test
import LinearAlgebra: norm
function test()
    A = 1500.0 * phun("m") # length  of loaded rectangle
    B = 2000.0 * phun("m") # length  of loaded rectangle
    C = 1000.0 * phun("m") # span of the plate

    Meshing = T3blockx
    # Select how find the mesh should be
    Refinement = 4
    nA, nB, nC = Refinement * 4, Refinement * 6, Refinement * 4
    xs = collect(linearspace(0.0, A, nA + 1))
    ys = collect(linearspace(0.0, B, nB + 1))
    fens, fes = Meshing(xs, ys)
    fens, fes = T3toT6(fens, fes)

    centroidpc = centroidparametric(fes)
    N = bfun(fes, centroidpc)
    NT = transpose(N)

    nf = NodalField(zeros(count(fens), 1))
    for i in eachindex(fens)
        x, y = fens.xyz[i, :]
        nf.values[i, :] .= sin(2 * x / A) * cos(6.5 * y / B)
    end
    File = "mnf2ef2maxa10-coarse-nf.vtk"
    VTK.vtkexportmesh(File, fens, fes; scalars = [("nf", nf.values)])
    try
        rm(File)
    catch
    end
    # @async run(`"paraview.exe" $File`)
    # try rm(File) catch end

    femm = FEMMBase(IntegDomain(fes, TriRule(3)))
    geom = NodalField(fens.xyz)
    ef = ElementalField(zeros(count(fes), 1))
    field_nodal_to_elem!(femm, geom, nf, ef; kind = :max)
    File = "mnf2ef2maxa10-coarse-ef.vtk"
    VTK.vtkexportmesh(File, fens, fes; scalars = [("ef", ef.values)])
    try
        rm(File)
    catch
    end
    @test true
end
end
using .mnf2ef2maxa10
mnf2ef2maxa10.test()


module mef2nf2maxa11
using FinEtools
using FinEtools.MeshSelectionModule: vselect
using FinEtools.MeshExportModule
using FinEtools.MeshExportModule.VTK
using FinEtools.FEMMBaseModule: field_elem_to_nodal!
using Test
import LinearAlgebra: norm
function test()
    A = 1500.0 * phun("m") # length  of loaded rectangle
    B = 2000.0 * phun("m") # length  of loaded rectangle
    C = 1000.0 * phun("m") # span of the plate

    Meshing = T3blockx
    # Select how find the mesh should be
    Refinement = 4
    nA, nB, nC = Refinement * 2, Refinement * 2, Refinement * 1
    xs = collect(linearspace(0.0, A, nA + 1))
    ys = collect(linearspace(0.0, B, nB + 1))
    fens, fes = Meshing(xs, ys)
    fens, fes = T3toT6(fens, fes)

    centroidpc = centroidparametric(fes)
    N = bfun(fes, centroidpc)
    NT = transpose(N)

    ef = ElementalField(zeros(count(fes), 1))
    for i in eachindex(fes)
        c = [k for k in fes.conn[i]]
        centroid = NT * fens.xyz[c, :]
        x, y = centroid
        ef.values[i, :] .= sin(2 * x / A) * cos(6.5 * y / B)
    end
    File = "mef2nf2maxa11-coarse-ef.vtk"
    VTK.vtkexportmesh(File, fens, fes; scalars = [("ef", ef.values)])
    try
        rm(File)
    catch
    end
    # @async run(`"paraview.exe" $File`)
    # try rm(File) catch end

    femm = FEMMBase(IntegDomain(fes, TriRule(3)))
    geom = NodalField(fens.xyz)
    nf = NodalField(zeros(count(fens), 1))
    field_elem_to_nodal!(femm, geom, ef, nf; kind = :max)
    File = "mef2nf2maxa11-coarse-nf.vtk"
    VTK.vtkexportmesh(File, fens, fes; scalars = [("nf", nf.values)])
    try
        rm(File)
    catch
    end
    @test true
end
end
using .mef2nf2maxa11
mef2nf2maxa11.test()

module mnf2ef2maxa11
using FinEtools
using FinEtools.MeshSelectionModule: vselect
using FinEtools.MeshExportModule
using FinEtools.MeshExportModule.VTK
using FinEtools.FEMMBaseModule: field_elem_to_nodal!
using FinEtools.FEMMBaseModule: field_nodal_to_elem!
using Test
import LinearAlgebra: norm
function test()
    A = 1500.0 * phun("m") # length  of loaded rectangle
    B = 2000.0 * phun("m") # length  of loaded rectangle
    C = 1000.0 * phun("m") # span of the plate

    Meshing = T3blockx
    # Select how find the mesh should be
    Refinement = 4
    nA, nB, nC = Refinement * 3, Refinement * 2, Refinement * 3
    xs = collect(linearspace(0.0, A, nA + 1))
    ys = collect(linearspace(0.0, B, nB + 1))
    fens, fes = Meshing(xs, ys)
    fens, fes = T3toT6(fens, fes)

    centroidpc = centroidparametric(fes)
    N = bfun(fes, centroidpc)
    NT = transpose(N)

    nf = NodalField(zeros(count(fens), 1))
    for i in eachindex(fens)
        x, y = fens.xyz[i, :]
        nf.values[i, :] .= sin(2 * x / A) * cos(6.5 * y / B)
    end
    File = "mnf2ef2maxa11-coarse-nf.vtk"
    VTK.vtkexportmesh(File, fens, fes; scalars = [("nf", nf.values)])
    try
        rm(File)
    catch
    end
    # @async run(`"paraview.exe" $File`)
    # try rm(File) catch end

    femm = FEMMBase(IntegDomain(fes, TriRule(3)))
    geom = NodalField(fens.xyz)
    ef = ElementalField(zeros(count(fes), 1))
    field_nodal_to_elem!(femm, geom, nf, ef; kind = :max)
    File = "mnf2ef2maxa11-coarse-ef.vtk"
    VTK.vtkexportmesh(File, fens, fes; scalars = [("ef", ef.values)])
    try
        rm(File)
    catch
    end
    @test true
end
end
using .mnf2ef2maxa11
mnf2ef2maxa11.test()

module mef2nf2maxa12
using FinEtools
using FinEtools.MeshSelectionModule: vselect
using FinEtools.MeshExportModule
using FinEtools.MeshExportModule.VTK
using FinEtools.FEMMBaseModule: field_elem_to_nodal!
using FinEtools.FEMMBaseModule: field_nodal_to_elem!
using Test
import LinearAlgebra: norm
function test()
    A = 1500.0 * phun("m") # length  of loaded rectangle
    B = 2000.0 * phun("m") # length  of loaded rectangle
    C = 1000.0 * phun("m") # span of the plate

    Meshing = T3blockx
    # Select how find the mesh should be
    Refinement = 4
    nA, nB, nC = Refinement * 5, Refinement * 5, Refinement * 1
    xs = collect(linearspace(0.0, A, nA + 1))
    ys = collect(linearspace(0.0, B, nB + 1))
    fens, fes = Meshing(xs, ys)
    # fens,fes = T3toT6(fens,fes)   

    centroidpc = centroidparametric(fes)
    N = bfun(fes, centroidpc)
    NT = transpose(N)

    ef = ElementalField(zeros(count(fes), 1))
    for i in eachindex(fes)
        c = [k for k in fes.conn[i]]
        centroid = NT * fens.xyz[c, :]
        x, y = centroid
        ef.values[i, :] .= sin(2 * x / A) * cos(6.5 * y / B)
    end
    File = "mef2nf2maxa12-coarse-ef.vtk"
    VTK.vtkexportmesh(File, fens, fes; scalars = [("ef", ef.values)])
    try
        rm(File)
    catch
    end
    # @async run(`"paraview.exe" $File`)
    # try rm(File) catch end

    femm = FEMMBase(IntegDomain(fes, TriRule(3)))
    geom = NodalField(fens.xyz)
    nf = NodalField(zeros(count(fens), 1))
    field_elem_to_nodal!(femm, geom, ef, nf; kind = :max)

    for sm in 1:3
        field_nodal_to_elem!(femm, geom, nf, ef; kind = :max)
        field_elem_to_nodal!(femm, geom, ef, nf; kind = :max)
    end

    File = "mef2nf2maxa12-coarse-nf-smoothed.vtk"
    VTK.vtkexportmesh(File, fens, fes; scalars = [("nf", nf.values)])
    try
        rm(File)
    catch
    end
    File = "mef2nf2maxa12-coarse-ef-smoothed.vtk"
    VTK.vtkexportmesh(File, fens, fes; scalars = [("ef", ef.values)])
    try
        rm(File)
    catch
    end
    @test true
end
end
using .mef2nf2maxa12
mef2nf2maxa12.test()

module mimptexptm1
using FinEtools
using FinEtools.MeshImportModule
using FinEtools.MeshExportModule.MESH
using Test
function test()
    output = MeshImportModule.import_NASTRAN(
        dirname(@__FILE__) * "/" * "cylinder.nas";
        allocationchunk = 13,
        expectfixedformat = true,
    )
    @test count(output["fens"]) == 1406
    @test count(output["fesets"][1]) == 829
    # show(fes.conn[count(fes), :])
    fens = output["fens"]
    fes = output["fesets"][1]

    File = "cylinder"
    MESH.write_MESH(File, fens, fes)

    output = MeshImportModule.import_MESH(File)
    @test count(output["fens"]) == 1406
    @test count(output["fesets"][1]) == 829
    try
        rm(File)
    catch
    end
    try
        rm("cylinder-conn.dat")
    catch
    end
    try
        rm("cylinder-xyz.dat")
    catch
    end
    try
        rm("cylinder.mesh")
    catch
    end
    true
end
end
using .mimptexptm1
mimptexptm1.test()

module mAbaqimportm1
using FinEtools
using FinEtools.MeshExportModule: VTK
using FinEtools.MeshImportModule
using Test
function test()
    output = MeshImportModule.import_ABAQUS(
        "./problem_LE10NAFEMS_MST10.inp";
        allocationchunk = 11,
    )

    fens, fes = output["fens"], output["fesets"][1]
    @test count(fes) == 60
    @test length(output["nsets"]["L3"]) == 7
end
end
using .mAbaqimportm1
mAbaqimportm1.test()


module mh5t01
using FinEtools
using FinEtools.MeshExportModule: H5MESH
using FinEtools.MeshImportModule: import_H5MESH
using LinearAlgebra: norm
using Test
function test()
    rho = 1.21 * 1e-9# mass density
    c = 345.0 * 1000# millimeters per second
    bulk = c^2 * rho
    Lx = 1900.0# length of the box, millimeters
    Ly = 800.0 # length of the box, millimeters

    fens, fes = Q4block(Lx, Ly, 40, 20) # 0Mesh
    setlabel!(fes, 3)

    H5MESH.write_H5MESH("block", fens, fes)

    output = import_H5MESH("block.h5mesh")
    #@show output
    @test count(output["fens"]) == count(fens)
    @test count(output["fesets"][1]) == count(fes)
    @test norm(output["fesets"][1].label - fes.label) == 0

    try
        rm("block.h5mesh")
    catch
    end
    true
end
end
using .mh5t01
mh5t01.test()


module mt3ext01
using FinEtools
using FinEtools.MeshTetrahedronModule: T4extrudeT3
using FinEtools.MeshExportModule: VTK
using LinearAlgebra: norm
using Test
function test()

    Lx = 1900.0# length of the box, millimeters
    Ly = 800.0 # length of the box, millimeters
    Lz = 1200.0 # length of the box, millimeters
    nlayers = 3

    surfens, surfes = T3block(Lx, Ly, 4, 2)
    setlabel!(surfes, 3)

    fens, fes = T4extrudeT3(
        surfens,
        surfes,
        nlayers,
        (X, layer) -> [X[1], X[2], layer * Lz / nlayers],
    )
    @test count(fens) == (nlayers + 1) * count(surfens)
    # @show fens, fes
    @test count(fes) == 3 * nlayers * count(surfes)
    # @show count(fens)
    goodnodes = true
    for i in eachindex(fes)
        for j in 1:4
            goodnodes = goodnodes && (1 <= fes.conn[i][j] <= count(fens))
        end
    end
    @test goodnodes == true
    File = "mesh.vtk"
    VTK.vtkexportmesh(File, fens, fes)
    try
        rm(File)
    catch
    end
    # @async run(`"paraview.exe" $File`)

    bfes = meshboundary(fes)
    File = "boundary_mesh.vtk"
    VTK.vtkexportmesh(File, fens, bfes)
    try
        rm(File)
    catch
    end
    # @async run(`"paraview.exe" $File`)

    true
end
end
using .mt3ext01
mt3ext01.test()


module mt3ext02
using FinEtools
using FinEtools.MeshTetrahedronModule: T4extrudeT3
using FinEtools.MeshExportModule: VTK
using FinEtools.MeshImportModule
using LinearAlgebra: norm
using Test
function test()

    radius = 30
    output = MeshImportModule.import_ABAQUS("./weird_sphere.inp"; allocationchunk = 11)
    surfens, surfes = output["fens"], output["fesets"][1]

    nlayers = 5

    function radially(X, layer)
        n = X ./ norm(X)
        X .= n .* (1.0 + layer / nlayers / 5) * radius
    end
    fens, fes = T4extrudeT3(surfens, surfes, nlayers, radially)
    @test count(fens) == (nlayers + 1) * count(surfens)

    @test count(fes) == 3 * nlayers * count(surfes)
    # @show count(fens)
    goodnodes = true
    for i in eachindex(fes)
        for j in 1:4
            goodnodes = goodnodes && (1 <= fes.conn[i][j] <= count(fens))
        end
    end
    @test goodnodes == true


    surbfes = meshboundary(surfes)

    bfes = meshboundary(fes)

    @test count(bfes) == 2 * count(surfes) + 2 * nlayers * count(surbfes)

    File = "mesh.vtk"
    VTK.vtkexportmesh(File, fens, fes)
    try
        rm(File)
    catch
    end
    # @async run(`"paraview.exe" $File`)

    File = "boundary_mesh.vtk"
    VTK.vtkexportmesh(File, fens, bfes)
    try
        rm(File)
    catch
    end
    # @async run(`"paraview.exe" $File`)


    true
end
end
using .mt3ext02
mt3ext02.test()


module mselecte1a1
using FinEtools
using Test
function test()
    W = 10.1
    L = 12.3
    t = 32.0
    nl = 2
    nt = 1
    nw = 1
    ref = 1
    tolerance = 0.1

    fens, fes = H8block(L, W, t, nl * ref, nw * ref, nt * ref)

    box = initbox!([], [L, W, t] .+ tolerance)
    cornerlist = selectelem(fens, fes, overlappingbox = box, inflate = 0.0)
    @test length(cornerlist) == 0

    box = initbox!([], [L, W, t] .+ tolerance)
    cornerlist = selectelem(fens, fes, overlappingbox = box)
    @test length(cornerlist) == 0

    box = initbox!([], [L, W, t] .+ tolerance)
    cornerlist = selectelem(fens, fes, overlappingbox = box, inflate = 0.1)
    @test length(cornerlist) == 1

    box = initbox!([], [L, W, t] .+ tolerance)
    cornerlist = selectelem(fens, fes, overlappingbox = box, inflate = 0.09)
    @test length(cornerlist) == 0

    # # # Traction on the opposite edge
    # boundaryfes  =   meshboundary(fes);
    # Toplist   = selectelem(fens,boundaryfes, facing = true,
    # direction = [-1.0 0.0 0.0], dotmin = 0.999);
    # # display(Toplist)
    # @test length(Toplist) == 49
    # Toplist   = selectelem(fens,boundaryfes, box = [L L -Inf Inf -Inf Inf],
    # allin = true, inflate  =  tolerance);
    # # display(Toplist)
    # @test length(Toplist) == 49
    # Toplist   = selectelem(fens,boundaryfes, box = [L L -Inf Inf -Inf Inf],
    # allin = false, inflate  =  tolerance);
    # # display(Toplist)
    # @test length(Toplist) == 77

    # File =  "b.vtk"
    # vtkexportmesh(File, subset(boundaryfes, Toplist).conn, fens.xyz,
    #     FinEtools.MeshExportModule.Q4)
    # File =  "a.vtk"
    # vtkexportmesh(File, fes.conn, fens.xyz,  FinEtools.MeshExportModule.H8)
    # @async run(`"paraview.exe" $File`)
end
end
using .mselecte1a1
mselecte1a1.test()


module mmAbaqusimport7
using FinEtools
using FinEtools.MeshImportModule
using FinEtools.MeshExportModule: VTKWrite
using FinEtools.MeshExportModule: VTK
using Test
function test()

    output = MeshImportModule.import_ABAQUS(
        "./barrelvault_s3r_fineirreg.inp";
        allocationchunk = 11,
    )
    fens, fes = output["fens"], output["fesets"][1]
    @test boundingbox(fens.xyz) == [229.813339, 300.0, 0.0, 192.836288, 0.0, 300.0]
    @test output["nsets"]["APEX"] == [
        1,
        2,
        3,
        4,
        5,
        6,
        7,
        8,
        9,
        10,
        11,
        12,
        13,
        14,
        15,
        16,
        17,
        18,
        19,
        20,
        21,
        22,
        23,
        24,
        25,
    ]
    File = "junk.vtk"
    VTK.vtkexportmesh(File, fens, fes)
    # @async run(`"paraview.exe" $File`)
    # try rm(File) catch end
end
end
using .mmAbaqusimport7
mmAbaqusimport7.test()

module mouter001
using FinEtools
using FinEtools.MeshImportModule
using FinEtools.MeshExportModule: VTK
using Statistics
using LinearAlgebra
using Test

function test()
    output = MeshImportModule.import_ABAQUS("./sphere_with_cavity.inp")
    fens, fes = output["fens"], output["fesets"][1]
    bdry_fes = meshboundary(fes)
    outerfes = outer_surface_of_solid(fens, bdry_fes)
    # File = "mesh.vtk"
    # result =  VTK.vtkexportmesh(File, fens, outerfes)
    @test count(outerfes) == 670

    output = MeshImportModule.import_ABAQUS("./cones_with_cavities.inp")
    fens, fes = output["fens"], output["fesets"][1]
    bdry_fes = meshboundary(fes)
    outerfes = outer_surface_of_solid(fens, bdry_fes)
    # File = "mesh.vtk"
    # result =  VTK.vtkexportmesh(File, fens, outerfes)
    @test count(outerfes) == 914

    output = MeshImportModule.import_ABAQUS("./nonconvex_with_cavities.inp")
    fens, fes = output["fens"], output["fesets"][1]
    bdry_fes = meshboundary(fes)
    outerfes = outer_surface_of_solid(fens, bdry_fes)
    # File = "mesh.vtk"
    # result =  VTK.vtkexportmesh(File, fens, outerfes)
    @test count(outerfes) == 852
    true
end
test()

end


module mAbaqusmelset1
using FinEtools
using FinEtools.MeshExportModule: VTK
using FinEtools.MeshImportModule
using Test
function test()

    output =
        MeshImportModule.import_ABAQUS("./cyl-w-voids-two-mats.inp"; allocationchunk = 11)
    fens, fes = output["fens"], output["fesets"][1]
    @test length(output["elsets"]["SET-1"]) == 266
    @test length(output["elsets"]["SET-2"]) == 260
    File = "cyl-w-voids-two-mats-elset-1.vtk"
    VTK.vtkexportmesh(File, fens, subset(fes, output["elsets"]["SET-1"]))
    File = "cyl-w-voids-two-mats-elset-2.vtk"
    VTK.vtkexportmesh(File, fens, subset(fes, output["elsets"]["SET-2"]))
    # @async run(`"paraview.exe" $File`)
    # try rm(File) catch end

end
end
using .mAbaqusmelset1
mAbaqusmelset1.test()


module mT4quartercylnm1
using FinEtools
using FinEtools.MeshExportModule: VTK
using Test
function test()
    Radius, Length, nperradius, nL = 3.0, 5.0, 3, 5
    fens, fes = T4quartercyln(Radius, Length, nperradius, nL)
    @test count(fes) == 270
    File = "mT4quartercylnm1-mesh.vtk"
    VTK.vtkexportmesh(File, fens, fes)
    rm(File)
    # @async run(`"paraview.exe" $File`)
    bfes = meshboundary(fes)
    @test count(bfes) == 156
    File = "mT4quartercylnm1-bmesh.vtk"
    VTK.vtkexportmesh(File, fens, fes)
    rm(File)
    # @async run(`"paraview.exe" $File`)
end
test()
end


module mT4cylinderm1
using FinEtools
using FinEtools.MeshExportModule: VTK
using Test
function test()
    Radius, Length, nperradius, nL = 3.0, 5.0, 3, 5
    fens, fes = T4quartercyln(Radius, Length, nperradius, nL)
    fens1, fes1 = mirrormesh(
        fens,
        fes,
        [0.0, -1.0, 0.0],
        [0.0, 0.0, 0.0],
        renumb = (c) -> c[[1, 3, 2, 4]],
    )
    meshes = Array{Tuple{FENodeSet,AbstractFESet},1}()
    push!(meshes, (fens, fes))
    push!(meshes, (fens1, fes1))
    fens, fesa = mergenmeshes(meshes, 0.0001)
    fes = cat(fesa[1], fesa[2])
    fens1, fes1 = mirrormesh(
        fens,
        fes,
        [-1.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
        renumb = (c) -> c[[1, 3, 2, 4]],
    )
    meshes = Array{Tuple{FENodeSet,AbstractFESet},1}()
    push!(meshes, (fens, fes))
    push!(meshes, (fens1, fes1))
    fens, fesa = mergenmeshes(meshes, 0.0001)
    fes = cat(fesa[1], fesa[2])
    @test count(fes) == 270 * 4
    File = "mT4cylinderm1-mesh.vtk"
    VTK.vtkexportmesh(File, fens, fes)
    rm(File)
    # @async run(`"paraview.exe" $File`)
    bfes = meshboundary(fes)

    File = "mT4cylinderm1-bmesh.vtk"
    VTK.vtkexportmesh(File, fens, fes)
    rm(File)
    # @async run(`"paraview.exe" $File`)
end
test()
end



module mT4cylinderm2a
using FinEtools
using FinEtools.MeshExportModule: VTK
using Test
function __volume(fens, fes)
    geom = NodalField(fens.xyz)
    femm = FEMMBase(IntegDomain(fes, TetRule(1)))
    return integratefunction(femm, geom, (x) -> 1.0)
end
function test()
    Radius, Length, nperradius, nL = 3.0, 5.0, 3, 5
    fens, fes = T4quartercyln(Radius, Length, nperradius, nL; orientation = :a)
    fens1, fes1 = mirrormesh(
        fens,
        fes,
        [0.0, -1.0, 0.0],
        [0.0, 0.0, 0.0],
        renumb = (c) -> c[[1, 3, 2, 4]],
    )
    meshes = Array{Tuple{FENodeSet,AbstractFESet},1}()
    push!(meshes, (fens, fes))
    push!(meshes, (fens1, fes1))
    fens, fesa = mergenmeshes(meshes, 0.0001)
    fes = cat(fesa[1], fesa[2])
    fens1, fes1 = mirrormesh(
        fens,
        fes,
        [-1.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
        renumb = (c) -> c[[1, 3, 2, 4]],
    )
    meshes = Array{Tuple{FENodeSet,AbstractFESet},1}()
    push!(meshes, (fens, fes))
    push!(meshes, (fens1, fes1))
    fens, fesa = mergenmeshes(meshes, 0.0001)
    fes = cat(fesa[1], fesa[2])
    @test count(fes) == 270 * 4
    File = "mT4cylinderm2a-mesh.vtk"
    VTK.vtkexportmesh(File, fens, fes)
    rm(File)
    # @async run(`"paraview.exe" $File`)
    bfes = meshboundary(fes)

    File = "mT4cylinderm2a-bmesh.vtk"
    VTK.vtkexportmesh(File, fens, fes)
    rm(File)
    # @async run(`"paraview.exe" $File`)
    V = pi * Radius^2 * Length
    @test (V - __volume(fens, fes)) / V < 0.02
end
test()
end

module mT4cylinderm2b
using FinEtools
using FinEtools.MeshExportModule: VTK
using Test
function __volume(fens, fes)
    geom = NodalField(fens.xyz)
    femm = FEMMBase(IntegDomain(fes, TetRule(1)))
    return integratefunction(femm, geom, (x) -> 1.0)
end
function test()
    Radius, Length, nperradius, nL = 3.0, 5.0, 3, 5
    fens, fes = T4quartercyln(Radius, Length, nperradius, nL; orientation = :b)
    fens1, fes1 = mirrormesh(
        fens,
        fes,
        [0.0, -1.0, 0.0],
        [0.0, 0.0, 0.0],
        renumb = (c) -> c[[1, 3, 2, 4]],
    )
    meshes = Array{Tuple{FENodeSet,AbstractFESet},1}()
    push!(meshes, (fens, fes))
    push!(meshes, (fens1, fes1))
    fens, fesa = mergenmeshes(meshes, 0.0001)
    fes = cat(fesa[1], fesa[2])
    fens1, fes1 = mirrormesh(
        fens,
        fes,
        [-1.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
        renumb = (c) -> c[[1, 3, 2, 4]],
    )
    meshes = Array{Tuple{FENodeSet,AbstractFESet},1}()
    push!(meshes, (fens, fes))
    push!(meshes, (fens1, fes1))
    fens, fesa = mergenmeshes(meshes, 0.0001)
    fes = cat(fesa[1], fesa[2])
    @test count(fes) == 270 * 4
    File = "mT4cylinderm2b-mesh.vtk"
    VTK.vtkexportmesh(File, fens, fes)
    rm(File)
    # @async run(`"paraview.exe" $File`)
    bfes = meshboundary(fes)

    File = "mT4cylinderm2b-bmesh.vtk"
    VTK.vtkexportmesh(File, fens, fes)
    rm(File)
    # @async run(`"paraview.exe" $File`)
    V = pi * Radius^2 * Length
    @test (V - __volume(fens, fes)) / V < 0.02
end
test()
end

module mT4cylinderm2ca
using FinEtools
using FinEtools.MeshExportModule: VTK
using Test
function __volume(fens, fes)
    geom = NodalField(fens.xyz)
    femm = FEMMBase(IntegDomain(fes, TetRule(1)))
    return integratefunction(femm, geom, (x) -> 1.0)
end
function test()
    Radius, Length, nperradius, nL = 3.0, 5.0, 3, 5
    fens, fes = T4quartercyln(Radius, Length, nperradius, nL; orientation = :ca)
    fens1, fes1 = mirrormesh(
        fens,
        fes,
        [0.0, -1.0, 0.0],
        [0.0, 0.0, 0.0],
        renumb = (c) -> c[[1, 3, 2, 4]],
    )
    meshes = Array{Tuple{FENodeSet,AbstractFESet},1}()
    push!(meshes, (fens, fes))
    push!(meshes, (fens1, fes1))
    fens, fesa = mergenmeshes(meshes, 0.0001)
    fes = cat(fesa[1], fesa[2])
    fens1, fes1 = mirrormesh(
        fens,
        fes,
        [-1.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
        renumb = (c) -> c[[1, 3, 2, 4]],
    )
    meshes = Array{Tuple{FENodeSet,AbstractFESet},1}()
    push!(meshes, (fens, fes))
    push!(meshes, (fens1, fes1))
    fens, fesa = mergenmeshes(meshes, 0.0001)
    fes = cat(fesa[1], fesa[2])
    @test count(fes) == 900
    File = "mT4cylinderm2ca-mesh.vtk"
    VTK.vtkexportmesh(File, fens, fes)
    rm(File)
    # @async run(`"paraview.exe" $File`)
    bfes = meshboundary(fes)

    File = "mT4cylinderm2ca-bmesh.vtk"
    VTK.vtkexportmesh(File, fens, fes)
    rm(File)
    # @async run(`"paraview.exe" $File`)
    V = pi * Radius^2 * Length
    @test (V - __volume(fens, fes)) / V < 0.02
end
test()
end

module mT4cylinderm2cb
using FinEtools
using FinEtools.MeshExportModule: VTK
using Test
function __volume(fens, fes)
    geom = NodalField(fens.xyz)
    femm = FEMMBase(IntegDomain(fes, TetRule(1)))
    return integratefunction(femm, geom, (x) -> 1.0)
end
function test()
    Radius, Length, nperradius, nL = 3.0, 5.0, 3, 5
    fens, fes = T4quartercyln(Radius, Length, nperradius, nL; orientation = :cb)
    fens1, fes1 = mirrormesh(
        fens,
        fes,
        [0.0, -1.0, 0.0],
        [0.0, 0.0, 0.0],
        renumb = (c) -> c[[1, 3, 2, 4]],
    )
    meshes = Array{Tuple{FENodeSet,AbstractFESet},1}()
    push!(meshes, (fens, fes))
    push!(meshes, (fens1, fes1))
    fens, fesa = mergenmeshes(meshes, 0.0001)
    fes = cat(fesa[1], fesa[2])
    fens1, fes1 = mirrormesh(
        fens,
        fes,
        [-1.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
        renumb = (c) -> c[[1, 3, 2, 4]],
    )
    meshes = Array{Tuple{FENodeSet,AbstractFESet},1}()
    push!(meshes, (fens, fes))
    push!(meshes, (fens1, fes1))
    fens, fesa = mergenmeshes(meshes, 0.0001)
    fes = cat(fesa[1], fesa[2])
    @test count(fes) == 900
    File = "mT4cylinderm2cb-mesh.vtk"
    VTK.vtkexportmesh(File, fens, fes)
    rm(File)
    # @async run(`"paraview.exe" $File`)
    bfes = meshboundary(fes)

    File = "mT4cylinderm2cb-bmesh.vtk"
    VTK.vtkexportmesh(File, fens, fes)
    rm(File)
    # @async run(`"paraview.exe" $File`)
    V = pi * Radius^2 * Length
    @test (V - __volume(fens, fes)) / V < 0.02
end
test()
end


module mT10cylinderm2b
using FinEtools
using FinEtools.MeshExportModule: VTK
using LinearAlgebra
using Test
function __volume(fens, fes)
    geom = NodalField(fens.xyz)
    femm = FEMMBase(IntegDomain(fes, TetRule(4)))
    return integratefunction(femm, geom, (x) -> 1.0)
end
function test()
    Radius, Length, nperradius, nL = 3.0, 5.0, 3, 5
    fens, fes = T4quartercyln(Radius, Length, nperradius, nL)
    fens, fes = T4toT10(fens, fes)
    bfes = meshboundary(fes)
    el = selectelem(fens, bfes, facing = true, direction = [1.0, 1.0, 0.0])
    cbfes = subset(bfes, el)
    for i in eachindex(cbfes)
        for k in cbfes.conn[i]
            fens.xyz[k, 1:2] = fens.xyz[k, 1:2] * Radius / norm(fens.xyz[k, 1:2])
        end
    end
    fens1, fes1 = mirrormesh(
        fens,
        fes,
        [0.0, -1.0, 0.0],
        [0.0, 0.0, 0.0],
        renumb = (c) -> c[[1, 3, 2, 4, 7, 6, 5, 8, 10, 9]],
    )
    meshes = Array{Tuple{FENodeSet,AbstractFESet},1}()
    push!(meshes, (fens, fes))
    push!(meshes, (fens1, fes1))
    fens, fesa = mergenmeshes(meshes, 0.0001)
    fes = cat(fesa[1], fesa[2])
    fens1, fes1 = mirrormesh(
        fens,
        fes,
        [-1.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
        renumb = (c) -> c[[1, 3, 2, 4, 7, 6, 5, 8, 10, 9]],
    )
    meshes = Array{Tuple{FENodeSet,AbstractFESet},1}()
    push!(meshes, (fens, fes))
    push!(meshes, (fens1, fes1))
    fens, fesa = mergenmeshes(meshes, 0.0001)
    fes = cat(fesa[1], fesa[2])
    @test count(fes) == 270 * 4
    File = "mT10cylinderm2b-mesh.vtk"
    VTK.vtkexportmesh(File, fens, fes)
    rm(File)
    # @async run(`"paraview.exe" $File`)
    bfes = meshboundary(fes)

    File = "mT10cylinderm2b-bmesh.vtk"
    VTK.vtkexportmesh(File, fens, fes)
    rm(File)
    # @async run(`"paraview.exe" $File`)
    V = pi * Radius^2 * Length
    @test (V - __volume(fens, fes)) / V < 0.002
end
test()
end


module mT4cylinderm4
using FinEtools
using FinEtools.MeshExportModule: VTK
using LinearAlgebra
using Test
function __volume(fens, fes)
    geom = NodalField(fens.xyz)
    femm = FEMMBase(IntegDomain(fes, TetRule(4)))
    return integratefunction(femm, geom, (x) -> 1.0)
end
function test()
    Radius, Length, nperradius, nL = 3.0, 5.0, 3, 5
    fens, fes = T4cylindern(Radius, Length, nperradius, nL)
    @test count(fes) == 270 * 4
    File = "mT4cylinderm4-mesh.vtk"
    VTK.vtkexportmesh(File, fens, fes)
    # rm(File)
    # @async run(`"paraview.exe" $File`)
    bfes = meshboundary(fes)

    File = "mT4cylinderm4-bmesh.vtk"
    VTK.vtkexportmesh(File, fens, fes)
    # rm(File)
    # @async run(`"paraview.exe" $File`)
    V = pi * Radius^2 * Length
    @test (V - __volume(fens, fes)) / V < 0.02
end
test()
end

module mT10cylinderm4
using FinEtools
using FinEtools.MeshExportModule: VTK
using LinearAlgebra
using Test
function __volume(fens, fes)
    geom = NodalField(fens.xyz)
    femm = FEMMBase(IntegDomain(fes, TetRule(4)))
    return integratefunction(femm, geom, (x) -> 1.0)
end
function test()
    Radius, Length, nperradius, nL = 3.0, 5.0, 3, 5
    fens, fes = T10cylindern(Radius, Length, nperradius, nL)
    @test count(fes) == 270 * 4
    File = "mT10cylinderm4-mesh.vtk"
    VTK.vtkexportmesh(File, fens, fes)
    # rm(File)
    # @async run(`"paraview.exe" $File`)
    bfes = meshboundary(fes)

    File = "mT10cylinderm4-bmesh.vtk"
    VTK.vtkexportmesh(File, fens, fes)
    # rm(File)
    # @async run(`"paraview.exe" $File`)
    V = pi * Radius^2 * Length
    @test (V - __volume(fens, fes)) / V < 0.002
end
test()
end

module mmmmaps1
using FinEtools
using Test
function test(n = 20)
    W, L, H = 3.5, 7.1, 9.3
    fens, fes = H8block(W, L, H, n, n, n)
    # println("Number of elements: $(count(fes))")
    
    geom = NodalField(fens.xyz)
    u = NodalField(zeros(size(fens.xyz, 1), 3)) # displacement field

    femm = FEMMBase(IntegDomain(fes, GaussRule(3, 2)))

    numberdofs!(u)

    # println("nalldofs(u) = $(nalldofs(u))").#
    n2e = FENodeToFEMap(fes.conn, count(fens))
    n2n = FENodeToNeighborsMap(n2e, fes.conn)
    
    found = true
    for i in eachindex(fes.conn)
        for k in fes.conn[i]
            for m in fes.conn[i]
                found = found && ((k in n2n.map[m]) && (m in n2n.map[k]))
            end
        end
    end
    @test found
end
test(13)
test(17)
test(18)
nothing
end

module mesh_test_3
using FinEtools
using FinEtools.MeshExportModule
using Statistics
using LinearAlgebra
using Test
function test()
    W = 25.0
    H = 50.0
    L = 50.0
    N = 11
    
    fens, fes = H8block(W, L, H, N, N, N)
      
    # vtkexportmesh("H8block.vtk", fes.conn, fens.xyz,  FinEtools.MeshExportModule.VTK.H8)
    
    n2e = FENodeToFEMap(fes.conn, count(fens))

    element_colors, unique_colors = element_coloring(fes, n2e)

    @test norm(unique(element_colors) - [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]) == 0

    # @show partitions = unique(element_colors)
    # for j in 1:length(partitions)
    #     sfes = subset(fes, findall(v -> v == partitions[j], element_colors))
    #     @show count(sfes)
    #     vtkexportmesh("H8block-c=$(partitions[j]).vtk", fens, sfes)
    # end
    
end
test()
end


module mesh_test_4
using FinEtools
using FinEtools.MeshExportModule
using Statistics
using LinearAlgebra
using Test
function test()
    W = 25.0
    H = 50.0
    L = 50.0
    N = 11

    fens, fes = T4block(W, L, H, N, N, N)

    # vtkexportmesh("T4block.vtk", fes.conn, fens.xyz,  FinEtools.MeshExportModule.VTK.T4)

    n2e = FENodeToFEMap(fes.conn, count(fens))

    element_colors, unique_colors = element_coloring(fes, n2e)

    @test norm(unique(element_colors) - [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36]) == 0

    # @show partitions = unique(element_colors)
    # for j in 1:length(partitions)
    #     sfes = subset(fes, findall(v -> v == partitions[j], element_colors))
    #     @show count(sfes)
    #     vtkexportmesh("T4block-c=$(partitions[j]).vtk", fens, sfes)
    # end

end
test()
end
