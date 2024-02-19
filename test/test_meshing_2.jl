
module momap2para3_2
using FinEtools
using FinEtools.MeshSelectionModule: vselect
using FinEtools.MeshExportModule
using FinEtools.MeshExportModule: VTKWrite
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
    File = "momap2para3_2-coarse.vtu"
    VTKWrite.vtkwrite(File, fensc, fesc; scalars = [("fc", fc.values)])
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
    File = "momap2para3_2-reference.vtu"
    VTKWrite.vtkwrite(File, fensf, fesf; scalars = [("referenceff", referenceff.values)])
    # @async run(`"paraview.exe" $File`)
    try
        rm(File)
    catch
    end

    ff = transferfield!(ff, fensf, fesf, fc, fensc, fesc, tolerance)
    File = "momap2para3_2-fine.vtu"
    VTKWrite.vtkwrite(File, fensf, fesf; scalars = [("ff", ff.values)])
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
using .momap2para3_2
momap2para3_2.test()

module momap2para3_3
using FinEtools
using FinEtools.MeshSelectionModule: vselect
using FinEtools.MeshExportModule
using FinEtools.MeshExportModule: VTKWrite
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
    File = "momap2para3_3-coarse.vtu"
    VTKWrite.vtkwrite(File, fensc, fesc; scalars = [("fc", fc.values)])
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
    File = "momap2para3_3-reference.vtu"
    VTKWrite.vtkwrite(File, fensf, fesf; vectors = [("geometry", deepcopy(fensf.xyz)),], scalars = [("referenceff", referenceff.values)])
    # @async run(`"paraview.exe" $File`)
    try
        # rm(File)
    catch
    end

end
end
using .momap2para3_3
momap2para3_3.test()

module mmmQ4blockneous2mm_2
using FinEtools
using FinEtools.MeshExportModule: VTKWrite
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
    File = "mesh.vtu"
    VTKWrite.vtkwrite(File, fens, fes)
    # rm(File)
    # @async run(`"paraview.exe" $File`)
end
end
using .mmmQ4blockneous2mm_2
mmmQ4blockneous2mm_2.test()

module miimportexportm_2
using FinEtools
using FinEtools.MeshImportModule
using FinEtools.MeshExportModule: VTKWrite
using Test
function test()
    output = MeshImportModule.import_NASTRAN(
        dirname(@__FILE__) * "/" * "Slot-coarser.nas";
        allocationchunk = 13,
    )
    # show(fes.conn[count(fes), :])
    File = "Slot-coarser.vtu"
    VTKWrite.vtkwrite(File, output["fens"], output["fesets"][1])
    rm(File)
    @test output["fesets"][1].conn[count(output["fesets"][1]), :] ==
          NTuple{10,Int64}[(143, 140, 144, 138, 361, 363, 176, 519, 781, 520)]
    # @async run(`"paraview.exe" $File`)
end
end
using .miimportexportm_2
miimportexportm_2.test()

module mmmLshapemmm_2
using FinEtools
using FinEtools.MeshExportModule: VTKWrite
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

    File = "L_shape.vtu"
    VTKWrite.vtkwrite(File, connasarray(fes), geom.values, VTKWrite.Q4)
    # @async run(`"paraview.exe" $File`)
    rm(File)
    VTKWrite.vtkwrite(File, fes.conn, geom.values, VTKWrite.Q4)
    # @async run(`"paraview.exe" $File`)
    rm(File)

    # println("Done")
    true
end
end
using .mmmLshapemmm_2
0mmmLshapemmm_2.test()

module mmsmoothingm1_2
using FinEtools
using FinEtools.MeshExportModule: VTKWrite
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

    File = "mesh_smoothing_before.vtu"
    VTKWrite.vtkwrite(File, fens, fes)
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

    File = "mesh_smoothing_after.vtu"
    VTKWrite.vtkwrite(File, connasarray(fes), geom.values, VTKWrite.T3)
    # @async run(`"paraview.exe" $File`)
    try
        rm(File)
    catch
    end

    # println("Done")
    true
end
end
using .mmsmoothingm1_2
mmsmoothingm1_2.test()


module mmAbaqusexport3_2
using FinEtools
using FinEtools.MeshImportModule
using FinEtools.MeshExportModule: VTKWrite
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

    File = "LE11NAFEMS_H8.vtu"
    VTKWrite.vtkwrite(File, fens, fes)
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
using .mmAbaqusexport3_2
mmAbaqusexport3_2.test()




module momap2para3_2a
using FinEtools
using FinEtools.MeshSelectionModule: vselect
using FinEtools.MeshExportModule
using FinEtools.MeshExportModule: VTKWrite
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
    File = "momap2para3_2a-coarse.vtu"
    VTKWrite.vtkwrite(File, fensc, fesc; scalars = [("fc", fc.values)])
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
    File = "momap2para3_2a-reference.vtu"
    VTKWrite.vtkwrite(File, fensf, fesf; scalars = [("referenceff", referenceff.values)])
    # @async run(`"paraview.exe" $File`)
    try
        rm(File)
    catch
    end

    ff = transferfield!(ff, fensf, fesf, fc, fensc, fesc, tolerance)
    File = "momap2para3_2a-fine.vtu"
    VTKWrite.vtkwrite(File, fensf, fesf; scalars = [("ff", ff.values)])
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
using .momap2para3_2a
momap2para3_2a.test()

module momap2para4_2
using FinEtools
using FinEtools.MeshSelectionModule: vselect
using FinEtools.MeshExportModule
using FinEtools.MeshExportModule: VTKWrite
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
    File = "momap2para3-coarse.vtu"
    VTKWrite.vtkwrite(File, fensc, fesc; scalars = [("fc", fc.values)])
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
    File = "momap2para3-reference.vtu"
    VTKWrite.vtkwrite(File, fensf, fesf; scalars = [("referenceff", referenceff.values)])
    # @async run(`"paraview.exe" $File`)
    try
        rm(File)
    catch
    end

    ff = transferfield!(ff, fensf, fesf, fc, fensc, fesc, tolerance)
    File = "momap2para3-fine.vtu"
    VTKWrite.vtkwrite(File, fensf, fesf; scalars = [("ff", ff.values)])
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
using .momap2para4_2
momap2para4_2.test()

module momap2para5_2
using FinEtools
using FinEtools.MeshSelectionModule: vselect
using FinEtools.MeshExportModule
using FinEtools.MeshExportModule: VTKWrite
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
    File = "momap2para3-coarse.vtu"
    VTKWrite.vtkwrite(File, fensc, fesc; scalars = [("fc", fc.values)])
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
    File = "momap2para3-reference.vtu"
    VTKWrite.vtkwrite(File, fensf, fesf; scalars = [("referenceff", referenceff.values)])
    # @async run(`"paraview.exe" $File`)
    try
        rm(File)
    catch
    end

    ff = transferfield!(ff, fensf, fesf, fc, fensc, fesc, tolerance)
    File = "momap2para3-fine.vtu"
    VTKWrite.vtkwrite(File, fensf, fesf; scalars = [("ff", ff.values)])
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
using .momap2para5_2
momap2para5_2.test()


module momap2para6_2
using FinEtools
using FinEtools.MeshSelectionModule: vselect
using FinEtools.MeshExportModule
using FinEtools.MeshExportModule: VTKWrite
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
    File = "momap2para3-coarse.vtu"
    VTKWrite.vtkwrite(File, fensc, fesc; scalars = [("fc", fc.values)])
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
    File = "momap2para3-reference.vtu"
    VTKWrite.vtkwrite(File, fensf, fesf; scalars = [("referenceff", referenceff.values)])
    # @async run(`"paraview.exe" $File`)
    try
        rm(File)
    catch
    end

    ff = transferfield!(ff, fensf, fesf, fc, fensc, fesc, tolerance)
    File = "momap2para3-fine.vtu"
    VTKWrite.vtkwrite(File, fensf, fesf; scalars = [("ff", ff.values)])
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
using .momap2para6_2
momap2para6_2.test()

module momap2para12_2
using FinEtools
using FinEtools.MeshSelectionModule: vselect
using FinEtools.MeshExportModule
using FinEtools.MeshExportModule: VTKWrite
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
    File = "momap2para12_2-coarse.vtu"
    VTKWrite.vtkwrite(File, fensc, fesc; scalars = [("fc", fc.values)])
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
    File = "momap2para12_2-reference.vtu"
    VTKWrite.vtkwrite(File, fensf, fesf; scalars = [("referenceff", referenceff.values)])
    # @async run(`"paraview.exe" $File`)
    try
        rm(File)
    catch
    end

    ff = transferfield!(ff, fensf, fesf, fc, fensc, fesc, tolerance)
    File = "momap2para12_2-fine.vtu"
    VTKWrite.vtkwrite(File, fensf, fesf; scalars = [("ff", ff.values)])
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
using .momap2para12_2
momap2para12_2.test()


module mt4refine1_2
using FinEtools
using FinEtools.MeshExportModule
using FinEtools.MeshExportModule: VTKWrite
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

    File = "Refine-T4-a.vtu"
    VTKWrite.vtkwrite(File, fens, bfes)
    rm(File)
    # @async run(`"paraview.exe" $File`)
end
end
using .mt4refine1_2
mt4refine1_2.test()



module mimportexportm1_2
using FinEtools
using FinEtools.MeshImportModule
using FinEtools.MeshExportModule
using FinEtools.MeshExportModule: VTKWrite
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
    File = "cylinder.vtu"
    VTKWrite.vtkwrite(File, output["fens"], output["fesets"][1])
    rm(File)
    #   @test output["fesets"][1].conn[count(output["fesets"][1]), :] == NTuple{10,Int64}[(143, 140, 144, 138, 361, 363, 176, 519, 781, 520)]
    # @async run(`"paraview.exe" $File`)
end
end
using .mimportexportm1_2
mimportexportm1_2.test()




module momap2para61_2
using FinEtools
using FinEtools.MeshSelectionModule: vselect
using FinEtools.MeshExportModule
using FinEtools.MeshExportModule: VTKWrite
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
    File = "momap2para61_2-coarse.vtu"
    VTKWrite.vtkwrite(File, fensc, fesc; scalars = [("fc", fc.values)])
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
    File = "momap2para61_2-reference.vtu"
    VTKWrite.vtkwrite(File, fensf, fesf; scalars = [("referenceff", referenceff.values)])
    # @async run(`"paraview.exe" $File`)
    try
        rm(File)
    catch
    end

    ff = transferfield!(ff, fensf, fesf, fc, fensc, fesc, tolerance)
    File = "momap2para61_2-fine.vtu"
    VTKWrite.vtkwrite(
        File,
        fensf,
        fesf;
        scalars = [("ff", ff.values), ("ffcopy", ff.values)],
    )
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
using .momap2para61_2
momap2para61_2.test()

module momap2para6378_2
using FinEtools
using FinEtools.MeshSelectionModule: vselect
using FinEtools.MeshExportModule
using FinEtools.MeshExportModule: VTKWrite
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
    File = "momap2para3-coarse.vtu"
    VTKWrite.vtkwrite(File, fensc, fesc; scalars = [("fc", fc.values)])
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
    File = "momap2para3-reference.vtu"
    VTKWrite.vtkwrite(File, fensf, fesf; scalars = [("referenceff", referenceff.values)])
    # @async run(`"paraview.exe" $File`)
    try
        rm(File)
    catch
    end

    ff = transferfield!(ff, fensf, fesf, fc, fensc, fesc, tolerance)
    File = "momap2para3-fine.vtu"
    VTKWrite.vtkwrite(
        File,
        fensf,
        fesf;
        scalars = [("ff", ff.values), ("ffcopy", ff.values)],
    )
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
using .momap2para6378_2
momap2para6378_2.test()

module mcircle3
using FinEtools
using FinEtools.MeshExportModule: VTK
using Test
function test()
    fens, fes = T3circlen(1.3, 5) # Mesh
    # File = "mesh.vtk"
    #   VTK.vtkexportmesh(File, fens, fes)
    l = selectnode(fens, box = [0.0, 1.3, 0.0, 0.0], inflate = 0.01)
    @test l == [1, 3, 7, 12, 16, 28, 32, 44, 49]
    true
end
end
using .mcircle3
mcircle3.test()

module mcircle3a
using FinEtools
using FinEtools.MeshExportModule: VTK
using Test
function test()
    r = 1.3
    m, n = 8, 4
    fens, fes = T3block(pi / 2, r, m, n) # Mesh
    for i in eachindex(fens)
        a = pi / 2 - fens.xyz[i, 1]
        r = fens.xyz[i, 2]
        fens.xyz[i, :] .= (r * cos(a), r * sin(a))
    end
    fens, newn = fusenodes(fens, fens, r / 1000)
    updateconn!(fes, newn)
    connected = findunconnnodes(fens, fes)
    fens, newn = compactnodes(fens, connected)
    @test count(fens) == (m + 1) * (n + 1) - (m + 1) + 1
    fes = renumberconn!(fes, newn)
    keep = fill(true, count(fes))
    ca = connasarray(fes)
    for j in eachindex(fes)
        keep[j] = (length(unique(ca[j, :])) == 3)
    end
    l = findall(x -> x == true, keep)
    fes = subset(fes, l)
    @test count(fes) == 2 * m * n - m
    # File = "mesh.vtk"
    # VTK.vtkexportmesh(File, fens, fes)
    bfes = meshboundary(fes)
    # File = "mesh-boundary.vtk"
    # VTK.vtkexportmesh(File, fens, bfes)
    @test count(bfes) == m + 2 * n
    true
end
end
using .mcircle3a
mcircle3a.test()

module mcircle3b
using FinEtools
using FinEtools.MeshExportModule: VTK
using Test
function test()
    r = 1.3
    nc, nr = 8, 4
    fens, fes = T3circleseg(pi / 2, r, nc, nr) # Mesh
    # File = "mesh.vtk"
    # VTK.vtkexportmesh(File, fens, fes)
    bfes = meshboundary(fes)
    # File = "mesh-boundary.vtk"
    # VTK.vtkexportmesh(File, fens, bfes)
    @test count(fens) == (nc + 1) * (nr + 1) - (nc + 1) + 1
    @test count(fes) == 2 * nc * nr - nc
    @test count(bfes) == nc + 2 * nr
    true
end
end
using .mcircle3b
mcircle3b.test()

module ml2extrude1
using FinEtools
using FinEtools.MeshExportModule: VTK
using Test
function test()
    X = [
        0.0 0.0
        0.0 -0.2
        0.3 -0.2
        0.3 -0.4
    ]
    l2fens = FENodeSet(X)
    conn = [1 2; 2 3; 3 4]
    l2fes = FESetL2(conn)
    # fens,fes = Q4extrudeL2(l2fens, l2fes, 5, (x, k) -> vcat(vec(x), (k-1)*1)); # Mesh
    ex = function (x, k)
        x .+= k * 0.1
        x
    end
    fens, fes = Q4extrudeL2(l2fens, l2fes, 5, ex) # Mesh
    # File = "mesh.vtk"
    # VTK.vtkexportmesh(File, fens, fes)
    # bfes = meshboundary(fes)
    # File = "mesh-boundary.vtk"
    # VTK.vtkexportmesh(File, fens, bfes)
    @test count(fens) == 24
    @test count(fes) == 15
    # @test count(bfes) == nc+2*nr
    true
end
end
using .ml2extrude1
ml2extrude1.test()


module ml2extrude2
using FinEtools
using FinEtools.MeshExportModule: VTK
using Test
function test()
    X = [
        0.0 0.0 0.0
        0.0 -0.2 0.0
        0.3 -0.2 0.0
        0.3 -0.4 0.0
    ]
    l2fens = FENodeSet(X)
    conn = [1 2; 2 3; 3 4]
    l2fes = FESetL2(conn)
    # fens,fes = Q4extrudeL2(l2fens, l2fes, 5, (x, k) -> vcat(vec(x), (k-1)*1)); # Mesh
    ex = function (x, k)
        x[3] += (k - 1) * 0.1
        x
    end
    fens, fes = Q4extrudeL2(l2fens, l2fes, 5, ex) # Mesh
    # File = "mesh.vtk"
    # VTK.vtkexportmesh(File, fens, fes)
    bfes = meshboundary(fes)
    # File = "mesh-boundary.vtk"
    # VTK.vtkexportmesh(File, fens, bfes)
    @test count(fens) == 24
    @test count(fes) == 15
    @test count(bfes) == 5 + 5 + 3 + 3
    true
end
end
using .ml2extrude2
ml2extrude2.test()

module miscellaneous_1g
using FinEtools
using Test
using FinEtools.MeshExportModule
function test()
    rho = 1.21 * 1e-9# mass density
    c = 345.0 * 1000# millimeters per second
    bulk = c^2 * rho
    Lx = 900.0# length of the box, millimeters
    Ly = 800.0 # length of the box, millimeters
    nx = 8
    ny = 4
    a = Lx / nx * 3 / 2
    b = Ly / ny * 3 / 2

    fens, fes = T3block(Lx, Ly, nx, ny) # Mesh
    fens = distortblock(fens, 1.0 / nx, -1.0 / ny)

    File = "mesh_1g.vtk"
    MeshExportModule.VTK.vtkexportmesh(File, fens, fes)

end
end
using .miscellaneous_1g
miscellaneous_1g.test()


module miscellaneous_2g
using FinEtools
using Test
using FinEtools.MeshExportModule
function test()
    rho = 1.21 * 1e-9# mass density
    c = 345.0 * 1000# millimeters per second
    bulk = c^2 * rho
    Lx = 900.0# length of the box, millimeters
    Ly = 800.0 # length of the box, millimeters
    nx = 8
    ny = 4
    a = Lx / nx * 3 / 2
    b = Ly / ny * 3 / 2

    fens, fes = Q8block(Lx, Ly, nx, ny) # Mesh
    fens = distortblock(fens, 1.0 / nx, -1.0 / ny)

    File = "mesh_2g.vtk"
    MeshExportModule.VTK.vtkexportmesh(File, fens, fes)

end
end
using .miscellaneous_2g
miscellaneous_2g.test()

module miscellaneous_3g
using FinEtools
using Test
using FinEtools.MeshExportModule
function test()
    rho = 1.21 * 1e-9# mass density
    c = 345.0 * 1000# millimeters per second
    bulk = c^2 * rho
    Lx = 900.0# length of the box, millimeters
    Ly = 800.0 # length of the box, millimeters
    nx = 8
    ny = 5
    a = Lx / nx * 1.2
    b = Ly / ny * 3 / 2

    fens, fes = Q4block(Lx, Ly, nx, ny) # Mesh
    fens = distortblock(fens, 1.0 / nx, -1.0 / ny)

    File = "mesh_3g.vtk"
    MeshExportModule.VTK.vtkexportmesh(File, fens, fes)

end
end
using .miscellaneous_3g
miscellaneous_3g.test()

module miscellaneous_4g
using FinEtools
using Test
using FinEtools.MeshExportModule
function test()
    rho = 1.21 * 1e-9# mass density
    c = 345.0 * 1000# millimeters per second
    bulk = c^2 * rho
    Lx = 900.0# length of the box, millimeters
    Ly = 800.0 # length of the box, millimeters
    nx = 8
    ny = 5
    a = Lx / nx * 1.2
    b = Ly / ny * 3 / 2

    fens, fes = Q4block(Lx, Ly, nx, ny) # Mesh
    fens = distortblock(fens, 1.0 / nx, 1.0 / ny)

    File = "mesh_4g.vtk"
    MeshExportModule.VTK.vtkexportmesh(File, fens, fes)

end
end
using .miscellaneous_4g
miscellaneous_4g.test()


module miscellaneous_5g
using FinEtools
using Test
using FinEtools.MeshExportModule
function test()
    rho = 1.21 * 1e-9# mass density
    c = 345.0 * 1000# millimeters per second
    bulk = c^2 * rho
    Lx = 900.0# length of the box, millimeters
    Ly = 800.0 # length of the box, millimeters
    nx = 8
    ny = 5
    a = Lx / nx * 1.2
    b = Ly / ny * 3 / 2

    fens, fes = distortblock(Q4block, Lx, Ly, nx, ny, 1.0 / nx, 1.0 / ny)

    File = "mesh_5g.vtk"
    MeshExportModule.VTK.vtkexportmesh(File, fens, fes)

end
end
using .miscellaneous_5g
miscellaneous_5g.test()


module miscellaneous_6g
using FinEtools
using Test
using FinEtools.MeshExportModule
function test()
    rho = 1.21 * 1e-9# mass density
    c = 345.0 * 1000# millimeters per second
    bulk = c^2 * rho
    Lx = 900.0# length of the box, millimeters
    Ly = 800.0 # length of the box, millimeters
    nx = 8
    ny = 5
    a = Lx / nx * 1.2
    b = Ly / ny * 3 / 2

    fens, fes = distortblock(Q4block, Lx, Ly, nx, ny, 1.0 / nx, 0.0 * Ly / ny)

    File = "mesh_6g.vtk"
    MeshExportModule.VTK.vtkexportmesh(File, fens, fes)

end
end
using .miscellaneous_6g
miscellaneous_6g.test()


module mt3rand1x
using FinEtools
using FinEtools.MeshExportModule: VTK
using LinearAlgebra
using Test
function test()
    Lx = 900.0# length of the box, millimeters
    Ly = 800.0 # length of the box, millimeters
    nx = 8
    ny = 11

    fens, fes = T3blockrand(Lx, Ly, nx, ny) # Mesh
    File = "mesh_c.vtk"
    VTK.vtkexportmesh(File, fens, fes)

    true
end
end
using .mt3rand1x
mt3rand1x.test()


module mt3rand2x
using FinEtools
using FinEtools.MeshExportModule: VTK
using LinearAlgebra
using Test
function test()
    Lx = 1900.0# length of the box, millimeters
    Ly = 800.0 # length of the box, millimeters
    nx = 38
    ny = 35

    fens, fes = T3blockrand(Lx, Ly, nx, ny) # Mesh
    File = "mesh_d.vtk"
    VTK.vtkexportmesh(File, fens, fes)

    true
end
end
using .mt3rand2x
mt3rand2x.test()


module mvtkcollection1
using FinEtools
using FinEtools.MeshSelectionModule: vselect
using FinEtools.MeshExportModule
using FinEtools.MeshExportModule: VTKWrite
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
    scalars = []
    times = []
    for t in 0.0:0.1:1.0
        for i in eachindex(fensc)
            x, y, z = fensc.xyz[i, :]
            fc.values[i, :] .=
                sin(5 * t) * sin(2 * x / A) * cos(6.5 * y / B) * sin(3 * z / C - 1.0)
        end
        push!(scalars, ("p", deepcopy(fc.values)))
        push!(times, t)
    end
    File = "mvtkcollection1"
    VTKWrite.vtkwritecollection(File, fensc, fesc, times; scalars = scalars)
    # @async run(`"paraview.exe" $File`)
    try
        rm(File)
    catch
    end
end
end
using .mvtkcollection1
mvtkcollection1.test()


module mtq8_integrate1
using FinEtools
using FinEtools.MeshExportModule
using FinEtools.MeshExportModule: VTKWrite
using Test
function test()
    xs = collect(linearspace(0.0, pi / 2, 5))
    ys = collect(linearspace(0.0, 10.0, 6) .^ 2)
    fens, fes = Q8blockx(xs, ys,)


    geom = NodalField(fens.xyz)
    femm = FEMMBase(IntegDomain(fes, GaussRule(2, 4)))
    V = integratefunction(femm, geom, (x) -> 1.0)
    @test V ≈ pi / 2 * 10.0^2
    nothing
end
test()
end


module mtq8_integrate1a
using FinEtools
using FinEtools.MeshExportModule
using FinEtools.MeshExportModule: VTKWrite
using Test
function test()
    xs = collect(linearspace(0.0, pi / 2, 5))
    ys = collect(linearspace(0.0, 10.0, 6))
    fens, fes = Q8blockx(xs, ys,)


    geom = NodalField(fens.xyz)
    femm = FEMMBase(IntegDomain(fes, GaussRule(2, 4)))
    V = integratefunction(femm, geom, (x) -> 1.0)
    @test V ≈ pi / 2 * 10.0
    nothing
end
test()
end

module mtq8_integrate2
using FinEtools
using FinEtools.MeshExportModule
using FinEtools.MeshExportModule: VTK
using Test
function test()
    xs = collect(linearspace(0.0, 5 * pi / 2, 6))
    ys = collect(linearspace(0.0, 10.0, 7))
    fens, fes = Q8blockx(xs, ys,)
    @test count(fens) == (5 + 1) * (6 + 1) + 5 * (6 + 1) + 6 * (5 + 1)
    @test count(fes) == 5 * 6

    geom = NodalField(fens.xyz)
    femm = FEMMBase(IntegDomain(fes, GaussRule(2, 4)))
    V = integratefunction(femm, geom, (x) -> 1.0)
    @test V ≈ 5 * pi / 2 * 10.0
    File = "mesh_d.vtk"
    VTK.vtkexportmesh(File, fens, fes)
    nothing
end
test()
end


module mtq9_integrate2
using FinEtools
using FinEtools.MeshExportModule
using FinEtools.MeshExportModule: VTK
using Test
function test()
    xs = collect(linearspace(0.0, 5 * pi / 2, 6))
    ys = collect(linearspace(0.0, 10.0, 7))
    fens, fes = Q9blockx(xs, ys,)
    @test count(fens) == (5 + 1) * (6 + 1) + 5 * (6 + 1) + 6 * (5 + 1) + 5 * 6
    @test count(fes) == 5 * 6

    geom = NodalField(fens.xyz)
    femm = FEMMBase(IntegDomain(fes, GaussRule(2, 4)))
    V = integratefunction(femm, geom, (x) -> 1.0)
    @test V ≈ 5 * pi / 2 * 10.0
    File = "mesh_d.vtk"
    VTK.vtkexportmesh(File, fens, fes)
    nothing
end
test()
end



module mtq9_integrate3
using FinEtools
using FinEtools.MeshExportModule
using FinEtools.MeshExportModule: VTK
using Test
function test()
    xs = collect(linearspace(0.0, 5 * pi / 2, 6))
    ys = collect(linearspace(0.0, 10.0, 7))
    fens, fes = Q9blockx(xs, ys,)
    @test count(fens) == (5 + 1) * (6 + 1) + 5 * (6 + 1) + 6 * (5 + 1) + 5 * 6
    @test count(fes) == 5 * 6

    bfes  = meshboundary(fes)

    geom = NodalField(fens.xyz)
    femm = FEMMBase(IntegDomain(bfes, GaussRule(1, 4)))
    L = integratefunction(femm, geom, (x) -> 1.0)
    @test L ≈ 2 * (5 * pi / 2 + 10.0)
    File = "mesh_d.vtk"
    VTK.vtkexportmesh(File, fens, bfes)
    nothing
end
test()
end




module mtq9_integrate4
using FinEtools
using FinEtools.MeshExportModule
using FinEtools.MeshExportModule: VTKWrite
using Test
function test()
    xs = collect(linearspace(0.0, 5 * pi / 2, 6))
    ys = collect(linearspace(0.0, 10.0, 7))
    fens, fes = Q9blockx(xs, ys,)
    fens.xyz = xyz3(fens)
    @test count(fens) == (5 + 1) * (6 + 1) + 5 * (6 + 1) + 6 * (5 + 1) + 5 * 6
    @test count(fes) == 5 * 6

    bfes  = meshboundary(fes)

    geom = NodalField(fens.xyz)
    femm = FEMMBase(IntegDomain(bfes, GaussRule(1, 4)))
    L = integratefunction(femm, geom, (x) -> 1.0)
    @test L ≈ 2 * (5 * pi / 2 + 10.0)
    File = "mesh_vec.vtk"
    v = deepcopy(fens.xyz)
    t = v[:, 1]; v[:, 1] = -v[:, 2]; v[:, 2] = t
    VTKWrite.vtkwrite(File, fens, fes; vectors = [("v", v), ])
    nothing
end
test()
end



module mtq9_integrate5
using FinEtools
using FinEtools.MeshExportModule
using FinEtools.MeshExportModule: VTKWrite
using Test
function test()
    xs = collect(linearspace(0.0, 5 * pi / 2, 6))
    ys = collect(linearspace(0.0, 10.0, 7))
    fens, fes = Q9blockx(xs, ys,)
    fens.xyz = xyz3(fens)
    @test count(fens) == (5 + 1) * (6 + 1) + 5 * (6 + 1) + 6 * (5 + 1) + 5 * 6
    @test count(fes) == 5 * 6

    bfes  = meshboundary(fes)

    geom = NodalField(fens.xyz)
    femm = FEMMBase(IntegDomain(bfes, GaussRule(1, 4)))
    L = integratefunction(femm, geom, (x) -> 1.0)
    @test L ≈ 2 * (5 * pi / 2 + 10.0)
    File = "mesh_vec.vtk"
    v = deepcopy(fens.xyz)
    t = v[:, 1]; v[:, 1] = -v[:, 2]; v[:, 2] = t
    VTKWrite.vtkwrite(File, fens, fes; vectors = [("v", v), ])
    nothing
end
test()
end



module mtt3q4_integrate5
using FinEtools
using FinEtools.MeshExportModule
using FinEtools.MeshExportModule: VTKWrite
using Test
function test()
    xs = collect(linearspace(0.0, 5 * pi / 2, 6))
    ys = collect(linearspace(0.0, 10.0, 7))
    fens, fes = T3blockx(xs, ys,)
    fens, fes = T3toQ4(fens, fes,)
    fens.xyz = xyz3(fens)



    bfes  = meshboundary(fes)

    geom = NodalField(fens.xyz)
    femm = FEMMBase(IntegDomain(bfes, GaussRule(1, 2)))
    L = integratefunction(femm, geom, (x) -> 1.0)
    @test L ≈ 2 * (5 * pi / 2 + 10.0)
    File = "mesh_vec_1.vtk"
    v = deepcopy(fens.xyz)
    t = v[:, 1]; v[:, 1] = -v[:, 2]; v[:, 2] = t
    VTKWrite.vtkwrite(File, fens, fes; vectors = [("v", v), ])
    nothing
end
test()
end


module misc_importexport_1
using FinEtools
using FinEtools.MeshExportModule.H5MESH: write_H5MESH
using FinEtools.MeshImportModule: import_H5MESH
using FinEtools.MeshExportModule: VTKWrite
using Test
function test()
    xs = collect(linearspace(0.0, 5 * pi / 2, 6))
    ys = collect(linearspace(0.0, 10.0, 7))

    for f in (T3blockx, Q4blockx, T6blockx, Q8blockx, Q9blockx)
        fens, fes = f(xs, ys,)
        fens.xyz = xyz3(fens)
        write_H5MESH("file", fens, fes)
        output = import_H5MESH("file")
        @test fens.xyz ≈ output["fens"].xyz
        @test connasarray(fes) ≈ connasarray(output["fesets"][1])
        File = "mesh_vec_1.vtk"
        v = deepcopy(fens.xyz)
        t = v[:, 1]; v[:, 1] = -v[:, 2]; v[:, 2] = t
        VTKWrite.vtkwrite(File, fens, fes; vectors = [("v", v), ])
    end
    nothing
end
test()
end


module misc_importexport_2
using FinEtools
using FinEtools.MeshExportModule.H5MESH: write_H5MESH
using FinEtools.MeshImportModule: import_H5MESH
using FinEtools.MeshExportModule: VTKWrite
using Test
function test()
    xs = collect(linearspace(0.0, 5 * pi / 2, 6))
    ys = collect(linearspace(0.0, 10.0, 7))
    zs = collect(linearspace(0.0, 10.0, 7))

    for f in (H8blockx, T4blockx, H20blockx, T10blockx)
        fens, fes = f(xs, ys, zs)
        fens.xyz = xyz3(fens)
        write_H5MESH("file", fens, fes)
        output = import_H5MESH("file")
        @test fens.xyz ≈ output["fens"].xyz
        @test connasarray(fes) ≈ connasarray(output["fesets"][1])
        File = "mesh_vec_1.vtk"
        v = deepcopy(fens.xyz)
        t = v[:, 1]; v[:, 1] = -v[:, 2]; v[:, 2] = t
        VTKWrite.vtkwrite(File, fens, fes; vectors = [("v", v), ])
    end
    nothing
end
test()
end


module misc_importexport_3
using FinEtools
using FinEtools.MeshExportModule.MESH: write_MESH
using FinEtools.MeshImportModule: import_MESH
using FinEtools.MeshExportModule: VTKWrite
using LinearAlgebra
using Test
function test()
    xs = collect(linearspace(0.0, 5 * pi / 2, 6))
    ys = collect(linearspace(0.0, 10.0, 7))

    for f in (T3blockx, Q4blockx, T6blockx, Q8blockx, Q9blockx)
        fens, fes = f(xs, ys,)
        fens.xyz = xyz3(fens)
        write_MESH("file", fens, fes)
        output = import_MESH("file")
        @test norm(fens.xyz - output["fens"].xyz) / norm(fens.xyz) < 1.0e-6
        @test connasarray(fes) ≈ connasarray(output["fesets"][1])
        File = "mesh_vec_1.vtk"
        v = deepcopy(fens.xyz)
        t = v[:, 1]; v[:, 1] = -v[:, 2]; v[:, 2] = t
        VTKWrite.vtkwrite(File, fens, fes; vectors = [("v", v), ])
    end
    nothing
end
test()
end


module misc_importexport_4
using FinEtools
using FinEtools.MeshExportModule.MESH: write_MESH
using FinEtools.MeshImportModule: import_MESH
using FinEtools.MeshExportModule: VTKWrite
using LinearAlgebra
using Test
function test()
    xs = collect(linearspace(0.0, 5 * pi / 2, 6))
    ys = collect(linearspace(0.0, 10.0, 7))
    zs = collect(linearspace(0.0, 10.0, 7))

    for f in (H8blockx, T4blockx, H20blockx, T10blockx)
        fens, fes = f(xs, ys, zs)
        fens.xyz = xyz3(fens)
        write_MESH("file", fens, fes)
        output = import_MESH("file")
        @test norm(fens.xyz - output["fens"].xyz) / norm(fens.xyz) < 1.0e-6
        @test connasarray(fes) ≈ connasarray(output["fesets"][1])
        File = "mesh_vec_1.vtk"
        v = deepcopy(fens.xyz)
        t = v[:, 1]; v[:, 1] = -v[:, 2]; v[:, 2] = t
        VTKWrite.vtkwrite(File, fens, fes; vectors = [("v", v), ])
    end
    nothing
end
test()
end

module misc_reorder_1
using FinEtools
using SymRCM
# using FinEtools.MeshExportModule: VTKWrite
using LinearAlgebra
using Test
function test()
    xs = collect(linearspace(0.0, 5 * pi / 2, 6))
    ys = collect(linearspace(0.0, 10.0, 7))
    zs = collect(linearspace(0.0, 10.0, 7))

    for f in (H8blockx, H20blockx )
    # for f in (H8blockx, T4blockx, H20blockx, T10blockx)
        fens, fes = f(xs, ys, zs)
        geom = NodalField(fens.xyz)
        femm = FEMMBase(IntegDomain(fes, GaussRule(3, 3)))
        V0 = integratefunction(femm, geom, (x) -> 1.0)

        C = connectionmatrix(FEMMBase(IntegDomain(fes, GaussRule(3, 3))), count(fens))
        ordering = symrcm(C)
        fens, fes = reordermesh(fens, fes, ordering)
        # File = "misc_reorder_1.vtk"
        #  VTKWrite.vtkwrite(File, fens, fes)
        geom = NodalField(fens.xyz)
        femm = FEMMBase(IntegDomain(fes, GaussRule(3, 3)))
        V = integratefunction(femm, geom, (x) -> 1.0)
        @test abs(V - V0) / V0 < 1.0e-6
    end

    
    # for f in (H8blockx, H20blockx )
    for f in (T4blockx, T10blockx)
        fens, fes = f(xs, ys, zs)
        geom = NodalField(fens.xyz)
        femm = FEMMBase(IntegDomain(fes, TetRule(4)))
        V0 = integratefunction(femm, geom, (x) -> 1.0)

        C = connectionmatrix(FEMMBase(IntegDomain(fes, TetRule(4))), count(fens))
        ordering = symrcm(C)
        fens, fes = reordermesh(fens, fes, ordering)
        # File = "misc_reorder_1.vtk"
        #  VTKWrite.vtkwrite(File, fens, fes)
        geom = NodalField(fens.xyz)
        femm = FEMMBase(IntegDomain(fes, TetRule(4)))
        V = integratefunction(femm, geom, (x) -> 1.0)
        @test abs(V - V0) / V0 < 1.0e-6
    end
    nothing
end
test()
end