
module mmassembly2ya1
using FinEtools
using Test
import LinearAlgebra: norm, cholesky
function test()
    m1 = [
        0.24406 0.599773 0.833404 0.0420141
        0.786024 0.00206713 0.995379 0.780298
        0.845816 0.198459 0.355149 0.224996
    ]
    m1 = m1' * m1
    i1 = vec([5 2 1 4])
    m2 = [
        0.146618 0.53471 0.614342 0.737833
        0.479719 0.41354 0.00760941 0.836455
        0.254868 0.476189 0.460794 0.00919633
        0.159064 0.261821 0.317078 0.77646
        0.643538 0.429817 0.59788 0.958909
    ]
    m2 = m2' * m2
    i2 = vec([2 3 1 5])
    testA = [
        2.85928 1.21875 0.891063 0.891614 2.56958 0.0 0.0
        1.21875 1.15515 0.716396 0.0714644 1.56825 0.0 0.0
        0.891063 0.716396 0.936979 0.0 1.36026 0.0 0.0
        0.891614 0.0714644 0.0 0.661253 0.813892 0.0 0.0
        2.56958 1.56825 1.36026 0.813892 4.15934 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0
    ]

    M = zeros(7, 7)
    M[i1, i1] .+= m1
    M[i2, i2] .+= m2
    # @show norm(M-testA)

    a = SysmatAssemblerSparse(0.0)
    startassembly!(a, 5*5*3, 7, 7)
    assemble!(a, m1, i1, i1)
    assemble!(a, m2, i2, i2)
    Au = makematrix!(a)
    @test maximum(abs.(testA - Au)) < 1.0e-5

    a = SysmatAssemblerSparseSymm(0.0)
    startassembly!(a, 5*5*3, 7, 7)
    assemble!(a, m1, i1, i1)
    assemble!(a, m2, i2, i2)
    A = makematrix!(a)
    @test maximum(abs.(A - Au)) < 1.0e-5

    @test maximum(abs.(testA - A)) < 1.0e-5
    @test abs.(maximum(A - transpose(A))) < 1.0e-5
end
end
using .mmassembly2ya1
mmassembly2ya1.test()

module mmassembly2ya2
using FinEtools
using Test
import LinearAlgebra: norm, cholesky, diag, diagm
function test()
    m1 = [
        0.24406 0.599773 0.833404 0.0420141
        0.786024 0.00206713 0.995379 0.780298
        0.845816 0.198459 0.355149 0.224996
    ]
    m1 = m1' * m1
    i1 = vec([5 2 1 4])
    m2 = [
        0.146618 0.53471 0.614342 0.737833
        0.479719 0.41354 0.00760941 0.836455
        0.254868 0.476189 0.460794 0.00919633
        0.159064 0.261821 0.317078 0.77646
        0.643538 0.429817 0.59788 0.958909
    ]
    m2 = m2' * m2
    i2 = vec([2 3 1 5])
    testA = [
        2.85928 1.21875 0.891063 0.891614 2.56958 0.0 0.0
        1.21875 1.15515 0.716396 0.0714644 1.56825 0.0 0.0
        0.891063 0.716396 0.936979 0.0 1.36026 0.0 0.0
        0.891614 0.0714644 0.0 0.661253 0.813892 0.0 0.0
        2.56958 1.56825 1.36026 0.813892 4.15934 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0 0.0 0.0 0.0
    ]

    a = SysmatAssemblerSparseDiag(0.0)
    startassembly!(a, 5*5*3, 7, 7)
    assemble!(a, m1, i1, i1)
    assemble!(a, m2, i2, i2)
    A = makematrix!(a)
    # @show norm(A-testA)

    @test maximum(abs.(diagm(0 => diag(testA)) - A)) < 1.0e-5
    true
end
end
using .mmassembly2ya2
mmassembly2ya2.test()


module mbas001
using FinEtools
using Test
function test()
    fens, fes = Q4block(1.3, 3.1, 3, 2) # Mesh
    @test manifdim(fes) == 2
    bfes = meshboundary(fes)
    @test manifdim(bfes) == 1
    true
end
end
using .mbas001
mbas001.test()

module mbas002
using FinEtools
using Test
function test()
    fens, fes = H8block(1.3, 3.1, 2.7, 3, 2, 4) # Mesh
    @test manifdim(fes) == 3
    bfes = meshboundary(fes)
    @test manifdim(bfes) == 2
    true
end
end
using .mbas002
mbas002.test()

module mbas003
using FinEtools
using Test
function test()
    fens, fes = T4block(1.3, 3.1, 2.7, 3, 2, 4) # Mesh
    @test manifdim(fes) == 3
    bfes = meshboundary(fes)
    @test manifdim(bfes) == 2
    true
end
end
using .mbas003
mbas003.test()


module mbas004
using FinEtools
using Test
function test()
    fens, fes = Q4block(1.3, 3.1, 3, 2) # Mesh
    @test nodesperelem(fes) == 4
    bfes = meshboundary(fes)
    @test nodesperelem(bfes) == 2
    true
end
end
using .mbas004
mbas004.test()

module mbas005
using FinEtools
using Test
function test()
    fens, fes = H8block(1.3, 3.1, 2.7, 3, 2, 4) # Mesh
    @test nodesperelem(fes) == 8
    bfes = meshboundary(fes)
    @test nodesperelem(bfes) == 4
    true
end
end
using .mbas005
mbas005.test()

module mbas006
using FinEtools
using Test
function test()
    fens, fes = T4block(1.3, 3.1, 2.7, 3, 2, 4) # Mesh
    @test nodesperelem(fes) == 4
    bfes = meshboundary(fes)
    @test nodesperelem(bfes) == 3
    true
end
end
using .mbas006
mbas006.test()


module mbas007
using FinEtools
using Test
function test()
    fens, fes = T10block(1.3, 3.1, 2.7, 3, 2, 4) # Mesh
    @test nodesperelem(fes) == 10
    bfes = meshboundary(fes)
    @test nodesperelem(bfes) == 6
    true
end
end
using .mbas007
mbas007.test()


module mbas008
using FinEtools
using Test
function test()
    fens, fes = L2block(1.37, 3) # Mesh
    @test nodesperelem(fes) == 2
    bfes = meshboundary(fes)
    @test nodesperelem(bfes) == 1
    true
end
end
using .mbas008
mbas008.test()


module mbas009
using FinEtools
using Test
function test()
    fens, fes = Q8block(1.3, 3.1, 3, 2) # Mesh
    @test nodesperelem(fes) == 8
    bfes = meshboundary(fes)
    @test nodesperelem(bfes) == 3
    true
end
end
using .mbas009
mbas009.test()

module mbas010
using FinEtools
using FinEtools.AlgoBaseModule: fieldnorm
using Test
function test()
    fens, fes = Q8block(1.3, 3.1, 13, 12) # Mesh
    geom = NodalField(fens.xyz)
    femm = FEMMBase(IntegDomain(fes, GaussRule(2, 3)))
    regions = [Dict{String,Any}("femm" => femm)]
    # x = NodalField(reshape(fens.xyz[:, 1], count(fens), 1))
    # targetfields = [x]
    c = NodalField(fill(1.0, count(fens), 1))
    targetfields = [c]
    modeldata = Dict{String,Any}(
        "fens" => fens,
        "regions" => regions,
        "targetfields" => targetfields,
        "geom" => geom,
        "elementsize" => 3.1 / 12,
    )
    @test abs(fieldnorm(modeldata) - sqrt(1.3 * 3.1)) <= 1.0e-5
    true
end
end
using .mbas010
mbas010.test()

module mbascs01
using FinEtools
using FinEtools.AlgoBaseModule: fieldnorm
using LinearAlgebra
using Test
function test()
    center = [0.0, 0.0, 0.0]
    ez = [0.0, 0.0, 1.0]
    v = rand(3)
    XYZ = reshape([0.2, 0.3, 0.4], 1, 3)
    tangents = rand(3, 3)
    rfcsmat = [
        0.5547001962252291 -0.8320502943378437 0.0
        0.8320502943378436 0.5547001962252293 0.0
        0.0 0.0 1.0
    ]
    function compute!(csmatout, XYZ, tangents, fe_label)
        # Cylindrical coordinate system
        xyz = XYZ[:] .- center
        xyz[3] = 0.0
        csmatout[:, 1] = xyz / norm(xyz)
        cross3!(v, ez, vec(csmatout[:, 1]))
        csmatout[:, 2] = v / norm(v)
        csmatout[:, 3] = ez
        return csmatout
    end
    csys = CSys(3, 3, compute!)
    updatecsmat!(csys, XYZ, tangents, 0)
    @test norm(csmat(csys) - rfcsmat) / norm(rfcsmat) <= 1.0e-6

    csys = CSys(3, 3, zero(Float32), compute!)
    updatecsmat!(csys, XYZ, tangents, 0)
    @test norm(csmat(csys) - rfcsmat) / norm(rfcsmat) <= 1.0e-6

    true
end
end
using .mbascs01
mbascs01.test()

module mbascs02
using FinEtools
using FinEtools.AlgoBaseModule: fieldnorm
using LinearAlgebra
using Test
function test()
    center = [0.0, 0.0, 0.0]
    ez = [0.0, 0.0, 1.0]
    v = rand(3)
    XYZ = reshape([0.2, 0.3, 0.4], 1, 3)
    tangents = rand(3, 3)
    rfcsmat = [
        0.5547001962252291 -0.8320502943378437 0.0
        0.8320502943378436 0.5547001962252293 0.0
        0.0 0.0 1.0
    ]
    # function compute!(csmatout::FFltMat, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)
    #     # Cylindrical coordinate system
    #     xyz = XYZ[:] .- center
    #     xyz[3] = 0.0
    #     csmatout[:, 1] = xyz / norm(xyz)
    #     cross3!(v, ez, vec(csmatout[:, 1]))
    #     csmatout[:, 2] = v / norm(v)
    #     csmatout[:, 3] = ez
    #     return csmatout
    # end
    csys = CSys(rfcsmat)
    updatecsmat!(csys, XYZ, tangents, 0)
    @test norm(csmat(csys) - rfcsmat) / norm(rfcsmat) <= 1.0e-6
    true
end
end
using .mbascs02
mbascs02.test()

module mbascs03
using FinEtools
using FinEtools.AlgoBaseModule: fieldnorm
using LinearAlgebra
using Test
function test()
    center = [0.0, 0.0, 0.0]
    ez = [0.0, 0.0, 1.0]
    v = rand(3)
    XYZ = reshape([0.2, 0.3, 0.4], 1, 3)
    tangents = rand(3, 3)
    rfcsmat = Matrix(1.0 * I, 3, 3)
    # function compute!(csmatout::FFltMat, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)
    #     # Cylindrical coordinate system
    #     xyz = XYZ[:] .- center
    #     xyz[3] = 0.0
    #     csmatout[:, 1] = xyz / norm(xyz)
    #     cross3!(v, ez, vec(csmatout[:, 1]))
    #     csmatout[:, 2] = v / norm(v)
    #     csmatout[:, 3] = ez
    #     return csmatout
    # end
    csys = CSys(3)
    updatecsmat!(csys, XYZ, tangents, 0)
    @test norm(csmat(csys) - rfcsmat) / norm(rfcsmat) <= 1.0e-6

    true
end
end
using .mbascs03
mbascs03.test()

module mbascs04
using FinEtools
using FinEtools.AlgoBaseModule: fieldnorm
using LinearAlgebra
using Test
function test()
    XYZ = reshape([0.2, 0.3, 0.4], 1, 3)
    rfcsmat = Matrix(1.0 * I, 3, 3)
    tangents = rand(3, 3)
    csys = CSys(3, 3)
    updatecsmat!(csys, XYZ, tangents, 0)
    @test norm(csmat(csys) - rfcsmat) / norm(rfcsmat) <= 1.0e-6
    true
end
end
using .mbascs04
mbascs04.test()

module mbascs05
using FinEtools
using FinEtools.AlgoBaseModule: fieldnorm
using LinearAlgebra
using Test
function test()
    XYZ = reshape([0.2, 0.3, 0.4], 1, 3)
    rfcsmat = reshape([0.37139067635410367; 0.5570860145311555; 0.7427813527082073], 3, 1)
    tangents = reshape([0.2, 0.3, 0.4], 3, 1)
    csys = CSys(3, 1)
    updatecsmat!(csys, XYZ, tangents, 0)
    # @show csmat(csys)
    @test norm(csmat(csys) - rfcsmat) / norm(rfcsmat) <= 1.0e-6
    true
end
end
using .mbascs05
mbascs05.test()

module mbascs06
using FinEtools
using FinEtools.AlgoBaseModule: fieldnorm
using LinearAlgebra
using Test
function test()
    XYZ = reshape([0.2, 0.3, 0.4], 1, 3)
    rfcsmat = reshape([0.37139067635410367; 0.5570860145311555; 0.7427813527082073], 3, 1)
    tangents = reshape([0.2 0.0; 0.0 0.5; 0.4 0.2], 3, 2)
    csys = CSys(3, 2)
    updatecsmat!(csys, XYZ, tangents, 0)
    # @show csmat(csys)
    n1 = cross(vec(tangents[:, 1]), vec(tangents[:, 2]))
    n1 = n1 / norm(n1)
    # @show n1
    n2 = cross(vec(csmat(csys)[:, 1]), vec(csmat(csys)[:, 2]))
    n2 = n2 / norm(n2)
    # @show n2
    @test norm(n1 - n2) / norm(n2) <= 1.0e-6
    true
end
end
using .mbascs06
mbascs06.test()

module mbascs06Alpha
using FinEtools
using FinEtools.MatrixUtilityModule:
    complete_lt!,
    add_mggt_ut_only!,
    add_gkgt_ut_only!,
    add_btdb_ut_only!,
    add_btsigma!,
    add_nnt_ut_only!
using LinearAlgebra
using Test
function test()
    let
        N = 8
        Ke::FFltMat, gradN::FFltMat, mult::FFlt = fill(0.0, N, N), rand(N, 2), 3.0
        add_mggt_ut_only!(Ke::FFltMat, gradN::FFltMat, mult::FFlt)
        complete_lt!(Ke::FFltMat)
        @test norm(Ke - (mult .* (gradN * gradN'))) / norm(Ke) <= 1.0e-9
    end


    let
        N = 8
        Ke::FFltMat, gradN::FFltMat, Jac_w::FFlt = fill(0.0, N, N), rand(N, 3), 0.33
        kappa_bar::FFltMat = rand(3, 3)
        kappa_bar = kappa_bar + kappa_bar'
        kappa_bargradNT::FFltMat = fill(0.0, 3, N)
        add_gkgt_ut_only!(
            Ke::FFltMat,
            gradN::FFltMat,
            Jac_w::FFlt,
            kappa_bar::FFltMat,
            kappa_bargradNT::FFltMat,
        )
        complete_lt!(Ke::FFltMat)
        @test norm(Ke - (Jac_w .* (gradN * kappa_bar * gradN'))) / norm(Ke) <= 1.0e-9
    end


    let
        N = 12
        Ke::FFltMat, B::FFltMat, Jac_w::FFlt = fill(0.0, N, N), rand(3, N), 0.33
        D::FFltMat = rand(3, 3)
        D = D + D'
        DB::FFltMat = fill(0.0, 3, N)
        add_btdb_ut_only!(Ke::FFltMat, B::FFltMat, Jac_w::FFlt, D::FFltMat, DB::FFltMat)
        complete_lt!(Ke::FFltMat)
        @test norm(Ke - (Jac_w .* (B' * D * B))) / norm(Ke) <= 1.0e-9
    end


    let
        N = 16
        Fe::FFltVec, B::FFltMat, coefficient::FFlt, sigma::FFltVec =
            fill(0.0, N), rand(3, N), 0.33, rand(3)
        D::FFltMat = rand(3, 3)
        D = D + D'
        add_btsigma!(Fe::FFltVec, B::FFltMat, coefficient::FFlt, sigma::FFltVec)
        @test norm(Fe - (coefficient * (B' * sigma))) / norm(Fe) <= 1.0e-9
    end


    let
        M = 21
        Ke::FFltMat, Nn::FFltMat, Jac_w_coeff::FFlt = fill(0.0, M, M), rand(M, 1), 0.533
        add_nnt_ut_only!(Ke::FMat{FFlt}, Nn::FFltMat, Jac_w_coeff::FFlt)
        complete_lt!(Ke::FFltMat)
        @test norm(Ke - ((Nn * Nn') * (Jac_w_coeff))) / norm(Ke) <= 1.0e-9
    end

    true
end
end
using .mbascs06Alpha
mbascs06Alpha.test()

module msymmx
using FinEtools
using FinEtools.MatrixUtilityModule: symmetrize!
using LinearAlgebra
using Test
function test()
    a = rand(6, 6)
    b = deepcopy(a)
    symmetrize!(a)
    @test norm(a - (b + b') / 2) / norm(a) <= 1.0e-9
    true
end
end
using .msymmx
msymmx.test()

module mgausr111
using FinEtools
using FinEtools.MeshExportModule: VTK
using LinearAlgebra
using Test
function test()
    fens, fes = H8spheren(3.1, 5)
    @test (count(fens), count(fes)) == (175, 108)
    # File = "mesh.vtk"
    # VTK.vtkexportmesh(File, fens, fes)
    geom = NodalField(fens.xyz)
    results = []
    for order in 1:10
        femm = FEMMBase(IntegDomain(fes, GaussRule(3, order)))
        V = integratefunction(femm, geom, (x) -> 1.0)
        push!(results, V)
    end
    @test norm(
        results - [
            15.06847962861544,
            15.137585913524491,
            15.137585913524477,
            15.137585913524472,
            15.137585913524491,
            15.13758591352459,
            15.137585913524601,
            15.137585913524363,
            15.137585913524278,
            15.137585913524612,
        ],
    ) <= 1.0e-9
    true
end
end
using .mgausr111
mgausr111.test()

module mtrapr111
using FinEtools
using FinEtools.MeshExportModule: VTK
using LinearAlgebra
using Test
function test()
    fens, fes = H8spheren(3.1, 5)
    @test (count(fens), count(fes)) == (175, 108)
    # File = "mesh.vtk"
    # VTK.vtkexportmesh(File, fens, fes)
    geom = NodalField(fens.xyz)
    femm = FEMMBase(IntegDomain(fes, TrapezoidalRule(3)))
    V = integratefunction(femm, geom, (x) -> 1.0)
    @test abs(V - 15.27579848334253) / V <= eps(V)
    true
end
end
using .mtrapr111
mtrapr111.test()

module mtrapr112
using FinEtools
using FinEtools.MeshExportModule: VTK
using LinearAlgebra
using Test
function test()
    rin::FFlt, rex::FFlt, nr::FInt, nc::FInt, Angl::FFlt, orientation::Symbol =
        100.0, 200.0, 3, 6, pi / 3, :a
    fens, fes =
        T3annulus(rin::FFlt, rex::FFlt, nr::FInt, nc::FInt, Angl::FFlt, orientation::Symbol)
    # File = "mesh.vtk"
    # VTK.vtkexportmesh(File, fens, fes)
    geom = NodalField(fens.xyz)
    femm = FEMMBase(IntegDomain(fes, TrapezoidalRule(2)))
    V1 = integratefunction(femm, geom, (x) -> 1.0)
    # 
    femm = FEMMBase(IntegDomain(fes, GaussRule(2, 2)))
    V2 = integratefunction(femm, geom, (x) -> 1.0)
    @test abs(V1 - V2) / V1 <= 1000 * eps(V1)
    true
end
end
using .mtrapr112
mtrapr112.test()


module mfld1a1
using FinEtools
using LinearAlgebra
using Test

cleandchi() = deepcopy(
    FinEtools.NodalFieldModule.NodalField{Float64, Int}(
        [
            0.0 0.0 0.0 -0.002084946198827584 0.0 0.0
            6.938893903907228e-18 1.0842021724855044e-19 -0.14078343955580946 -0.001474279595600101 -2.710505431213761e-20 1.0842021724855044e-19
            2.7755575615628914e-17 0.0 -0.19909784957735863 -8.605854744103691e-19 0.0 -5.421010862427522e-20
            2.7755575615628914e-17 1.3552527156068805e-20 -0.14078343955580952 0.001474279595600101 0.0 0.0
            0.0 0.0 0.0 0.0020849461988275853 0.0 0.0
        ],
        [21 22 23 1 24 25; 2 3 4 5 6 7; 8 9 10 11 12 13; 14 15 16 17 18 19; 26 27 28 20 29 30],
        Bool[1 1 1 0 1 1; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 1 1 1 0 1 1],
        20,
    ),
)
cleandchiv() = deepcopy([
    -0.002084946198827584,
    6.938893903907228e-18,
    1.0842021724855044e-19,
    -0.14078343955580946,
    -0.001474279595600101,
    -2.710505431213761e-20,
    1.0842021724855044e-19,
    2.7755575615628914e-17,
    0.0,
    -0.19909784957735863,
    -8.605854744103691e-19,
    0.0,
    -5.421010862427522e-20,
    2.7755575615628914e-17,
    1.3552527156068805e-20,
    -0.14078343955580952,
    0.001474279595600101,
    0.0,
    0.0,
    0.0020849461988275853,
])
function test()
    dchi = cleandchi()
    v = gathersysvec(dchi)
    @test norm(v - cleandchiv()) / norm(cleandchiv()) <= 1.0e-6
    v[:] .= 0.0
    gathersysvec!(dchi, v)
    @test norm(v - cleandchiv()) / norm(cleandchiv()) <= 1.0e-6
    elv = fill(0.0, 12)
    gathervalues_asvec!(dchi, elv, [1, 2])
    @test norm(
        elv - [
            0.0,
            0.0,
            0.0,
            -0.002084946198827584,
            0.0,
            0.0,
            6.938893903907228e-18,
            1.0842021724855044e-19,
            -0.14078343955580946,
            -0.001474279595600101,
            -2.710505431213761e-20,
            1.0842021724855044e-19,
        ],
    ) / norm(elv) <= 1.0e-6
    gatherfixedvalues_asvec!(dchi, elv, [1, 2])
    norm(elv) <= 1.0e-15
    @test FinEtools.FieldModule.anyfixedvaluenz(dchi, [1, 2]) == false
    edn = fill(0, 12)
    gatherdofnums!(dchi, edn, [1, 2])
    @test norm(edn - vec([21 22 23 1 24 25; 2 3 4 5 6 7]')) == 0
    dchi = cleandchi()
    numberdofs!(dchi)
    @test norm(dchi.dofnums - cleandchi().dofnums) == 0

    dchi = cleandchi()
    setebc!(dchi, 1, true, 1, -0.013)
    @test dchi.is_fixed[1, 1] == true
    @test dchi.values[1, 1] == -0.013
    setebc!(dchi, 5, true, 4, 0.61)
    @test dchi.is_fixed[5, 4] == true
    @test dchi.values[5, 4] == 0.61

    dchi = cleandchi()
    setebc!(dchi, [1, 5], true, 2, [-0.013, 2.0])
    @test dchi.is_fixed[1, 2] == true
    @test dchi.values[1, 2] == -0.013
    @test dchi.is_fixed[5, 2] == true
    @test dchi.values[5, 2] == 2.0

    dchi = cleandchi()
    setebc!(dchi, [1, 5], true, 2, -10.0)
    @test dchi.is_fixed[1, 2] == true
    @test dchi.values[1, 2] == -10.0
    @test dchi.is_fixed[5, 2] == true
    @test dchi.values[5, 2] == -10.0

    dchi = cleandchi()
    setebc!(dchi, [1, 5], true, 2)
    @test dchi.is_fixed[1, 2] == true
    @test dchi.values[1, 2] == 0.0
    @test dchi.is_fixed[5, 2] == true
    @test dchi.values[5, 2] == 0.0

    dchi = cleandchi()
    setebc!(dchi, [1, 5], 2, [-10.0, 10.0])
    @test dchi.is_fixed[1, 2] == true
    @test dchi.values[1, 2] == -10.0
    @test dchi.is_fixed[5, 2] == true
    @test dchi.values[5, 2] == 10.0

    dchi = cleandchi()
    setebc!(dchi, [1, 5], 2)
    @test dchi.is_fixed[1, 2] == true
    @test dchi.values[1, 2] == 0.0
    @test dchi.is_fixed[5, 2] == true
    @test dchi.values[5, 2] == 0.0

    dchi = cleandchi()
    setebc!(dchi, [1, 5], true, [2, 3], -10.0)
    @test dchi.is_fixed[1, 2] == true
    @test dchi.values[1, 2] == -10.0
    @test dchi.is_fixed[5, 3] == true
    @test dchi.values[5, 3] == -10.0

    dchi = cleandchi()
    setebc!(dchi, [1, 5], true, [2, 3])
    @test dchi.is_fixed[1, 2] == true
    @test dchi.values[1, 2] == 0.0
    @test dchi.is_fixed[5, 3] == true
    @test dchi.values[5, 3] == 0.0

    dchi = cleandchi()
    setebc!(dchi, [1, 5])
    for i in [1, 5]
        @test any([!dchi.is_fixed[i, idx] for idx in 1:ndofs(dchi)]) == false
        @test any([dchi.values[i, idx] != 0.0 for idx in 1:ndofs(dchi)]) == false
    end

    dchi = cleandchi()
    setebc!(dchi, 3)
    for i in [3]
        @test any([!dchi.is_fixed[i, idx] for idx in 1:ndofs(dchi)]) == false
        @test any([dchi.values[i, idx] != 0.0 for idx in 1:ndofs(dchi)]) == false
    end

    dchi = cleandchi()
    setebc!(dchi)
    for i in 1:nents(dchi)
        @test any([dchi.is_fixed[i, idx] for idx in 1:ndofs(dchi)]) == false
        @test any([dchi.values[i, idx] != 0.0 for idx in 1:ndofs(dchi)]) == false
    end

    dchi = cleandchi()
    v = gathersysvec(dchi)
    dchi.values[:] .= 0.0
    scattersysvec!(dchi, v)
    v = gathersysvec(dchi)
    @test norm(v - cleandchiv()) / norm(cleandchiv()) <= 1.0e-6
    v[:] .= 0.0
    gathersysvec!(dchi, v)
    @test norm(v - cleandchiv()) / norm(cleandchiv()) <= 1.0e-6

    dchi = cleandchi()
    v = gathersysvec(dchi)
    dchi.values[:] .= 0.0
    incrscattersysvec!(dchi, v)
    incrscattersysvec!(dchi, v)
    incrscattersysvec!(dchi, v)
    v = gathersysvec(dchi)
    @test norm(v - 3 .* cleandchiv()) / norm(cleandchiv()) <= 1.0e-6

    dchi = cleandchi()
    setebc!(dchi)
    numberdofs!(dchi)
    dofnums, prescribedvalues = prescribeddofs(cleandchi(), dchi)
    @test norm(dofnums - [1, 2, 3, 5, 6, 25, 26, 27, 29, 30]) == 0
    @test norm(prescribedvalues - [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]) <=
          1.0e-9
    true
end
end
using .mfld1a1
mfld1a1.test()

module mvass1
using FinEtools
using LinearAlgebra
using Test
function test()
    N = 30
    v = fill(0.0, N)
    a = SysvecAssembler()
    startassembly!(a, N)
    vec, dofnums = rand(3), [1, 7, 2]
    assemble!(a, vec, dofnums)
    for (i, p) in zip(dofnums, vec)
        v[i] += p
    end
    vec, dofnums = rand(7), [29, 15, 1, 7, 3, 6, 2]
    assemble!(a, vec, dofnums)
    for (i, p) in zip(dofnums, vec)
        v[i] += p
    end
    w = makevector!(a)
    @test norm(v - w) / norm(w) <= 1.0e-9
end
end
using .mvass1
mvass1.test()

module mvass2
using FinEtools
using LinearAlgebra
using Test
function test()
    N = 30
    v = fill(0.0, N)
    a = SysvecAssembler()
    startassembly!(a, N)
    vec, dofnums = rand(3), [1, 7, 2]
    assemble!(a, vec, dofnums)
    for (i, p) in zip(dofnums, vec)
        v[i] += p
    end
    vec, dofnums = rand(7), [29, 6, 1, 7, 3, 5, 2]
    assemble!(a, vec, dofnums)
    for (i, p) in zip(dofnums, vec)
        v[i] += p
    end
    w = makevector!(a)
    @test norm(v - w) / norm(w) <= 1.0e-9
end
end
using .mvass2
mvass2.test()


module mmhrzass1
using FinEtools
using LinearAlgebra
using Test
function test()
    N = 30
    elem_mat_dim = 8
    elem_mat_nmatrices = 3
    elmat = rand(elem_mat_dim, elem_mat_dim)
    elmat = (elmat + elmat')
    Mref = fill(0.0, N, N)
    add!(Mref, e, d) = begin
        em2 = sum(sum(e, dims = 1))
        dem2 = sum(diag(e))
        ffactor = em2 / dem2
        for i in eachindex(d)
            Mref[d[i], d[i]] += ffactor * e[i, i]
        end
        Mref
    end
    a = SysmatAssemblerSparseHRZLumpingSymm()
    startassembly!(a, elem_mat_dim * elem_mat_nmatrices, N, N)
    dofnums = [10, 29, 15, 1, 7, 3, 6, 2]
    assemble!(a, elmat, dofnums, dofnums)
    add!(Mref, elmat, dofnums)
    dofnums = [15, 1, 7, 11, 29, 3, 8, 4]
    assemble!(a, elmat, dofnums, dofnums)
    add!(Mref, elmat, dofnums)
    dofnums = reshape([25, 1, 17, 13, 30, 13, 18, 14], elem_mat_dim, 1)
    assemble!(a, elmat, dofnums, dofnums)
    add!(Mref, elmat, dofnums)
    M = makematrix!(a)
    s = 0.0
    for i in 1:N
        s += abs(M[i, i] - diag(Mref)[i])
    end
    @test s / norm(elmat) <= 1.0e-9
end
end
using .mmhrzass1
mmhrzass1.test()

module mgram1
using FinEtools
using LinearAlgebra
using Test
function test()
    W = 1.1
    L = 12.0
    t = 0.32
    nl, nt, nw = 2, 3, 4

    fens, fes = H8block(L, W, t, nl, nw, nt)
    geom = NodalField(fens.xyz)
    psi = NodalField(fill(1.0, count(fens), 1))
    numberdofs!(psi)

    femm = FEMMBase(IntegDomain(fes, GaussRule(3, 2)))
    v = gathersysvec(psi)
    G = innerproduct(femm, geom, psi)
    # @show v' * G * v
    @test abs(v' * G * v - (W*L*t)) / (W*L*t) <= 1.0e-5
    G = bilform_dot(femm, geom, psi, DataCache(LinearAlgebra.I(1)))
    @test abs(v' * G * v - (W*L*t)) / (W*L*t) <= 1.0e-5
    true
end
end
using .mgram1
mgram1.test()

module mgram1a
using FinEtools
using LinearAlgebra
using Test
function test()
    W = 1.1
    L = 12.0
    t = 0.32
    nl, nt, nw = 80, 40, 40


    fens, fes = H8block(L, W, t, nl, nw, nt)
    geom = NodalField(fens.xyz)
    psi = NodalField(fill(1.0, count(fens), 3))
    numberdofs!(psi)

    # @show nl * nt * nw * (nodesperelem(fes) * ndofs(psi))^2

    femm = FEMMBase(IntegDomain(fes, GaussRule(3, 2)))
    v = gathersysvec(psi)

    @time G1 = bilform_dot(femm, geom, psi, DataCache(LinearAlgebra.I(3)))
    # @show v' * G1 * v, 3 * (W*L*t)
    @test abs(v' * G1 * v - 3 * (W*L*t)) / (W*L*t) <= 1.0e-5
    true
end
end
using .mgram1a
mgram1a.test()
mgram1a.test()

module mdistributedl1
using FinEtools
using Test
function test()
    W = 1.1
    L = 12.0
    t = 4.32
    nl, nt, nw = 2, 3, 4

    fens, fes = H8block(L, W, t, nl, nw, nt)
    geom = NodalField(fens.xyz)
    psi = NodalField(fill(1.0, count(fens), 1))
    numberdofs!(psi)

    femm = FEMMBase(IntegDomain(fes, GaussRule(3, 2)))
    fi = ForceIntensity([11.0])
    F = distribloads(femm, geom, psi, fi, 3)
    @test abs(sum(F) - (*).(L, W, t, 11.0)) / 667 <= 1.0e-5
    fi = ForceIntensity(11.0)
    F = distribloads(femm, geom, psi, fi, 3)
    @test abs(sum(F) - (*).(L, W, t, 11.0)) / 667 <= 1.0e-5

    psi = NodalField(fill(1.0 + 1.0im, count(fens), 1))
    numberdofs!(psi)

    femm = FEMMBase(IntegDomain(fes, GaussRule(3, 2)))
    fi = ForceIntensity([11.0im])
    F = distribloads(femm, geom, psi, fi, 3)
    @test abs(sum(F) - (*).(L, W, t, 11.0im)) / 667 <= 1.0e-5
    fi = ForceIntensity(11.0im)
    F = distribloads(femm, geom, psi, fi, 3)
    @test abs(sum(F) - (*).(L, W, t, 11.0im)) / 667 <= 1.0e-5
    true
end
end
using .mdistributedl1
mdistributedl1.test()

module m1djac1
using FinEtools
using Test
function test()
    Lx = 1900.0# length of the box, millimeters
    Ly = 800.0 # length of the box, millimeters
    th = 21.0

    fens, fes = Q4block(Lx, Ly, 3, 2) # Mesh
    bfes = meshboundary(fes)
    bfemm = FEMMBase(IntegDomain(bfes, GaussRule(1, 2), th))

    geom = NodalField(fens.xyz)
    circumarea = integratefunction(bfemm, geom, (x) -> 1.0; m = 2)
    @test abs(circumarea - 2 * (Lx + Ly) * th) <= 1.0e-9
    true
end
end
using .m1djac1
m1djac1.test()

module m1djac2
using FinEtools
using Test
function test()
    Lx = 1900.0# length of the box, millimeters
    Ly = 800.0 # length of the box, millimeters
    th = 21.0

    fens, fes = Q4block(Lx, Ly, 3, 2) # Mesh
    bfes = meshboundary(fes)
    bfemm = FEMMBase(IntegDomain(bfes, GaussRule(1, 2)))

    geom = NodalField(fens.xyz)
    circumference = integratefunction(bfemm, geom, (x) -> 1.0; m = 1)
    @test abs(circumference - 2 * (Lx + Ly)) <= 1.0e-9
    true
end
end
using .m1djac2
m1djac2.test()

module m1djac3
using FinEtools
using Test
function test()
    Lx = 1900.0# length of the box, millimeters
    Ly = 800.0 # length of the box, millimeters
    th = 21.0

    fens, fes = Q4block(Lx, Ly, 3, 2) # Mesh
    bfes = meshboundary(fes)
    bfemm = FEMMBase(IntegDomain(bfes, GaussRule(1, 2), 7.0))

    geom = NodalField(fens.xyz)
    circumference = integratefunction(bfemm, geom, (x) -> 1.0; m = 3)
    @test abs(circumference / 7 - 2 * (Lx + Ly)) <= 1.0e-9
    true
end
end
using .m1djac3
m1djac3.test()

module m2djac3
using FinEtools
using Test
function test()
    Lx = 1900.0# length of the box, millimeters
    Ly = 800.0 # length of the box, millimeters
    th = 7.0

    fens, fes = Q4block(Lx, Ly, 3, 2) # Mesh
    bfes = meshboundary(fes)
    femm = FEMMBase(IntegDomain(fes, GaussRule(2, 2), th))

    geom = NodalField(fens.xyz)
    volume = integratefunction(femm, geom, (x) -> 1.0; m = 3)
    @test abs(volume - (Lx * Ly) * th) <= 1.0e-7
    true
end
end
using .m2djac3
m2djac3.test()

module mgrad2mm
using FinEtools
using FinEtools.FESetModule: gradN!
using LinearAlgebra
using Test
function test()
    rho = 1.21 * 1e-9# mass density
    c = 345.0 * 1000# millimeters per second
    bulk = c^2 * rho
    Lx = 1900.0# length of the box, millimeters
    Ly = 800.0 # length of the box, millimeters

    fens, fes = Q4block(Lx, Ly, 3, 2) # Mesh
    # show(fes.conn)
    gradN, gradNparams, redJ =
        fill(0.0, nodesperelem(fes), 2), bfundpar(fes, [0.57, 0.57]), Matrix(1.0 * I, 2, 2)

    gradN!(fes, gradN::FFltMat, gradNparams::FFltMat, redJ::FFltMat)
    @test norm(
        gradN - [
            -0.10750000000000001 -0.10750000000000001
            0.10750000000000001 -0.39249999999999996
            0.39249999999999996 0.39249999999999996
            -0.39249999999999996 0.10750000000000001
        ],
    ) <= 1.0e-6
end
end
using .mgrad2mm
mgrad2mm.test()

module mgrad3mm
using FinEtools
using FinEtools.FESetModule: gradN!
using LinearAlgebra
using Test
function test()
    rho = 1.21 * 1e-9# mass density
    c = 345.0 * 1000# millimeters per second
    bulk = c^2 * rho
    Lx = 1900.0# length of the box, millimeters
    Ly = 800.0 # length of the box, millimeters
    Lz = 800.0 # length of the box, millimeters

    fens, fes = H8block(Lx, Ly, Lz, 3, 2, 4) # Mesh
    # show(fes.conn)
    gradN, gradNparams, redJ = fill(0.0, nodesperelem(fes), 3),
    bfundpar(fes, [0.57, 0.57, -0.57]),
    Matrix(1.0 * I, 3, 3)

    gradN!(fes, gradN::FFltMat, gradNparams::FFltMat, redJ::FFltMat)
    # @show gradN
    @test norm(
        gradN - [
            -0.0843875 -0.0843875 -0.023112500000000005
            0.0843875 -0.30811249999999996 -0.0843875
            0.30811249999999996 0.30811249999999996 -0.30811249999999996
            -0.30811249999999996 0.0843875 -0.0843875
            -0.023112500000000005 -0.023112500000000005 0.023112500000000005
            0.023112500000000005 -0.0843875 0.0843875
            0.0843875 0.0843875 0.30811249999999996
            -0.0843875 0.023112500000000005 0.0843875
        ],
    ) <= 1.0e-6
end
end
using .mgrad3mm
mgrad3mm.test()

module mgffip1
using FinEtools
using LinearAlgebra
using FinEtools.MeshExportModule: VTK
using FinEtools.AlgoBaseModule: evalconvergencestudy
import FinEtools.FEMMBaseModule: inspectintegpoints
using FinEtools.MatrixUtilityModule: locjac!
using Test

Lx = 1900.0# length of the box, millimeters
Ly = 800.0 # length of the box, millimeters
nx, ny = 10, 9
th = 7.0

f(x) = cos(0.93 * pi * x[1] / Lx) + sin(1.7 * pi * x[2] / Ly)
g(x) = -cos(0.93 * pi * x[1] / Lx) + 2 * sin(1.7 * pi * x[2] / Ly)

function test()

    fens, fes = Q4block(Lx, Ly, nx, ny) # Mesh
    femm = FEMMBase(IntegDomain(fes, GaussRule(2, 2), th))
    psi = GeneralField(fill(0.0, count(fens), 2))
    for i in eachindex(fens)
        psi.values[i, 1] = f(fens.xyz[i, :])
        psi.values[i, 2] = g(fens.xyz[i, :])
    end
    numberdofs!(psi)
    @test nfreedofs(psi) == 220

    psi = GeneralField(fill(0.0, count(fens)))
    for i in eachindex(fens)
        psi.values[i, 1] = f(fens.xyz[i, :])
    end
    numberdofs!(psi)
    @test nfreedofs(psi) == 110
    true
end
end
using .mgffip1
mgffip1.test()

module msurfn1
using FinEtools
using LinearAlgebra
using FinEtools.MeshExportModule: VTK
using Test

function test()
    ndimensions = 2
    tangents, fe_label = reshape([1.0; 1.0], 2, 1), 0
    n = SurfaceNormal(ndimensions::FInt)
    normal = updatenormal!(n, [0.0 0.0 0.0], tangents::FFltMat, fe_label::FInt)
    @test norm(normal - [0.7071067811865475, -0.7071067811865475]) <= 1.0e-5

    ndimensions = 3
    tangents, fe_label = reshape([1.0 0.0; 1.0 0.0; 0.0 1.0], 3, 2), 0
    n = SurfaceNormal(ndimensions::FInt)
    normal = updatenormal!(n, [0.0 0.0 0.0], tangents::FFltMat, fe_label::FInt)
    @test norm(normal - [0.7071067811865475, -0.7071067811865475, 0.0]) <= 1.0e-5
    true
end
end
using .msurfn1
msurfn1.test()

module mnsimpr1
using FinEtools
using FinEtools.MeshExportModule
using Test
function test()
    xs = collect(linearspace(0.0, pi / 2, 5))
    fens, fes = L2blockx(xs)
    @test count(fes) == 4

    geom = NodalField(fens.xyz)
    femm = FEMMBase(IntegDomain(fes, NodalSimplexRule(3)))
    V1 = integratefunction(femm, geom, (x) -> 1.0)
    femm = FEMMBase(IntegDomain(fes, SimplexRule(3, 4)))
    V2 = integratefunction(femm, geom, (x) -> 1.0)
    @test abs(V1 - V2) / V1 <= 1.0e-5
    true
end
end
using .mnsimpr1
mnsimpr1.test()

module mnsimpr2
using FinEtools
using FinEtools.MeshExportModule
using Test
function test()
    xs = collect(linearspace(0.0, pi / 2, 5))
    ys = collect(linearspace(0.0, 1.0, 6) .^ 2)
    fens, fes = T3blockx(xs, ys, :a)
    for i in eachindex(fens)
        a, y = fens.xyz[i, :]
        fens.xyz[i, 1] = sin(a) * (y + 0.5)
        fens.xyz[i, 2] = cos(a) * (y + 0.5)
    end
    @test count(fes) == 40

    geom = NodalField(fens.xyz)
    femm = FEMMBase(IntegDomain(fes, NodalSimplexRule(3)))
    V1 = integratefunction(femm, geom, (x) -> 1.0)
    femm = FEMMBase(IntegDomain(fes, SimplexRule(3, 4)))
    V2 = integratefunction(femm, geom, (x) -> 1.0)
    @test abs(V1 - V2) / V1 <= 1.0e-5
    true
end
end
using .mnsimpr2
mnsimpr2.test()

module mnsimpr3
using FinEtools
using FinEtools.MeshExportModule
using Test
function test()
    xs = collect(linearspace(0.0, pi / 2, 5))
    ys = collect(linearspace(0.0, 1.0, 6) .^ 2)
    zs = collect(linearspace(0.0, 1.0, 3))
    fens, fes = T4blockx(xs, ys, zs, :a)
    for i in eachindex(fens)
        a, y, z = fens.xyz[i, :]
        fens.xyz[i, 1] = sin(a) * (y + 0.5)
        fens.xyz[i, 2] = cos(a) * (y + 0.5)
        fens.xyz[i, 3] = z
    end
    @test count(fes) == 240

    geom = NodalField(fens.xyz)
    femm = FEMMBase(IntegDomain(fes, NodalSimplexRule(3)))
    V1 = integratefunction(femm, geom, (x) -> 1.0)
    femm = FEMMBase(IntegDomain(fes, SimplexRule(3, 4)))
    V2 = integratefunction(femm, geom, (x) -> 1.0)
    @test abs(V1 - V2) / V1 <= 1.0e-5
    true
end
end
using .mnsimpr3
mnsimpr3.test()

module mntpr1
using FinEtools
using FinEtools.MeshExportModule
using Test
function test()
    xs = collect(linearspace(0.0, pi / 2, 5))
    fens, fes = L2blockx(xs)
    @test count(fes) == 4

    geom = NodalField(fens.xyz)
    femm = FEMMBase(IntegDomain(fes, NodalTensorProductRule(1)))
    V1 = integratefunction(femm, geom, (x) -> 1.0)
    femm = FEMMBase(IntegDomain(fes, TrapezoidalRule(1)))
    V2 = integratefunction(femm, geom, (x) -> 1.0)
    @test abs(V1 - V2) / V1 <= 1.0e-5
    true
end
end
using .mntpr1
mntpr1.test()

module mntpr2
using FinEtools
using FinEtools.MeshExportModule
using Test
function test()
    xs = collect(linearspace(0.0, pi / 2, 5))
    ys = collect(linearspace(0.0, 1.0, 6) .^ 2)
    fens, fes = Q4blockx(xs, ys)
    for i in eachindex(fens)
        a, y = fens.xyz[i, :]
        fens.xyz[i, 1] = sin(a) * (y + 0.5)
        fens.xyz[i, 2] = cos(a) * (y + 0.5)
    end
    @test count(fes) == 20

    geom = NodalField(fens.xyz)
    femm = FEMMBase(IntegDomain(fes, NodalTensorProductRule(2)))
    V1 = integratefunction(femm, geom, (x) -> 1.0)
    femm = FEMMBase(IntegDomain(fes, TrapezoidalRule(2)))
    V2 = integratefunction(femm, geom, (x) -> 1.0)
    @test abs(V1 - V2) / V1 <= 1.0e-5
    true
end
end
using .mntpr2
mntpr2.test()

module mntpr3
using FinEtools
using FinEtools.MeshExportModule
using Test
function test()
    xs = collect(linearspace(0.0, pi / 2, 5))
    ys = collect(linearspace(0.0, 1.0, 6) .^ 2)
    zs = collect(linearspace(0.0, 1.0, 3))
    fens, fes = H8blockx(xs, ys, zs)
    for i in eachindex(fens)
        a, y, z = fens.xyz[i, :]
        fens.xyz[i, 1] = sin(a) * (y + 0.5)
        fens.xyz[i, 2] = cos(a) * (y + 0.5)
        fens.xyz[i, 3] = z
    end
    @test count(fes) == 40

    geom = NodalField(fens.xyz)
    femm = FEMMBase(IntegDomain(fes, NodalTensorProductRule(3)))
    V1 = integratefunction(femm, geom, (x) -> 1.0)
    femm = FEMMBase(IntegDomain(fes, TrapezoidalRule(3)))
    V2 = integratefunction(femm, geom, (x) -> 1.0)
    @test abs(V1 - V2) / V1 <= 1.0e-5
    true
end
end
using .mntpr3
mntpr3.test()

module mr3m1
using FinEtools
using LinearAlgebra
using Test
function test()
    a = rand(3)
    R = rotmat3(a)
    R1 = fill(0.0, 3, 3)
    rotmat3!(R1, a)
    @test norm(R - R1) <= 1.0e-9
    @test norm(R * R1' - I) <= 1.0e-9
end
end
using .mr3m1
mr3m1.test()

module MMATD1
using FinEtools
using LinearAlgebra
using Test
mutable struct Mat <: AbstractMat
    mass_density::Float64
end

function test()
    m = Mat(133.0)
    @test massdensity(m) == 133.0
end
end
using .MMATD1
MMATD1.test()

module mmassembl2
using FinEtools
using Test
import LinearAlgebra: norm, cholesky
function test()
    refa = zeros(7,7)
    a = SysmatAssemblerSparse(0.0)
    startassembly!(a, 5*5*3, 7, 7)
    m = [
        0.24406 0.599773 0.833404 0.0420141
        0.786024 0.00206713 0.995379 0.780298
        0.845816 0.198459 0.355149 0.224996
    ]
    assemble!(a, m, [1 7 5], [5 2 1 4])
    refa[vec([1 7 5]), vec([5 2 1 4])] += m
    m = [
        0.146618 0.53471 0.614342 0.737833
        0.479719 0.41354 0.00760941 0.836455
        0.254868 0.476189 0.460794 0.00919633
        0.159064 0.261821 0.317078 0.77646
        0.643538 0.429817 0.59788 0.958909
    ]
    assemble!(a, m, [2 3 1 4 5], [6 7 3 4])
    refa[vec([2 3 1 4 5]), vec([6 7 3 4])] += m
    A = makematrix!(a)
    # @show Matrix(A)
    @test abs.(
        maximum(
            refa - A,
        )
    ) < 1.0e-5
    # @test abs(maximum(T_i)-1380.5883006341187) < 1.0e-3
end
end
using .mmassembl2
mmassembl2.test()


module mbas010x1
using FinEtools
using FinEtools.AlgoBaseModule: fieldnorm
using Test
function test()
    a, b, na, nb = 1.3, 3.1, 3, 2
    fens, fes = Q8block(a, b, na, nb) # Mesh
    geom = NodalField(fens.xyz)
    femm = FEMMBase(IntegDomain(fes, GaussRule(2, 3), true))

    V = integratefunction(femm, geom, (x) -> 1.0; m = 3)
    @test abs(V - pi * a^2 * b) / V <= eps(V)
    true
end
end
using .mbas010x1
mbas010x1.test()


module mbas010x2
using FinEtools
using FinEtools.AlgoBaseModule: fieldnorm
using Test
function test()
    a, b, na, nb = 1.3, 3.1, 3, 2
    fens, fes = Q8block(a, b, na, nb) # Mesh
    geom = NodalField(fens.xyz)
    bfes = meshboundary(fes)
    femm = FEMMBase(IntegDomain(bfes, GaussRule(1, 3), true))

    S = integratefunction(femm, geom, (x) -> 1.0; m = 2)
    @test abs(S - (2 * pi * a * b + 2 * pi * a^2)) / S <= eps(S)
    true
end
end
using .mbas010x2
mbas010x2.test()


module mbas010x3
using FinEtools
using FinEtools.AlgoBaseModule: fieldnorm
using Test
function test()
    a, b, na, nb = 1.3, 3.1, 3, 2
    fens, fes = T3block(a, b, na, nb) # Mesh
    geom = NodalField(fens.xyz)
    femm = FEMMBase(IntegDomain(fes, NodalSimplexRule(2), true))

    V = integratefunction(femm, geom, (x) -> 1.0; m = 3)
    @test abs(V - pi * a^2 * b) / V <= eps(V)
    true
end
end
using .mbas010x3
mbas010x3.test()

module mbas010x4
using FinEtools
using Test
function test()
    a, b, na, nb = 1.3, 3.1, 3, 2
    fens, fes = T3block(a, b, na, nb) # Mesh
    geom = NodalField(fens.xyz)
    bfes = meshboundary(fes)
    femm = FEMMBase(IntegDomain(bfes, NodalSimplexRule(1), true))

    S = integratefunction(femm, geom, (x) -> 1.0; m = 2)
    @test abs(S - (2 * pi * a * b + 2 * pi * a^2)) / S <= eps(S)
    true
end
end
using .mbas010x4
mbas010x4.test()

module massnomat1
using SparseArrays: sparse, spzeros
using FinEtools
using FinEtools.MeshImportModule
using FinEtools.MeshExportModule.MESH
using Test
function test()
    a = SysmatAssemblerSparse(0.0, true)
    startassembly!(a, 5*5*3, 7, 7)
    m = [
        0.24406 0.599773 0.833404 0.0420141
        0.786024 0.00206713 0.995379 0.780298
        0.845816 0.198459 0.355149 0.224996
    ]
    assemble!(a, m, [1 7 5], [5 2 1 4])
    m = [
        0.146618 0.53471 0.614342 0.737833
        0.479719 0.41354 0.00760941 0.836455
        0.254868 0.476189 0.460794 0.00919633
        0.159064 0.261821 0.317078 0.77646
        0.643538 0.429817 0.59788 0.958909
    ]
    assemble!(a, m, [2 3 1 7 5], [6 7 3 4])
    A = makematrix!(a)
    @test A == spzeros(7, 7)
    a.nomatrixresult = false
    A = makematrix!(a)
    @test A[1, 1] ≈ 0.833404
    @test A[5, 1] ≈ 0.355149
    @test A[7, 6] ≈ 0.159064
    @test A[3, 7] ≈ 0.41354
    @test A[7, 7] ≈ 0.261821
    true
end
end
using .massnomat1
massnomat1.test()


module mwithnodes4
using FinEtools
using Test
function test()
    fens, fes = Q4block(1.3, 3.1, 3, 2) # Mesh
    l = selectelem(fens, fes, withnodes = [3, 4, 7, 8, 11, 12])
    @test l == [5, 6]
    true
end
end
using .mwithnodes4
mwithnodes4.test()


module mutil_th001

using Test
using LinearAlgebra
using FinEtools
using FinEtools.AssemblyModule
using ChunkSplitters

function _test()
    # @show Base.Threads.nthreads()
    m = [
    0.24406   0.599773    0.833404  0.0420141
    0.786024  0.00206713  0.995379  0.780298
    0.845816  0.198459    0.355149  0.224996
    ]
    m1 = m'*m; r1 = [5 2 1 4]; c1 = r1
    m = [
    0.146618  0.53471   0.614342    0.737833
    0.479719  0.41354   0.00760941  0.836455
    0.254868  0.476189  0.460794    0.00919633
    0.159064  0.261821  0.317078    0.77646
    0.643538  0.429817  0.59788     0.958909
    ]
    m2 = m'*m; r2 = [2 3 1 5]; c2 = r2
    m = [
    -0.146618  0.53471   0.614342    0.737833
    0.479719  0.41354   0.00760941  -0.836455
    0.254868  -0.476189  0.460794    0.00919633
    0.159064  0.261821  0.317078    0.77646
    ]
    m3 = m; r3 = [2 3 1 5]; c3 = [7 6 2 4]
    m4 = Matrix(m'); r4 = [7 6 2 4]; c4 = [2 3 1 5];
    assembly_line = [
    (m1, r1, c1),
    (m1, r1, c1),
    (m2, r2, c2),
    (m2, r2, c2),
    (m2, r2, c2),
    (m2, r2, c2),
    (m2, r2, c2),
    (m2, r2, c2),
    (m2, r2, c2),
    (m2, r2, c2),
    (m2, r2, c2),
    (m2, r2, c2),
    (m2, r2, c2),
    (m2, r2, c2),
    (m3, r3, c3),
    (m3, r3, c3),
    (m3, r3, c3),
    (m3, r3, c3),
    (m3, r3, c3),
    (m3, r3, c3),
    (m4, r4, c4),
    (m4, r4, c4),
    (m4, r4, c4),
    (m4, r4, c4),
    (m4, r4, c4),
    (m4, r4, c4),
    (m1, r1, c1),
    (m1, r1, c1),
    (m1, r1, c1),
    (m1, r1, c1),
    (m1, r1, c1),
    (m1, r1, c1),
    (m2, r2, c2),
    ]
    elem_mat_nrows = 5
    elem_mat_ncols = 5
    elem_mat_nmatrices = length(assembly_line)
    ndofs_row = 7
    ndofs_col = 7
    a = SysmatAssemblerSparse(0.0)
    startassembly!(a, elem_mat_nrows*elem_mat_ncols*elem_mat_nmatrices, ndofs_row, ndofs_col)
    for i in eachindex(assembly_line)
        assemble!(a, assembly_line[i]...)
    end
    refA = makematrix!(a)

    # @show a

    nth = Base.Threads.nthreads()

    # Serial execution
    a = SysmatAssemblerSparse(0.0)
    startassembly!(a, elem_mat_nrows*elem_mat_ncols*elem_mat_nmatrices, ndofs_row, ndofs_col)
    ntasks = Base.Threads.nthreads()
    istart = 1; iend = 0;
    for ch in chunks(1:length(assembly_line), ntasks)
        # @show ch
        buffer_length = 5 * 5 * length(ch[1])
        iend = iend + buffer_length
        matbuffer = view(a.matbuffer, istart:iend)
        rowbuffer = view(a.rowbuffer, istart:iend)
        colbuffer = view(a.colbuffer, istart:iend)
        # @show length(colbuffer), istart, iend, buffer_length
        buffer_pointer = 1
        a1 = SysmatAssemblerSparse(buffer_length, matbuffer, rowbuffer, colbuffer, buffer_pointer, ndofs_row, ndofs_col, true, false)
        for i in ch[1]
            assemble!(a1, assembly_line[i]...)
        end
        makematrix!(a1)
        istart = iend + 1
    end
    a.buffer_pointer = iend
    A = makematrix!(a)
    @test norm(A - refA) / norm(refA) < 1.0e-9

    # Parallel execution
    a = SysmatAssemblerSparse(0.0)
    startassembly!(a, elem_mat_nrows*elem_mat_ncols*elem_mat_nmatrices, ndofs_row, ndofs_col)
    ntasks = Base.Threads.nthreads()
    istart = 1; iend = 0;
    Threads.@sync begin
        for ch in chunks(1:length(assembly_line), ntasks)
            # @show ch[2], ch[1]
            buffer_length = 5 * 5 * length(ch[1])
            iend = iend + buffer_length
            matbuffer = view(a.matbuffer, istart:iend)
            rowbuffer = view(a.rowbuffer, istart:iend)
            colbuffer = view(a.colbuffer, istart:iend)
            matbuffer .= 0.0
            rowbuffer .= 1
            colbuffer .= 1
            buffer_pointer = 1
            a1 = SysmatAssemblerSparse(buffer_length, matbuffer, rowbuffer, colbuffer, buffer_pointer, ndofs_row, ndofs_col, true, false)
            Threads.@spawn let r =  $ch[1]
                for i in r
                    assemble!(a1, assembly_line[i]...)
                end
                makematrix!(a1)
            end
            istart = iend + 1
        end
    end
    a.buffer_pointer = iend
    A = makematrix!(a)
    @test norm(A - refA) / norm(refA) < 1.0e-9
    return true
end

_test()

end # module

module mutil_th002

using Test
using LinearAlgebra
using FinEtools
using FinEtools.AssemblyModule
using ChunkSplitters
using Random

function _test()
    Random.seed!(1234);
        # @show Base.Threads.nthreads()
    m = [
    0.24406   0.599773    0.833404  0.0420141
    0.786024  0.00206713  0.995379  0.780298
    0.845816  0.198459    0.355149  0.224996
    ]
    m1 = m'*m;
    m = [
    0.146618  -0.53471   0.614342    0.737833
    0.479719  0.41354   -0.00760941  0.836455
    0.254868  0.476189  0.460794    0.00919633
    0.159064  0.261821  0.317078    0.77646
    -0.643538  0.429817  0.59788     0.958909
    ]
    m2 = m'*m;
    m = [
    -0.146618  0.53471   0.614342    0.737833
    0.479719  0.41354   0.00760941  -0.836455
    0.479719  0.41354   0.00760941  -0.836455
    0.254868  -0.476189  0.460794    0.00919633
    0.159064  0.261821  0.317078    0.77646
    0.159064  0.261821  0.317078    0.77646
    ]
    m3 = m;
    m4 = Matrix(m');
    ms = [m1, m2, m3, m4, ]
    elem_mat_nrows = maximum(size.(ms, 1))
    elem_mat_ncols = maximum(size.(ms, 2))
    N = 56567
    p = randperm(N)
    assembly_line = []
    for i in 1:N
        for k in eachindex(ms)
            if length(p) < sum(size(m))
                p = randperm(N)
            end
            m = ms[k]
            r = [popfirst!(p) for _ in 1:size(m, 1)]
            c = [popfirst!(p) for _ in 1:size(m, 2)]
            push!(assembly_line, (m, r, c))
        end
    end

    elem_mat_nmatrices = length(assembly_line)
    ndofs_row = N
    ndofs_col = N
    # @info "Serial execution"
    start = time()
    a = SysmatAssemblerSparse(0.0)
    startassembly!(a, elem_mat_nrows*elem_mat_ncols*elem_mat_nmatrices, ndofs_row, ndofs_col)
    for i in eachindex(assembly_line)
        assemble!(a, assembly_line[i]...)
    end
    refA = makematrix!(a)
    # @show time() - start

    # @show a

    nth = Base.Threads.nthreads()

    # Serial execution
    # @info "Serial chunked execution"
    start = time()
    a = SysmatAssemblerSparse(0.0)
    startassembly!(a, elem_mat_nrows*elem_mat_ncols*elem_mat_nmatrices, ndofs_row, ndofs_col)
    ntasks = Base.Threads.nthreads()
    istart = 1; iend = 0;
    for ch in chunks(1:length(assembly_line), ntasks)
        # @show ch[2], ch[1]
        buffer_length = elem_mat_nrows * elem_mat_ncols * length(ch[1])
        iend = iend + buffer_length
        matbuffer = view(a.matbuffer, istart:iend)
        rowbuffer = view(a.rowbuffer, istart:iend)
        colbuffer = view(a.colbuffer, istart:iend)
        buffer_pointer = 1
        a1 = SysmatAssemblerSparse(buffer_length, matbuffer, rowbuffer, colbuffer, buffer_pointer, ndofs_row, ndofs_col, true, false)
        for i in ch[1]
            assemble!(a1, assembly_line[i]...)
        end
        makematrix!(a1)
        # @show "done $(ch[2])"
        istart = iend + 1
    end
    a.buffer_pointer = iend
    A = makematrix!(a)
    # @show time() - start
    @test norm(A - refA) / norm(refA) < 1.0e-9

    # Parallel execution
    # @info "Parallel execution"
    start = time()
    a = SysmatAssemblerSparse(0.0)
    startassembly!(a, elem_mat_nrows*elem_mat_ncols*elem_mat_nmatrices, ndofs_row, ndofs_col)
    ntasks = Base.Threads.nthreads()
    istart = 1; iend = 0;
    Threads.@sync begin
        for ch in chunks(1:length(assembly_line), ntasks)
            # @show ch[2], ch[1]
            buffer_length = elem_mat_nrows * elem_mat_ncols * length(ch[1])
            iend = iend + buffer_length
            matbuffer = view(a.matbuffer, istart:iend)
            rowbuffer = view(a.rowbuffer, istart:iend)
            colbuffer = view(a.colbuffer, istart:iend)
            buffer_pointer = 1
            Threads.@spawn let r =  $ch[1]
                a1 = SysmatAssemblerSparse(buffer_length, $matbuffer, $rowbuffer, $colbuffer, $buffer_pointer, ndofs_row, ndofs_col, true, false)
                # @show ch[2], r
                for i in r
                    assemble!(a1, assembly_line[i]...)
                end
                makematrix!(a1)
            end
            # @show "done $(ch[2])"
            istart = iend + 1
        end
    end
    a.buffer_pointer = iend
    A = makematrix!(a)
    # @show time() - start
    @test norm(A - refA) / norm(refA) < 1.0e-9
    return true
end

_test()

end # module

module mutil_th003

using Test
using LinearAlgebra
using FinEtools
using FinEtools.AssemblyModule
using ChunkSplitters
using Random

function _test()
    Random.seed!(1234);
        # @show Base.Threads.nthreads()
    m = [
    0.24406   0.599773    0.833404  0.0420141
    0.786024  0.00206713  0.995379  0.780298
    0.845816  0.198459    0.355149  0.224996
    ]
    m1 = m'*m;
    m = [
    0.146618  -0.53471   0.614342    0.737833
    0.479719  0.41354   -0.00760941  0.836455
    0.254868  0.476189  0.460794    0.00919633
    0.159064  0.261821  0.317078    0.77646
    -0.643538  0.429817  0.59788     0.958909
    ]
    m2 = m'*m;
    m = [
    -0.146618  0.53471   0.614342    0.737833
    0.479719  0.41354   0.00760941  -0.836455
    0.479719  0.41354   0.00760941  -0.836455
    0.254868  -0.476189  0.460794    0.00919633
    0.159064  0.261821  0.317078    0.77646
    0.159064  0.261821  0.317078    0.77646
    ]
    m3 = m;
    m4 = Matrix(m');
    ms = [m1, m2, m3, m4, ]
    elem_mat_nrows = maximum(size.(ms, 1))
    elem_mat_ncols = maximum(size.(ms, 2))
    N = 56567
    p = randperm(N)
    assembly_line = []
    for i in 1:N
        for k in eachindex(ms)
            if length(p) < sum(size(m))
                p = randperm(N)
            end
            m = ms[k]
            r = [popfirst!(p) for _ in 1:size(m, 1)]
            c = [popfirst!(p) for _ in 1:size(m, 2)]
            push!(assembly_line, (m, r, c))
        end
    end

    elem_mat_nmatrices = length(assembly_line)
    ndofs_row = N
    ndofs_col = N
    # @info "Serial execution"
    start = time()
    a = SysmatAssemblerSparse(0.0)
    startassembly!(a, elem_mat_nrows*elem_mat_ncols*elem_mat_nmatrices, ndofs_row, ndofs_col)
    for i in eachindex(assembly_line)
        assemble!(a, assembly_line[i]...)
    end
    refA = makematrix!(a)
    # @show time() - start

    # @show a

    nth = Base.Threads.nthreads()

    # Serial execution
    # @info "Serial chunked execution"
    start = time()
    a = SysmatAssemblerSparse(0.0)
    startassembly!(a, elem_mat_nrows*elem_mat_ncols*elem_mat_nmatrices, ndofs_row, ndofs_col)
    ntasks = Base.Threads.nthreads()
    istart = 1; iend = 0;
    for ch in chunks(1:length(assembly_line), ntasks)
        # @show ch[2], ch[1]
        buffer_length = elem_mat_nrows * elem_mat_ncols * length(ch[1])
        iend = iend + buffer_length
        matbuffer = view(a.matbuffer, istart:iend)
        rowbuffer = view(a.rowbuffer, istart:iend)
        colbuffer = view(a.colbuffer, istart:iend)
        matbuffer .= 0.0
        rowbuffer .= 1
        colbuffer .= 1
        buffer_pointer = 1
        a1 = SysmatAssemblerSparse(buffer_length, matbuffer, rowbuffer, colbuffer, buffer_pointer, ndofs_row, ndofs_col, true, false)
        for i in ch[1]
            assemble!(a1, assembly_line[i]...)
        end
        # @show "done $(ch[2]), $(istart), $(iend)"
        istart = iend + 1
    end
    a.buffer_pointer = iend
    A = makematrix!(a)
    # @show time() - start
    @test norm(A - refA) / norm(refA) < 1.0e-9

    function _update_buffer_range(elem_mat_nrows, elem_mat_ncols, range, iend)
        buffer_length = elem_mat_nrows * elem_mat_ncols * length(range)
        istart = iend + 1
        iend = iend + buffer_length
        buffer_range = istart:iend
        return buffer_range, iend
    end

    function _task_local_assembler(a, buffer_range)
        buffer_length = maximum(buffer_range) - minimum(buffer_range) + 1
        matbuffer = view(a.matbuffer, buffer_range)
        rowbuffer = view(a.rowbuffer, buffer_range)
        colbuffer = view(a.colbuffer, buffer_range)
        buffer_pointer = 1
        matbuffer .= 0.0
        rowbuffer .= 1
        colbuffer .= 1
        a1 = SysmatAssemblerSparse(buffer_length, matbuffer, rowbuffer, colbuffer, buffer_pointer, ndofs_row, ndofs_col, true, false)
    end

    # Parallel execution
    # @info "Parallel execution"
    start = time()
    a = SysmatAssemblerSparse(0.0)
    startassembly!(a, elem_mat_nrows*elem_mat_ncols*elem_mat_nmatrices, ndofs_row, ndofs_col)
    ntasks = Base.Threads.nthreads()
    iend = 0;
    Threads.@sync begin
        for ch in chunks(1:length(assembly_line), ntasks)
            # @show ch[2], ch[1]
            buffer_range, iend = _update_buffer_range(elem_mat_nrows, elem_mat_ncols, ch[1], iend)
            Threads.@spawn let r =  $ch[1]
                a1 = _task_local_assembler(a, $buffer_range)
                for i in r
                    assemble!(a1, assembly_line[i]...)
                end
            end
            # @show "done $(ch[2]), $(iend)"
        end
    end
    a.buffer_pointer = iend
    A = makematrix!(a)
    # @show time() - start
    @test norm(A - refA) / norm(refA) < 1.0e-9
    return true
end

_test()

end # module

module mbasicass001

using Test
using LinearAlgebra
using FinEtools
using FinEtools.AssemblyModule
using SparseArrays
using Random

function _test()
    Random.seed!(1234);
        # @show Base.Threads.nthreads()
    m = [
    0.24406   0.599773    0.833404  0.0420141
    0.786024  0.00206713  0.995379  0.780298
    0.845816  0.198459    0.355149  0.224996
    ]
    m1 = m'*m;
    m5 = m*m';
    m = [
    0.146618  -0.53471   0.614342    0.737833
    0.479719  0.41354   -0.00760941  0.836455
    0.254868  0.476189  0.460794    0.00919633
    0.159064  0.261821  0.317078    0.77646
    -0.643538  0.429817  0.59788     0.958909
    ]
    m2 = m'*m;
    m = [
    -0.146618  0.53471   0.614342    0.737833
    0.479719  0.41354   0.00760941  -0.836455
    0.479719  0.41354   0.00760941  -0.836455
    0.254868  -0.476189  0.460794    0.00919633
    0.159064  0.261821  0.317078    0.77646
    0.159064  0.261821  0.317078    0.77646
    ]
    m3 = m'*m;
    m4 = Matrix(m*m');
    ms = [m1, m2, m3, m4, m5, ]
    # Testing symmetric matrices
    for i in eachindex(ms)
        m = ms[i][1]
        @assert norm(m - m') == 0
    end
    elem_mat_nrows = maximum(size.(ms, 1))
    elem_mat_ncols = maximum(size.(ms, 2))
    N = 133
    p = randperm(N)
    assembly_line = []
    for i in 1:N
        for k in eachindex(ms)
            m = ms[k]
            if length(p) <= sum(size(m))
                p = randperm(N)
            end
            r = [popfirst!(p) for _ in 1:size(m, 1)]
            c = [popfirst!(p) for _ in 1:size(m, 2)]
            # Symmetric matrices: the rows are the same as columns
            push!(assembly_line, (m, r, r))
            push!(assembly_line, (m, c, c))
        end
    end

    elem_mat_nmatrices = length(assembly_line)
    ndofs_row = N
    ndofs_col = N

    refA = spzeros(ndofs_row, ndofs_col)
    for i in eachindex(assembly_line)
        m, r, c = assembly_line[i]
        @assert length(r) == length(c)
        refA[r, c] += m
    end

    a = SysmatAssemblerSparse(0.0)
    startassembly!(a, elem_mat_nrows*elem_mat_ncols*elem_mat_nmatrices, ndofs_row, ndofs_col)
    for i in eachindex(assembly_line)
        assemble!(a, assembly_line[i]...)
    end
    A = makematrix!(a)
    # @show time() - start
    @test norm(A - refA) / norm(refA) < 1.0e-9


    a = SysmatAssemblerSparseSymm(0.0)
    startassembly!(a, elem_mat_nrows*elem_mat_ncols*elem_mat_nmatrices, ndofs_row, ndofs_col)
    for i in eachindex(assembly_line)
        assemble!(a, assembly_line[i]...)
    end
    A = makematrix!(a)

    # @show time() - start
    @test norm(A - refA) / norm(refA) < 1.0e-9

    return true
end

_test()

end # module

module mbasicass002
using Test
using LinearAlgebra
using FinEtools
using FinEtools.AssemblyModule
using SparseArrays
using Random

function _test()
    Random.seed!(1234);
        # @show Base.Threads.nthreads()
    m = [
    0.24406   0.599773    0.833404  0.0420141
    0.786024  0.00206713  0.995379  0.780298
    0.845816  0.198459    0.355149  0.224996
    ]
    m1 = m'*m;
    m5 = m*m';
    m = [
    0.146618  -0.53471   0.614342    0.737833
    0.479719  0.41354   -0.00760941  0.836455
    0.254868  0.476189  0.460794    0.00919633
    0.159064  0.261821  0.317078    0.77646
    -0.643538  0.429817  0.59788     0.958909
    ]
    m2 = m'*m;
    m = [
    -0.146618  0.53471   0.614342    0.737833
    0.479719  0.41354   0.00760941  -0.836455
    0.479719  0.41354   0.00760941  -0.836455
    0.254868  -0.476189  0.460794    0.00919633
    0.159064  0.261821  0.317078    0.77646
    0.159064  0.261821  0.317078    0.77646
    ]
    m3 = m'*m;
    m4 = Matrix(m*m');
    ms = [m1, m2, m3, m4, m5, ]
    # Testing symmetric matrices
    for i in eachindex(ms)
        m = ms[i][1]
        @assert norm(m - m') == 0
    end
    elem_mat_nrows = maximum(size.(ms, 1))
    elem_mat_ncols = maximum(size.(ms, 2))
    N = 1333
    p = randperm(N)
    assembly_line = []
    for i in 1:N
        for k in eachindex(ms)
            m = ms[k]
            if length(p) <= sum(size(m))
                p = randperm(N)
            end
            r = [popfirst!(p) for _ in 1:size(m, 1)]
            c = [popfirst!(p) for _ in 1:size(m, 2)]
            # Symmetric matrices: the rows are the same as columns
            push!(assembly_line, (m, r, r))
            push!(assembly_line, (m, c, c))
        end
    end

    elem_mat_nmatrices = length(assembly_line)
    ndofs_row = N
    ndofs_col = N

    refA = spzeros(ndofs_row, ndofs_col)
    for i in eachindex(assembly_line)
        m, r, c = assembly_line[i]
        @assert length(r) == length(c)
        refA[r, c] += m
    end

    a = SysmatAssemblerSparse(0.0)
    # We are testing resizing of buffers by under sizing initially
    startassembly!(a, elem_mat_nrows*elem_mat_ncols*elem_mat_nmatrices, ndofs_row, ndofs_col)
    for i in eachindex(assembly_line)
        assemble!(a, assembly_line[i]...)
    end
    A = makematrix!(a)
    # @show time() - start
    @test norm(A - refA) / norm(refA) < 1.0e-9


    a = SysmatAssemblerSparseSymm(0.0)
    # We are testing resizing of buffers by under sizing initially
    startassembly!(a, elem_mat_nrows*elem_mat_ncols*elem_mat_nmatrices, ndofs_row, ndofs_col)
    for i in eachindex(assembly_line)
        assemble!(a, assembly_line[i]...)
    end
    A = makematrix!(a)

    # @show time() - start
    @test norm(A - refA) / norm(refA) < 1.0e-9

    return true
end

_test()

end # module

module mbasicass003

using Test
using LinearAlgebra
using FinEtools
using FinEtools.AssemblyModule
using SparseArrays
using Random

function _test()
    Random.seed!(1234);
        # @show Base.Threads.nthreads()
    m = [
    0.24406   0.599773    0.833404  0.0420141
    0.786024  0.00206713  0.995379  0.780298
    0.845816  0.198459    0.355149  0.224996
    ]
    m1 = m'*m;
    m5 = m*m';
    m = [
    0.146618  -0.53471   0.614342    0.737833
    0.479719  0.41354   -0.00760941  0.836455
    0.254868  0.476189  0.460794    0.00919633
    0.159064  0.261821  0.317078    0.77646
    -0.643538  0.429817  0.59788     0.958909
    ]
    m2 = m'*m;
    m = [
    -0.146618  0.53471   0.614342    0.737833
    0.479719  0.41354   0.00760941  -0.836455
    0.479719  0.41354   0.00760941  -0.836455
    0.254868  -0.476189  0.460794    0.00919633
    0.159064  0.261821  0.317078    0.77646
    0.159064  0.261821  0.317078    0.77646
    ]
    m3 = m'*m;
    m4 = Matrix(m*m');
    ms = [m1, m2, m3, m4, m5, ]
    # Testing symmetric matrices
    for i in eachindex(ms)
        m = ms[i][1]
        @assert norm(m - m') == 0
    end
    elem_mat_nrows = maximum(size.(ms, 1))
    elem_mat_ncols = maximum(size.(ms, 2))
    N = 1333
    p = randperm(N)
    assembly_line = []
    for i in 1:N
        for k in eachindex(ms)
            m = ms[k]
            if length(p) <= sum(size(m))
                p = randperm(N)
            end
            r = [popfirst!(p) for _ in 1:size(m, 1)]
            c = [popfirst!(p) for _ in 1:size(m, 2)]
            # Symmetric matrices: the rows are the same as columns
            push!(assembly_line, (m, r, r))
            push!(assembly_line, (m, c, c))
        end
    end

    elem_mat_nmatrices = length(assembly_line)
    ndofs_row = N
    ndofs_col = N

    refA = spzeros(ndofs_row, ndofs_col)
    for i in eachindex(assembly_line)
        m, r, c = assembly_line[i]
        @assert length(r) == length(c)
        refA[r, c] += m
    end

    a = SysmatAssemblerSparse(0.0)
    startassembly!(a, elem_mat_nrows*elem_mat_ncols*110, ndofs_row, ndofs_col)
    for i in eachindex(assembly_line)
        assemble!(a, assembly_line[i]...)
    end
    A = makematrix!(a)
    # @show time() - start
    @test norm(A - refA) / norm(refA) < 1.0e-9
    # Here we test that we can start assembly again
    startassembly!(a, elem_mat_nrows*elem_mat_ncols*110, ndofs_row, ndofs_col)
    for i in eachindex(assembly_line)
        assemble!(a, assembly_line[i]...)
    end
    A = makematrix!(a)
    # @show time() - start
    @test norm(A - refA) / norm(refA) < 1.0e-9

    a = SysmatAssemblerSparseSymm(0.0)
    startassembly!(a, elem_mat_nrows*elem_mat_ncols*elem_mat_nmatrices, ndofs_row, ndofs_col)
    for i in eachindex(assembly_line)
        assemble!(a, assembly_line[i]...)
    end
    A = makematrix!(a)
    # @show time() - start
    @test norm(A - refA) / norm(refA) < 1.0e-9
    # Here we test that we can start assembly again
    startassembly!(a, elem_mat_nrows*elem_mat_ncols*elem_mat_nmatrices, ndofs_row, ndofs_col)
    for i in eachindex(assembly_line)
        assemble!(a, assembly_line[i]...)
    end
    A = makematrix!(a)
    # @show time() - start
    @test norm(A - refA) / norm(refA) < 1.0e-9

    return true
end

_test()

end # module


module testing_cache_1
using Test
using LinearAlgebra
using FinEtools

XYZ, tangents, fe_label = (reshape([0.0, 0.0], 1, 2), [1.0 0.0; 0.0 1.0], 1)
nentries = 3
function fillcache!(cacheout::Vector{CT},
        XYZ::VecOrMat{T}, tangents::Matrix{T},
        fe_label) where {CT, T}
    cacheout .= 13
    return cacheout
end
c = DataCache(zeros(nentries), fillcache!)
v = c(XYZ, tangents, fe_label)
@test v == [13, 13, 13]
end


module testing_cache_2
using Test
using LinearAlgebra
using FinEtools

XYZ, tangents, fe_label = (reshape([0.0, 0.0], 1, 2), [1.0 0.0; 0.0 1.0], 1)
nentries = 3
function fillcache!(cacheout::Vector{CT},
        XYZ::VecOrMat{T}, tangents::Matrix{T},
        fe_label) where {CT, T}
    cacheout .= 13
    return cacheout
end
c = DataCache(zeros(Float32, nentries), fillcache!)
v = c(XYZ, tangents, fe_label)
@test v == [13, 13, 13]
end

module testing_cache_3
using Test
using LinearAlgebra
using FinEtools

XYZ, tangents, fe_label = (reshape([0.0, 0.0], 1, 2), [1.0 0.0; 0.0 1.0], 1)
nentries = 3
function fillcache!(cacheout::Array{CT, N},
        XYZ::VecOrMat{T}, tangents::Matrix{T},
        fe_label) where {CT, N, T}
    cacheout .= 13
    return cacheout
end
data = rand(3, 3)
c = DataCache(data)
v = c(XYZ, tangents, fe_label)
@test v == data
end


module testing_cache_4
using Test
using LinearAlgebra
using FinEtools

XYZ, tangents, fe_label = (reshape([0.0, 0.0], 1, 2), [1.0 0.0; 0.0 1.0], 1)
nentries = 3
function fillcache!(cacheout::Array{CT, N},
        XYZ::VecOrMat{T}, tangents::Matrix{T},
        fe_label) where {CT, N, T}
    cacheout .= 13
    return cacheout
end
data = rand(3, 3)
c = DataCache(data, fillcache!)
v = c(XYZ, tangents, fe_label)
data .= 13
@test v == data
end

module testing_cache_5
using Test
using LinearAlgebra
using FinEtools

XYZ, tangents, fe_label = (reshape([0.0, 0.0], 1, 2), [1.0 0.0; 0.0 1.0], 1)
nentries = 3
function fillcache!(cacheout::Array{CT, N},
        XYZ::VecOrMat{T}, tangents::Matrix{T},
        fe_label) where {CT, N, T}
    cacheout .= 13
    return cacheout
end
data = rand(3, 3)
c = DataCache(data, fillcache!)
v = c(XYZ, tangents, fe_label)
data .= 13
@test v == data
end

module testing_cache_6
using Test
using LinearAlgebra
using FinEtools

XYZ, tangents, fe_label = (reshape([0.0, 0.0], 1, 2), [1.0 0.0; 0.0 1.0], 1)
nentries = 3
function fillcache!(cacheout::Array{CT, N},
        XYZ::VecOrMat{T}, tangents::Matrix{T},
        fe_label) where {CT, N, T}
    cacheout .= fe_label
    return cacheout
end
for t in (Float32, Float64, ComplexF32, ComplexF64)
    data = rand(t, 3, 3)
    c = DataCache(data, fillcache!)
    v = c(XYZ, tangents, fe_label)
    data .= fe_label
    @test v == data
end
end


module testing_cache_7
using Test
using LinearAlgebra
using FinEtools

function fillcache!(cacheout::Array{CT, N},
        XYZ::VecOrMat{T}, tangents::Matrix{T},
        fe_label) where {CT, N, T}
    cacheout .= LinearAlgebra.I(3)
    return cacheout
end
c = DataCache(zeros(Float32, 3, 3), fillcache!)
function f(c)
    XYZ, tangents, fe_label = (reshape([0.0, 0.0], 1, 2), [1.0 0.0; 0.0 1.0], 1)
    data = c(XYZ, tangents, fe_label)
end
@test f(c) == LinearAlgebra.I(3)
end

module testing_cache_9
using Test
using LinearAlgebra
using FinEtools

function fillcache!(cacheout::Array{CT, N},
        XYZ::VecOrMat{T}, tangents::Matrix{T},
        fe_label, time) where {CT, N, T}
    if time > 1.0
        cacheout .= LinearAlgebra.I(3)
    else
        cacheout .= 0
    end
    return cacheout
end
t = Ref(0.0)

c = DataCache(zeros(Float32, 3, 3), (cacheout, XYZ, tangents, fe_label) -> fillcache!(cacheout, XYZ, tangents, fe_label, t[]))
function f(c)
    XYZ, tangents, fe_label = (reshape([0.0, 0.0], 1, 2), [1.0 0.0; 0.0 1.0], 1)
    t[] = 0.0
    data = c(XYZ, tangents, fe_label)
    @test data == zeros(Float32, 3, 3)
    t[] = 2.1
    data = c(XYZ, tangents, fe_label)
end
@test f(c) == LinearAlgebra.I(3)
end


module testing_cache_10
using Test
using LinearAlgebra
using FinEtools

XYZ, tangents, _ = (reshape([0.0, 0.0], 1, 2), [1.0 0.0; 0.0 1.0], 1)

function fillcache!(cacheout::D,
        XYZ::VecOrMat{T}, tangents::Matrix{T},
        fe_label) where {D, T}
    cacheout = D(fe_label)
    return cacheout
end
for t in (Float32, Float64, ComplexF32, ComplexF64)
    data = t(-42)
     c = DataCache(data, fillcache!)
    for fe_label in 1:5
         v = c(XYZ, tangents, fe_label)
        @test v == t(fe_label)
    end
end
end


module mmacousticcouplingpanelsm1
using FinEtools
using FinEtools.SurfaceNormalModule: SurfaceNormal, updatenormal!
using LinearAlgebra
using Test

function __computenormal!(normalout::FFltVec, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt, time = 0.0)
    fill!(normalout, 0.0)
    # We are assuming a surface element here!
    if (size(tangents,1) == 3) && (size(tangents,2) == 2)# surface in three dimensions
        normalout[:] .= cross(vec(tangents[:,1]), vec(tangents[:,2]));# outer normal to the surface
    else
        error("No definition of normal vector");
    end
    nn = norm(normalout);
    if  nn != 0.0 # otherwise return an unnormalized normal
        normalout ./= nn
    end
    return normalout
end

function test()
    n = SurfaceNormal(3, __computenormal!)
    XYZ, tangents, l = (reshape([0.0, 0.0, 0.0], 1, 3), [1.0 0.0; 0.0 1.0; 0.0 0.0], 1)
    updatenormal!(n, XYZ, tangents, l)
    true

end
test()
end
