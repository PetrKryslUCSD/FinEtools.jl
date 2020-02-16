
module mmiscellaneous1mmm
using FinEtools
using Test
function test()
  rho=1.21*1e-9;# mass density
  c =345.0*1000;# millimeters per second
  bulk= c^2*rho;
  Lx=1900.0;# length of the box, millimeters
  Ly=800.0; # length of the box, millimeters

  fens,fes = Q4block(Lx,Ly,3,2); # Mesh
  # show(fes.conn)
  length(fes.conn)

  bfes = meshboundary(fes)
  @test bfes.conn == Tuple{Int64,Int64}[(1, 2), (5, 1), (2, 3), (3, 4), (4, 8), (9, 5), (8, 12), (10, 9), (11, 10), (12, 11)]
end
end
using .mmiscellaneous1mmm
mmiscellaneous1mmm.test()

module mrichmmm
using FinEtools
using FinEtools.AlgoBaseModule
using Test
import LinearAlgebra: norm
function test()
  xs = [93.0734, 92.8633, 92.7252]
    hs = [0.1000, 0.0500, 0.0250]
  solnestim, beta, c, residual = AlgoBaseModule.richextrapol(xs, hs)
  # # println("$((solnestim, beta, c))")
  @test norm([solnestim, beta, c] - [92.46031652777476, 0.6053628424093497, -2.471055221256022]) < 1.e-5

  sol = [231.7, 239.1, 244.8]
  h = [(4.0/5.)^i for i in 0:1:2]
  solnestim, beta, c, residual = AlgoBaseModule.richextrapol(sol, h)
  # # println("$((solnestim, beta, c))")
  @test norm([solnestim, beta, c] - [263.91176470588067, 1.1697126080157385, 32.21176470588068]) < 1.e-5

end
end
using .mrichmmm
mrichmmm.test()

module mmmeasurementm1
using FinEtools
using Test
function test()
    W = 1.1;
    L = 12.;
    t =  0.32;
    nl, nt, nw = 2, 3, 4;

    fens,fes  = H8block(L,W,t, nl,nw,nt)
    geom  =  NodalField(fens.xyz)

    femm  =  FEMMBase(IntegDomain(fes, GaussRule(3, 2)))
    V = integratefunction(femm, geom, (x) ->  1.0)
    @test abs(V - W*L*t)/V < 1.0e-5
end
end
using .mmmeasurementm1
mmmeasurementm1.test()

module mmmeasurementm2
using FinEtools
using Test
function test()
    W = 1.1;
    L = 12.;
    t =  0.32;
    nl, nt, nw = 2, 3, 4;

    fens,fes  = H8block(L,W,t, nl,nw,nt)
    geom  =  NodalField(fens.xyz)
    bfes = meshboundary(fes)
    femm  =  FEMMBase(IntegDomain(bfes, GaussRule(2, 2)))
    S = integratefunction(femm, geom, (x) ->  1.0)
    @test abs(S - 2*(W*L + L*t + W*t))/S < 1.0e-5
end
end
using .mmmeasurementm2
mmmeasurementm2.test()


module mmmeasurementm3
using FinEtools
using Test
function test()
    W = 1.1;
    L = 12.;
    t =  0.32;
    nl, nt, nw = 2, 3, 4;

    fens,fes  = Q4block(L,W,nl,nw)
    geom  =  NodalField(fens.xyz)
    bfes = meshboundary(fes)
    femm  =  FEMMBase(IntegDomain(bfes, GaussRule(1, 2)))
    S = integratefunction(femm, geom, (x) ->  1.0)
    @test abs(S - 2*(W + L))/S < 1.0e-5
end
end
using .mmmeasurementm3
mmmeasurementm3.test()


module mmmeasurementm4
using FinEtools
using Test
function test()
    W = 1.1;
    L = 12.;
    t =  0.32;
    nl, nt, nw = 2, 3, 4;

    fens,fes  = L2block(L,nl)
    geom  =  NodalField(fens.xyz)
    bfes = meshboundary(fes)
    femm  =  FEMMBase(IntegDomain(bfes, PointRule()))
    S = integratefunction(femm, geom, (x) ->  1.0)
    @test abs(S - 2)/S < 1.0e-5
end
end
using .mmmeasurementm4
mmmeasurementm4.test()


module mmmeasurementm5
using FinEtools
using Test
function test()
    W = 1.1;
    L = 12.;
    t =  0.32;
    nl, nt, nw = 2, 3, 4;

    fens,fes  = Q4block(L,W,nl,nw)
    geom  =  NodalField(fens.xyz)

    nfes = FESetP1(reshape(collect(1:count(fens)), count(fens), 1))
    femm  =  FEMMBase(IntegDomain(nfes, PointRule()))
    S = integratefunction(femm, geom, (x) ->  1.0)
    @test abs(S - count(fens))/S < 1.0e-5

end
end
using .mmmeasurementm5
mmmeasurementm5.test()


module mmmeasurementm6
using FinEtools
using Test
function test()
    W = 1.1;
    L = 12.;
    t =  0.32;
    nl, nt, nw = 2, 3, 4;

    fens,fes  = H8block(L,W,t, nl,nw,nt)
    geom  =  NodalField(fens.xyz)

    nfes = FESetP1(reshape(collect(1:count(fens)), count(fens), 1))
    femm  =  FEMMBase(IntegDomain(nfes, PointRule()))
    S = integratefunction(femm, geom, (x) ->  1.0)
    @test abs(S - count(fens))/S < 1.0e-5

end
end
using .mmmeasurementm6
mmmeasurementm6.test()


module mmmeasurementm7
using FinEtools
using Test
function test()
    W = 1.1;
    L = 12.;
    t =  0.32;
    nl, nt, nw = 2, 3, 4;

    fens,fes  = L2block(L,nl)
    geom  =  NodalField(fens.xyz)
    bfes = FESetP1(reshape([nl+1], 1, 1))
    axisymmetric = true
    femm  =  FEMMBase(IntegDomain(bfes, PointRule(), axisymmetric))
    S = integratefunction(femm, geom, (x) ->  1.0, 1)
    # # println(" Length  of the circle = $(S)")
    @test abs(S - 2*pi*L)/S < 1.0e-5
end
end
using .mmmeasurementm7
mmmeasurementm7.test()

module mmmeasurementm8
using FinEtools
using Test
function test()
    W = 1.1;
    L = 12.;
    t =  0.32;
    nl, nt, nw = 2, 3, 4;

    fens,fes  = L2block(L,nl)
    geom  =  NodalField(fens.xyz)
    bfes = FESetP1(reshape([nl+1], 1, 1))

    femm  =  FEMMBase(IntegDomain(bfes, PointRule(), t))
    S = integratefunction(femm, geom, (x) ->  1.0, 1)
    # # println("Length  of the boundary curve = $(S)")
    @test abs(S - t)/S < 1.0e-5
end
end
using .mmmeasurementm8
mmmeasurementm8.test()

module mmmeasurementm9
using FinEtools
using Test
function test()
    W = 1.1;
    L = 12.;
    t =  0.32;
    nl, nt, nw = 2, 3, 4;

    fens,fes  = L2block(L,nl)
    geom  =  NodalField(fens.xyz)
    bfes = FESetP1(reshape([nl+1], 1, 1))

    femm  =  FEMMBase(IntegDomain(bfes, PointRule(), t*W))
    S = integratefunction(femm, geom, (x) ->  1.0, 2)
    # # println("Length  of the boundary curve = $(S)")
    @test abs(S - t*W)/S < 1.0e-5
end
end
using .mmmeasurementm9
mmmeasurementm9.test()

module mmmeasurementm10
using FinEtools
using Test
function test()
    W = 1.1;
    L = 12.;
    t =  0.32;
    nl, nt, nw = 2, 3, 4;

    fens,fes  = L2block(L,nl)
    geom  =  NodalField(fens.xyz)
    bfes = FESetP1(reshape([nl+1], 1, 1))
axisymmetric = true
    femm  =  FEMMBase(IntegDomain(bfes, PointRule(), axisymmetric, W))
    S = integratefunction(femm, geom, (x) ->  1.0, 2)
    # # println("Length  of the boundary curve = $(S)")
    @test abs(S - 2*pi*L*W)/S < 1.0e-5
end
end
using .mmmeasurementm10
mmmeasurementm10.test()

module mmmeasurementm11
using FinEtools
using Test
function test()
    W = 1.1;
    L = 12.;
    t =  0.32;
    nl, nt, nw = 2, 3, 4;

    fens,fes  = L2block(L,nl)
    geom  =  NodalField(fens.xyz)
    bfes = FESetP1(reshape([nl+1], 1, 1))
axisymmetric = true
    femm  =  FEMMBase(IntegDomain(bfes, PointRule(), axisymmetric, W*t))
    S = integratefunction(femm, geom, (x) ->  1.0, 3)
    # # println("Length  of the boundary curve = $(S)")
    @test abs(S - 2*pi*L*W*t)/S < 1.0e-5
end
end
using .mmmeasurementm11
mmmeasurementm11.test()

module mmmeasurementm12
using FinEtools
using Test
function test()
    W = 1.1;
    L = 12.;
    t =  0.32;
    nl, nt, nw = 2, 3, 4;

    fens,fes  = L2block(L,nl)
    geom  =  NodalField(fens.xyz)
    bfes = FESetP1(reshape([nl+1], 1, 1))
axisymmetric = false
    femm  =  FEMMBase(IntegDomain(bfes, PointRule(), axisymmetric, L*W*t))
    S = integratefunction(femm, geom, (x) ->  1.0, 3)
    # # println("Length  of the boundary curve = $(S)")
    @test abs(S - L*W*t)/S < 1.0e-5
end
end
using .mmmeasurementm12
mmmeasurementm12.test()

module mmpartitioning1m
using FinEtools
using Test
import LinearAlgebra: norm
function test()
    a = 10. # radius of the hole
    nC = 20
    nR = 4
    fens,fes = Q4annulus(a, 1.5*a, nR, nC, 1.9*pi)

    npartitions = 8
    partitioning = nodepartitioning(fens, npartitions)
    partitionnumbers = unique(partitioning)
    @test norm(sort(partitionnumbers) - sort([1
    3
    2
    5
    6
    7
    8
    4])) < 1.e-5
end
end
using .mmpartitioning1m
mmpartitioning1m.test()


module mmpartitioning2m
using FinEtools
using Test
import LinearAlgebra: norm
function test()
    H = 100. # strip width
    a = 10. # radius of the hole
    L = 200. # length of the strip
    nL = 15
    nH = 10
    nR = 50
    fens,fes = Q4elliphole(a, a, L/2, H/2, nL, nH, nR)
@test count(fes) == 1250
    npartitions = 4
    partitioning = nodepartitioning(fens, npartitions)
    partitionnumbers = unique(partitioning)
    @test norm(sort(partitionnumbers) - sort([1
    3
    4
    2])) < 1.e-5
end
end
using .mmpartitioning2m
mmpartitioning2m.test()

module mmpartitioning3m
using FinEtools
using Test
import LinearAlgebra: norm
function test()
    H = 30. # strip width
    R = 10. # radius of the hole
    L = 20. # length of the strip
    nL = 15
    nH = 10
    nR = 5
    fens,fes = H8block(L, H, R, nL, nH, nR)

    npartitions = 16
    partitioning = nodepartitioning(fens, npartitions)
    partitionnumbers = unique(partitioning)
    @test norm(sort(partitionnumbers) - sort(1:npartitions)) < 1.e-5

    # for gp = partitionnumbers
    #   groupnodes = findall(k -> k == gp, partitioning)
    #   File =  "partition-nodes-$(gp).vtk"
    #   vtkexportmesh(File, fens, FESetP1(reshape(groupnodes, length(groupnodes), 1)))
    # end
    # File =  "partition-mesh.vtk"
    # vtkexportmesh(File, fens, fes)
    # @async run(`"paraview.exe" $File`)
end
end
using .mmpartitioning3m
mmpartitioning3m.test()

module mmboxm1
using FinEtools
using Test
import LinearAlgebra: norm
function test()
    a = [0.431345 0.611088 0.913161;
    0.913581 0.459229 0.82186;
    0.999429 0.965389 0.571139;
    0.390146 0.711732 0.302495;
    0.873037 0.248077 0.51149;
    0.999928 0.832524 0.93455]
    b1 = boundingbox(a)
    @test norm(b1 - [0.390146, 0.999928, 0.248077, 0.965389, 0.302495, 0.93455]) < 1.0e-4
    b2 = updatebox!(b1, a)
    @test norm(b1 - b2) < 1.0e-4
    b2 = updatebox!([], a)
    @test norm(b1 - b2) < 1.0e-4
    c = [-1.0, 3.0, -0.5]
    b3 = updatebox!(b1, c)
    # # println("$(b3)")
    @test norm(b3 - [-1.0, 0.999928, 0.248077, 3.0, -0.5, 0.93455]) < 1.0e-4
    x = [0.25 1.1 -0.3]
    @test inbox(b3, x)
    @test inbox(b3, c)
    @test inbox(b3, a[2, :])
    b4 = boundingbox(-a)
    # # println("$(b3)")
    # # println("$(b4)")
    # # println("$(boxesoverlap(b3, b4))")
    @test !boxesoverlap(b3, b4)
    b5 = updatebox!(b3, [0.0 -0.4 0.0])
    # # println("$(b5)")
    # # println("$(boxesoverlap(b5, b4))")
    @test boxesoverlap(b5, b4)
end
end
using .mmboxm1
mmboxm1.test()


module mmMeasurement_1
using FinEtools
using Test
function test()
    W = 1.1;
    L = 12.;
    t =  3.32;
    nl, nt, nw = 5, 3, 4;

    for or in [:a :b :ca :cb]
        fens,fes  = T4block(L,W,t, nl,nw,nt, or)
        geom  =  NodalField(fens.xyz)

        femm  =  FEMMBase(IntegDomain(fes, TetRule(5)))
        V = integratefunction(femm, geom, (x) ->  1.0)
        @test abs(V - W*L*t)/V < 1.0e-5
    end

end
end
using .mmMeasurement_1
mmMeasurement_1.test()

module mphunm2
using FinEtools
using Test
function test()
    for btu in [:SEC :MIN :HR :DY :YR :WK]
        # @show btu
        t = 0.333*phun(base_time_units = btu, "s")
        v = 2.0*phun(base_time_units = btu, "m/s")
        # display(v*t)
        @test abs(v*t - 0.666) < 1.0e-6
    end
end
end
using .mphunm2
mphunm2.test()

module mphunm3
using FinEtools
using Test
function test()
    for sou in [:SI :US :IMPERIAL :CGS :SIMM]
        t = 0.333*phun(system_of_units = sou, "s")
        v = 2.0*phun(system_of_units = sou, "m/s")
        # display(v*t/phun(system_of_units = sou, "m"))
        @test abs(v*t/phun(system_of_units = sou, "m") - 0.666) < 1.0e-6
    end
end
end
using .mphunm3
mphunm3.test()


module mxmeasurementm1
using FinEtools
using Test
function test()
    W = 1.1;
    L = 12.;
    t =  0.32;
    nl, nt, nw = 2, 3, 4;

    fens,fes  = H27block(L,W,t, nl,nw,nt)
    geom  =  NodalField(fens.xyz)

    femm  =  FEMMBase(IntegDomain(fes, GaussRule(3, 3)))
    V = integratefunction(femm, geom, (x) ->  1.0)
    @test abs(V - W*L*t)/V < 1.0e-5
end
end
using .mxmeasurementm1
mxmeasurementm1.test()

module mxmeasurementm2
using FinEtools
using Test
function test()
    W = 1.1;
    L = 12.;
    t =  0.32;
    nl, nt, nw = 2, 3, 4;

    fens,fes  = H27block(L,W,t, nl,nw,nt)
    geom  =  NodalField(fens.xyz)

bfes = meshboundary(fes)
    femm  =  FEMMBase(IntegDomain(bfes, GaussRule(2, 4)))
    V = integratefunction(femm, geom, (x) ->  1.0)
    @test abs(V - 2*(W*L+L*t+W*t))/V < 1.0e-5
end
end
using .mxmeasurementm2
mxmeasurementm2.test()

module mxmeasurementm3
using FinEtools
using FinEtools.MeshExportModule
using Test
function test()
    W = 4.1;
    L = 12.;
    t =  5.32;
    a = 0.3
    nl, nt, nw = 6, 8, 9;

    fens,fes  = H20block(L,W,t, nl,nw,nt)
    for ixxxx = 1:count(fens)
        x,y,z = fens.xyz[ixxxx, :]
        fens.xyz[ixxxx, :] = [x+a*sin(y) y+a*sin(z) z+a*sin(x)]
    end
    geom  =  NodalField(fens.xyz)
    # File = "mesh.vtk"
    # MeshExportModule.vtkexportmesh(File, fens, fes)

    bfes = meshboundary(fes)
    femm  =  FEMMBase(IntegDomain(bfes, GaussRule(2, 4)))
    S20 = integratefunction(femm, geom, (x) ->  1.0)

    fens,fes  = H27block(L,W,t, nl,nw,nt)
    for ixxxx = 1:count(fens)
        x,y,z = fens.xyz[ixxxx, :]
        fens.xyz[ixxxx, :] = [x+a*sin(y) y+a*sin(z) z+a*sin(x)]
    end
    geom  =  NodalField(fens.xyz)
    # File = "mesh.vtk"
    # MeshExportModule.vtkexportmesh(File, fens, fes)

    bfes = meshboundary(fes)
    femm  =  FEMMBase(IntegDomain(bfes, GaussRule(2, 4)))
    S27 = integratefunction(femm, geom, (x) ->  1.0)
    # # println("$((S20-S27)/S20)")

    @test abs((S20-S27)/S20) < 1.0e-5
end
end
using .mxmeasurementm3
mxmeasurementm3.test()

module mxmeasurementm4
using FinEtools
using FinEtools.MeshExportModule
using Test
function test()
    W = 4.1;
    L = 12.;
    t =  5.32;
    a = 0.3
    nl, nt, nw = 6, 8, 9;

    fens,fes  = H20block(L,W,t, nl,nw,nt)
    for ixxxx = 1:count(fens)
        x,y,z = fens.xyz[ixxxx, :]
        fens.xyz[ixxxx, :] = [x+a*sin(y) y+a*sin(z) z+a*sin(x)]
    end
    geom  =  NodalField(fens.xyz)
    # File = "mesh.vtk"
    # MeshExportModule.vtkexportmesh(File, fens, fes)

    femm  =  FEMMBase(IntegDomain(fes, GaussRule(3, 4)))
    V20 = integratefunction(femm, geom, (x) ->  1.0)

    fens,fes  = H27block(L,W,t, nl,nw,nt)
    for ixxxx = 1:count(fens)
        x,y,z = fens.xyz[ixxxx, :]
        fens.xyz[ixxxx, :] = [x+a*sin(y) y+a*sin(z) z+a*sin(x)]
    end
    geom  =  NodalField(fens.xyz)
    # File = "mesh.vtk"
    # MeshExportModule.vtkexportmesh(File, fens, fes)

    femm  =  FEMMBase(IntegDomain(fes, GaussRule(3, 4)))
    V27 = integratefunction(femm, geom, (x) ->  1.0)
    # # println("$((S20-S27)/S20)")

    @test abs((V20-V27)/V20) < 1.0e-5
end
end
using .mxmeasurementm4
mxmeasurementm4.test()

module mmmiscellaneous2
using FinEtools
using FinEtools.MeshExportModule
using Test
function test()
  rho=1.21*1e-9;# mass density
  c =345.0*1000;# millimeters per second
  bulk= c^2*rho;
  Lx=1900.0;# length of the box, millimeters
  Ly=800.0; # length of the box, millimeters

  fens,fes = Q4block(Lx,Ly,3,2); # Mesh
  X = xyz3(fens)
  fens.xyz = deepcopy(X)
  # show(fes.conn)
  # File = "mesh.vtk"
  # MeshExportModule.vtkexportmesh(File, fens, fes)

  bfes = meshboundary(fes)
  @test bfes.conn == Tuple{Int64,Int64}[(1, 2), (5, 1), (2, 3), (3, 4), (4, 8), (9, 5), (8, 12), (10, 9), (11, 10), (12, 11)]
end
end
using .mmmiscellaneous2
mmmiscellaneous2.test()

module mmmiscellaneous3
using FinEtools
using FinEtools.MeshExportModule
using Test
function test()
  rho=1.21*1e-9;# mass density
  c =345.0*1000;# millimeters per second
  bulk= c^2*rho;
  Lx=1900.0;# length of the box, millimeters
  Ly=800.0; # length of the box, millimeters

  fens,fes = L2block(Lx,4); # Mesh
  @test size(fens.xyz,2) == 1
  X = xyz3(fens)
  fens.xyz = deepcopy(X)
  @test size(fens.xyz,2) == 3
  X = xyz3(fens)
  fens.xyz = deepcopy(X)
  @test size(fens.xyz,2) == 3

  # show(fes.conn)
  # File = "mesh.vtk"
  # MeshExportModule.vtkexportmesh(File, fens, fes)

  bfes = meshboundary(fes)
  # show(fes.conn)
  @test bfes.conn == Tuple{Int64}[(1,), (5,)]
end
end
using .mmmiscellaneous3
mmmiscellaneous3.test()

module mmmfieldmm1
using FinEtools
using Test
import LinearAlgebra: norm
function test()
    W = 4.1;
    L = 12.;
    t =  5.32;
    a = 0.3
    nl, nt, nw = 6, 8, 9;

    fens,fes  = H8block(L,W,t, nl,nw,nt)
    geom  =  NodalField(fens.xyz)
    u = deepcopy(geom)
    setebc!(u)
    copyto!(u, geom)
    @test norm(u.values - geom.values) < 1.0e-5
    wipe!(u)
    numberdofs!(u)
    @test norm(u.dofnums) > 1.0e-3
    wipe!(u)
    @test norm(u.dofnums) == 0
    setebc!(u, [1, 3])
    numberdofs!(u)
    @test norm(u.dofnums[1,:]) == 0
    @test norm(u.dofnums[3,:]) == 0
    @test norm(u.dofnums[2,:]) > 0.0
end
end
using .mmmfieldmm1
mmmfieldmm1.test()

module mmcrossm
using FinEtools
using FinEtools.RotationUtilModule
using Test
import LinearAlgebra: norm, cross
function test()
    a = vec([0.1102, -0.369506, -0.0167305])
    b = vec([0.0824301, -0.137487, 0.351721])
    c = cross(a, b)
    # # println("$(c)")
    A = zeros(3, 3)
    skewmat!(A, a)
    d = A * b
    # # println("$(d)")
    @test norm(c - d) < 1.0e-6
    e = cross(a, b)
    @test norm(c - e) < 1.0e-6
    f = zeros(3)
    cross3!(f, a, b)
    @test norm(c - f) < 1.0e-6

    a = vec([0.1102, -0.135])
    b = vec([-0.137487, 0.351721])
    c = cross2(a, b)
    # # println("$(c)")
    @test norm(c - 0.0201989092) < 1.0e-6
end
end
using .mmcrossm
mmcrossm.test()

module mgen_iso_csmat1
using FinEtools
using FinEtools.FESetModule
using FinEtools.CSysModule
using Test
import LinearAlgebra: norm
import Statistics: mean
function test()
    L = 2.0
    nl = 1

    fens,fes  = L2block(L, nl)
    fens.xyz = xyz3(fens)
    fens.xyz[2, 3] += L/2
    csmatout = zeros(FFlt, 3, 1)
    gradNparams = FESetModule.bfundpar(fes, vec([0.0]));
    J = transpose(fens.xyz) * gradNparams
    CSysModule.gen_iso_csmat!(csmatout, mean(fens.xyz, dims = 1), J, 0)
    # # println("$(csmatout)")
    # # println("$(norm(vec(csmatout)))")
    @test norm(csmatout - [0.894427; 0.0; 0.447214]) < 1.0e-5
end
end
using .mgen_iso_csmat1
mgen_iso_csmat1.test()

module mgen_iso_csmat2
using FinEtools
using FinEtools.FESetModule
using FinEtools.CSysModule
using FinEtools.MeshExportModule
using Test
import LinearAlgebra: norm
import Statistics: mean
function test()
    L = 2.0
    nl = 1

    fens,fes  = Q4block(L, 2*L, nl, nl)
    fens.xyz = xyz3(fens)
    fens.xyz[2, 3] += L/2
    File = "mesh.vtk"
    MeshExportModule.VTK.vtkexportmesh(File, fens, fes)
    csmatout = zeros(FFlt, 3, 2)
    gradNparams = FESetModule.bfundpar(fes, vec([0.0 0.0]));
    J = zeros(3, 2)
    for i = 1:length(fes.conn[1])
        J += reshape(fens.xyz[fes.conn[1][i], :], 3, 1) * reshape(gradNparams[i, :], 1, 2)
    end
    # # println("J = $(J)")
    @test norm(J - [1.0 0.0; 0.0 2.0; 0.25 -0.25]) < 1.0e-5
    CSysModule.gen_iso_csmat!(csmatout, mean(fens.xyz, dims = 1), J, 0)
    # # println("csmatout = $(csmatout)")
    @test norm(csmatout - [0.970143 0.0291979; 0.0 0.992727; 0.242536 -0.116791]) < 1.0e-5
    try rm(File); catch end
end
end
using .mgen_iso_csmat2
mgen_iso_csmat2.test()

module mmfieldx1
using FinEtools
using Test
function test()
    f = NodalField(zeros(5, 1))
    setebc!(f, [3,4], true, 1, 7.0)
    # display(f)
    applyebc!(f)
    dest = zeros(2,1)
    # length([1,4])
    gatherfixedvalues_asmat!(f, dest, [1,4])
    # display(dest)
    @test (dest[1,1] == 0.0) && (dest[2,1] == 7.0)

    f = NodalField(zeros(5, 1))
    setebc!(f, [3,4], true, 1)
    # display(f)
    applyebc!(f)
    dest = zeros(2,1)
    # length([1,4])
    gatherfixedvalues_asmat!(f, dest, [1,4])
    # display(dest)
    @test (dest[1,1] == 0.0) && (dest[2,1] == 0.0)


    f = NodalField(zeros(5, 1))
    setebc!(f, [3,4], 1, 8.2)
    # display(f)
    applyebc!(f)
    dest = zeros(2,1)
    # length([1,4])
    gatherfixedvalues_asmat!(f, dest, [1,4])
    # display(dest)
    @test (dest[1,1] == 0.0) && (dest[2,1] == 8.2)
end
end
using .mmfieldx1
mmfieldx1.test()


module mmconnection1
using FinEtools
using Test
function test()
    h = 0.05*phun("M");
    l = 10*h;
    nh = 3; nl  = 4; nc = 5;

    fens,fes  = H8block(h,l,2.0*pi,nh,nl,nc)
    femm = FEMMBase(IntegDomain(fes, GaussRule(3, 2)))
    C = connectionmatrix(femm, count(fens))
    Degree = [length(findall(x->x!=0, C[j,:])) for j in 1:size(C, 1)]
    # # println("Maximum degree  = $(maximum(Degree))")
    # # println("Minimum degree  = $(minimum(Degree))")
    @test maximum(Degree) == 27
    @test minimum(Degree) == 8
end
end
using .mmconnection1
mmconnection1.test()

module mmquadrature3
using FinEtools
using Test
function test()
    gr = GaussRule(1,1)
    @test gr.param_coords[1] == 0.0
    @test_throws ErrorException gr = GaussRule(1,15)
    @test_throws AssertionError gr = GaussRule(4,1)

    @test_throws ErrorException tr = TriRule(8)

    @test_throws ErrorException tr = TetRule(3)
end
end
using .mmquadrature3
mmquadrature3.test()

module mmmatchingmm
using FinEtools.AlgoBaseModule: FDataDict, dcheck!
using Test
function test()
    d = FDataDict("postp"=>true, "preprocessing"=>nothing)
    recognized_keys = ["postprocessing", "something",  "else"]
    notmatched = dcheck!(d, recognized_keys)
    # display(notmatched)
    @test !isempty(notmatched)
end
end
using .mmmatchingmm
mmmatchingmm.test()


module mboundary13
using FinEtools
using FinEtools.MeshExportModule
using Test
function test()
    xs = collect(linearspace( 1.0, 3.0, 4))
    ys = collect(linearspace(-1.0, 5.0, 5))
    fens, fes = T3blockx(xs, ys, :a)
    @test count(fes) == 4*3*2

    boundaryfes  =   meshboundary(fes);
    ytl   = selectelem(fens, boundaryfes, facing = true,
        direction = [0.0 +1.0], dotmin = 0.999);
    @test length(ytl) == 3

    geom  =  NodalField(fens.xyz)

    femm  =  FEMMBase(IntegDomain(subset(boundaryfes, ytl), GaussRule(1, 2)))
    L = integratefunction(femm, geom, (x) ->  1.0)
    @test abs(L - 2.0)/L < 1.0e-5

    # File = "playground.vtk"
    # MeshExportModule.vtkexportmesh(File, fens, fes)
    # @async run(`"paraview.exe" $File`)
    # try rm(File) catch end

end
end
using .mboundary13
mboundary13.test()


module mboundary14
using FinEtools
using FinEtools.MeshExportModule
using Test
function test()
    xs = collect(linearspace( 1.0, 3.0, 4))
    ys = collect(linearspace(-1.0, 5.0, 5))
    fens, fes = Q4blockx(xs, ys)
    fens, fes = Q4toQ8(fens, fes)
    @test count(fes) == 4*3

    boundaryfes  =   meshboundary(fes);
    ytl   = selectelem(fens, boundaryfes, facing = true,
        direction = [0.0 +1.0], dotmin = 0.999);
    @test length(ytl) == 3

    geom  =  NodalField(fens.xyz)

    femm  =  FEMMBase(IntegDomain(subset(boundaryfes, ytl), GaussRule(1, 3)))
    L = integratefunction(femm, geom, (x) ->  1.0)
    @test abs(L - 2.0)/L < 1.0e-5

    # File = "playground.vtk"
    # MeshExportModule.vtkexportmesh(File, fens, fes)
    # @async run(`"paraview.exe" $File`)
    # try rm(File) catch end

end
end
using .mboundary14
mboundary14.test()


module mboundary15
using FinEtools
using FinEtools.MeshExportModule
using Test
function test()
    xs = collect(linearspace( 1.0, 3.0, 4))
    ys = collect(linearspace(-1.0, 5.0, 5))
    zs = collect(linearspace(+2.0, 5.0, 7))
    fens, fes = T4blockx(xs, ys, zs, :a)
    @test count(fes) == 432

    boundaryfes  =   meshboundary(fes);
    ytl   = selectelem(fens, boundaryfes, facing = true,
        direction = [0.0 +1.0 0.0], dotmin = 0.999);
    @test length(ytl) == 36

    geom  =  NodalField(fens.xyz)

    femm  =  FEMMBase(IntegDomain(subset(boundaryfes, ytl), SimplexRule(2, 1)))
    S = integratefunction(femm, geom, (x) ->  1.0)
    @test abs(S - 6.0)/S < 1.0e-5


    femm  =  FEMMBase(IntegDomain(fes, SimplexRule(3, 1)))
    V = integratefunction(femm, geom, (x) ->  1.0)
    @test abs(V - 36.0)/V < 1.0e-5

    # File = "playground.vtk"
    # MeshExportModule.vtkexportmesh(File, fens, fes)
    # @async run(`"paraview.exe" $File`)
    # try rm(File) catch end

end
end
using .mboundary15
mboundary15.test()

module mmmCSVm
using FinEtools
using FinEtools.MeshExportModule
using Test
function test()
    @test savecsv("a", a = rand(3), b = rand(3))
    rm( "a.csv")
end
end
using .mmmCSVm
mmmCSVm.test()


module mxmeasurementm3a1
using FinEtools
using FinEtools.MeshExportModule
using Test
function test()
    W = 4.1;
    L = 12.;
    t =  5.32;
    a = 0.4
    nl, nt, nw = 6, 8, 9;

    fens,fes  = H20block(L,W,t, nl,nw,nt)
    # # println("Mesh: $(count(fes))")
    for ixxxx = 1:count(fens)
        x,y,z = fens.xyz[ixxxx, :]
        fens.xyz[ixxxx, :] = [x+a*sin(y) y+x/10*a*sin(z) z+y*a*sin(x)]
    end
    geom  =  NodalField(fens.xyz)
    # File = "mesh.vtk"
    # MeshExportModule.vtkexportmesh(File, fens, fes)
    # @async run(`"paraview.exe" $File`)

    femm  =  FEMMBase(IntegDomain(fes, GaussRule(3, 4)))
    V20 = integratefunction(femm, geom, (x) ->  1.0)

    subregion1list = selectelem(fens, fes, box = [0.0 L/2 -Inf Inf -Inf Inf], inflate = t/1000)
    subregion2list = setdiff(1:count(fes), subregion1list)
    # # println("Sub mesh 1: $(length(subregion1list))")
    # # println("Sub mesh 2: $(length(subregion2list))")

    fes1 = subset(fes, subregion1list)
    connected1 = findunconnnodes(fens, fes1);
    fens1, new_numbering1 = compactnodes(fens, connected1);
    fes1 = renumberconn!(fes1, new_numbering1);
    present = findall(x -> x > 0, new_numbering1)
    geom1  =  NodalField(fens.xyz[present, :])

    fes2 = subset(fes, subregion2list)
    connected2 = findunconnnodes(fens, fes2);
    fens2, new_numbering2 = compactnodes(fens, connected2);
    fes2 = renumberconn!(fes2, new_numbering2);
    present = findall(x -> x > 0, new_numbering2)
    geom2  =  NodalField(fens.xyz[present, :])

    femm1  =  FEMMBase(IntegDomain(fes1, GaussRule(3, 4)))
    V20p = integratefunction(femm1, geom1, (x) ->  1.0)
    femm2  =  FEMMBase(IntegDomain(fes2, GaussRule(3, 4)))
    V20p += integratefunction(femm2, geom2, (x) ->  1.0)
    # # println("V20p = $(V20p)")
    # # println("V20 = $(V20)")

    @test abs(V20 - V20p)/V20 < 1.0e-6

end
end
using .mxmeasurementm3a1
mxmeasurementm3a1.test()

module mboxintersection_1
using FinEtools
using Test
import LinearAlgebra: norm
function test()
    a = [ 0.042525  0.455813  0.528458
    0.580612  0.933498  0.929843
    0.99648   0.800709  0.00175703
    0.433793  0.119944  0.966154
    0.793678  0.693062  0.919114]
    b1 = boundingbox(a)
    # println("b1 = $(b1)")
    a = [ 0.86714   0.311569   0.780585
    0.415177  0.60264    0.906292
    0.114056  0.0389293  0.733558
    0.657139  0.156761   0.83009
    0.890426  0.310158   0.516064]
    b2 = boundingbox(a)
    # println("b2 = $(b2)")
    b = intersectboxes(b1, b2)
    # println("b = $(b)")
    @test norm(b - [0.114056, 0.890426, 0.119944, 0.60264, 0.516064, 0.906292]) < 1.0e-4

    b1 = [0.042525, 0.49648, 0.119944, 0.933498, 0.00175703, 0.966154]
    b2 = [0.514056, 0.890426, 0.0389293, 0.60264, 0.516064, 0.906292]
    # println("b1 = $(b1)")
    # println("b2 = $(b2)")
    b = intersectboxes(b1, b2)
    # println("b = $(b)")
    @test length(b) == 0

    b1 = [0.042525, 0.69648, 0.119944, 0.933498, 0.00175703, 0.966154]
    b2 = [0.514056, 0.890426, 0.0389293, 0.060264, 0.516064, 0.906292]
    # println("b1 = $(b1)")
    # println("b2 = $(b2)")
    b = intersectboxes(b1, b2)
    # println("b = $(b)")
    @test length(b) == 0

    b1 = [0.042525, 0.69648, 0.119944, 0.933498]
    b2 = [0.514056, 0.890426, 0.0389293, 0.060264]
    # println("b1 = $(b1)")
    # println("b2 = $(b2)")
    b = intersectboxes(b1, b2)
    # println("b = $(b)")
    @test length(b) == 0

    b1 = [0.042525, 0.69648, 0.119944, 0.933498]
    b2 = [0.514056, 0.890426, 0.0389293, 0.160264]
    # println("b1 = $(b1)")
    # println("b2 = $(b2)")
    b = intersectboxes(b1, b2)
    # # println("b = $(b)")
    @test norm(b - [0.514056, 0.69648, 0.119944, 0.160264]) < 1.0e-4
end
end
using .mboxintersection_1
mboxintersection_1.test()

module mmmsearcconnectedelements
using FinEtools
using Test
function test()
    rin =  1.0;#internal radius
    rex =  2.0;#external radius
    nr = 20; nc = 180;
    Angle = 2*pi;
    tolerance = min(rin/nr,  rin/nc/2/pi)/10000;


    fens, fes  =  Q4annulus(rin, rex, nr, nc, Angle)
    fens, fes  =  mergenodes(fens,  fes,  tolerance);

    anode = [13, 3, 1961, 61, 43]
    global connectedcount = 0
    for i = 1:count(fes)
        for m = 1:length(anode)
            if in(anode[m], fes.conn[i])
                global connectedcount = connectedcount + 1
                break
            end
        end
    end
    # println("connectedcount = $(connectedcount)")

    connectedele = connectedelems(fes, anode, count(fens))
    # println("connectedele = $(connectedele)")

    @test connectedcount == length(connectedele)
end
end
using .mmmsearcconnectedelements
mmmsearcconnectedelements.test()

module mconjugategradient1
using FinEtools
using Test
using SparseArrays
import LinearAlgebra: norm, dot, lu, diff, cross
import FinEtools.AlgoBaseModule: conjugategradient

function test()
    # >> A = gallery('lehmer',6),
    A =[  1.0000    0.5000    0.3333    0.2500    0.2000    0.1667
        0.5000    1.0000    0.6667    0.5000    0.4000    0.3333
        0.3333    0.6667    1.0000    0.7500    0.6000    0.5000
        0.2500    0.5000    0.7500    1.0000    0.8000    0.6667
        0.2000    0.4000    0.6000    0.8000    1.0000    0.8333
        0.1667    0.3333    0.5000    0.6667    0.8333    1.0000]
    b =[0.8147
        0.9058
        0.1270
        0.9134
        0.6324
        0.0975]
    x0 = fill(zero(FFlt), 6)
    maxiter = 20
    x = conjugategradient(A, b, x0, maxiter)
end
end
using .mconjugategradient1
mconjugategradient1.test()

module mmbisecty
using FinEtools
using FinEtools.AlgoBaseModule: bisect
using Test
function test()
    xs = [231.7, 239.1, 244.8]
    hs = [(4.0/5.)^i for i in range(0, length = 3)]
    approxerrors = diff(xs)
    fun = y ->  approxerrors[1] / (hs[1]^y - hs[2]^y) - approxerrors[2] / (hs[2]^y - hs[3]^y)

    beta = bisect(fun, 0.001, 2.0, 0.00001, 0.00001)
    beta = (beta[1] + beta[2]) / 2.0
    @assert abs(beta - 1.1697) < 1.0e-4
    gamm1 = hs[2]^beta / hs[1]^beta
    gamm2 = hs[3]^beta / hs[2]^beta
    estimtrueerror1 = approxerrors[1] ./ (1.0 - gamm1)
    estimtrueerror2 = approxerrors[2] ./ (1.0 - gamm2)
    # println("xs[1] + estimtrueerror1 = $(xs[1] + estimtrueerror1)")
    # println("xs[2] + estimtrueerror2 = $(xs[2] + estimtrueerror2)")
    true
end
end
using .mmbisecty
mmbisecty.test()

module mmh2libexporttri
using FinEtools
import FinEtools.MeshExportModule.H2Lib: h2libexporttri
using Test
function test()
triangles = [4 3 5
3 1 2
5 2 6
6 2 7
7 8 11
8 10 11
9 10 8
2 9 8
1 9 2
3 2 5]
vertices = [0.943089 0.528795; 0.704865 0.101195; 0.980527 0.225325; 0.126678 0.323197; 0.969361 0.376956; 0.339188 0.306131; 0.322906 0.244326; 0.327952 0.769903; 0.820541 0.766872; 0.919947 0.390552; 0.958515 0.533267]
h2libexporttri("sample.tri", triangles, vertices)
try rm("sample.tri"); catch end
end
end
using .mmh2libexporttri
mmh2libexporttri.test()


module mmqtrapm1
using FinEtools
using Test
using FinEtools.AlgoBaseModule: qtrap, qcovariance, qvariance
function test()
    ps = range(0, stop=1.0, length = 1000)
    xs = sin.(40.0 * ps)
    ys = sin.(39.0 * ps)
    # println("qcovariance(ps, xs, ys) = $(qcovariance(ps, xs, ys))")
    @test abs(qcovariance(ps, xs, ys) - 0.422761407481519) < 1.0e-3
    true
end
end
using .mmqtrapm1
mmqtrapm1.test()

module mmqtrapm2
using FinEtools
using Test
using FinEtools.AlgoBaseModule: qtrap, qcovariance, qvariance
function test()
    ps = range(0, stop=1.0, length = 1000)
    xs = sin.(40.0 * ps)
    # @show qvariance(ps, xs)
    # @show cov(xs)
    # println("qvariance(ps, xs) = $(qvariance(ps, xs))")
    @test abs(qvariance(ps, xs) - 0.504472271593515) < 1.0e-3
    true
end
end
using .mmqtrapm2
mmqtrapm2.test()

module mmqtrapm3
using FinEtools
using Test
using FinEtools.AlgoBaseModule: qtrap, qcovariance, qvariance
using StatsBase
import Statistics: mean, cov
function test()
    ps = vec([0.0201161  0.0567767  0.123692  0.141182  0.196189  0.311646  0.463708   0.798094  0.832338  0.875213  0.880719  0.947033  0.993938  0.99884]) #sort(rand(20))
    xs = sin.(40.0 * ps)
    ys = sin.(39.0 * ps)
    # @show qcovariance(ps, xs, ys)
    # println("qcovariance(ps, xs, ys) = $(qcovariance(ps, xs, ys))")
    @test abs(qcovariance(ps, xs, ys) - 0.273726042904099) < 1.0e-6
    # @show cov(xs, ys)
    # println("cov(xs, ys) = $(cov(xs, ys))")
    @test abs(cov(xs, ys) - 0.3917048440575396) < 1.0e-6
    true
end
end
using .mmqtrapm3
mmqtrapm3.test()


module mfixupdecimal1
using FinEtools
using Test
using FinEtools.MeshImportModule: _fixupdecimal
function test()
    s = _fixupdecimal("-10.0-1")
    @test parse(Float64, s) == -1.0
    s = _fixupdecimal("-10.033333-1")
    @test parse(Float64, s) == -1.0033333
    s = _fixupdecimal("-100.33333-002")
    @test parse(Float64, s) == -1.0033333
    s = _fixupdecimal("-1.0033333+0")
    @test parse(Float64, s) == -1.0033333
    s = _fixupdecimal("-1.0033333+000")
    @test parse(Float64, s) == -1.0033333
    s = _fixupdecimal(" +1.0033333+000 ")
    @test parse(Float64, s) == +1.0033333
    true
end
end
using .mfixupdecimal1
mfixupdecimal1.test()

module mPartitioning6mmmmmm
using FinEtools
using LinearAlgebra: norm
using Test
function test()
    H = 30. # strip width
    L = 20. # length of the strip
    nL = 15
    nH = 10
    fens,fes = Q4block(L, H, nL, nH)
    reg1 = selectelem(fens, fes; box = [0.0 L/3 0.0 H], inflate = L/nL/1000)
    reg2 = selectelem(fens, fes; box = [L/3 L 0.0 H/2], inflate = L/nL/1000)
    reg3 = selectelem(fens, fes; box = [L/3 L H/2 H], inflate = L/nL/1000)
    fesarr = vec([subset(fes, reg1) subset(fes, reg2) subset(fes, reg3)])

    # Partitioning of all the nodes
    partitioning = nodepartitioning(fens, fesarr, vec([2 4 2]))
   @test norm(partitioning .- [2, 2, 2, 2, 2, 2, 6, 6, 6, 6, 6, 5, 5, 5, 5, 5, 2, 2, 2, 2, 2, 2, 6, 6, 6, 6, 6, 5, 5, 5, 5, 5, 2, 2, 2, 2, 2, 2, 6, 6,
   6, 6, 6, 5, 5, 5, 5, 5, 2, 2, 2, 2, 2, 2, 4, 4, 4, 4, 4, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 4, 4, 4, 4, 4, 3, 3, 3, 3, 3, 1, 1, 1, 2, 2, 2, 4, 4, 4, 4, 4, 3, 3, 3, 3, 3, 1, 1, 1, 1, 1, 1, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 1, 1, 1, 1, 1, 1, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 1, 1, 1, 1, 1, 1, 7, 7, 7, 7, 7, 8, 8, 8, 8, 8, 1, 1, 1, 1, 1, 1, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 1, 1, 1, 1, 1, 1, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7]) == 0
    partitionnumbers = unique(partitioning)

    # Visualize partitioning
    # for gp = partitionnumbers
    #   groupnodes = findall(k -> k == gp, partitioning)
    #   @show groupnodes
    # #   File =  "partition-nodes-$(gp).vtk"
    # #   vtkexportmesh(File, fens, FESetP1(reshape(groupnodes, length(groupnodes), 1)))
    # end
    # for i = 1:length(fesarr)
    #     File =  "mesh-$(i).vtk"
    #     vtkexportmesh(File, fens, fesarr[i])
    # end
    # File =  "mesh-1.vtk"
    # @async run(`"paraview.exe" $File`)
end
end
using .mPartitioning6mmmmmm
mPartitioning6mmmmmm.test()

module mmsparsem1
using FinEtools
using SparseArrays
using Test
function test()
    N = 10
    a = vec(rand(N))
    i = convert(Vector{Int64}, round.(rand(N) * N) .+ 1)
    j = convert(Vector{Int64}, round.(rand(N) * N) .+ 1)
    # println("i = $(i)")
    # println("j = $(j)")
    A = sparse(i, j, a, N+1, N+1)
    @test typeof(A * A) <: AbstractSparseMatrix
    true
end
end
using .mmsparsem1
mmsparsem1.test()

module mmMeasurement_1a
using FinEtools
using Test
function test()
    W = 1.1;
    L = 12.;
    t =  3.32;
    nl, nt, nw = 5, 3, 4;

    # println("New segmentation fault?")
    for orientation in [:a :b :ca :cb]
        fens,fes  = T4block(L,W,t, nl,nw,nt, orientation)
        geom  =  NodalField(fens.xyz)

        femm  =  FEMMBase(IntegDomain(fes, NodalSimplexRule(3)))
        V = integratefunction(femm, geom, (x) ->  1.0)
        @test abs(V - W*L*t)/V < 1.0e-5
    end

end
end
using .mmMeasurement_1a
mmMeasurement_1a.test()

module mmMeasurement_2a
using FinEtools
using Test
function test()
    W = 1.1;
    L = 12.;
    t =  3.32;
    nl, nt, nw = 5, 3, 4;

    # println("New segmentation fault?")
    for orientation in [:a :b :ca :cb]
        fens,fes  = T4block(L,W,t, nl,nw,nt, orientation)
        geom  =  NodalField(fens.xyz)

        femm  =  FEMMBase(IntegDomain(fes, NodalSimplexRule(3)))
        V = integratefunction(femm, geom, (x) ->  1.0)
        @test abs(V - W*L*t)/V < 1.0e-5
    end

end
end
using .mmMeasurement_2a
mmMeasurement_2a.test()

module mmAdjugate
using FinEtools
using Test
using LinearAlgebra
import FinEtools.MatrixUtilityModule: adjugate3!
function test()
    A = rand(3, 3)
    B = rand(3, 3)
    adjugate3!(B, A)
    @test norm(det(A)*inv(A)-B) < 1.0e-6 * norm(A)
end
end
using .mmAdjugate
mmAdjugate.test()

module mxmeasurementm3b1
using FinEtools
using FinEtools.MeshExportModule
using Test
function test()
    W = 4.1;
    L = 12.;
    t =  5.32;
    a = 0.4
    nl, nt, nw = 6, 8, 9;

    fens,fes  = H8block(L,W,t, nl,nw,nt)
    # # println("Mesh: $(count(fes))")
    for ixxxx = 1:count(fens)
        x,y,z = fens.xyz[ixxxx, :]
        fens.xyz[ixxxx, :] = [x+a*sin(y) y+x/10*a*sin(z) z+y*a*sin(x)]
    end
    geom  =  NodalField(fens.xyz)
    # File = "mesh.vtk"
    # MeshExportModule.vtkexportmesh(File, fens, fes)
    # @async run(`"paraview.exe" $File`)

    femm  =  FEMMBase(IntegDomain(fes, TrapezoidalRule(3)))
    V8 = integratefunction(femm, geom, (x) ->  1.0)

    subregion1list = selectelem(fens, fes, box = [0.0 L/2 -Inf Inf -Inf Inf], inflate = t/1000)
    subregion2list = setdiff(1:count(fes), subregion1list)
    # # println("Sub mesh 1: $(length(subregion1list))")
    # # println("Sub mesh 2: $(length(subregion2list))")

    fes1 = subset(fes, subregion1list)
    connected1 = findunconnnodes(fens, fes1);
    fens1, new_numbering1 = compactnodes(fens, connected1);
    fes1 = renumberconn!(fes1, new_numbering1);
    present = findall(x -> x > 0, new_numbering1)
    geom1  =  NodalField(fens.xyz[present, :])

    fes2 = subset(fes, subregion2list)
    connected2 = findunconnnodes(fens, fes2);
    fens2, new_numbering2 = compactnodes(fens, connected2);
    fes2 = renumberconn!(fes2, new_numbering2);
    present = findall(x -> x > 0, new_numbering2)
    geom2  =  NodalField(fens.xyz[present, :])

    femm1  =  FEMMBase(IntegDomain(fes1, GaussRule(3, 4)))
    V8p = integratefunction(femm1, geom1, (x) ->  1.0)
    femm2  =  FEMMBase(IntegDomain(fes2, GaussRule(3, 4)))
    V8p += integratefunction(femm2, geom2, (x) ->  1.0)
    # # println("V20p = $(V20p)")
    # # println("V20 = $(V20)")

    @test abs(V8 - V8p)/V8 < 1.0e-4

end
end
using .mxmeasurementm3b1
mxmeasurementm3b1.test()

module mxmeasurementm3c1
using FinEtools
using FinEtools.MeshExportModule
using Test
function test()
    L = 12.1331;
    nl = 9;

    fens,fes  = L2block(L, nl)
    geom  =  NodalField(fens.xyz)
    # File = "mesh.vtk"
    # MeshExportModule.vtkexportmesh(File, fens, fes)
    # @async run(`"paraview.exe" $File`)

    femm  =  FEMMBase(IntegDomain(fes, TrapezoidalRule(1)))
    La = integratefunction(femm, geom, (x) ->  1.0)

    @test abs(L - La)/L < 1.0e-9

end
end
using .mxmeasurementm3c1
mxmeasurementm3c1.test()

module mmmBisectm1
using FinEtools
using FinEtools.AlgoBaseModule: bisect
using Test
function test()
    fun = x -> 1.0 - x^2
    xl = 0.0
    xu = 1.0
    tolx = eps(1.0)
    tolf = eps(1.0)
    @test_throws AssertionError bracket = bisect(fun, xl, xu, tolx, tolf)

    fun = x -> (x - 0.5)^3
    xl = -0.5
    xu = 1.5
    tolx = eps(1.0)
    tolf = eps(1.0)
    bracket = bisect(fun, xl, xu, tolx, tolf)
    @test bracket == (0.5, 0.5)
end
end
using .mmmBisectm1
mmmBisectm1.test()

module mmmmtest_trirules
using FinEtools
using Test
import LinearAlgebra: cholesky
function test()
    for NPTS = [1, 3, 4, 6, 7, 9, 12, 13]
        r = TriRule(NPTS)
        # @show NPTS
        @test abs(sum(r.weights) - 0.5) <1.e-5
    end
    true
end
end
using .mmmmtest_trirules
mmmmtest_trirules.test()

module mtestRichardson1
using FinEtools
using Test
import LinearAlgebra: norm
using FinEtools.AlgoBaseModule: richextrapol, richextrapoluniform
function test()
    # Tests of Richardson extrapolation for sequence with uniform refinement factor
    testsets = [ (alpha = 0.5, beta = 2.0, c = 0.3, truesol = 1.0),
           (alpha = 0.85, beta = 1.2, c = -0.3, truesol = 10.0),
           (alpha = 0.85, beta = 1.2, c = -0.3, truesol = -10.0),
           (alpha = 0.25, beta = 1.2, c = -1.3, truesol = -10.0),
           (alpha = 0.66, beta = 1.02, c = -1.3, truesol = -10.0e4),
           (alpha = 0.96, beta = 0.3, c = -1.3, truesol = -1.0e4),
           ]
    for t in testsets
        parameters = [1.0, t.alpha, t.alpha^2]
        model = p -> (t.truesol) - (t.c)*p^(t.beta)
        solns = [model(p) for p in parameters]
        solnestim, betaestim, cestim, residual = richextrapol(solns, parameters)
        @test abs(t.truesol - solnestim) < 1.0e-6 * abs(t.truesol)
        @test abs(t.beta - betaestim) < 1.0e-5 * t.beta
        @test abs(t.c - cestim) < 1.0e-5 * abs(t.c)
    end
    for t in testsets
        parameters = [1.0, t.alpha, t.alpha^2]
        model = p -> (t.truesol) - (t.c)*p^(t.beta)
        solns = [model(p) for p in parameters]
        solnestim, betaestim, cestim, residual = richextrapoluniform(solns, parameters)
        @test abs(t.truesol - solnestim) < 1.0e-6 * abs(t.truesol)
        @test abs(t.beta - betaestim) < 1.0e-4 * t.beta
        @test abs(t.c - cestim) < 1.0e-4 * abs(t.c)
    end
    true
end
end
using .mtestRichardson1
mtestRichardson1.test()

module mtestRichardson2
using FinEtools
using Test
import LinearAlgebra: norm
using FinEtools.AlgoBaseModule: richextrapol
function test()
    # Tests of Richardson extrapolation for sequence with Non-uniform refinement
    testsets = [ (parameters = [0.5, 0.3, 0.1], beta = 2.0, c = 0.3, truesol = 1.0),
           (parameters = [0.5, 0.4, 0.1], beta = 1.2, c = -0.3, truesol = 10.0),
           (parameters = [0.54, 0.14, 0.1], beta = 1.2, c = -0.3, truesol = 10.0),
           (parameters = [0.95, 0.41, 0.1], beta = 1.2, c = -0.3, truesol = 10.0),
           (parameters = [0.5, 0.43, 0.1], beta = 1.2, c = -0.3, truesol = 10.0),
           (parameters = [0.5, 0.3, 0.2], beta = 1.2, c = -0.3, truesol = -10.0),
           (parameters = [0.5, 0.3, 0.01], beta = 1.2, c = -1.3, truesol = -10.0),
           (parameters = [0.5, 0.49, 0.01], beta = 1.2, c = -1.3, truesol = -10.0),
           (parameters = [0.5, 0.03, 0.01], beta = 1.2, c = -1.3, truesol = -10.0),
           (parameters = [0.5, 0.13, 0.1], beta = 1.02, c = -1.3, truesol = -10.0e4),
           (parameters = [0.35, 0.3, 0.1], beta = 0.3, c = -1.3, truesol = -1.0e4),
           ]
    for t in testsets
      parameters = t.parameters
      model = p -> (t.truesol) - (t.c)*p^(t.beta)
      solns = [model(p) for p in parameters]
      solnestim, betaestim, cestim, residual = richextrapol(solns, parameters)
      @test abs(t.truesol - solnestim) < 1.0e-6 * abs(t.truesol)
      @test abs(t.beta - betaestim) < 1.0e-5 * t.beta
      @test abs(t.c - cestim) < 1.0e-5 * abs(t.c)
    end
    true
end
end
using .mtestRichardson2
mtestRichardson2.test()

module mmassembly
using FinEtools
using Test
import LinearAlgebra: norm, cholesky
function test()
	a = SysmatAssemblerSparse(0.0)
	startassembly!(a, 5, 5, 3, 7, 7)
	m = [0.24406   0.599773    0.833404  0.0420141
	0.786024  0.00206713  0.995379  0.780298
	0.845816  0.198459    0.355149  0.224996]
	assemble!(a, m, [1 7 5], [5 2 1 4])
	m = [0.146618  0.53471   0.614342    0.737833
 0.479719  0.41354   0.00760941  0.836455
 0.254868  0.476189  0.460794    0.00919633
 0.159064  0.261821  0.317078    0.77646
 0.643538  0.429817  0.59788     0.958909]
	assemble!(a, m, [2 3 1 7 5], [6 7 3 4])
	A = makematrix!(a)
	@test abs.(maximum([0.833404 0.599773 0.460794 0.0512104 0.24406 0.254868 0.476189; 0.0 0.0 0.614342 0.737833 0.0 0.146618 0.53471; 0.0 0.0 0.00760941 0.836455 0.0 0.479719 0.41354; 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.355149 0.198459 0.59788 1.1839 0.845816 0.643538 0.429817; 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.995379 0.00206713 0.317078 1.55676 0.786024 0.159064 0.261821]  - A)) < 1.0e-5
	# @test abs(maximum(T_i)-1380.5883006341187) < 1.0e-3
end
end
using .mmassembly
mmassembly.test()

module mmassembly2
using FinEtools
using Test
import LinearAlgebra: norm, cholesky
function test()
	a = SysmatAssemblerSparseSymm(0.0)
		startassembly!(a, 5, 5, 3, 7, 7)
		m = [0.24406   0.599773    0.833404  0.0420141
			0.786024  0.00206713  0.995379  0.780298
			0.845816  0.198459    0.355149  0.224996]
		assemble!(a, m'*m, [5 2 1 4], [5 2 1 4])
		m = [0.146618  0.53471   0.614342    0.737833
			 0.479719  0.41354   0.00760941  0.836455
			 0.254868  0.476189  0.460794    0.00919633
			 0.159064  0.261821  0.317078    0.77646
			 0.643538  0.429817  0.59788     0.958909]
		assemble!(a, m'*m, [2 3 1 5], [2 3 1 5])
		A = makematrix!(a)
	@test abs.(maximum([ 2.85928   1.21875    0.891063  0.891614   2.56958   0.0  0.0
 1.21875   1.15515    0.716396  0.0714644  1.56825   0.0  0.0
 0.891063  0.716396   0.936979  0.0        1.36026   0.0  0.0
 0.891614  0.0714644  0.0       0.661253   0.813892  0.0  0.0
 2.56958   1.56825    1.36026   0.813892   4.15934   0.0  0.0
 0.0       0.0        0.0       0.0        0.0       0.0  0.0
 0.0       0.0        0.0       0.0        0.0       0.0  0.0 ]  - A)) < 1.0e-5
	@test abs.(maximum(A - transpose(A))) < 1.0e-5
end
end
using .mmassembly2
mmassembly2.test()

module vectorcachetest1
using FinEtools
using Test
function test()
   XYZ = reshape([0.0, 0.0], 2, 1)
   tangents = reshape([0.0, 1.0], 2, 1)
   fe_label = 0
   c = VectorCache([10.0])
   updateretrieve!(c, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)
   @test c.cache == [10.0]
end
end
using .vectorcachetest1
vectorcachetest1.test()

module vectorcachetest2
using FinEtools
using Test
function test()
   XYZ = reshape([0.0, 0.0], 2, 1)
   tangents = reshape([0.0, 1.0], 2, 1)
   fe_label = 0
   setvector!(v, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt; time::FFlt = 0.0) = begin
      v .= [10.0]
   end
   c = VectorCache(FFlt, 1, setvector!)
   updateretrieve!(c, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)
   @test c.cache == [10.0]
end
end
using .vectorcachetest2
vectorcachetest2.test()

module vectorcachetest3
using FinEtools
using Test
function test()
   XYZ = reshape([0.0, 0.0], 2, 1)
   tangents = reshape([0.0, 1.0], 2, 1)
   fe_label = 0
   setvector!(v, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt; time::FFlt = 0.0) = begin
      if time < 5.0
         v .= [10.0]
      else
         v .= [0.0]
      end
      return v
   end
   c = VectorCache(FFlt, 1, setvector!, 0.0)
   updateretrieve!(c, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)
   @test c.cache == [10.0]
   settime!(c, 6.0)
   updateretrieve!(c, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)
   @test c.cache == [0.0]
end
end
using .vectorcachetest3
vectorcachetest3.test()

module vectorcachetest4
using FinEtools
using Test
function test()
   XYZ = reshape([0.0, 0.0], 2, 1)
   tangents = reshape([0.0, 1.0], 2, 1)
   fe_label = 0
   setvector!(v, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt) = begin
      v .= [10.0]
      return v
   end
   c = VectorCache(FFlt, 1, setvector!)
   updateretrieve!(c, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)
   @test c.cache == [10.0]
end
end
using .vectorcachetest4
vectorcachetest4.test()

module forceintensitytest1
using FinEtools
using Test
function test()
   XYZ = reshape([0.0, 0.0], 2, 1)
   tangents = reshape([0.0, 1.0], 2, 1)
   fe_label = 0
   vector = [10.0]
   fi = ForceIntensity(vector)
   v = updateforce!(fi, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)
   @test v == [10.0]
end
end
using .forceintensitytest1
forceintensitytest1.test()

module forceintensitytest2
using FinEtools
using Test
function test()
   XYZ = reshape([0.0, 0.0], 2, 1)
   tangents = reshape([0.0, 1.0], 2, 1)
   fe_label = 0
   setvector!(v, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt; time::FFlt = 0.0) = begin
      v .= [10.0]
   end
   vector = [10.0]
   fi = ForceIntensity(FFlt, length(vector), setvector!)
   v = updateforce!(fi, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)
   @test v == [10.0]
end
end
using .forceintensitytest2
forceintensitytest2.test()

module forceintensitytest3
using FinEtools
using Test
function test()
   XYZ = reshape([0.0, 0.0], 2, 1)
   tangents = reshape([0.0, 1.0], 2, 1)
   fe_label = 0
   setvector!(v, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt; time::FFlt = 0.0) = begin
      if time < 5.0
         v .= [10.0]
      else
         v .= [0.0]
      end
      return v
   end
   vector = [10.0]
   fi = ForceIntensity(FFlt, length(vector), setvector!, 0.0)
   v = updateforce!(fi, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)
   @test v == [10.0]
   settime!(fi, 6.0)
   v = updateforce!(fi, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)
   @test v == [0.0]
end
end
using .forceintensitytest3
forceintensitytest3.test()

module forceintensitytest4
using FinEtools
using Test
function test()
   XYZ = reshape([0.0, 0.0], 2, 1)
   tangents = reshape([0.0, 1.0], 2, 1)
   fe_label = 0
   setvector!(v, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt) = begin
      v .= [10.0]
      return v
   end
   vector = [10.0]
   fi = ForceIntensity(FFlt, length(vector), setvector!)
   v = updateforce!(fi, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)
   @test v == [10.0]
   settime!(fi, 6.0)
   v = updateforce!(fi, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)
   @test v == [10.0]
end
end
using .forceintensitytest4
forceintensitytest4.test()

module surfacenormaltest1
using FinEtools
using Test
function test()
   XYZ = reshape([0.0, 0.0], 2, 1)
   tangents = reshape([0.0, 1.0], 2, 1)
   fe_label = 0
   vector = [10.0, -3.0]
   fi = SurfaceNormal(vector)
   v = updatenormal!(fi, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)
   @test v == [10.0, -3.0]
end
end
using .surfacenormaltest1
surfacenormaltest1.test()

module surfacenormaltest2
using FinEtools
using Test
function test()
   XYZ = reshape([0.0, 0.0], 2, 1)
   tangents = reshape([0.0, 1.0], 2, 1)
   fe_label = 0
   setvector!(v, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt; time::FFlt = 0.0) = begin
      v .= [10.0, -3.13]
   end
   fi = SurfaceNormal(2, setvector!)
   v = updatenormal!(fi, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)
   @test v == [10.0, -3.13]
end
end
using .surfacenormaltest2
surfacenormaltest2.test()

module surfacenormaltest3
using FinEtools
using Test
function test()
   XYZ = reshape([0.0, 0.0], 2, 1)
   tangents = reshape([-1.0, 1.0], 2, 1)
   fe_label = 0
   fi = SurfaceNormal(2)
   v = updatenormal!(fi, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)
   @test v == [0.7071067811865475, 0.7071067811865475]
end
end
using .surfacenormaltest3
surfacenormaltest3.test()

module mxmeasure1
using FinEtools
using Test
function test()
    W = 4.1;
    L = 12.;
    t =  3.32;
    nl, nt, nw = 2, 3, 4;

    fens,fes  = H27block(L,W,t, nl,nw,nt)
    geom  =  NodalField(fens.xyz)

    femm  =  FEMMBase(IntegDomain(fes, GaussRule(3, 3)))

    # Test the calculation of the volume
    V = integratefunction(femm, geom, (x) ->  1.0)
    @test abs(V - W*L*t)/V < 1.0e-5

    # Test the calculation of the center of gravity
    Sx = integratefunction(femm, geom, (x) ->  x[1])
    @test Sx / V ≈ L / 2

    # Test the calculation of the moments of inertia
    # The block is translated to be centered at the origin
    fens.xyz[:, 1] .-= L / 2
    fens.xyz[:, 2] .-= W / 2
    fens.xyz[:, 3] .-= t / 2
    geom  =  NodalField(fens.xyz)
    Sx = integratefunction(femm, geom, (x) ->  x[1])
    @test abs(Sx / V) / L <= eps(1.0)
    Sy = integratefunction(femm, geom, (x) ->  x[2])
    @test abs(Sy / V) / W <= eps(1.0)
    Sz = integratefunction(femm, geom, (x) ->  x[3])
    @test abs(Sz / V) / t <= eps(1.0)
    Ixx = integratefunction(femm, geom, (x) ->  x[2]^2 + x[3]^2)
    @test abs(Ixx - V * (W^2 + t^2) / 12) / V <= eps(V)
    Iyy = integratefunction(femm, geom, (x) ->  x[1]^2 + x[3]^2)
    @test abs(Iyy - V * (L^2 + t^2) / 12) / V <= eps(V)
    Ixy = integratefunction(femm, geom, (x) ->  (-x[1] * x[2]))
    @test abs(Ixy) / V <= eps(V)
    true
end
end
using .mxmeasure1
mxmeasure1.test()

module maddbts1
using FinEtools
using FinEtools.MatrixUtilityModule: add_btsigma!
using LinearAlgebra: norm
using Test
function test()
	for i in 1:10
		B = rand(6, 30)
		s = rand(6)
		c = rand()
		F1 = B' * (c * s)
		F2 = fill(0.0, size(B, 2))
		add_btsigma!(F2, B, c, s)
		@test norm(F1 - F2) / norm(F1) < 1.0e-6
	end
	true
end
end
using .maddbts1
maddbts1.test()

module mxmatmul3a1
using FinEtools
using FinEtools.MatrixUtilityModule: mulCAB!, mulCAtB!, mulCABt!
using LinearAlgebra: norm
using Test

function test()
	for i in 1:10
		A = rand(3, 3)
		B = rand(3, 3)
		C = rand(3, 3)
		@test norm(mulCAB!(Val(3), C, A, B) .- A*B) <= 1.0e-6 * norm(C)
	end
    for i in 1:10
    	A = rand(3, 3)
    	B = rand(3, 3)
    	C = rand(3, 3)
    	@test norm(mulCABt!(Val(3), C, A, B) .- A*B') <= 1.0e-6 * norm(C)
    end
    for i in 1:10
    	A = rand(3, 3)
    	B = rand(3, 3)
    	C = rand(3, 3)
    	@test norm(mulCAtB!(Val(3), C, A, B) .- A'*B) <= 1.0e-6 * norm(C)
    end
    true
end
end
using .mxmatmul3a1
mxmatmul3a1.test()


module mxmatmul3a2
using FinEtools
using FinEtools.MatrixUtilityModule: mulCAB!, mulCAtB!, mulCABt!
using LinearAlgebra: norm
using Test

function test()
	for i in 1:10
		A = rand(3, 5)
		B = rand(5, 3)
		C = A * B
		@test norm(mulCAB!(C, A, B) .- A*B) <= 1.0e-6 * norm(C)
	end
    for i in 1:10
    	A = rand(8, 3)
    	B = rand(8, 3)
    	C = A' * B
    	@test norm(mulCAtB!(C, A, B) .- A'*B) <= 1.0e-6 * norm(C)
    end
    for i in 1:10
    	A = rand(8, 3)
    	B = rand(8, 3)
    	C = A * B'
    	@test norm(mulCABt!(C, A, B) .- A*B') <= 1.0e-6 * norm(C)
    end
    true
end
end
using .mxmatmul3a2
mxmatmul3a2.test()

module mxmatmul3a3
using FinEtools
using FinEtools.MatrixUtilityModule: detC
using LinearAlgebra: norm, det
using Test

function test()
	for i in 1:10
		C = rand(3, 3)
		@test abs(detC(Val(3), C) - det(C)) <= 1.0e-6 * norm(C)
	end
    true
end
end
using .mxmatmul3a3
mxmatmul3a3.test()

module mrotationmatrixs1
using FinEtools
using LinearAlgebra: norm, I
# using BenchmarkTools
using Test
function test()
    for i in 1:10
        Rmout = rand(3, 3)
        a = rand(3)
        # Rmout = [0.056575544213388396 0.8235145288079604 0.4566335160952417; 0.6335215149040396 0.8852762577618685 0.803417762628938; 0.5940995203626245 0.35601280936101065 0.6825575362008556]                                                                                                 
        # a = [0.11011521831193871, 0.5097478695998647, 0.8139760429477749]  
        rotmat3!(Rmout, a)
        # Rmout = [0.5736166081143759 -0.6670428981551146 0.47541325067375245; 0.7189364978350996 0.6881251218379588 0.0980516638109532; -0.3925484670406451 
        # 0.2855478746485778 0.8742814834523946]      
        @test norm(Rmout'*Rmout - 1.0*I) <= 1.0e-6
        # @btime rotmat3!($Rmout, $a)
    end
    true
end
end
using .mrotationmatrixs1
mrotationmatrixs1.test()

module mrotationmatrixs2
using FinEtools
using LinearAlgebra: norm, I
# using BenchmarkTools
using Test
function rotmat3original!(Rmout, a) 
    na = norm(a);
    thetatilde = zeros(3,3);
    skewmat!(thetatilde, a);
    a = a/na;
    ca = cos(na);
    sa = sin(na);
    aa = a*a';
    copyto!(Rmout, ca * (1.0*I-aa) + sa/na*thetatilde + aa);
    return Rmout
end
function test()
    for i in 1:10
        Rmout = rand(3, 3)
        a = rand(3)
        # Rmout = [0.056575544213388396 0.8235145288079604 0.4566335160952417; 0.6335215149040396 0.8852762577618685 0.803417762628938; 0.5940995203626245 0.35601280936101065 0.6825575362008556]                                                                                                 
        # a = [0.11011521831193871, 0.5097478695998647, 0.8139760429477749]  
        rotmat3!(Rmout, a)
        Rmout2 = rand(3, 3)
        rotmat3original!(Rmout2, a) 
        # Rmout = [0.5736166081143759 -0.6670428981551146 0.47541325067375245; 0.7189364978350996 0.6881251218379588 0.0980516638109532; -0.3925484670406451 
        # 0.2855478746485778 0.8742814834523946]      
        @test norm(Rmout-Rmout2) <= 1.0e-6
        # @btime rotmat3!($Rmout, $a)
    end
    true
end
end
using .mrotationmatrixs2
mrotationmatrixs2.test()

