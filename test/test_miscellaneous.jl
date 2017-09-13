module mmiscellaneous1mmm
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
using mmiscellaneous1mmm
mmiscellaneous1mmm.test()

module mstressconversionm
using FinEtools
using Base.Test
function test()
  symmtens(N) = begin t=rand(N, N); t = (t+t')/2.0; end
  t = symmtens(2)
  v = zeros(3)
  strain2x2tto3v!(v, t)
  to = zeros(2, 2)
  strain3vto2x2t!(to, v)
  @test norm(t-to) < eps(1.0)

  t = symmtens(3)
  v = zeros(6)
  strain3x3tto6v!(v, t)
  to = zeros(3, 3)
  strain6vto3x3t!(to, v)
  @test norm(t-to) < eps(1.0)

  t = symmtens(2)
  v = zeros(3)
  stress2x2to3v!(v, t)
  to = zeros(2, 2)
  stress3vto2x2t!(to, v)
  @test norm(t-to) < eps(1.0)

  v = vec([1. 2. 3.])
  t = zeros(3, 3)
  stress3vto3x3t!(t, v)
  to = [1. 3. 0; 3. 2. 0; 0 0 0]
  @test norm(t-to) < eps(1.0)

  v = vec([1. 2 3 4])
  t = zeros(3, 3)
  stress4vto3x3t!(t, v)
  to = [1. 3 0; 3 2 0; 0 0 4]
  @test norm(t-to) < eps(1.0)

  v = rand(6)
  t = zeros(3, 3)
  stress6vto3x3t!(t, v)
  vo = zeros(6)
  stress3x3tto6v!(vo, t)
  @test norm(v-vo) < eps(1.0)

  v = rand(9)
  t = zeros(3, 3)
  strain9vto3x3t!(t, v)
  t = (t + t')/2.0 # symmetrize
  strain3x3tto9v!(v, t)
  v6 = zeros(6)
  strain9vto6v!(v6, v)
  v9 = zeros(9)
  strain6vto9v!(v9, v6)
  @test norm(v-v9) < eps(1.0)

  v = vec([1. 2 3 4 4 5 5 6 6])
  v6 = zeros(6)
  stress9vto6v!(v6, v)
  v9 = zeros(9)
  stress6vto9v!(v9, v6)
  @test norm(v-v9) < eps(1.0)

end
end
using mstressconversionm
mstressconversionm.test()

module mmtwistedeexportmm
using FinEtools
using Base.Test
using FinEtools.MeshExportModule
function test()
  E = 0.29e8;
  nu = 0.22;
  W = 1.1;
  L = 12.;
  t =  0.32;
  nl = 2; nt = 1; nw = 1; ref = 3;
  p =   1/W/t;
  #  Loading in the Z direction
  loadv = [0;0;p]; dir = 3; uex = 0.005424534868469; # Harder: 5.424e-3;
  #   Loading in the Y direction
  #loadv = [0;p;0]; dir = 2; uex = 0.001753248285256; # Harder: 1.754e-3;
  tolerance  = t/1000;

  fens,fes  = H8block(L,W,t, nl*ref,nw*ref,nt*ref)

  # Reshape into a twisted beam shape
  for i = 1:count(fens)
    a = fens.xyz[i,1]/L*(pi/2); y = fens.xyz[i,2]-(W/2); z = fens.xyz[i,3]-(t/2);
    fens.xyz[i,:] = [fens.xyz[i,1],y*cos(a)-z*sin(a),y*sin(a)+z*cos(a)];
  end

  # Clamped end of the beam
  l1  = selectnode(fens; box = [0 0 -100*W 100*W -100*W 100*W], inflate  =  tolerance)
  e1 = FDataDict("node_list"=>l1, "component"=>1, "displacement"=>0.0)
  e2 = FDataDict("node_list"=>l1, "component"=>2, "displacement"=>0.0)
  e3 = FDataDict("node_list"=>l1, "component"=>3, "displacement"=>0.0)

  # Traction on the opposite edge
  boundaryfes  =   meshboundary(fes);
  Toplist   = selectelem(fens,boundaryfes, box =  [L L -100*W 100*W -100*W 100*W], inflate =   tolerance);
  el1femm  = FEMMBase(GeoD(subset(boundaryfes,Toplist), GaussRule(2, 2)))
  flux1 = FDataDict("femm"=>el1femm, "traction_vector"=>loadv)


  # Make the region
  MR = DeforModelRed3D
  material = MatDeforElastIso(MR, 00.0, E, nu, 0.0)
  region1 = FDataDict("femm"=>FEMMDeforLinearMSH8(MR, GeoD(fes, GaussRule(3,2)),
            material))

  # Make model data
  modeldata =  FDataDict(
  "fens"=> fens, "regions"=>  [region1],
  "essential_bcs"=>[e1, e2, e3], "traction_bcs"=>  [flux1])


  AE = AbaqusExporter("twisted_beam");
  HEADING(AE, "Twisted beam example");
  PART(AE, "part1");
  END_PART(AE);
  ASSEMBLY(AE, "ASSEM1");
  INSTANCE(AE, "INSTNC1", "PART1");
  NODE(AE, fens.xyz);
  ELEMENT(AE, "c3d8rh", "AllElements", 1, region1["femm"].geod.fes.conn)
  ELEMENT(AE, "SFM3D4", "TractionElements",
    1+count(region1["femm"].geod.fes), flux1["femm"].geod.fes.conn)
  NSET_NSET(AE, "l1", l1)
  ORIENTATION(AE, "GlobalOrientation", vec([1. 0 0]), vec([0 1. 0]));
  SOLID_SECTION(AE, "elasticity", "GlobalOrientation", "AllElements", "Hourglass");
  SURFACE_SECTION(AE, "TractionElements")
  END_INSTANCE(AE);
  END_ASSEMBLY(AE);
  MATERIAL(AE, "elasticity")
  ELASTIC(AE, E, nu)
  SECTION_CONTROLS(AE, "section1", "HOURGLASS=ENHANCED")
  STEP_PERTURBATION_STATIC(AE)
  BOUNDARY(AE, "ASSEM1.INSTNC1.l1", 1)
  BOUNDARY(AE, "ASSEM1.INSTNC1.l1", 2)
  BOUNDARY(AE, "ASSEM1.INSTNC1.l1", 3)
  DLOAD(AE, "ASSEM1.INSTNC1.TractionElements", vec(flux1["traction_vector"]))
  END_STEP(AE)
  close(AE)
  nlines = 0
  open("twisted_beam.inp") do f
    s = readlines(f)
    nlines = length(s)
  end
  @test nlines == 223
  rm("twisted_beam.inp")

  true
end
end
using mmtwistedeexportmm
mmtwistedeexportmm.test()


module mmtwistedeexport2mm
using FinEtools
using Base.Test
using FinEtools.MeshExportModule
function test()
  E = 0.29e8;
  nu = 0.22;
  W = 1.1;
  L = 12.;
  t =  0.32;
  nl = 2; nt = 1; nw = 1; ref = 3;
  p =   1/W/t;
  #  Loading in the Z direction
  loadv = [0;0;p]; dir = 3; uex = 0.005424534868469; # Harder: 5.424e-3;
  #   Loading in the Y direction
  #loadv = [0;p;0]; dir = 2; uex = 0.001753248285256; # Harder: 1.754e-3;
  tolerance  = t/1000;

  fens,fes  = H8block(L,W,t, nl*ref,nw*ref,nt*ref)

  # Reshape into a twisted beam shape
  for i = 1:count(fens)
    a = fens.xyz[i,1]/L*(pi/2); y = fens.xyz[i,2]-(W/2); z = fens.xyz[i,3]-(t/2);
    fens.xyz[i,:] = [fens.xyz[i,1],y*cos(a)-z*sin(a),y*sin(a)+z*cos(a)];
  end

  # Clamped end of the beam
  l1  = selectnode(fens; box = [0 0 -100*W 100*W -100*W 100*W], inflate  =  tolerance)
  e1 = FDataDict("node_list"=>l1, "component"=>1, "displacement"=>0.0)
  e2 = FDataDict("node_list"=>l1, "component"=>2, "displacement"=>0.0)
  e3 = FDataDict("node_list"=>l1, "component"=>3, "displacement"=>0.0)

  # Traction on the opposite edge
  boundaryfes  =   meshboundary(fes);
  Toplist   = selectelem(fens,boundaryfes, box =  [L L -100*W 100*W -100*W 100*W], inflate =   tolerance);
  el1femm  = FEMMBase(GeoD(subset(boundaryfes,Toplist), GaussRule(2, 2)))
  flux1 = FDataDict("femm"=>el1femm, "traction_vector"=>loadv)


  # Make the region
  MR = DeforModelRed3D
  material = MatDeforElastIso(MR, 00.0, E, nu, 0.0)
  region1 = FDataDict("femm"=>FEMMDeforLinearMSH8(MR, GeoD(fes, GaussRule(3,2)),
            material))

  # Make model data
  modeldata =  FDataDict(
  "fens"=> fens, "regions"=>  [region1],
  "essential_bcs"=>[e1, e2, e3], "traction_bcs"=>  [flux1])


  AE = AbaqusExporter("twisted_beam");
  HEADING(AE, "Twisted beam example");
  PART(AE, "part1");
  END_PART(AE);
  ASSEMBLY(AE, "ASSEM1");
  INSTANCE(AE, "INSTNC1", "PART1");
  NODE(AE, fens.xyz);
  ELEMENT(AE, "c3d8rh", "AllElements", 1, region1["femm"].geod.fes.conn)
  ELEMENT(AE, "SFM3D4", "TractionElements",
    1+count(region1["femm"].geod.fes), flux1["femm"].geod.fes.conn)
  NSET_NSET(AE, "l1", l1)
  ORIENTATION(AE, "GlobalOrientation", vec([1. 0 0]), vec([0 1. 0]));
  SOLID_SECTION(AE, "elasticity", "GlobalOrientation", "AllElements", "Hourglass");
  SURFACE_SECTION(AE, "TractionElements")
  END_INSTANCE(AE);
  END_ASSEMBLY(AE);
  MATERIAL(AE, "elasticity")
  ELASTIC(AE, E, nu)
  SECTION_CONTROLS(AE, "section1", "HOURGLASS=ENHANCED")
  STEP_PERTURBATION_STATIC(AE)
  BOUNDARY(AE, "ASSEM1.INSTNC1.l1", 1, 0.0)
  BOUNDARY(AE, "ASSEM1.INSTNC1.l1", 2, 0.0)
  BOUNDARY(AE, "ASSEM1.INSTNC1.l1", 3, 0.0)
  DLOAD(AE, "ASSEM1.INSTNC1.TractionElements", vec(flux1["traction_vector"]))
  END_STEP(AE)
  close(AE)
  nlines = 0
  open("twisted_beam.inp") do f
    s = readlines(f)
    nlines = length(s)
  end
  @test nlines == 223
  rm("twisted_beam.inp")

  true
end
end
using mmtwistedeexport2mm
mmtwistedeexport2mm.test()

module mrichmmm
using FinEtools
using FinEtools.AlgoBaseModule
using Base.Test
function test()
  xs = [93.0734, 92.8633, 92.7252]
    hs = [0.1000, 0.0500, 0.0250]
  solnestim, beta, c, residual = richextrapol(xs, hs)
  # println("$((solnestim, beta, c))")
  @test norm([solnestim, beta, c] - [92.46031652777476, 0.6053628424093497, -2.471055221256022]) < 1.e-5

  sol = [231.7, 239.1, 244.8]
  h = [(4.0/5.)^i for i in 0:1:2]
  solnestim, beta, c, residual = richextrapol(sol, h)
  # println("$((solnestim, beta, c))")
  @test norm([solnestim, beta, c] - [263.91176470588067, 1.1697126080157385, 32.21176470588068]) < 1.e-5

end
end
using mrichmmm
mrichmmm.test()

module mmmeasurementm1
using FinEtools
using Base.Test
function test()
    W = 1.1;
    L = 12.;
    t =  0.32;
    nl, nt, nw = 2, 3, 4;

    fens,fes  = H8block(L,W,t, nl,nw,nt)
    geom  =  NodalField(fens.xyz)

    femm  =  FEMMBase(GeoD(fes, GaussRule(3, 2)))
    V = integratefunction(femm, geom, (x) ->  1.0)
    @test abs(V - W*L*t)/V < 1.0e-5
end
end
using mmmeasurementm1
mmmeasurementm1.test()

module mmmeasurementm2
using FinEtools
using Base.Test
function test()
    W = 1.1;
    L = 12.;
    t =  0.32;
    nl, nt, nw = 2, 3, 4;

    fens,fes  = H8block(L,W,t, nl,nw,nt)
    geom  =  NodalField(fens.xyz)
    bfes = meshboundary(fes)
    femm  =  FEMMBase(GeoD(bfes, GaussRule(2, 2)))
    S = integratefunction(femm, geom, (x) ->  1.0)
    @test abs(S - 2*(W*L + L*t + W*t))/S < 1.0e-5
end
end
using mmmeasurementm2
mmmeasurementm2.test()


module mmmeasurementm3
using FinEtools
using Base.Test
function test()
    W = 1.1;
    L = 12.;
    t =  0.32;
    nl, nt, nw = 2, 3, 4;

    fens,fes  = Q4block(L,W,nl,nw)
    geom  =  NodalField(fens.xyz)
    bfes = meshboundary(fes)
    femm  =  FEMMBase(GeoD(bfes, GaussRule(1, 2)))
    S = integratefunction(femm, geom, (x) ->  1.0)
    @test abs(S - 2*(W + L))/S < 1.0e-5
end
end
using mmmeasurementm3
mmmeasurementm3.test()


module mmmeasurementm4
using FinEtools
using Base.Test
function test()
    W = 1.1;
    L = 12.;
    t =  0.32;
    nl, nt, nw = 2, 3, 4;

    fens,fes  = L2block(L,nl)
    geom  =  NodalField(fens.xyz)
    bfes = meshboundary(fes)
    femm  =  FEMMBase(GeoD(bfes, PointRule()))
    S = integratefunction(femm, geom, (x) ->  1.0)
    @test abs(S - 2)/S < 1.0e-5
end
end
using mmmeasurementm4
mmmeasurementm4.test()


module mmmeasurementm5
using FinEtools
using Base.Test
function test()
    W = 1.1;
    L = 12.;
    t =  0.32;
    nl, nt, nw = 2, 3, 4;

    fens,fes  = Q4block(L,W,nl,nw)
    geom  =  NodalField(fens.xyz)

    nfes = FESetP1(reshape(collect(1:count(fens)), count(fens), 1))
    femm  =  FEMMBase(GeoD(nfes, PointRule()))
    S = integratefunction(femm, geom, (x) ->  1.0)
    @test abs(S - count(fens))/S < 1.0e-5

end
end
using mmmeasurementm5
mmmeasurementm5.test()


module mmmeasurementm6
using FinEtools
using Base.Test
function test()
    W = 1.1;
    L = 12.;
    t =  0.32;
    nl, nt, nw = 2, 3, 4;

    fens,fes  = H8block(L,W,t, nl,nw,nt)
    geom  =  NodalField(fens.xyz)

    nfes = FESetP1(reshape(collect(1:count(fens)), count(fens), 1))
    femm  =  FEMMBase(GeoD(nfes, PointRule()))
    S = integratefunction(femm, geom, (x) ->  1.0)
    @test abs(S - count(fens))/S < 1.0e-5

end
end
using mmmeasurementm6
mmmeasurementm6.test()


module mmmeasurementm7
using FinEtools
using Base.Test
function test()
    W = 1.1;
    L = 12.;
    t =  0.32;
    nl, nt, nw = 2, 3, 4;

    fens,fes  = L2block(L,nl)
    geom  =  NodalField(fens.xyz)
    bfes = FESetP1(reshape([nl+1], 1, 1))
    axisymmetric = true
    femm  =  FEMMBase(GeoD(bfes, PointRule(), axisymmetric))
    S = integratefunction(femm, geom, (x) ->  1.0, 1)
    # println(" Length  of the circle = $(S)")
    @test abs(S - 2*pi*L)/S < 1.0e-5
end
end
using mmmeasurementm7
mmmeasurementm7.test()

module mmmeasurementm8
using FinEtools
using Base.Test
function test()
    W = 1.1;
    L = 12.;
    t =  0.32;
    nl, nt, nw = 2, 3, 4;

    fens,fes  = L2block(L,nl)
    geom  =  NodalField(fens.xyz)
    bfes = FESetP1(reshape([nl+1], 1, 1))

    femm  =  FEMMBase(GeoD(bfes, PointRule(), t))
    S = integratefunction(femm, geom, (x) ->  1.0, 1)
    # println("Length  of the boundary curve = $(S)")
    @test abs(S - t)/S < 1.0e-5
end
end
using mmmeasurementm8
mmmeasurementm8.test()

module mmmeasurementm9
using FinEtools
using Base.Test
function test()
    W = 1.1;
    L = 12.;
    t =  0.32;
    nl, nt, nw = 2, 3, 4;

    fens,fes  = L2block(L,nl)
    geom  =  NodalField(fens.xyz)
    bfes = FESetP1(reshape([nl+1], 1, 1))

    femm  =  FEMMBase(GeoD(bfes, PointRule(), t*W))
    S = integratefunction(femm, geom, (x) ->  1.0, 2)
    # println("Length  of the boundary curve = $(S)")
    @test abs(S - t*W)/S < 1.0e-5
end
end
using mmmeasurementm9
mmmeasurementm9.test()

module mmmeasurementm10
using FinEtools
using Base.Test
function test()
    W = 1.1;
    L = 12.;
    t =  0.32;
    nl, nt, nw = 2, 3, 4;

    fens,fes  = L2block(L,nl)
    geom  =  NodalField(fens.xyz)
    bfes = FESetP1(reshape([nl+1], 1, 1))
axisymmetric = true
    femm  =  FEMMBase(GeoD(bfes, PointRule(), axisymmetric, W))
    S = integratefunction(femm, geom, (x) ->  1.0, 2)
    # println("Length  of the boundary curve = $(S)")
    @test abs(S - 2*pi*L*W)/S < 1.0e-5
end
end
using mmmeasurementm10
mmmeasurementm10.test()

module mmmeasurementm11
using FinEtools
using Base.Test
function test()
    W = 1.1;
    L = 12.;
    t =  0.32;
    nl, nt, nw = 2, 3, 4;

    fens,fes  = L2block(L,nl)
    geom  =  NodalField(fens.xyz)
    bfes = FESetP1(reshape([nl+1], 1, 1))
axisymmetric = true
    femm  =  FEMMBase(GeoD(bfes, PointRule(), axisymmetric, W*t))
    S = integratefunction(femm, geom, (x) ->  1.0, 3)
    # println("Length  of the boundary curve = $(S)")
    @test abs(S - 2*pi*L*W*t)/S < 1.0e-5
end
end
using mmmeasurementm11
mmmeasurementm11.test()

module mmmeasurementm12
using FinEtools
using Base.Test
function test()
    W = 1.1;
    L = 12.;
    t =  0.32;
    nl, nt, nw = 2, 3, 4;

    fens,fes  = L2block(L,nl)
    geom  =  NodalField(fens.xyz)
    bfes = FESetP1(reshape([nl+1], 1, 1))
axisymmetric = false
    femm  =  FEMMBase(GeoD(bfes, PointRule(), axisymmetric, L*W*t))
    S = integratefunction(femm, geom, (x) ->  1.0, 3)
    # println("Length  of the boundary curve = $(S)")
    @test abs(S - L*W*t)/S < 1.0e-5
end
end
using mmmeasurementm12
mmmeasurementm12.test()

module mmpartitioning1m
using FinEtools
using Base.Test
function test()
    a = 10. # radius of the hole
    nC = 20
    nR = 4
    fens,fes = Q4annulus(a, 1.5*a, nR, nC, 2*pi)

    npartitions = 8
    partitioning = nodepartitioning(fens, npartitions)
    partitionnumbers = unique(partitioning)
    @test norm(partitionnumbers - [1
    3
    2
    5
    6
    7
    8
    4]) < 1.e-5
end
end
using mmpartitioning1m
mmpartitioning1m.test()


module mmpartitioning2m
using FinEtools
using Base.Test
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
    @test norm(partitionnumbers - [1
    3
    4
    2]) < 1.e-5
end
end
using mmpartitioning2m
mmpartitioning2m.test()

module mmboxm1
using FinEtools
using Base.Test
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
    # println("$(b3)")
    @test norm(b3 - [-1.0, 0.999928, 0.248077, 3.0, -0.5, 0.93455]) < 1.0e-4
    x = [0.25 1.1 -0.3]
    @test inbox(b3, x)
    @test inbox(b3, c)
    @test inbox(b3, a[2, :])
    b4 = boundingbox(-a)
    # println("$(b3)")
    # println("$(b4)")
    # println("$(boxesoverlap(b3, b4))")
    @test !boxesoverlap(b3, b4)
    b5 = updatebox!(b3, [0.0 -0.4 0.0])
    # println("$(b5)")
    # println("$(boxesoverlap(b5, b4))")
    @test boxesoverlap(b5, b4)
end
end
using mmboxm1
mmboxm1.test()


module mmMeasurement_1
using FinEtools
using Base.Test
function test()
    W = 1.1;
    L = 12.;
    t =  3.32;
    nl, nt, nw = 5, 3, 4;

    for or in [:a :b :ca :cb]
        fens,fes  = T4block(L,W,t, nl,nw,nt, or)
        geom  =  NodalField(fens.xyz)

        femm  =  FEMMBase(GeoD(fes, TetRule(5)))
        V = integratefunction(femm, geom, (x) ->  1.0)
        @test abs(V - W*L*t)/V < 1.0e-5
    end

end
end
using mmMeasurement_1
mmMeasurement_1.test()

module mphunm2
using FinEtools
using Base.Test
function test()
    for btu in [:SEC :MIN :HR :DY :YR :WK]
        t = 0.333*phun(base_time_units = btu, "s")
        v = 2.0*phun(base_time_units = btu, "m/s")
        # display(v*t)
        @test abs(v*t - 0.666) < 1.0e-6
    end
end
end
using mphunm2
mphunm2.test()

module mphunm3
using FinEtools
using Base.Test
function test()
    for sou in [:SI :US :IMPERIAL :CGS :SIMM]
        t = 0.333*phun(system_of_units = sou, "s")
        v = 2.0*phun(system_of_units = sou, "m/s")
        # display(v*t/phun(system_of_units = sou, "m"))
        @test abs(v*t/phun(system_of_units = sou, "m") - 0.666) < 1.0e-6
    end
end
end
using mphunm3
mphunm3.test()


module mxmeasurementm1
using FinEtools
using Base.Test
function test()
    W = 1.1;
    L = 12.;
    t =  0.32;
    nl, nt, nw = 2, 3, 4;

    fens,fes  = H27block(L,W,t, nl,nw,nt)
    geom  =  NodalField(fens.xyz)

    femm  =  FEMMBase(GeoD(fes, GaussRule(3, 3)))
    V = integratefunction(femm, geom, (x) ->  1.0)
    @test abs(V - W*L*t)/V < 1.0e-5
end
end
using mxmeasurementm1
mxmeasurementm1.test()

module mxmeasurementm2
using FinEtools
using Base.Test
function test()
    W = 1.1;
    L = 12.;
    t =  0.32;
    nl, nt, nw = 2, 3, 4;

    fens,fes  = H27block(L,W,t, nl,nw,nt)
    geom  =  NodalField(fens.xyz)

bfes = meshboundary(fes)
    femm  =  FEMMBase(GeoD(bfes, GaussRule(2, 4)))
    V = integratefunction(femm, geom, (x) ->  1.0)
    @test abs(V - 2*(W*L+L*t+W*t))/V < 1.0e-5
end
end
using mxmeasurementm2
mxmeasurementm2.test()

module mxmeasurementm3
using FinEtools
using FinEtools.MeshExportModule
using Base.Test
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
    femm  =  FEMMBase(GeoD(bfes, GaussRule(2, 4)))
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
    femm  =  FEMMBase(GeoD(bfes, GaussRule(2, 4)))
    S27 = integratefunction(femm, geom, (x) ->  1.0)
    # println("$((S20-S27)/S20)")

    @test abs((S20-S27)/S20) < 1.0e-5
end
end
using mxmeasurementm3
mxmeasurementm3.test()

module mxmeasurementm4
using FinEtools
using FinEtools.MeshExportModule
using Base.Test
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

    femm  =  FEMMBase(GeoD(fes, GaussRule(3, 4)))
    V20 = integratefunction(femm, geom, (x) ->  1.0)

    fens,fes  = H27block(L,W,t, nl,nw,nt)
    for ixxxx = 1:count(fens)
        x,y,z = fens.xyz[ixxxx, :]
        fens.xyz[ixxxx, :] = [x+a*sin(y) y+a*sin(z) z+a*sin(x)]
    end
    geom  =  NodalField(fens.xyz)
    # File = "mesh.vtk"
    # MeshExportModule.vtkexportmesh(File, fens, fes)

    femm  =  FEMMBase(GeoD(fes, GaussRule(3, 4)))
    V27 = integratefunction(femm, geom, (x) ->  1.0)
    # println("$((S20-S27)/S20)")

    @test abs((V20-V27)/V20) < 1.0e-5
end
end
using mxmeasurementm4
mxmeasurementm4.test()

module mmmiscellaneous2
using FinEtools
using FinEtools.MeshExportModule
using Base.Test
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
  @test bfes.conn == [1 2; 5 1; 2 3; 3 4; 4 8; 9 5; 8 12; 10 9; 11 10; 12 11]
end
end
using mmmiscellaneous2
mmmiscellaneous2.test()

module mmmiscellaneous3
using FinEtools
using FinEtools.MeshExportModule
using Base.Test
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
  @test vec(bfes.conn) == vec([1 5])
end
end
using mmmiscellaneous3
mmmiscellaneous3.test()

module mmmfieldmm1
using FinEtools
using Base.Test
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
    copy!(u, geom)
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
using mmmfieldmm1
mmmfieldmm1.test()

module mmcrossm
using FinEtools
using FinEtools.RotationUtilModule
using Base.Test
function test()
    a = vec([0.1102, -0.369506, -0.0167305])
    b = vec([0.0824301, -0.137487, 0.351721])
    c = cross(a, b)
    # println("$(c)")
    A = zeros(3, 3)
    skewmat!(A, a)
    d = A * b
    # println("$(d)")
    @test norm(c - d) < 1.0e-6
    e = RotationUtilModule.cross3(a, b)
    @test norm(c - e) < 1.0e-6
    f = zeros(3)
    RotationUtilModule.cross3!(f, a, b)
    @test norm(c - f) < 1.0e-6

    a = vec([0.1102, -0.135])
    b = vec([-0.137487, 0.351721])
    c = RotationUtilModule.cross2(a, b)
    # println("$(c)")
    @test norm(c - 0.0201989092) < 1.0e-6
end
end
using mmcrossm
mmcrossm.test()

module mstresscomponentmap
using FinEtools
using Base.Test
function test()
    MR = DeforModelRed1D
    @test stresscomponentmap(MR)[:x] == 1
    MR = DeforModelRed2DAxisymm
    @test stresscomponentmap(MR)[:x] == 1
    @test stresscomponentmap(MR)[:zz] == 3
end
end
using mstresscomponentmap
mstresscomponentmap.test()

module mgen_iso_csmat1
using FinEtools
using FinEtools.FESetModule
using FinEtools.CSysModule
using Base.Test
function test()
    L = 2.0
    nl = 1

    fens,fes  = L2block(L, nl)
    fens.xyz = xyz3(fens)
    fens.xyz[2, 3] += L/2
    csmatout = zeros(FFlt, 3, 1)
    gradNparams = FESetModule.bfundpar(fes, vec([0.0]));
    J = transpose(fens.xyz) * gradNparams
    CSysModule.gen_iso_csmat!(csmatout, mean(fens.xyz, 1), J, 0)
    # println("$(csmatout)")
    # println("$(norm(vec(csmatout)))")
    @test norm(csmatout - [0.894427; 0.0; 0.447214]) < 1.0e-5
end
end
using mgen_iso_csmat1
mgen_iso_csmat1.test()

module mgen_iso_csmat2
using FinEtools
using FinEtools.FESetModule
using FinEtools.CSysModule
using FinEtools.MeshExportModule
using Base.Test
function test()
    L = 2.0
    nl = 1

    fens,fes  = Q4block(L, 2*L, nl, nl)
    fens.xyz = xyz3(fens)
    fens.xyz[2, 3] += L/2
    File = "mesh.vtk"
    MeshExportModule.vtkexportmesh(File, fens, fes)
    csmatout = zeros(FFlt, 3, 2)
    gradNparams = FESetModule.bfundpar(fes, vec([0.0 0.0]));
    J = transpose(fens.xyz[fes.conn[1, :], :]) * gradNparams
    # println("J = $(J)")
    @test norm(J - [1.0 0.0; 0.0 2.0; 0.25 -0.25]) < 1.0e-5
    CSysModule.gen_iso_csmat!(csmatout, mean(fens.xyz, 1), J, 0)
    # println("csmatout = $(csmatout)")
    @test norm(csmatout - [0.970143 0.0291979; 0.0 0.992727; 0.242536 -0.116791]) < 1.0e-5
end
end
using mgen_iso_csmat2
mgen_iso_csmat2.test()

module mmfieldx1
using FinEtools
using Base.Test
function test()
    f = NodalField(zeros(5, 1))
    setebc!(f, [3,4], true, 1;        val=7.0)
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
    setebc!(f, [3,4], 1; val=8.2)
    # display(f)
    applyebc!(f)
    dest = zeros(2,1)
    # length([1,4])
    gatherfixedvalues_asmat!(f, dest, [1,4])
    # display(dest)
    @test (dest[1,1] == 0.0) && (dest[2,1] == 8.2)
end
end
using mmfieldx1
mmfieldx1.test()


module mmconnection1
using FinEtools
using Base.Test
function test()
    h = 0.05*phun("M");
    l = 10*h;
    nh = 3; nl  = 4; nc = 5;

    fens,fes  = H8block(h,l,2.0*pi,nh,nl,nc)
    femm = FEMMBase(GeoD(fes, GaussRule(3, 2)))
    C = connectionmatrix(femm, count(fens))
    Degree = [length(find(C[j,:])) for j in 1:size(C, 1)]
    # println("Maximum degree  = $(maximum(Degree))")
    # println("Minimum degree  = $(minimum(Degree))")
    @test maximum(Degree) == 27
    @test minimum(Degree) == 8
end
end
using mmconnection1
mmconnection1.test()

module mmquadrature3
using FinEtools
using Base.Test
function test()
    gr = GaussRule(1,1)
    @test gr.param_coords[1] == 0.0
    @test_throws ErrorException gr = GaussRule(1,5)
    @test_throws AssertionError gr = GaussRule(4,1)

    @test_throws ErrorException tr = TriRule(4)

    @test_throws ErrorException tr = TetRule(3)
end
end
using mmquadrature3
mmquadrature3.test()

module mmmatchingmm
using FinEtools.AlgoBaseModule: FDataDict, dcheck!
using Base.Test
function test()
    d = FDataDict("postp"=>true, "preprocessing"=>nothing)
    recognized_keys = ["postprocessing", "something",  "else"]
    notmatched = dcheck!(d, recognized_keys)
    # display(notmatched)
    @test !isempty(notmatched)
end
end
using mmmatchingmm
mmmatchingmm.test()


module mboundary13
using FinEtools
using FinEtools.MeshExportModule
using Base.Test
function test()
    xs = collect(linspace( 1.0, 3.0, 4))
    ys = collect(linspace(-1.0, 5.0, 5))
    fens, fes = T3blockx(xs, ys, :a)
    @test count(fes) == 4*3*2

    boundaryfes  =   meshboundary(fes);
    ytl   = selectelem(fens, boundaryfes, facing = true,
        direction = [0.0 +1.0], dotmin = 0.999);
    @test length(ytl) == 3

    geom  =  NodalField(fens.xyz)

    femm  =  FEMMBase(GeoD(subset(boundaryfes, ytl), GaussRule(1, 2)))
    L = integratefunction(femm, geom, (x) ->  1.0)
    @test abs(L - 2.0)/L < 1.0e-5

    # File = "playground.vtk"
    # MeshExportModule.vtkexportmesh(File, fens, fes)
    # @async run(`"paraview.exe" $File`)
    # try rm(File) catch end

end
end
using mboundary13
mboundary13.test()


module mboundary14
using FinEtools
using FinEtools.MeshExportModule
using Base.Test
function test()
    xs = collect(linspace( 1.0, 3.0, 4))
    ys = collect(linspace(-1.0, 5.0, 5))
    fens, fes = Q4blockx(xs, ys)
    fens, fes = Q4toQ8(fens, fes)
    @test count(fes) == 4*3

    boundaryfes  =   meshboundary(fes);
    ytl   = selectelem(fens, boundaryfes, facing = true,
        direction = [0.0 +1.0], dotmin = 0.999);
    @test length(ytl) == 3

    geom  =  NodalField(fens.xyz)

    femm  =  FEMMBase(GeoD(subset(boundaryfes, ytl), GaussRule(1, 3)))
    L = integratefunction(femm, geom, (x) ->  1.0)
    @test abs(L - 2.0)/L < 1.0e-5

    # File = "playground.vtk"
    # MeshExportModule.vtkexportmesh(File, fens, fes)
    # @async run(`"paraview.exe" $File`)
    # try rm(File) catch end

end
end
using mboundary14
mboundary14.test()


module mboundary15
using FinEtools
using FinEtools.MeshExportModule
using Base.Test
function test()
    xs = collect(linspace( 1.0, 3.0, 4))
    ys = collect(linspace(-1.0, 5.0, 5))
    zs = collect(linspace(+2.0, 5.0, 7))
    fens, fes = T4blockx(xs, ys, zs, :a)
    @test count(fes) == 432

    boundaryfes  =   meshboundary(fes);
    ytl   = selectelem(fens, boundaryfes, facing = true,
        direction = [0.0 +1.0 0.0], dotmin = 0.999);
    @test length(ytl) == 36

    geom  =  NodalField(fens.xyz)

    femm  =  FEMMBase(GeoD(subset(boundaryfes, ytl), SimplexRule(2, 1)))
    S = integratefunction(femm, geom, (x) ->  1.0)
    @test abs(S - 6.0)/S < 1.0e-5


    femm  =  FEMMBase(GeoD(fes, SimplexRule(3, 1)))
    V = integratefunction(femm, geom, (x) ->  1.0)
    @test abs(V - 36.0)/V < 1.0e-5

    # File = "playground.vtk"
    # MeshExportModule.vtkexportmesh(File, fens, fes)
    # @async run(`"paraview.exe" $File`)
    # try rm(File) catch end

end
end
using mboundary15
mboundary15.test()
