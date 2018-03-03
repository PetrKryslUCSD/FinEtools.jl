module Meyer_Piening_examples
using FinEtools
using FinEtools.AlgoDeforLinearModule
using FinEtools.MeshUtilModule
function Meyer_Piening_sandwich()
    println("""
    Meyer-Piening sandwich plate
    """)
    
    # Reference results from:
    # [1] Application of the Elasticity Solution
    # to Linear Sandwich Beam, Plate
    # and Shell Analyses
    # H.-R. MEYER -PIENING
    # Journal of SANDWICH STRUCTURES AND MATERIALS , Vol. 6—July 2004
    
    # Assessment of the refined sinus plate finite element:
    # Free edge effect and Meyer-Piening sandwich test
    # P. Vidal, O. Polit, M. D'Ottavio, E. Valot
    # http://dx.doi.org/10.1016/j.finel.2014.08.004
    #
    # The second study deals with a benchmark problem proposed
    # by Meyer-Piening [14]. It involves a simply -supported rectangular
    # sandwich plate submitted to a localized pressure applied on an
    # area of 5x20 mm. The geometry of the sandwich structure is
    # given in Fig.13. Due to the symmetry, only one quarter of the plate
    # is meshed. The faces have different thicknesses: h = 0.5 mm
    # (bottom face), h = 0.1 mm (top face). The thickness of the core
    # is h =  11.4 mm. The material properties are given in Table 3.
    # Note that this benchmark involves strong heterogeneities (very
    # different geometric and constitutive properties between core and
    # face) and local stress gradient due to the localized pressure load.
    #
    # [14] H.-R. Meyer-Piening, Experiences with exact linear sandwich beam and plate
    # analyses regarding bending, instability and frequency investigations, in:
    # Proceedings of the Fifth International Conference On Sandwich Constructions,
    # September 5–7, vol. I, Zurich, Switzerland, 2000, pp. 37–48.
    
    
    
    t0 = time()
    # Orthotropic material for the SKIN
    E1s = 70000.0*phun("MPa")
    E2s = 71000.0*phun("MPa")
    E3s = 69000.0*phun("MPa")
    nu12s = nu13s = nu23s = 0.3
    G12s = G13s = G23s = 26000.0*phun("MPa")
    CTE1 =  CTE2 =  CTE3 = 0.0
    # Orthotropic material for the CORE
    E1c = 3.0*phun("MPa")
    E2c = 3.0*phun("MPa")
    E3c = 2.8*phun("MPa")
    nu12c = nu13c = nu23c = 0.25
    G12c = G13c = G23c = 1.0*phun("MPa")
    CTE1 =  CTE2 =  CTE3 = 0.0
    
    Lx = 5.0*phun("mm") # length  of loaded rectangle
    Ly = 20.0*phun("mm") # length  of loaded rectangle
    Sx = 100.0*phun("mm") # span of the plate
    Sy = 200.0*phun("mm") # span of the plate
    
    # Here we define the layout and the thicknesses of the layers.
    angles = vec([0.0 0.0 0.0]);
    ts = vec([0.5  11.4  0.1])*phun("mm"); # layer thicknesses
    TH = sum(ts); # total thickness of the plate
    
    tolerance = 0.0001*TH
    
    # The line load is in the negative Z direction.
    q0 = 1*phun("MPa"); #    line load
    
    # Reference deflection under the load is
    wtopref = -3.789*phun("mm"); # From [1]
    wbottomref = -2.16*phun("mm"); # Not given in [1]; guessed from the figure
    
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
    
    fens,fes = H8layeredplatex(xs, ys, ts, nts)
    
    
    # This is the material  model
    MR = DeforModelRed3D
    skinmaterial = MatDeforElastOrtho(MR,
    0.0, E1s, E2s, E3s,
    nu12s, nu13s, nu23s,
    G12s, G13s, G23s,
    CTE1, CTE2, CTE3)
    corematerial = MatDeforElastOrtho(MR,
    0.0, E1c, E2c, E3c,
    nu12c, nu13c, nu23c,
    G12c, G13c, G23c,
    CTE1, CTE2, CTE3)
    
    # The material coordinate system function is defined as:
    function updatecs!(csmatout::FFltMat, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)
        rotmat3!(csmatout, angles[fe_label]/180.0*pi* [0.0; 0.0; 1.0]);
    end
    
    # The vvolume integrals are evaluated using this rule
    gr = GaussRule(3, 2)
    
    # We will create two regions, one for the skin,
    # and one for the core.
    rls = selectelem(fens, fes, label = 1)
    botskinregion = FDataDict("femm"=>FEMMDeforLinearMSH8(MR, IntegData(subset(fes, rls), gr), CSys(3, 3, updatecs!), skinmaterial))
    rls = selectelem(fens, fes, label = 3)
    topskinregion = FDataDict("femm"=>FEMMDeforLinearMSH8(MR, IntegData(subset(fes, rls), gr), CSys(3, 3, updatecs!), skinmaterial))
    rlc = selectelem(fens, fes, label = 2)
    coreregion = FDataDict("femm"=>FEMMDeforLinearMSH8(MR, IntegData(subset(fes, rlc), gr), CSys(3, 3, updatecs!), corematerial))
    
    # File =  "Meyer_Piening_sandwich-r1.vtk"
    # vtkexportmesh(File, skinregion["femm"].integdata.fes.conn, fens.xyz, FinEtools.MeshExportModule.H8)
    # # @async run(`"paraview.exe" $File`)
    # File =  "Meyer_Piening_sandwich-r2.vtk"
    # vtkexportmesh(File, coreregion["femm"].integdata.fes.conn, fens.xyz, FinEtools.MeshExportModule.H8)
    # @async run(`"paraview.exe" $File`)
    
    # The essential boundary conditions are applied on the symmetry planes.
    # First the plane X=0;...
    lx0 = selectnode(fens, box=[0.0 0.0 -Inf Inf -Inf Inf], inflate=tolerance)
    ex0 = FDataDict( "displacement"=>  0.0, "component"=> 1, "node_list"=>lx0 )
    # ... and then the plane Y=0.
    ly0 = selectnode(fens, box=[-Inf Inf 0.0 0.0 -Inf Inf], inflate=tolerance)
    ey0 = FDataDict( "displacement"=>  0.0, "component"=> 2, "node_list"=>ly0 )
    # The transverse displacement is fixed around the circumference.
    lz0 = vcat(selectnode(fens, box=[Sx/2 Sx/2 -Inf Inf -Inf Inf], inflate=tolerance),
    selectnode(fens, box=[-Inf Inf Sy/2 Sy/2 -Inf Inf], inflate=tolerance))
    ez0 = FDataDict( "displacement"=>  0.0, "component"=> 3, "node_list"=>lz0 )
    
    # The traction boundary condition is applied  along rectangle in the middle of the plate.
    bfes = meshboundary(fes)
    # From  the entire boundary we select those quadrilaterals that lie on the plane
    # Z = thickness
    tl = selectelem(fens, bfes, box = [0.0 Lx/2 0 Ly/2 TH TH], inflate=tolerance)
    Trac = FDataDict("traction_vector"=>vec([0.0; 0.0; -q0]),
    "femm"=>FEMMBase(IntegData(subset(bfes, tl), GaussRule(2, 2))))
    
    modeldata = FDataDict("fens"=>fens,
    "regions"=>[botskinregion, coreregion, topskinregion],
    "essential_bcs"=>[ex0, ey0, ez0],
    "traction_bcs"=> [Trac]
    )
    modeldata = AlgoDeforLinearModule.linearstatics(modeldata)
    
    modeldata["postprocessing"] = FDataDict("file"=>"Meyer_Piening_sandwich")
    modeldata = AlgoDeforLinearModule.exportdeformation(modeldata)
    
    u = modeldata["u"]
    geom = modeldata["geom"]
    
    # The results of the displacement and stresses will be reported at
    # nodes located at the appropriate points.
    nbottomcenter = selectnode(fens, box=[0.0 0.0 0.0 0.0 0.0 0.0], inflate=tolerance)
    ntopcenter = selectnode(fens, box=[0.0 0.0 0.0 0.0 TH TH], inflate=tolerance)
    ncenterline = selectnode(fens, box=[0.0 0.0 0.0 0.0 0.0 TH], inflate=tolerance)
    nintertop = selectnode(fens, box=[-Inf Inf 0.0 0.0 sum(ts[1:2]) sum(ts[1:2])], inflate=tolerance)
    ninterbot = selectnode(fens, box=[-Inf Inf 0.0 0.0 sum(ts[1:1]) sum(ts[1:1])], inflate=tolerance)
    
    zclo = sortperm(vec(geom.values[ncenterline, 3]))
    centerz = geom.values[ncenterline[zclo], 3]
    xclotop = sortperm(vec(geom.values[nintertop, 1]))
    topx = geom.values[nintertop[xclotop], 1]
    xclobot = sortperm(vec(geom.values[ninterbot, 1]))
    botx = geom.values[ninterbot[xclobot], 1]
    
    conninbotskin = intersect(connectednodes(botskinregion["femm"].integdata.fes), ncenterline)
    connincore = intersect(connectednodes(coreregion["femm"].integdata.fes), ncenterline)
    connintopskin = intersect(connectednodes(topskinregion["femm"].integdata.fes), ncenterline)
    inbotskin = [n in conninbotskin for n in ncenterline]
    incore = [n in connincore for n in ncenterline]
    intopskin = [n in connintopskin for n in ncenterline]
    
    println("")
    println("Top Center deflection: $(u.values[ntopcenter, 3]/phun("mm")) [mm]")
    println("Bottom Center deflection: $(u.values[nbottomcenter, 3]/phun("mm")) [mm]")
    
    # # extrap = :extrapmean
    # extrap = :extraptrend
    # nodevalmeth = :averaging
    extrap = :default
    nodevalmeth = :invdistance
    
    # Normal stress in the X direction
    modeldata["postprocessing"] = FDataDict("file"=>"Meyer_Piening_sandwich-sx",
    "quantity"=>:Cauchy, "component"=>1, "outputcsys"=>CSys(3),
    "nodevalmethod"=>nodevalmeth, "reportat"=>extrap)
    modeldata = AlgoDeforLinearModule.exportstress(modeldata)
    s = modeldata["postprocessing"]["exported"][1]["field"]
    sxbot = s.values[ncenterline[zclo], 1]
    s = modeldata["postprocessing"]["exported"][2]["field"]
    sxcore = s.values[ncenterline[zclo], 1]
    s = modeldata["postprocessing"]["exported"][3]["field"]
    sxtop = s.values[ncenterline[zclo], 1]
    
    # The graph data needs to be collected by going through each layer separately.
    # Some quantities may be discontinuous between layers.
    zs = vcat(  [z for (j,z) in enumerate(centerz) if inbotskin[j]],
    [z for (j,z) in enumerate(centerz) if incore[j]],
    [z for (j,z) in enumerate(centerz) if intopskin[j]]
    )
    sxs = vcat( [sxbot[j] for (j,z) in enumerate(centerz) if inbotskin[j]],
    [sxcore[j] for (j,z) in enumerate(centerz) if incore[j]],
    [sxtop[j] for (j,z) in enumerate(centerz) if intopskin[j]]
    )
    
    File = "Meyer_Piening_sandwich-sx-$(extrap).CSV"
    savecsv(File, zs=vec(zs)/phun("mm"), sx=vec(sxs)/phun("MPa"))
    
    # @async run(`"paraview.exe" $File`)
    
    # Inter laminar stress between the skin and the core
    modeldata["postprocessing"] = FDataDict("file"=>"Meyer_Piening_sandwich-sxz",
    "quantity"=>:Cauchy, "component"=>5, "outputcsys"=>CSys(3),
    "nodevalmethod"=>nodevalmeth, "reportat"=>extrap)
    modeldata = AlgoDeforLinearModule.exportstress(modeldata)
    s = modeldata["postprocessing"]["exported"][1]["field"]
    sxzskinbot = s.values[ninterbot[xclobot], 1]
    s = modeldata["postprocessing"]["exported"][2]["field"]
    sxzcoretop = s.values[nintertop[xclotop], 1]
    sxzcorebot = s.values[ninterbot[xclobot], 1]
    s = modeldata["postprocessing"]["exported"][3]["field"]
    sxzskintop = s.values[nintertop[xclotop], 1]
    
    File = "Meyer_Piening_sandwich-sxz-$(extrap).CSV"
    savecsv(File, xstop=vec(topx[xclotop])/phun("mm"), sxzskintop=vec(sxzskintop[xclotop])/phun("MPa"), sxzcoretop=vec(sxzcoretop[xclotop])/phun("MPa"),thexsbot=vec(botx[xclobot])/phun("mm"), sxzskinbot=vec(sxzskinbot[xclobot])/phun("MPa"), sxzcorebot=vec(sxzcorebot[xclobot])/phun("MPa"))
    
    println("Done")
    true
    
end # Meyer_Piening_sandwich


function Meyer_Piening_sandwich_H20()
    println("""
    Meyer-Piening sandwich plate, serendipity H20
    """)
    
    # Reference results from:
    # [1] Application of the Elasticity Solution
    # to Linear Sandwich Beam, Plate
    # and Shell Analyses
    # H.-R. MEYER -PIENING
    # Journal of SANDWICH STRUCTURES AND MATERIALS , Vol. 6—July 2004
    
    # Assessment of the refined sinus plate finite element:
    # Free edge effect and Meyer-Piening sandwich test
    # P. Vidal, O. Polit, M. D'Ottavio, E. Valot
    # http://dx.doi.org/10.1016/j.finel.2014.08.004
    #
    # The second study deals with a benchmark problem proposed
    # by Meyer-Piening [14]. It involves a simply -supported rectangular
    # sandwich plate submitted to a localized pressure applied on an
    # area of 5x20 mm. The geometry of the sandwich structure is
    # given in Fig.13. Due to the symmetry, only one quarter of the plate
    # is meshed. The faces have different thicknesses: h = 0.5 mm
    # (bottom face), h = 0.1 mm (top face). The thickness of the core
    # is h =  11.4 mm. The material properties are given in Table 3.
    # Note that this benchmark involves strong heterogeneities (very
    # different geometric and constitutive properties between core and
    # face) and local stress gradient due to the localized pressure load.
    #
    # [14] H.-R. Meyer-Piening, Experiences with exact linear sandwich beam and plate
    # analyses regarding bending, instability and frequency investigations, in:
    # Proceedings of the Fifth International Conference On Sandwich Constructions,
    # September 5–7, vol. I, Zurich, Switzerland, 2000, pp. 37–48.
    
    filebase = "Meyer-Piening-sandwich-H20"
    
    t0 = time()
    # Orthotropic material for the SKIN
    E1s = 70000.0*phun("MPa")
    E2s = 71000.0*phun("MPa")
    E3s = 69000.0*phun("MPa")
    nu12s = nu13s = nu23s = 0.3
    G12s = G13s = G23s = 26000.0*phun("MPa")
    CTE1 =  CTE2 =  CTE3 = 0.0
    # Orthotropic material for the CORE
    E1c = 3.0*phun("MPa")
    E2c = 3.0*phun("MPa")
    E3c = 2.8*phun("MPa")
    nu12c = nu13c = nu23c = 0.25
    G12c = G13c = G23c = 1.0*phun("MPa")
    CTE1 =  CTE2 =  CTE3 = 0.0
    
    Lx = 5.0*phun("mm") # length  of loaded rectangle
    Ly = 20.0*phun("mm") # length  of loaded rectangle
    Sx = 100.0*phun("mm") # span of the plate
    Sy = 200.0*phun("mm") # span of the plate
    
    # Here we define the layout and the thicknesses of the layers.
    angles = vec([0.0 0.0 0.0]);
    ts = vec([0.5  11.4  0.1])*phun("mm"); # layer thicknesses
    TH = sum(ts); # total thickness of the plate
    
    tolerance = 0.0001*TH
    
    # The line load is in the negative Z direction.
    q0 = 1*phun("MPa"); #    line load
    
    # Reference deflection under the load is
    wtopref = -3.789*phun("mm"); # From [1]
    wbottomref = -2.16*phun("mm"); # Not given in [1]; guessed from the figure
    
    # Select how find the mesh should be
    Refinement = 3
    nL = Refinement * 1;
    nSx = nL + Refinement * 4;
    nSy = 2 * nSx;
    
    # Each layer is modeled with a single element.
    nts= Refinement * [1, 2, 1];# number of elements per layer
    strength = 1.5
    sp = (a, b, n) -> MeshUtilModule.gradedspace(a, b, n, strength)
    sp = (a, b, n) -> linearspace(a, b, n)
    xs = unique(vcat(reverse(collect(sp(Lx/2, 0.0, nL+1))),
    collect(sp(Lx/2, Sx/2, nSx-nL+1))))
    ys = unique(vcat(reverse(collect(MeshUtilModule.gradedspace(Ly/2, 0.0, nL+1))),
    collect(sp(Ly/2, Sy/2, nSy-nL+1))))
    
    fens,fes = H8layeredplatex(xs, ys, ts, nts)
    fens,fes = H8toH20(fens,fes)
    
    # This is the material  model
    MR = DeforModelRed3D
    skinmaterial = MatDeforElastOrtho(MR,
    0.0, E1s, E2s, E3s,
    nu12s, nu13s, nu23s,
    G12s, G13s, G23s,
    CTE1, CTE2, CTE3)
    corematerial = MatDeforElastOrtho(MR,
    0.0, E1c, E2c, E3c,
    nu12c, nu13c, nu23c,
    G12c, G13c, G23c,
    CTE1, CTE2, CTE3)
    
    # The material coordinate system function is defined as:
    function updatecs!(csmatout::FFltMat, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)
        rotmat3!(csmatout, angles[fe_label]/180.0*pi* [0.0; 0.0; 1.0]);
    end
    
    # The volume integrals are evaluated using this rule
    gr = GaussRule(3, 3)
    
    # We will create two regions, one for the skin,
    # and one for the core.
    rls = selectelem(fens, fes, label = 1)
    botskinregion = FDataDict("femm"=>FEMMDeforLinear(MR, IntegData(subset(fes, rls), gr), CSys(3, 3, updatecs!), skinmaterial))
    rls = selectelem(fens, fes, label = 3)
    topskinregion = FDataDict("femm"=>FEMMDeforLinear(MR, IntegData(subset(fes, rls), gr), CSys(3, 3, updatecs!), skinmaterial))
    rlc = selectelem(fens, fes, label = 2)
    coreregion = FDataDict("femm"=>FEMMDeforLinear(MR, IntegData(subset(fes, rlc), gr), CSys(3, 3, updatecs!), corematerial))
    
    # File =  "Meyer_Piening_sandwich-r1.vtk"
    # vtkexportmesh(File, skinregion["femm"].integdata.fes.conn, fens.xyz, FinEtools.MeshExportModule.H8)
    # # @async run(`"paraview.exe" $File`)
    # File =  "Meyer_Piening_sandwich-r2.vtk"
    # vtkexportmesh(File, coreregion["femm"].integdata.fes.conn, fens.xyz, FinEtools.MeshExportModule.H8)
    # @async run(`"paraview.exe" $File`)
    
    # The essential boundary conditions are applied on the symmetry planes.
    # First the plane X=0;...
    lx0 = selectnode(fens, box=[0.0 0.0 -Inf Inf -Inf Inf], inflate=tolerance)
    ex0 = FDataDict( "displacement"=>  0.0, "component"=> 1, "node_list"=>lx0 )
    # ... and then the plane Y=0.
    ly0 = selectnode(fens, box=[-Inf Inf 0.0 0.0 -Inf Inf], inflate=tolerance)
    ey0 = FDataDict( "displacement"=>  0.0, "component"=> 2, "node_list"=>ly0 )
    # The transverse displacement is fixed around the circumference.
    lz0 = vcat(selectnode(fens, box=[Sx/2 Sx/2 -Inf Inf -Inf Inf], inflate=tolerance),
    selectnode(fens, box=[-Inf Inf Sy/2 Sy/2 -Inf Inf], inflate=tolerance))
    ez0 = FDataDict( "displacement"=>  0.0, "component"=> 3, "node_list"=>lz0 )
    
    # The traction boundary condition is applied  along rectangle in the middle of the plate.
    bfes = meshboundary(fes)
    # From  the entire boundary we select those quadrilaterals that lie on the plane
    # Z = thickness
    tl = selectelem(fens, bfes, box = [0.0 Lx/2 0 Ly/2 TH TH], inflate=tolerance)
    Trac = FDataDict("traction_vector"=>vec([0.0; 0.0; -q0]), "femm"=>FEMMBase(IntegData(subset(bfes, tl), GaussRule(2, 3))))
    
    modeldata = FDataDict("fens"=>fens,
    "regions"=>[botskinregion, coreregion, topskinregion],
    "essential_bcs"=>[ex0, ey0, ez0],
    "traction_bcs"=> [Trac]
    )
    modeldata = AlgoDeforLinearModule.linearstatics(modeldata)
    
    modeldata["postprocessing"] = FDataDict("file"=>filebase * "-u")
    modeldata = AlgoDeforLinearModule.exportdeformation(modeldata)
    
    u = modeldata["u"]
    geom = modeldata["geom"]
    
    # The results of the displacement and stresses will be reported at
    # nodes located at the appropriate points.
    nbottomcenter = selectnode(fens, box=[0.0 0.0 0.0 0.0 0.0 0.0], inflate=tolerance)
    ntopcenter = selectnode(fens, box=[0.0 0.0 0.0 0.0 TH TH], inflate=tolerance)
    ncenterline = selectnode(fens, box=[0.0 0.0 0.0 0.0 0.0 TH], inflate=tolerance)
    nintertop = selectnode(fens, box=[-Inf Inf 0.0 0.0 sum(ts[1:2]) sum(ts[1:2])], inflate=tolerance)
    ninterbot = selectnode(fens, box=[-Inf Inf 0.0 0.0 sum(ts[1:1]) sum(ts[1:1])], inflate=tolerance)
    
    zclo = sortperm(vec(geom.values[ncenterline, 3]))
    ncenterline = ncenterline[zclo]
    centerz = geom.values[ncenterline, 3]
    zclo = nothing
    
    xclotop = sortperm(vec(geom.values[nintertop, 1]))
    nintertop = nintertop[xclotop]
    topx = geom.values[nintertop, 1]
    xclobot = sortperm(vec(geom.values[ninterbot, 1]))
    ninterbot = ninterbot[xclobot]
    botx = geom.values[ninterbot, 1]
    xclotop = xclobot = nothing
    
    conninbotskin = intersect(connectednodes(botskinregion["femm"].integdata.fes), ncenterline)
    connincore = intersect(connectednodes(coreregion["femm"].integdata.fes), ncenterline)
    connintopskin = intersect(connectednodes(topskinregion["femm"].integdata.fes), ncenterline)
    inbotskin = [n in conninbotskin for n in ncenterline]
    incore = [n in connincore for n in ncenterline]
    intopskin = [n in connintopskin for n in ncenterline]
    
    println("")
    println("Top Center deflection: $(u.values[ntopcenter, 3]/phun("mm")) [mm]")
    println("Bottom Center deflection: $(u.values[nbottomcenter, 3]/phun("mm")) [mm]")
    
    # # extrap = :extrapmean
    # extrap = :extraptrend
    # nodevalmeth = :averaging
    extrap = :default
    nodevalmeth = :invdistance
    
    # Normal stress in the X direction
    modeldata["postprocessing"] = FDataDict("file"=>filebase * "-sx",
    "quantity"=>:Cauchy, "component"=>1, "outputcsys"=>CSys(3),
    "nodevalmethod"=>nodevalmeth, "reportat"=>extrap)
    modeldata = AlgoDeforLinearModule.exportstress(modeldata)
    s = modeldata["postprocessing"]["exported"][1]["field"]
    sxbot = s.values[ncenterline, 1]
    s = modeldata["postprocessing"]["exported"][2]["field"]
    sxcore = s.values[ncenterline, 1]
    s = modeldata["postprocessing"]["exported"][3]["field"]
    sxtop = s.values[ncenterline, 1]
    
    # The graph data needs to be collected by going through each layer separately.
    # Some quantities may be discontinuous between layers.
    zs = vcat(  [z for (j,z) in enumerate(centerz) if inbotskin[j]],
    [z for (j,z) in enumerate(centerz) if incore[j]],
    [z for (j,z) in enumerate(centerz) if intopskin[j]]
    )
    sxs = vcat( [sxbot[j] for (j,z) in enumerate(centerz) if inbotskin[j]],
    [sxcore[j] for (j,z) in enumerate(centerz) if incore[j]],
    [sxtop[j] for (j,z) in enumerate(centerz) if intopskin[j]]
    )
    
    File = filebase * "-sx-$(extrap).CSV"
    savecsv(File, zs=vec(zs)/phun("mm"), sx=vec(sxs)/phun("MPa"))
    
    # @async run(`"paraview.exe" $File`)
    
    # Inter laminar stress between the skin and the core
    modeldata["postprocessing"] = FDataDict("file"=>filebase * "-sxz",
    "quantity"=>:Cauchy, "component"=>5, "outputcsys"=>CSys(3),
    "nodevalmethod"=>nodevalmeth, "reportat"=>extrap)
    modeldata = AlgoDeforLinearModule.exportstress(modeldata)
    s = modeldata["postprocessing"]["exported"][1]["field"]
    sxzskinbot = s.values[ninterbot, 1]
    s = modeldata["postprocessing"]["exported"][2]["field"]
    sxzcoretop = s.values[nintertop, 1]
    sxzcorebot = s.values[ninterbot, 1]
    s = modeldata["postprocessing"]["exported"][3]["field"]
    sxzskintop = s.values[nintertop, 1]
    
    File = filebase * "-sxz-$(extrap).CSV"
    savecsv(File, xstop=vec(topx)/phun("mm"), sxzskintop=vec(sxzskintop)/phun("MPa"), sxzcoretop=vec(sxzcoretop)/phun("MPa"), xsbot=vec(botx)/phun("mm"), sxzskinbot=vec(sxzskinbot)/phun("MPa"), sxzcorebot=vec(sxzcorebot)/phun("MPa"))
    
    @async run(`"paraview.exe" $File`)
    
    println("Done")
    true
    
end # Meyer_Piening_sandwich_H20


function Meyer_Piening_sandwich_H8()
    println("""
    Meyer-Piening sandwich plate, plain-vanilla H8
    """)
    
    # Reference results from:
    # [1] Application of the Elasticity Solution
    # to Linear Sandwich Beam, Plate
    # and Shell Analyses
    # H.-R. MEYER -PIENING
    # Journal of SANDWICH STRUCTURES AND MATERIALS , Vol. 6—July 2004
    
    # Assessment of the refined sinus plate finite element:
    # Free edge effect and Meyer-Piening sandwich test
    # P. Vidal, O. Polit, M. D'Ottavio, E. Valot
    # http://dx.doi.org/10.1016/j.finel.2014.08.004
    #
    # The second study deals with a benchmark problem proposed
    # by Meyer-Piening [14]. It involves a simply -supported rectangular
    # sandwich plate submitted to a localized pressure applied on an
    # area of 5x20 mm. The geometry of the sandwich structure is
    # given in Fig.13. Due to the symmetry, only one quarter of the plate
    # is meshed. The faces have different thicknesses: h = 0.5 mm
    # (bottom face), h = 0.1 mm (top face). The thickness of the core
    # is h =  11.4 mm. The material properties are given in Table 3.
    # Note that this benchmark involves strong heterogeneities (very
    # different geometric and constitutive properties between core and
    # face) and local stress gradient due to the localized pressure load.
    #
    # [14] H.-R. Meyer-Piening, Experiences with exact linear sandwich beam and plate
    # analyses regarding bending, instability and frequency investigations, in:
    # Proceedings of the Fifth International Conference On Sandwich Constructions,
    # September 5–7, vol. I, Zurich, Switzerland, 2000, pp. 37–48.
    
    
    
    t0 = time()
    # Orthotropic material for the SKIN
    E1s = 70000.0*phun("MPa")
    E2s = 71000.0*phun("MPa")
    E3s = 69000.0*phun("MPa")
    nu12s = nu13s = nu23s = 0.3
    G12s = G13s = G23s = 26000.0*phun("MPa")
    CTE1 =  CTE2 =  CTE3 = 0.0
    # Orthotropic material for the CORE
    E1c = 3.0*phun("MPa")
    E2c = 3.0*phun("MPa")
    E3c = 2.8*phun("MPa")
    nu12c = nu13c = nu23c = 0.25
    G12c = G13c = G23c = 1.0*phun("MPa")
    CTE1 =  CTE2 =  CTE3 = 0.0
    
    Lx = 5.0*phun("mm") # length  of loaded rectangle
    Ly = 20.0*phun("mm") # length  of loaded rectangle
    Sx = 100.0*phun("mm") # span of the plate
    Sy = 200.0*phun("mm") # span of the plate
    
    # Here we define the layout and the thicknesses of the layers.
    angles = vec([0.0 0.0 0.0]);
    ts = vec([0.5  11.4  0.1])*phun("mm"); # layer thicknesses
    TH = sum(ts); # total thickness of the plate
    
    tolerance = 0.0001*TH
    
    # The line load is in the negative Z direction.
    q0 = 1*phun("MPa"); #    line load
    
    # Reference deflection under the load is
    wtopref = -3.789*phun("mm"); # From [1]
    wbottomref = -2.16*phun("mm"); # Not given in [1]; guessed from the figure
    
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
    
    fens,fes = H8layeredplatex(xs, ys, ts, nts)
    
    
    # This is the material  model
    MR = DeforModelRed3D
    skinmaterial = MatDeforElastOrtho(MR,
    0.0, E1s, E2s, E3s,
    nu12s, nu13s, nu23s,
    G12s, G13s, G23s,
    CTE1, CTE2, CTE3)
    corematerial = MatDeforElastOrtho(MR,
    0.0, E1c, E2c, E3c,
    nu12c, nu13c, nu23c,
    G12c, G13c, G23c,
    CTE1, CTE2, CTE3)
    
    # The material coordinate system function is defined as:
    function updatecs!(csmatout::FFltMat, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)
        rotmat3!(csmatout, angles[fe_label]/180.0*pi* [0.0; 0.0; 1.0]);
    end
    
    # The vvolume integrals are evaluated using this rule
    gr = GaussRule(3, 2)
    
    # We will create two regions, one for the skin,
    # and one for the core.
    rls = selectelem(fens, fes, label = 1)
    botskinregion = FDataDict("femm"=>FEMMDeforLinear(MR, IntegData(subset(fes, rls), gr), CSys(3, 3, updatecs!), skinmaterial))
    rls = selectelem(fens, fes, label = 3)
    topskinregion = FDataDict("femm"=>FEMMDeforLinear(MR, IntegData(subset(fes, rls), gr), CSys(3, 3, updatecs!), skinmaterial))
    rlc = selectelem(fens, fes, label = 2)
    coreregion = FDataDict("femm"=>FEMMDeforLinear(MR, IntegData(subset(fes, rlc), gr), CSys(3, 3, updatecs!), corematerial))
    
    # File =  "Meyer_Piening_sandwich-r1.vtk"
    # vtkexportmesh(File, skinregion["femm"].integdata.fes.conn, fens.xyz, FinEtools.MeshExportModule.H8)
    # # @async run(`"paraview.exe" $File`)
    # File =  "Meyer_Piening_sandwich-r2.vtk"
    # vtkexportmesh(File, coreregion["femm"].integdata.fes.conn, fens.xyz, FinEtools.MeshExportModule.H8)
    # @async run(`"paraview.exe" $File`)
    
    # The essential boundary conditions are applied on the symmetry planes.
    # First the plane X=0;...
    lx0 = selectnode(fens, box=[0.0 0.0 -Inf Inf -Inf Inf], inflate=tolerance)
    ex0 = FDataDict( "displacement"=>  0.0, "component"=> 1, "node_list"=>lx0 )
    # ... and then the plane Y=0.
    ly0 = selectnode(fens, box=[-Inf Inf 0.0 0.0 -Inf Inf], inflate=tolerance)
    ey0 = FDataDict( "displacement"=>  0.0, "component"=> 2, "node_list"=>ly0 )
    # The transverse displacement is fixed around the circumference.
    lz0 = vcat(selectnode(fens, box=[Sx/2 Sx/2 -Inf Inf -Inf Inf], inflate=tolerance),
    selectnode(fens, box=[-Inf Inf Sy/2 Sy/2 -Inf Inf], inflate=tolerance))
    ez0 = FDataDict( "displacement"=>  0.0, "component"=> 3, "node_list"=>lz0 )
    
    # The traction boundary condition is applied  along rectangle in the middle of the plate.
    bfes = meshboundary(fes)
    # From  the entire boundary we select those quadrilaterals that lie on the plane
    # Z = thickness
    tl = selectelem(fens, bfes, box = [0.0 Lx/2 0 Ly/2 TH TH], inflate=tolerance)
    Trac = FDataDict("traction_vector"=>vec([0.0; 0.0; -q0]),
    "femm"=>FEMMBase(IntegData(subset(bfes, tl), GaussRule(2, 2))))
    
    modeldata = FDataDict("fens"=>fens,
    "regions"=>[botskinregion, coreregion, topskinregion],
    "essential_bcs"=>[ex0, ey0, ez0],
    "traction_bcs"=> [Trac]
    )
    modeldata = AlgoDeforLinearModule.linearstatics(modeldata)
    
    modeldata["postprocessing"] = FDataDict("file"=>"Meyer_Piening_sandwich")
    modeldata = AlgoDeforLinearModule.exportdeformation(modeldata)
    
    u = modeldata["u"]
    geom = modeldata["geom"]
    
    # The results of the displacement and stresses will be reported at
    # nodes located at the appropriate points.
    nbottomcenter = selectnode(fens, box=[0.0 0.0 0.0 0.0 0.0 0.0], inflate=tolerance)
    ntopcenter = selectnode(fens, box=[0.0 0.0 0.0 0.0 TH TH], inflate=tolerance)
    ncenterline = selectnode(fens, box=[0.0 0.0 0.0 0.0 0.0 TH], inflate=tolerance)
    nintertop = selectnode(fens, box=[-Inf Inf 0.0 0.0 sum(ts[1:2]) sum(ts[1:2])], inflate=tolerance)
    ninterbot = selectnode(fens, box=[-Inf Inf 0.0 0.0 sum(ts[1:1]) sum(ts[1:1])], inflate=tolerance)
    
    zclo = sortperm(vec(geom.values[ncenterline, 3]))
    centerz = geom.values[ncenterline[zclo], 3]
    xclotop = sortperm(vec(geom.values[nintertop, 1]))
    topx = geom.values[nintertop[xclotop], 1]
    xclobot = sortperm(vec(geom.values[ninterbot, 1]))
    botx = geom.values[ninterbot[xclobot], 1]
    
    conninbotskin = intersect(connectednodes(botskinregion["femm"].integdata.fes), ncenterline)
    connincore = intersect(connectednodes(coreregion["femm"].integdata.fes), ncenterline)
    connintopskin = intersect(connectednodes(topskinregion["femm"].integdata.fes), ncenterline)
    inbotskin = [n in conninbotskin for n in ncenterline]
    incore = [n in connincore for n in ncenterline]
    intopskin = [n in connintopskin for n in ncenterline]
    
    println("")
    println("Top Center deflection: $(u.values[ntopcenter, 3]/phun("mm")) [mm]")
    println("Bottom Center deflection: $(u.values[nbottomcenter, 3]/phun("mm")) [mm]")
    
    # # extrap = :extrapmean
    # extrap = :extraptrend
    # nodevalmeth = :averaging
    extrap = :default
    nodevalmeth = :invdistance
    
    # Normal stress in the X direction
    modeldata["postprocessing"] = FDataDict("file"=>"Meyer_Piening_sandwich-sx",
    "quantity"=>:Cauchy, "component"=>1, "outputcsys"=>CSys(3),
    "nodevalmethod"=>nodevalmeth, "reportat"=>extrap)
    modeldata = AlgoDeforLinearModule.exportstress(modeldata)
    s = modeldata["postprocessing"]["exported"][1]["field"]
    sxbot = s.values[ncenterline[zclo], 1]
    s = modeldata["postprocessing"]["exported"][2]["field"]
    sxcore = s.values[ncenterline[zclo], 1]
    s = modeldata["postprocessing"]["exported"][3]["field"]
    sxtop = s.values[ncenterline[zclo], 1]
    
    # The graph data needs to be collected by going through each layer separately.
    # Some quantities may be discontinuous between layers.
    zs = vcat(  [z for (j,z) in enumerate(centerz) if inbotskin[j]],
    [z for (j,z) in enumerate(centerz) if incore[j]],
    [z for (j,z) in enumerate(centerz) if intopskin[j]]
    )
    sxs = vcat( [sxbot[j] for (j,z) in enumerate(centerz) if inbotskin[j]],
    [sxcore[j] for (j,z) in enumerate(centerz) if incore[j]],
    [sxtop[j] for (j,z) in enumerate(centerz) if intopskin[j]]
    )
     
    File = "Meyer_Piening_sandwich-sx-$(extrap).CSV"
    savecsv(File, zs=vec(zs)/phun("mm"), sx=vec(sxs)/phun("MPa"))
    
    # @async run(`"paraview.exe" $File`)
    
    # Inter laminar stress between the skin and the core
    modeldata["postprocessing"] = FDataDict("file"=>"Meyer_Piening_sandwich-sxz",
    "quantity"=>:Cauchy, "component"=>5, "outputcsys"=>CSys(3),
    "nodevalmethod"=>nodevalmeth, "reportat"=>extrap)
    modeldata = AlgoDeforLinearModule.exportstress(modeldata)
    s = modeldata["postprocessing"]["exported"][1]["field"]
    sxzskinbot = s.values[ninterbot[xclobot], 1]
    s = modeldata["postprocessing"]["exported"][2]["field"]
    sxzcoretop = s.values[nintertop[xclotop], 1]
    sxzcorebot = s.values[ninterbot[xclobot], 1]
    s = modeldata["postprocessing"]["exported"][3]["field"]
    sxzskintop = s.values[nintertop[xclotop], 1]
    
    File = "Meyer_Piening_sandwich-sxz-$(extrap).CSV"
    savecsv(File, xstop=vec(topx[xclotop])/phun("mm"), sxzskintop=vec(sxzskintop[xclotop])/phun("MPa"), sxzcoretop=vec(sxzcoretop[xclotop])/phun("MPa"), xsbot=vec(botx[xclobot])/phun("mm"), sxzskinbot=vec(sxzskinbot[xclobot])/phun("MPa"), sxzcorebot=vec(sxzcorebot[xclobot])/phun("MPa"))
    
    @async run(`"paraview.exe" $File`)
    
    println("Done")
    true
    
end # Meyer_Piening_sandwich_H8


function Meyer_Piening_sandwich_MSH8()
    println("""
    Meyer-Piening sandwich plate: mean-strain hexahedron
    """)
    
    # Reference results from:
    # [1] Application of the Elasticity Solution
    # to Linear Sandwich Beam, Plate
    # and Shell Analyses
    # H.-R. MEYER -PIENING
    # Journal of SANDWICH STRUCTURES AND MATERIALS , Vol. 6—July 2004
    
    # Assessment of the refined sinus plate finite element:
    # Free edge effect and Meyer-Piening sandwich test
    # P. Vidal, O. Polit, M. D'Ottavio, E. Valot
    # http://dx.doi.org/10.1016/j.finel.2014.08.004
    #
    # The second study deals with a benchmark problem proposed
    # by Meyer-Piening [14]. It involves a simply -supported rectangular
    # sandwich plate submitted to a localized pressure applied on an
    # area of 5x20 mm. The geometry of the sandwich structure is
    # given in Fig.13. Due to the symmetry, only one quarter of the plate
    # is meshed. The faces have different thicknesses: h = 0.5 mm
    # (bottom face), h = 0.1 mm (top face). The thickness of the core
    # is h =  11.4 mm. The material properties are given in Table 3.
    # Note that this benchmark involves strong heterogeneities (very
    # different geometric and constitutive properties between core and
    # face) and local stress gradient due to the localized pressure load.
    #
    # [14] H.-R. Meyer-Piening, Experiences with exact linear sandwich beam and plate
    # analyses regarding bending, instability and frequency investigations, in:
    # Proceedings of the Fifth International Conference On Sandwich Constructions,
    # September 5–7, vol. I, Zurich, Switzerland, 2000, pp. 37–48.
    
    
    filebase = "Meyer-Piening-sandwich-MSH8"
    
    t0 = time()
    # Orthotropic material for the SKIN
    E1s = 70000.0*phun("MPa")
    E2s = 71000.0*phun("MPa")
    E3s = 69000.0*phun("MPa")
    nu12s = nu13s = nu23s = 0.3
    G12s = G13s = G23s = 26000.0*phun("MPa")
    CTE1 =  CTE2 =  CTE3 = 0.0
    # Orthotropic material for the CORE
    E1c = 3.0*phun("MPa")
    E2c = 3.0*phun("MPa")
    E3c = 2.8*phun("MPa")
    nu12c = nu13c = nu23c = 0.25
    G12c = G13c = G23c = 1.0*phun("MPa")
    CTE1 =  CTE2 =  CTE3 = 0.0
    
    Lx = 5.0*phun("mm") # length  of loaded rectangle
    Ly = 20.0*phun("mm") # length  of loaded rectangle
    Sx = 100.0*phun("mm") # span of the plate
    Sy = 200.0*phun("mm") # span of the plate
    
    # Here we define the layout and the thicknesses of the layers.
    angles = vec([0.0 0.0 0.0]);
    ts = vec([0.5  11.4  0.1])*phun("mm"); # layer thicknesses
    TH = sum(ts); # total thickness of the plate
    
    tolerance = 0.0001*TH
    
    # The line load is in the negative Z direction.
    q0 = 1*phun("MPa"); #    line load
    
    # Reference deflection under the load is
    wtopref = -3.789*phun("mm"); # From [1]
    wbottomref = -2.16*phun("mm"); # Not given in [1]; guessed from the figure
    
    # Select how find the mesh should be
    Refinement = 7
    nL = Refinement * 1;
    nSx = nL + Refinement * 4;
    nSy = 2 * nSx;
    
    # Each layer is modeled with a single element.
    nts= Refinement * [1, 2, 1];# number of elements per layer
    strength = 1.5
    sp = (a, b, n) -> MeshUtilModule.gradedspace(a, b, n, strength)
    sp = (a, b, n) -> linearspace(a, b, n)
    xs = unique(vcat(reverse(collect(sp(Lx/2, 0.0, nL+1))),
    collect(sp(Lx/2, Sx/2, nSx-nL+1))))
    ys = unique(vcat(reverse(collect(MeshUtilModule.gradedspace(Ly/2, 0.0, nL+1))),
    collect(sp(Ly/2, Sy/2, nSy-nL+1))))
    
    fens,fes = H8layeredplatex(xs, ys, ts, nts)
    
    
    # This is the material  model
    MR = DeforModelRed3D
    skinmaterial = MatDeforElastOrtho(MR,
    0.0, E1s, E2s, E3s,
    nu12s, nu13s, nu23s,
    G12s, G13s, G23s,
    CTE1, CTE2, CTE3)
    corematerial = MatDeforElastOrtho(MR,
    0.0, E1c, E2c, E3c,
    nu12c, nu13c, nu23c,
    G12c, G13c, G23c,
    CTE1, CTE2, CTE3)
    
    # The material coordinate system function is defined as:
    function updatecs!(csmatout::FFltMat, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)
        rotmat3!(csmatout, angles[fe_label]/180.0*pi* [0.0; 0.0; 1.0]);
    end
    
    # The vvolume integrals are evaluated using this rule
    gr = GaussRule(3, 2)
    
    # We will create two regions, one for the skin,
    # and one for the core.
    rls = selectelem(fens, fes, label = 1)
    botskinregion = FDataDict("femm"=>FEMMDeforLinearMSH8(MR, IntegData(subset(fes, rls), gr), CSys(3, 3, updatecs!), skinmaterial))
    rls = selectelem(fens, fes, label = 3)
    topskinregion = FDataDict("femm"=>FEMMDeforLinearMSH8(MR, IntegData(subset(fes, rls), gr), CSys(3, 3, updatecs!), skinmaterial))
    rlc = selectelem(fens, fes, label = 2)
    coreregion = FDataDict("femm"=>FEMMDeforLinearMSH8(MR, IntegData(subset(fes, rlc), gr), CSys(3, 3, updatecs!), corematerial))
    
    # File =  "Meyer_Piening_sandwich-r1.vtk"
    # vtkexportmesh(File, skinregion["femm"].integdata.fes.conn, fens.xyz, FinEtools.MeshExportModule.H8)
    # # @async run(`"paraview.exe" $File`)
    # File =  "Meyer_Piening_sandwich-r2.vtk"
    # vtkexportmesh(File, coreregion["femm"].integdata.fes.conn, fens.xyz, FinEtools.MeshExportModule.H8)
    # @async run(`"paraview.exe" $File`)
    
    # The essential boundary conditions are applied on the symmetry planes.
    # First the plane X=0;...
    lx0 = selectnode(fens, box=[0.0 0.0 -Inf Inf -Inf Inf], inflate=tolerance)
    ex0 = FDataDict( "displacement"=>  0.0, "component"=> 1, "node_list"=>lx0 )
    # ... and then the plane Y=0.
    ly0 = selectnode(fens, box=[-Inf Inf 0.0 0.0 -Inf Inf], inflate=tolerance)
    ey0 = FDataDict( "displacement"=>  0.0, "component"=> 2, "node_list"=>ly0 )
    # The transverse displacement is fixed around the circumference.
    lz0 = vcat(selectnode(fens, box=[Sx/2 Sx/2 -Inf Inf -Inf Inf], inflate=tolerance),
    selectnode(fens, box=[-Inf Inf Sy/2 Sy/2 -Inf Inf], inflate=tolerance))
    ez0 = FDataDict( "displacement"=>  0.0, "component"=> 3, "node_list"=>lz0 )
    
    # The traction boundary condition is applied  along rectangle in the middle of the plate.
    bfes = meshboundary(fes)
    # From  the entire boundary we select those quadrilaterals that lie on the plane
    # Z = thickness
    tl = selectelem(fens, bfes, box = [0.0 Lx/2 0 Ly/2 TH TH], inflate=tolerance)
    Trac = FDataDict("traction_vector"=>vec([0.0; 0.0; -q0]),
    "femm"=>FEMMBase(IntegData(subset(bfes, tl), GaussRule(2, 2))))
    
    modeldata = FDataDict("fens"=>fens,
    "regions"=>[botskinregion, coreregion, topskinregion],
    "essential_bcs"=>[ex0, ey0, ez0],
    "traction_bcs"=> [Trac]
    )
    modeldata = AlgoDeforLinearModule.linearstatics(modeldata)
    
    modeldata["postprocessing"] = FDataDict("file"=>filebase * "-u")
    modeldata = AlgoDeforLinearModule.exportdeformation(modeldata)
    
    u = modeldata["u"]
    geom = modeldata["geom"]
    
    # The results of the displacement and stresses will be reported at
    # nodes located at the appropriate points.
    nbottomcenter = selectnode(fens, box=[0.0 0.0 0.0 0.0 0.0 0.0], inflate=tolerance)
    ntopcenter = selectnode(fens, box=[0.0 0.0 0.0 0.0 TH TH], inflate=tolerance)
    ncenterline = selectnode(fens, box=[0.0 0.0 0.0 0.0 0.0 TH], inflate=tolerance)
    nintertop = selectnode(fens, box=[-Inf Inf 0.0 0.0 sum(ts[1:2]) sum(ts[1:2])], inflate=tolerance)
    ninterbot = selectnode(fens, box=[-Inf Inf 0.0 0.0 sum(ts[1:1]) sum(ts[1:1])], inflate=tolerance)
    
    zclo = sortperm(vec(geom.values[ncenterline, 3]))
    ncenterline = ncenterline[zclo]
    centerz = geom.values[ncenterline, 3]
    zclo = nothing
    
    xclotop = sortperm(vec(geom.values[nintertop, 1]))
    nintertop = nintertop[xclotop]
    topx = geom.values[nintertop, 1]
    xclobot = sortperm(vec(geom.values[ninterbot, 1]))
    ninterbot = ninterbot[xclobot]
    botx = geom.values[ninterbot, 1]
    xclotop = xclobot = nothing
    
    conninbotskin = intersect(connectednodes(botskinregion["femm"].integdata.fes), ncenterline)
    connincore = intersect(connectednodes(coreregion["femm"].integdata.fes), ncenterline)
    connintopskin = intersect(connectednodes(topskinregion["femm"].integdata.fes), ncenterline)
    inbotskin = [n in conninbotskin for n in ncenterline]
    incore = [n in connincore for n in ncenterline]
    intopskin = [n in connintopskin for n in ncenterline]
    
    println("")
    println("Top Center deflection: $(u.values[ntopcenter, 3]/phun("mm")) [mm]")
    println("Bottom Center deflection: $(u.values[nbottomcenter, 3]/phun("mm")) [mm]")
    
    # # extrap = :extrapmean
    extrap = :extraptrend
    nodevalmeth = :averaging
    # extrap = :default
    # nodevalmeth = :invdistance
    
    # Normal stress in the X direction
    modeldata["postprocessing"] = FDataDict("file"=>filebase * "-sx",
    "quantity"=>:Cauchy, "component"=>1, "outputcsys"=>CSys(3),
    "nodevalmethod"=>nodevalmeth, "reportat"=>extrap)
    modeldata = AlgoDeforLinearModule.exportstress(modeldata)
    s = modeldata["postprocessing"]["exported"][1]["field"]
    sxbot = s.values[ncenterline, 1]
    s = modeldata["postprocessing"]["exported"][2]["field"]
    sxcore = s.values[ncenterline, 1]
    s = modeldata["postprocessing"]["exported"][3]["field"]
    sxtop = s.values[ncenterline, 1]
    
    # The graph data needs to be collected by going through each layer separately.
    # Some quantities may be discontinuous between layers.
    zs = vcat(  [z for (j,z) in enumerate(centerz) if inbotskin[j]],
    [z for (j,z) in enumerate(centerz) if incore[j]],
    [z for (j,z) in enumerate(centerz) if intopskin[j]]
    )
    sxs = vcat( [sxbot[j] for (j,z) in enumerate(centerz) if inbotskin[j]],
    [sxcore[j] for (j,z) in enumerate(centerz) if incore[j]],
    [sxtop[j] for (j,z) in enumerate(centerz) if intopskin[j]]
    )
    
    File = filebase * "-sx-$(extrap).CSV"
    savecsv(File, zs=vec(zs)/phun("mm"), sx=vec(sxs)/phun("MPa"))
    
    # @async run(`"paraview.exe" $File`)
    
    # Inter laminar stress between the skin and the core
    modeldata["postprocessing"] = FDataDict("file"=>filebase * "-sxz",
    "quantity"=>:Cauchy, "component"=>5, "outputcsys"=>CSys(3),
    "nodevalmethod"=>nodevalmeth, "reportat"=>extrap)
    modeldata = AlgoDeforLinearModule.exportstress(modeldata)
    s = modeldata["postprocessing"]["exported"][1]["field"]
    sxzskinbot = s.values[ninterbot, 1]
    s = modeldata["postprocessing"]["exported"][2]["field"]
    sxzcoretop = s.values[nintertop, 1]
    sxzcorebot = s.values[ninterbot, 1]
    s = modeldata["postprocessing"]["exported"][3]["field"]
    sxzskintop = s.values[nintertop, 1]
    
    File = filebase * "-sxz-$(extrap).CSV"
    savecsv(File, xstop=vec(topx)/phun("mm"), sxzskintop=vec(sxzskintop)/phun("MPa"), sxzcoretop=vec(sxzcoretop)/phun("MPa"), xsbot=vec(botx)/phun("mm"), sxzskinbot=vec(sxzskinbot)/phun("MPa"), sxzcorebot=vec(sxzcorebot)/phun("MPa"))
    
    @async run(`"paraview.exe" $File`)
    
    println("Done")
    true
    
end # Meyer_Piening_sandwich_MSH8


function Meyer_Piening_sandwich_MST10()
    println("""
    Meyer-Piening sandwich plate, mean-strain MST10
    """)
    
    # Reference results from:
    # [1] Application of the Elasticity Solution
    # to Linear Sandwich Beam, Plate
    # and Shell Analyses
    # H.-R. MEYER -PIENING
    # Journal of SANDWICH STRUCTURES AND MATERIALS , Vol. 6—July 2004
    
    # Assessment of the refined sinus plate finite element:
    # Free edge effect and Meyer-Piening sandwich test
    # P. Vidal, O. Polit, M. D'Ottavio, E. Valot
    # http://dx.doi.org/10.1016/j.finel.2014.08.004
    #
    # The second study deals with a benchmark problem proposed
    # by Meyer-Piening [14]. It involves a simply -supported rectangular
    # sandwich plate submitted to a localized pressure applied on an
    # area of 5x20 mm. The geometry of the sandwich structure is
    # given in Fig.13. Due to the symmetry, only one quarter of the plate
    # is meshed. The faces have different thicknesses: h = 0.5 mm
    # (bottom face), h = 0.1 mm (top face). The thickness of the core
    # is h =  11.4 mm. The material properties are given in Table 3.
    # Note that this benchmark involves strong heterogeneities (very
    # different geometric and constitutive properties between core and
    # face) and local stress gradient due to the localized pressure load.
    #
    # [14] H.-R. Meyer-Piening, Experiences with exact linear sandwich beam and plate
    # analyses regarding bending, instability and frequency investigations, in:
    # Proceedings of the Fifth International Conference On Sandwich Constructions,
    # September 5–7, vol. I, Zurich, Switzerland, 2000, pp. 37–48.
    
    filebase = "Meyer-Piening-sandwich-MST10"
    
    t0 = time()
    # Orthotropic material for the SKIN
    E1s = 70000.0*phun("MPa")
    E2s = 71000.0*phun("MPa")
    E3s = 69000.0*phun("MPa")
    nu12s = nu13s = nu23s = 0.3
    G12s = G13s = G23s = 26000.0*phun("MPa")
    CTE1 =  CTE2 =  CTE3 = 0.0
    # Orthotropic material for the CORE
    E1c = 3.0*phun("MPa")
    E2c = 3.0*phun("MPa")
    E3c = 2.8*phun("MPa")
    nu12c = nu13c = nu23c = 0.25
    G12c = G13c = G23c = 1.0*phun("MPa")
    CTE1 =  CTE2 =  CTE3 = 0.0
    
    Lx = 5.0*phun("mm") # length  of loaded rectangle
    Ly = 20.0*phun("mm") # length  of loaded rectangle
    Sx = 100.0*phun("mm") # span of the plate
    Sy = 200.0*phun("mm") # span of the plate
    
    # Here we define the layout and the thicknesses of the layers.
    angles = vec([0.0 0.0 0.0]);
    ts = vec([0.5  11.4  0.1])*phun("mm"); # layer thicknesses
    TH = sum(ts); # total thickness of the plate
    
    tolerance = 0.0001*TH
    
    # The line load is in the negative Z direction.
    q0 = 1*phun("MPa"); #    line load
    
    # Reference deflection under the load is
    wtopref = -3.789*phun("mm"); # From [1]
    wbottomref = -2.16*phun("mm"); # Not given in [1]; guessed from the figure
    
    # Select how find the mesh should be
    Refinement = 3
    nL = Refinement * 1;
    nSx = nL + Refinement * 4;
    nSy = 2 * nSx;
    
    # Each layer is modeled with a single element.
    nts= Refinement * [1, 2, 1];# number of elements per layer
    strength = 1.5
    sp = (a, b, n) -> MeshUtilModule.gradedspace(a, b, n, strength)
    sp = (a, b, n) -> linearspace(a, b, n)
    xs = unique(vcat(reverse(collect(sp(Lx/2, 0.0, nL+1))),
    collect(sp(Lx/2, Sx/2, nSx-nL+1))))
    ys = unique(vcat(reverse(collect(MeshUtilModule.gradedspace(Ly/2, 0.0, nL+1))),
    collect(sp(Ly/2, Sy/2, nSy-nL+1))))
    
    fens,fes = T10layeredplatex(xs, ys, ts, nts)
    
    # This is the material  model
    MR = DeforModelRed3D
    skinmaterial = MatDeforElastOrtho(MR,
    0.0, E1s, E2s, E3s,
    nu12s, nu13s, nu23s,
    G12s, G13s, G23s,
    CTE1, CTE2, CTE3)
    corematerial = MatDeforElastOrtho(MR,
    0.0, E1c, E2c, E3c,
    nu12c, nu13c, nu23c,
    G12c, G13c, G23c,
    CTE1, CTE2, CTE3)
    
    # The material coordinate system function is defined as:
    function updatecs!(csmatout::FFltMat, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)
        rotmat3!(csmatout, angles[fe_label]/180.0*pi* [0.0; 0.0; 1.0]);
    end
    
    # The volume integrals are evaluated using this rule
    gr = SimplexRule(3, 4)
    
    # We will create two regions, one for the skin,
    # and one for the core.
    rls = selectelem(fens, fes, label = 1)
    botskinregion = FDataDict("femm"=>FEMMDeforLinearMST10(MR, IntegData(subset(fes, rls), gr), CSys(3, 3, updatecs!), skinmaterial))
    rls = selectelem(fens, fes, label = 3)
    topskinregion = FDataDict("femm"=>FEMMDeforLinearMST10(MR, IntegData(subset(fes, rls), gr), CSys(3, 3, updatecs!), skinmaterial))
    rlc = selectelem(fens, fes, label = 2)
    coreregion = FDataDict("femm"=>FEMMDeforLinearMST10(MR, IntegData(subset(fes, rlc), gr), CSys(3, 3, updatecs!), corematerial))
    
    # File =  "Meyer_Piening_sandwich-r1.vtk"
    # vtkexportmesh(File, botskinregion["femm"].integdata.fes.conn, fens.xyz,
    #     FinEtools.MeshExportModule.T10)
    # # @async run(`"paraview.exe" $File`)
    # File =  "Meyer_Piening_sandwich-r2.vtk"
    # vtkexportmesh(File, coreregion["femm"].integdata.fes.conn, fens.xyz,
    #     FinEtools.MeshExportModule.T10)
    # @async run(`"paraview.exe" $File`)
    
    # The essential boundary conditions are applied on the symmetry planes.
    # First the plane X=0;...
    lx0 = selectnode(fens, box=[0.0 0.0 -Inf Inf -Inf Inf], inflate=tolerance)
    ex0 = FDataDict( "displacement"=>  0.0, "component"=> 1, "node_list"=>lx0 )
    # ... and then the plane Y=0.
    ly0 = selectnode(fens, box=[-Inf Inf 0.0 0.0 -Inf Inf], inflate=tolerance)
    ey0 = FDataDict( "displacement"=>  0.0, "component"=> 2, "node_list"=>ly0 )
    # The transverse displacement is fixed around the circumference.
    lz0 = vcat(selectnode(fens, box=[Sx/2 Sx/2 -Inf Inf -Inf Inf], inflate=tolerance),
    selectnode(fens, box=[-Inf Inf Sy/2 Sy/2 -Inf Inf], inflate=tolerance))
    ez0 = FDataDict( "displacement"=>  0.0, "component"=> 3, "node_list"=>lz0 )
    
    # The traction boundary condition is applied  along rectangle in the middle of the plate.
    bfes = meshboundary(fes)
    # From  the entire boundary we select those quadrilaterals that lie on the plane
    # Z = thickness
    tl = selectelem(fens, bfes, box = [0.0 Lx/2 0 Ly/2 TH TH], inflate=tolerance)
    Trac = FDataDict("traction_vector"=>vec([0.0; 0.0; -q0]),
    "femm"=>FEMMBase(IntegData(subset(bfes, tl), SimplexRule(2, 3))))
    
    modeldata = FDataDict("fens"=>fens,
    "regions"=>[botskinregion, coreregion, topskinregion],
    "essential_bcs"=>[ex0, ey0, ez0],
    "traction_bcs"=> [Trac]
    )
    modeldata = AlgoDeforLinearModule.linearstatics(modeldata)
    
    modeldata["postprocessing"] = FDataDict("file"=>filebase * "-u")
    modeldata = AlgoDeforLinearModule.exportdeformation(modeldata)
    
    u = modeldata["u"]
    geom = modeldata["geom"]
    
    # The results of the displacement and stresses will be reported at
    # nodes located at the appropriate points.
    nbottomcenter = selectnode(fens, box=[0.0 0.0 0.0 0.0 0.0 0.0], inflate=tolerance)
    ntopcenter = selectnode(fens, box=[0.0 0.0 0.0 0.0 TH TH], inflate=tolerance)
    ncenterline = selectnode(fens, box=[0.0 0.0 0.0 0.0 0.0 TH], inflate=tolerance)
    nintertop = selectnode(fens, box=[-Inf Inf 0.0 0.0 sum(ts[1:2]) sum(ts[1:2])], inflate=tolerance)
    ninterbot = selectnode(fens, box=[-Inf Inf 0.0 0.0 sum(ts[1:1]) sum(ts[1:1])], inflate=tolerance)
    
    zclo = sortperm(vec(geom.values[ncenterline, 3]))
    ncenterline = ncenterline[zclo]
    centerz = geom.values[ncenterline, 3]
    zclo = nothing
    
    xclotop = sortperm(vec(geom.values[nintertop, 1]))
    nintertop = nintertop[xclotop]
    topx = geom.values[nintertop, 1]
    xclobot = sortperm(vec(geom.values[ninterbot, 1]))
    ninterbot = ninterbot[xclobot]
    botx = geom.values[ninterbot, 1]
    xclotop = xclobot = nothing
    
    conninbotskin = intersect(connectednodes(botskinregion["femm"].integdata.fes), ncenterline)
    connincore = intersect(connectednodes(coreregion["femm"].integdata.fes), ncenterline)
    connintopskin = intersect(connectednodes(topskinregion["femm"].integdata.fes), ncenterline)
    inbotskin = [n in conninbotskin for n in ncenterline]
    incore = [n in connincore for n in ncenterline]
    intopskin = [n in connintopskin for n in ncenterline]
    
    println("")
    println("Top Center deflection: $(u.values[ntopcenter, 3]/phun("mm")) [mm]")
    println("Bottom Center deflection: $(u.values[nbottomcenter, 3]/phun("mm")) [mm]")
    
    # # extrap = :extrapmean
    extrap = :extraptrend
    nodevalmeth = :averaging
    # extrap = :default
    # nodevalmeth = :invdistance
    
    # Normal stress in the X direction
    modeldata["postprocessing"] = FDataDict("file"=>filebase * "-sx",
    "quantity"=>:Cauchy, "component"=>1, "outputcsys"=>CSys(3),
    "nodevalmethod"=>nodevalmeth, "reportat"=>extrap)
    modeldata = AlgoDeforLinearModule.exportstress(modeldata)
    s = modeldata["postprocessing"]["exported"][1]["field"]
    sxbot = s.values[ncenterline, 1]
    s = modeldata["postprocessing"]["exported"][2]["field"]
    sxcore = s.values[ncenterline, 1]
    s = modeldata["postprocessing"]["exported"][3]["field"]
    sxtop = s.values[ncenterline, 1]
    
    # The graph data needs to be collected by going through each layer separately.
    # Some quantities may be discontinuous between layers.
    zs = vcat(  [z for (j,z) in enumerate(centerz) if inbotskin[j]],
    [z for (j,z) in enumerate(centerz) if incore[j]],
    [z for (j,z) in enumerate(centerz) if intopskin[j]]
    )
    sxs = vcat( [sxbot[j] for (j,z) in enumerate(centerz) if inbotskin[j]],
    [sxcore[j] for (j,z) in enumerate(centerz) if incore[j]],
    [sxtop[j] for (j,z) in enumerate(centerz) if intopskin[j]]
    )
    
    File = filebase * "-sx-$(extrap).CSV"
    savecsv(File, zs=vec(zs)/phun("mm"), sx=vec(sxs)/phun("MPa"))
    
    # @async run(`"paraview.exe" $File`)
    
    # Inter laminar stress between the skin and the core
    modeldata["postprocessing"] = FDataDict("file"=>filebase * "-sxz",
    "quantity"=>:Cauchy, "component"=>5, "outputcsys"=>CSys(3),
    "nodevalmethod"=>nodevalmeth, "reportat"=>extrap)
    modeldata = AlgoDeforLinearModule.exportstress(modeldata)
    s = modeldata["postprocessing"]["exported"][1]["field"]
    sxzskinbot = s.values[ninterbot, 1]
    s = modeldata["postprocessing"]["exported"][2]["field"]
    sxzcoretop = s.values[nintertop, 1]
    sxzcorebot = s.values[ninterbot, 1]
    s = modeldata["postprocessing"]["exported"][3]["field"]
    sxzskintop = s.values[nintertop, 1]
    
    File = filebase * "-sxz-$(extrap).CSV"
    savecsv(File, xstop=vec(topx)/phun("mm"),
    sxzskintop=vec(sxzskintop)/phun("MPa"),
    sxzcoretop=vec(sxzcoretop)/phun("MPa"),
    xsbot=vec(botx)/phun("mm"),
    sxzskinbot=vec(sxzskinbot)/phun("MPa"),
    sxzcorebot=vec(sxzcorebot)/phun("MPa"))
    
    @async run(`"paraview.exe" $File`)
    
    println("Done")
    true
    
end # Meyer_Piening_sandwich_MST10


function Meyer_Piening_sandwich_MST10_timing()
    println("""
    Meyer-Piening sandwich plate, mean-strain MST10
    """)
    
    # Reference results from:
    # [1] Application of the Elasticity Solution
    # to Linear Sandwich Beam, Plate
    # and Shell Analyses
    # H.-R. MEYER -PIENING
    # Journal of SANDWICH STRUCTURES AND MATERIALS , Vol. 6—July 2004
    
    # Assessment of the refined sinus plate finite element:
    # Free edge effect and Meyer-Piening sandwich test
    # P. Vidal, O. Polit, M. D'Ottavio, E. Valot
    # http://dx.doi.org/10.1016/j.finel.2014.08.004
    #
    # The second study deals with a benchmark problem proposed
    # by Meyer-Piening [14]. It involves a simply -supported rectangular
    # sandwich plate submitted to a localized pressure applied on an
    # area of 5x20 mm. The geometry of the sandwich structure is
    # given in Fig.13. Due to the symmetry, only one quarter of the plate
    # is meshed. The faces have different thicknesses: h = 0.5 mm
    # (bottom face), h = 0.1 mm (top face). The thickness of the core
    # is h =  11.4 mm. The material properties are given in Table 3.
    # Note that this benchmark involves strong heterogeneities (very
    # different geometric and constitutive properties between core and
    # face) and local stress gradient due to the localized pressure load.
    #
    # [14] H.-R. Meyer-Piening, Experiences with exact linear sandwich beam and plate
    # analyses regarding bending, instability and frequency investigations, in:
    # Proceedings of the Fifth International Conference On Sandwich Constructions,
    # September 5–7, vol. I, Zurich, Switzerland, 2000, pp. 37–48.
    
    filebase = "Meyer-Piening-sandwich-MST10"
    
    t0 = time()
    # Orthotropic material for the SKIN
    E1s = 70000.0*phun("MPa")
    E2s = 71000.0*phun("MPa")
    E3s = 69000.0*phun("MPa")
    nu12s = nu13s = nu23s = 0.3
    G12s = G13s = G23s = 26000.0*phun("MPa")
    CTE1 =  CTE2 =  CTE3 = 0.0
    # Orthotropic material for the CORE
    E1c = 3.0*phun("MPa")
    E2c = 3.0*phun("MPa")
    E3c = 2.8*phun("MPa")
    nu12c = nu13c = nu23c = 0.25
    G12c = G13c = G23c = 1.0*phun("MPa")
    CTE1 =  CTE2 =  CTE3 = 0.0
    
    Lx = 5.0*phun("mm") # length  of loaded rectangle
    Ly = 20.0*phun("mm") # length  of loaded rectangle
    Sx = 100.0*phun("mm") # span of the plate
    Sy = 200.0*phun("mm") # span of the plate
    
    # Here we define the layout and the thicknesses of the layers.
    angles = vec([0.0 0.0 0.0]);
    ts = vec([0.5  11.4  0.1])*phun("mm"); # layer thicknesses
    TH = sum(ts); # total thickness of the plate
    
    tolerance = 0.0001*TH
    
    # The line load is in the negative Z direction.
    q0 = 1*phun("MPa"); #    line load
    
    # Reference deflection under the load is
    wtopref = -3.789*phun("mm"); # From [1]
    wbottomref = -2.16*phun("mm"); # Not given in [1]; guessed from the figure
    
    # Select how find the mesh should be
    Refinement = 3
    nL = Refinement * 1;
    nSx = nL + Refinement * 4;
    nSy = 2 * nSx;
    
    # Each layer is modeled with a single element.
    nts= Refinement * [1, 2, 1];# number of elements per layer
    strength = 1.5
    sp = (a, b, n) -> MeshUtilModule.gradedspace(a, b, n, strength)
    sp = (a, b, n) -> linearspace(a, b, n)
    xs = unique(vcat(reverse(collect(sp(Lx/2, 0.0, nL+1))),
    collect(sp(Lx/2, Sx/2, nSx-nL+1))))
    ys = unique(vcat(reverse(collect(MeshUtilModule.gradedspace(Ly/2, 0.0, nL+1))),
    collect(sp(Ly/2, Sy/2, nSy-nL+1))))
    
    fens,fes = T10layeredplatex(xs, ys, ts, nts)
    
    # This is the material  model
    MR = DeforModelRed3D
    skinmaterial = MatDeforElastOrtho(MR,
    0.0, E1s, E2s, E3s,
    nu12s, nu13s, nu23s,
    G12s, G13s, G23s,
    CTE1, CTE2, CTE3)
    corematerial = MatDeforElastOrtho(MR,
    0.0, E1c, E2c, E3c,
    nu12c, nu13c, nu23c,
    G12c, G13c, G23c,
    CTE1, CTE2, CTE3)
    
    # The material coordinate system function is defined as:
    function updatecs!(csmatout::FFltMat, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)
        rotmat3!(csmatout, angles[fe_label]/180.0*pi* [0.0; 0.0; 1.0]);
    end
    
    # The volume integrals are evaluated using this rule
    gr = SimplexRule(3, 4)
    
    # We will create two regions, one for the skin,
    # and one for the core.
    rls = selectelem(fens, fes, label = 1)
    botskinregion = FDataDict("femm"=>FEMMDeforLinearMST10(MR, IntegData(subset(fes, rls), gr), CSys(3, 3, updatecs!), skinmaterial))
    rls = selectelem(fens, fes, label = 3)
    topskinregion = FDataDict("femm"=>FEMMDeforLinearMST10(MR, IntegData(subset(fes, rls), gr), CSys(3, 3, updatecs!), skinmaterial))
    rlc = selectelem(fens, fes, label = 2)
    coreregion = FDataDict("femm"=>FEMMDeforLinearMST10(MR, IntegData(subset(fes, rlc), gr), CSys(3, 3, updatecs!), corematerial))
    
    # File =  "Meyer_Piening_sandwich-r1.vtk"
    # vtkexportmesh(File, botskinregion["femm"].integdata.fes.conn, fens.xyz,
    #     FinEtools.MeshExportModule.T10)
    # # @async run(`"paraview.exe" $File`)
    # File =  "Meyer_Piening_sandwich-r2.vtk"
    # vtkexportmesh(File, coreregion["femm"].integdata.fes.conn, fens.xyz,
    #     FinEtools.MeshExportModule.T10)
    # @async run(`"paraview.exe" $File`)
    
    # The essential boundary conditions are applied on the symmetry planes.
    # First the plane X=0;...
    lx0 = selectnode(fens, box=[0.0 0.0 -Inf Inf -Inf Inf], inflate=tolerance)
    ex0 = FDataDict( "displacement"=>  0.0, "component"=> 1, "node_list"=>lx0 )
    # ... and then the plane Y=0.
    ly0 = selectnode(fens, box=[-Inf Inf 0.0 0.0 -Inf Inf], inflate=tolerance)
    ey0 = FDataDict( "displacement"=>  0.0, "component"=> 2, "node_list"=>ly0 )
    # The transverse displacement is fixed around the circumference.
    lz0 = vcat(selectnode(fens, box=[Sx/2 Sx/2 -Inf Inf -Inf Inf], inflate=tolerance),
    selectnode(fens, box=[-Inf Inf Sy/2 Sy/2 -Inf Inf], inflate=tolerance))
    ez0 = FDataDict( "displacement"=>  0.0, "component"=> 3, "node_list"=>lz0 )
    
    # The traction boundary condition is applied  along rectangle in the middle of the plate.
    bfes = meshboundary(fes)
    # From  the entire boundary we select those quadrilaterals that lie on the plane
    # Z = thickness
    tl = selectelem(fens, bfes, box = [0.0 Lx/2 0 Ly/2 TH TH], inflate=tolerance)
    Trac = FDataDict("traction_vector"=>vec([0.0; 0.0; -q0]),
    "femm"=>FEMMBase(IntegData(subset(bfes, tl), SimplexRule(2, 3))))
    
    modeldata = FDataDict("fens"=>fens,
    "regions"=>[botskinregion, coreregion, topskinregion],
    "essential_bcs"=>[ex0, ey0, ez0],
    "traction_bcs"=> [Trac]
    )
    t0 = time()
    modeldata = AlgoDeforLinearModule.linearstatics(modeldata)
    tstiffness = modeldata["timing"]["stiffness"]
    tsolution = modeldata["timing"]["solution"]
    println("count(fes) = $(count(fes))")
    println("Timing: Assembly $(tstiffness), Solution $(tsolution)")
    
    modeldata["postprocessing"] = FDataDict("file"=>filebase * "-u")
    modeldata = AlgoDeforLinearModule.exportdeformation(modeldata)
    
    u = modeldata["u"]
    geom = modeldata["geom"]
    
    # The results of the displacement and stresses will be reported at
    # nodes located at the appropriate points.
    nbottomcenter = selectnode(fens, box=[0.0 0.0 0.0 0.0 0.0 0.0], inflate=tolerance)
    ntopcenter = selectnode(fens, box=[0.0 0.0 0.0 0.0 TH TH], inflate=tolerance)
    # ncenterline = selectnode(fens, box=[0.0 0.0 0.0 0.0 0.0 TH], inflate=tolerance)
    # nintertop = selectnode(fens, box=[-Inf Inf 0.0 0.0 sum(ts[1:2]) sum(ts[1:2])], inflate=tolerance)
    # ninterbot = selectnode(fens, box=[-Inf Inf 0.0 0.0 sum(ts[1:1]) sum(ts[1:1])], inflate=tolerance)
    #
    # zclo = sortperm(vec(geom.values[ncenterline, 3]))
    # ncenterline = ncenterline[zclo]
    # centerz = geom.values[ncenterline, 3]
    # zclo = nothing
    #
    # xclotop = sortperm(vec(geom.values[nintertop, 1]))
    # nintertop = nintertop[xclotop]
    # topx = geom.values[nintertop, 1]
    # xclobot = sortperm(vec(geom.values[ninterbot, 1]))
    # ninterbot = ninterbot[xclobot]
    # botx = geom.values[ninterbot, 1]
    # xclotop = xclobot = nothing
    #
    # conninbotskin = intersect(connectednodes(botskinregion["femm"].integdata.fes), ncenterline)
    # connincore = intersect(connectednodes(coreregion["femm"].integdata.fes), ncenterline)
    # connintopskin = intersect(connectednodes(topskinregion["femm"].integdata.fes), ncenterline)
    # inbotskin = [n in conninbotskin for n in ncenterline]
    # incore = [n in connincore for n in ncenterline]
    # intopskin = [n in connintopskin for n in ncenterline]
    
    println("")
    println("Top Center deflection: $(u.values[ntopcenter, 3]/phun("mm")) [mm]")
    println("Bottom Center deflection: $(u.values[nbottomcenter, 3]/phun("mm")) [mm]")
    
    # # # extrap = :extrapmean
    # extrap = :extraptrend
    # nodevalmeth = :averaging
    # # extrap = :default
    # # nodevalmeth = :invdistance
    #
    # # Normal stress in the X direction
    # modeldata["postprocessing"] = FDataDict("file"=>filebase * "-sx",
    #     "quantity"=>:Cauchy, "component"=>1, "outputcsys"=>CSys(3),
    #      "nodevalmethod"=>nodevalmeth, "reportat"=>extrap)
    # modeldata = AlgoDeforLinearModule.exportstress(modeldata)
    # s = modeldata["postprocessing"]["exported"][1]["field"]
    # sxbot = s.values[ncenterline, 1]
    # s = modeldata["postprocessing"]["exported"][2]["field"]
    # sxcore = s.values[ncenterline, 1]
    # s = modeldata["postprocessing"]["exported"][3]["field"]
    # sxtop = s.values[ncenterline, 1]
    #
    # # The graph data needs to be collected by going through each layer separately.
    # # Some quantities may be discontinuous between layers.
    # zs = vcat(  [z for (j,z) in enumerate(centerz) if inbotskin[j]],
    #             [z for (j,z) in enumerate(centerz) if incore[j]],
    #             [z for (j,z) in enumerate(centerz) if intopskin[j]]
    #             )
    # sxs = vcat( [sxbot[j] for (j,z) in enumerate(centerz) if inbotskin[j]],
    #             [sxcore[j] for (j,z) in enumerate(centerz) if incore[j]],
    #             [sxtop[j] for (j,z) in enumerate(centerz) if intopskin[j]]
    #             )
    
    # File = filebase * "-sx-$(extrap).CSV"
    # savecsv(File, zs=vec(zs)/phun("mm"), sx=vec(sxs)/phun("MPa"))
    #
    # # @async run(`"paraview.exe" $File`)
    #
    # # Inter laminar stress between the skin and the core
    # modeldata["postprocessing"] = FDataDict("file"=>filebase * "-sxz",
    #     "quantity"=>:Cauchy, "component"=>5, "outputcsys"=>CSys(3),
    #      "nodevalmethod"=>nodevalmeth, "reportat"=>extrap)
    # modeldata = AlgoDeforLinearModule.exportstress(modeldata)
    # s = modeldata["postprocessing"]["exported"][1]["field"]
    # sxzskinbot = s.values[ninterbot, 1]
    # s = modeldata["postprocessing"]["exported"][2]["field"]
    # sxzcoretop = s.values[nintertop, 1]
    # sxzcorebot = s.values[ninterbot, 1]
    # s = modeldata["postprocessing"]["exported"][3]["field"]
    # sxzskintop = s.values[nintertop, 1]
    
    # File = filebase * "-sxz-$(extrap).CSV"
    # savecsv(File, xstop=vec(topx)/phun("mm"),
    #     sxzskintop=vec(sxzskintop)/phun("MPa"),
    #     sxzcoretop=vec(sxzcoretop)/phun("MPa"),
    #     xsbot=vec(botx)/phun("mm"),
    #     sxzskinbot=vec(sxzskinbot)/phun("MPa"),
    #     sxzcorebot=vec(sxzcorebot)/phun("MPa"))
    #
    # @async run(`"paraview.exe" $File`)
    
    println("count(fes) = $(count(fes))")
    println("Timing: $( time() - t0 )")
    true
    
end # Meyer_Piening_sandwich_MST10_timing


function Meyer_Piening_sandwich_T10_timing()
    println("""
    Meyer-Piening sandwich plate, mean-strain T10
    """)
    
    # Reference results from:
    # [1] Application of the Elasticity Solution
    # to Linear Sandwich Beam, Plate
    # and Shell Analyses
    # H.-R. MEYER -PIENING
    # Journal of SANDWICH STRUCTURES AND MATERIALS , Vol. 6—July 2004
    
    # Assessment of the refined sinus plate finite element:
    # Free edge effect and Meyer-Piening sandwich test
    # P. Vidal, O. Polit, M. D'Ottavio, E. Valot
    # http://dx.doi.org/10.1016/j.finel.2014.08.004
    #
    # The second study deals with a benchmark problem proposed
    # by Meyer-Piening [14]. It involves a simply -supported rectangular
    # sandwich plate submitted to a localized pressure applied on an
    # area of 5x20 mm. The geometry of the sandwich structure is
    # given in Fig.13. Due to the symmetry, only one quarter of the plate
    # is meshed. The faces have different thicknesses: h = 0.5 mm
    # (bottom face), h = 0.1 mm (top face). The thickness of the core
    # is h =  11.4 mm. The material properties are given in Table 3.
    # Note that this benchmark involves strong heterogeneities (very
    # different geometric and constitutive properties between core and
    # face) and local stress gradient due to the localized pressure load.
    #
    # [14] H.-R. Meyer-Piening, Experiences with exact linear sandwich beam and plate
    # analyses regarding bending, instability and frequency investigations, in:
    # Proceedings of the Fifth International Conference On Sandwich Constructions,
    # September 5–7, vol. I, Zurich, Switzerland, 2000, pp. 37–48.
    
    filebase = "Meyer-Piening-sandwich-T10"
    
    
    # Orthotropic material for the SKIN
    E1s = 70000.0*phun("MPa")
    E2s = 71000.0*phun("MPa")
    E3s = 69000.0*phun("MPa")
    nu12s = nu13s = nu23s = 0.3
    G12s = G13s = G23s = 26000.0*phun("MPa")
    CTE1 =  CTE2 =  CTE3 = 0.0
    # Orthotropic material for the CORE
    E1c = 3.0*phun("MPa")
    E2c = 3.0*phun("MPa")
    E3c = 2.8*phun("MPa")
    nu12c = nu13c = nu23c = 0.25
    G12c = G13c = G23c = 1.0*phun("MPa")
    CTE1 =  CTE2 =  CTE3 = 0.0
    
    Lx = 5.0*phun("mm") # length  of loaded rectangle
    Ly = 20.0*phun("mm") # length  of loaded rectangle
    Sx = 100.0*phun("mm") # span of the plate
    Sy = 200.0*phun("mm") # span of the plate
    
    # Here we define the layout and the thicknesses of the layers.
    angles = vec([0.0 0.0 0.0]);
    ts = vec([0.5  11.4  0.1])*phun("mm"); # layer thicknesses
    TH = sum(ts); # total thickness of the plate
    
    tolerance = 0.0001*TH
    
    # The line load is in the negative Z direction.
    q0 = 1*phun("MPa"); #    line load
    
    # Reference deflection under the load is
    wtopref = -3.789*phun("mm"); # From [1]
    wbottomref = -2.16*phun("mm"); # Not given in [1]; guessed from the figure
    
    # Select how find the mesh should be
    Refinement = 3
    nL = Refinement * 1;
    nSx = nL + Refinement * 4;
    nSy = 2 * nSx;
    
    # Each layer is modeled with a single element.
    nts= Refinement * [1, 2, 1];# number of elements per layer
    strength = 1.5
    sp = (a, b, n) -> MeshUtilModule.gradedspace(a, b, n, strength)
    sp = (a, b, n) -> linearspace(a, b, n)
    xs = unique(vcat(reverse(collect(sp(Lx/2, 0.0, nL+1))),
    collect(sp(Lx/2, Sx/2, nSx-nL+1))))
    ys = unique(vcat(reverse(collect(MeshUtilModule.gradedspace(Ly/2, 0.0, nL+1))),
    collect(sp(Ly/2, Sy/2, nSy-nL+1))))
    
    fens,fes = T10layeredplatex(xs, ys, ts, nts)
    
    # This is the material  model
    MR = DeforModelRed3D
    skinmaterial = MatDeforElastOrtho(MR,
    0.0, E1s, E2s, E3s,
    nu12s, nu13s, nu23s,
    G12s, G13s, G23s,
    CTE1, CTE2, CTE3)
    corematerial = MatDeforElastOrtho(MR,
    0.0, E1c, E2c, E3c,
    nu12c, nu13c, nu23c,
    G12c, G13c, G23c,
    CTE1, CTE2, CTE3)
    
    # The material coordinate system function is defined as:
    function updatecs!(csmatout::FFltMat, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)
        rotmat3!(csmatout, angles[fe_label]/180.0*pi* [0.0; 0.0; 1.0]);
    end
    
    # The volume integrals are evaluated using this rule
    gr = SimplexRule(3, 4)
    
    # We will create two regions, one for the skin,
    # and one for the core.
    rls = selectelem(fens, fes, label = 1)
    botskinregion = FDataDict("femm"=>FEMMDeforLinear(MR, IntegData(subset(fes, rls), gr), CSys(3, 3, updatecs!), skinmaterial))
    rls = selectelem(fens, fes, label = 3)
    topskinregion = FDataDict("femm"=>FEMMDeforLinear(MR, IntegData(subset(fes, rls), gr), CSys(3, 3, updatecs!), skinmaterial))
    rlc = selectelem(fens, fes, label = 2)
    coreregion = FDataDict("femm"=>FEMMDeforLinear(MR, IntegData(subset(fes, rlc), gr), CSys(3, 3, updatecs!), corematerial))
    
    # File =  "Meyer_Piening_sandwich-r1.vtk"
    # vtkexportmesh(File, botskinregion["femm"].integdata.fes.conn, fens.xyz,
    #     FinEtools.MeshExportModule.T10)
    # # @async run(`"paraview.exe" $File`)
    # File =  "Meyer_Piening_sandwich-r2.vtk"
    # vtkexportmesh(File, coreregion["femm"].integdata.fes.conn, fens.xyz,
    #     FinEtools.MeshExportModule.T10)
    # @async run(`"paraview.exe" $File`)
    
    # The essential boundary conditions are applied on the symmetry planes.
    # First the plane X=0;...
    lx0 = selectnode(fens, box=[0.0 0.0 -Inf Inf -Inf Inf], inflate=tolerance)
    ex0 = FDataDict( "displacement"=>  0.0, "component"=> 1, "node_list"=>lx0 )
    # ... and then the plane Y=0.
    ly0 = selectnode(fens, box=[-Inf Inf 0.0 0.0 -Inf Inf], inflate=tolerance)
    ey0 = FDataDict( "displacement"=>  0.0, "component"=> 2, "node_list"=>ly0 )
    # The transverse displacement is fixed around the circumference.
    lz0 = vcat(selectnode(fens, box=[Sx/2 Sx/2 -Inf Inf -Inf Inf], inflate=tolerance),
    selectnode(fens, box=[-Inf Inf Sy/2 Sy/2 -Inf Inf], inflate=tolerance))
    ez0 = FDataDict( "displacement"=>  0.0, "component"=> 3, "node_list"=>lz0 )
    
    # The traction boundary condition is applied  along rectangle in the middle of the plate.
    bfes = meshboundary(fes)
    # From  the entire boundary we select those quadrilaterals that lie on the plane
    # Z = thickness
    tl = selectelem(fens, bfes, box = [0.0 Lx/2 0 Ly/2 TH TH], inflate=tolerance)
    Trac = FDataDict("traction_vector"=>vec([0.0; 0.0; -q0]),
    "femm"=>FEMMBase(IntegData(subset(bfes, tl), SimplexRule(2, 3))))
    
    modeldata = FDataDict("fens"=>fens,
    "regions"=>[botskinregion, coreregion, topskinregion],
    "essential_bcs"=>[ex0, ey0, ez0],
    "traction_bcs"=> [Trac]
    )
    t0 = time()
    modeldata = AlgoDeforLinearModule.linearstatics(modeldata)
    tstiffness = modeldata["timing"]["stiffness"]
    tsolution = modeldata["timing"]["solution"]
    println("count(fes) = $(count(fes))")
    println("Timing: Assembly $(tstiffness), Solution $(tsolution)")
    
    modeldata["postprocessing"] = FDataDict("file"=>filebase * "-u")
    modeldata = AlgoDeforLinearModule.exportdeformation(modeldata)
    
    u = modeldata["u"]
    geom = modeldata["geom"]
    
    # The results of the displacement and stresses will be reported at
    # nodes located at the appropriate points.
    nbottomcenter = selectnode(fens, box=[0.0 0.0 0.0 0.0 0.0 0.0], inflate=tolerance)
    ntopcenter = selectnode(fens, box=[0.0 0.0 0.0 0.0 TH TH], inflate=tolerance)
    # ncenterline = selectnode(fens, box=[0.0 0.0 0.0 0.0 0.0 TH], inflate=tolerance)
    # nintertop = selectnode(fens, box=[-Inf Inf 0.0 0.0 sum(ts[1:2]) sum(ts[1:2])], inflate=tolerance)
    # ninterbot = selectnode(fens, box=[-Inf Inf 0.0 0.0 sum(ts[1:1]) sum(ts[1:1])], inflate=tolerance)
    #
    # zclo = sortperm(vec(geom.values[ncenterline, 3]))
    # ncenterline = ncenterline[zclo]
    # centerz = geom.values[ncenterline, 3]
    # zclo = nothing
    #
    # xclotop = sortperm(vec(geom.values[nintertop, 1]))
    # nintertop = nintertop[xclotop]
    # topx = geom.values[nintertop, 1]
    # xclobot = sortperm(vec(geom.values[ninterbot, 1]))
    # ninterbot = ninterbot[xclobot]
    # botx = geom.values[ninterbot, 1]
    # xclotop = xclobot = nothing
    #
    # conninbotskin = intersect(connectednodes(botskinregion["femm"].integdata.fes), ncenterline)
    # connincore = intersect(connectednodes(coreregion["femm"].integdata.fes), ncenterline)
    # connintopskin = intersect(connectednodes(topskinregion["femm"].integdata.fes), ncenterline)
    # inbotskin = [n in conninbotskin for n in ncenterline]
    # incore = [n in connincore for n in ncenterline]
    # intopskin = [n in connintopskin for n in ncenterline]
    
    println("")
    println("Top Center deflection: $(u.values[ntopcenter, 3]/phun("mm")) [mm]")
    println("Bottom Center deflection: $(u.values[nbottomcenter, 3]/phun("mm")) [mm]")
    
    # # # extrap = :extrapmean
    # extrap = :extraptrend
    # nodevalmeth = :averaging
    # # extrap = :default
    # # nodevalmeth = :invdistance
    #
    # # Normal stress in the X direction
    # modeldata["postprocessing"] = FDataDict("file"=>filebase * "-sx",
    #     "quantity"=>:Cauchy, "component"=>1, "outputcsys"=>CSys(3),
    #      "nodevalmethod"=>nodevalmeth, "reportat"=>extrap)
    # modeldata = AlgoDeforLinearModule.exportstress(modeldata)
    # s = modeldata["postprocessing"]["exported"][1]["field"]
    # sxbot = s.values[ncenterline, 1]
    # s = modeldata["postprocessing"]["exported"][2]["field"]
    # sxcore = s.values[ncenterline, 1]
    # s = modeldata["postprocessing"]["exported"][3]["field"]
    # sxtop = s.values[ncenterline, 1]
    #
    # # The graph data needs to be collected by going through each layer separately.
    # # Some quantities may be discontinuous between layers.
    # zs = vcat(  [z for (j,z) in enumerate(centerz) if inbotskin[j]],
    #             [z for (j,z) in enumerate(centerz) if incore[j]],
    #             [z for (j,z) in enumerate(centerz) if intopskin[j]]
    #             )
    # sxs = vcat( [sxbot[j] for (j,z) in enumerate(centerz) if inbotskin[j]],
    #             [sxcore[j] for (j,z) in enumerate(centerz) if incore[j]],
    #             [sxtop[j] for (j,z) in enumerate(centerz) if intopskin[j]]
    #             )
    
    # File = filebase * "-sx-$(extrap).CSV"
    # savecsv(File, zs=vec(zs)/phun("mm"), sx=vec(sxs)/phun("MPa"))
    #
    # # @async run(`"paraview.exe" $File`)
    #
    # # Inter laminar stress between the skin and the core
    # modeldata["postprocessing"] = FDataDict("file"=>filebase * "-sxz",
    #     "quantity"=>:Cauchy, "component"=>5, "outputcsys"=>CSys(3),
    #      "nodevalmethod"=>nodevalmeth, "reportat"=>extrap)
    # modeldata = AlgoDeforLinearModule.exportstress(modeldata)
    # s = modeldata["postprocessing"]["exported"][1]["field"]
    # sxzskinbot = s.values[ninterbot, 1]
    # s = modeldata["postprocessing"]["exported"][2]["field"]
    # sxzcoretop = s.values[nintertop, 1]
    # sxzcorebot = s.values[ninterbot, 1]
    # s = modeldata["postprocessing"]["exported"][3]["field"]
    # sxzskintop = s.values[nintertop, 1]
    
    # File = filebase * "-sxz-$(extrap).CSV"
    # savecsv(File, xstop=vec(topx)/phun("mm"),
    #     sxzskintop=vec(sxzskintop)/phun("MPa"),
    #     sxzcoretop=vec(sxzcoretop)/phun("MPa"),
    #     xsbot=vec(botx)/phun("mm"),
    #     sxzskinbot=vec(sxzskinbot)/phun("MPa"),
    #     sxzcorebot=vec(sxzcorebot)/phun("MPa"))
    #
    # @async run(`"paraview.exe" $File`)
    
    
    true
    
end # Meyer_Piening_sandwich_T10_timing

function allrun()
    println("#####################################################") 
    println("# Meyer_Piening_sandwich ")
    Meyer_Piening_sandwich()
    println("#####################################################") 
    println("# Meyer_Piening_sandwich_H20 ")
    Meyer_Piening_sandwich_H20()
    println("#####################################################") 
    println("# Meyer_Piening_sandwich_H8 ")
    Meyer_Piening_sandwich_H8()
    println("#####################################################") 
    println("# Meyer_Piening_sandwich_MSH8 ")
    Meyer_Piening_sandwich_MSH8()
    println("#####################################################") 
    println("# Meyer_Piening_sandwich_MST10 ")
    Meyer_Piening_sandwich_MST10()
    println("#####################################################") 
    println("# Meyer_Piening_sandwich_MST10_timing ")
    Meyer_Piening_sandwich_MST10_timing()
    println("#####################################################") 
    println("# Meyer_Piening_sandwich_T10_timing ")
    Meyer_Piening_sandwich_T10_timing()
    return true
end # function allrun

end # module Meyer_Piening_examples
