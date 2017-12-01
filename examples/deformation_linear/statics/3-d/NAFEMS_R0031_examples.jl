module NAFEMS_R0031_examples
using FinEtools
using FinEtools.AlgoDeforLinearModule


function NAFEMS_R0031_1()
    println("""
    Laminated Strip Under Three-Point Bending
    """)
    
    # Determine the central transverse displacement in a simply-supported seven
    # layer symmetric strip with a central line load. A 0/90/0/90/0/90/0
    # material lay-up is specified with the center ply being four times as
    # thick as the others.
    # Reference: NAFEMS Report R0031, Test No.1, 17-Dec-1998.
    
    # Because of the symmetries of the geometry and load, only the
    # first-quadrant   (in XY) quarter of the plate is modeled.
    
    # The coordinate system is centered at point E (at the difference with
    # respect to the original benchmark definition).  The  load is applied
    # along a curve passing through point C. The simple support is applied
    # along the curve passing through point B.
    
    # We realize the simple supports along the lines  A, B and the line load at
    # point C  are illegal from the point of view of convergence.  No
    # convergence can be hoped for as the stress underneath the load and above
    # the simple supports  is infinite in the limit (these locations are stress
    # singularities).   However, for relatively coarse meshes the results away
    # from the singularities are still meaningful.
    
    # The target quantities are displacement at the bottom surface at point E,
    # the tensile axial stress at the same point,  and of the transverse shear
    # stress at point D  in between the bottommost two layers (See figure 1).
    
    t0 = time()
    # Orthotropic material parameters of the material of the layers
    E1s = 100.0*phun("GPa")
    E2s = E3s = 5.0*phun("GPa")
    nu12s = 0.4
    nu13s = 0.3
    nu23s = 0.3
    G12s = 3.0*phun("GPa")
    G13s = G23s = 2.0*phun("GPa")
    CTE1 = 3.0e-6
    CTE2 = 2.0e-5
    CTE3 = 2.0e-5
    
    AB = 30.0*phun("mm") # span between simple supports
    CD = 4.0*phun("mm") # distance between the point D and the center
    OH = 10.0*phun("mm") # overhang
    W = 10.0*phun("mm") # width of the strip
    
    # Here we define the layout and the thicknesses of the layers.
    angles = vec([0 90 0 90 0 90 0]);
    ts = vec([0.1  0.1  0.1  0.4  0.1  0.1  0.1])*phun("mm"); # layer thicknesses
    TH = sum(ts); # total thickness of the plate
    
    tolerance = 0.0001*TH
    
    # The line load is in the negative Z direction.
    q0 = 10*phun("N/mm"); #    line load
    
    # Reference deflection under the load is
    wEref = -1.06*phun("mm");
    
    # The reference tensile stress at the bottom of the lowest layer is
    sigma11Eref = 684*phun("MPa");
    
    # Because we model the first-quadrant quarter of the plate using
    # coordinate axes centered  at the point E  the shear at the point D is
    # positive instead of negative as in the benchmark where the coordinate
    # system is located at the outer corner of the strip.
    sigma13Dref=4.1*phun("MPa");
    
    Refinement = 10
    # We select 8 elements spanwise and 2 elements widthwise.  The overhang
    # of the plate is given one element.
    nL = Refinement * 4; nO = Refinement * 1; nW = Refinement * 1;
    
    # Each layer is modeled with a single element.
    nts= Refinement * ones(Int, length(angles));# number of elements per layer
    
    xs = unique(vcat(collect(linspace(0,AB/2,nL+1)), [CD],
    collect(linspace(AB/2,AB/2+OH,nO+1))))
    xs = xs[sortperm(xs)]
    ys = collect(linspace(0,W/2,nW+1));
    
    fens,fes = H8layeredplatex(xs, ys, ts, nts)
    
    
    # This is the material  model
    MR = DeforModelRed3D
    material = MatDeforElastOrtho(MR,
    0.0, E1s, E2s, E3s,
    nu12s, nu13s, nu23s,
    G12s, G13s, G23s,
    CTE1, CTE2, CTE3)
    
    # The material coordinate system function is defined as:
    function updatecs!(csmatout::FFltMat, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)
        rotmat3!(csmatout, angles[fe_label]/180.0*pi* [0.0; 0.0; 1.0]);
    end
    
    # The vvolume integrals are evaluated using this rule
    gr = GaussRule(3, 2)
    
    # We will create two regions, one for the layers with 0°  orientation,
    # and one for the layers with 90° orientation.
    rl1 = vcat(selectelem(fens, fes, label = 1), selectelem(fens, fes, label = 3),
    selectelem(fens, fes, label = 5), selectelem(fens, fes, label = 7))
    rl2 = vcat(selectelem(fens, fes, label = 2), selectelem(fens, fes, label = 4),
    selectelem(fens, fes, label = 6))
    region1 = FDataDict("femm"=>FEMMDeforLinearMSH8(MR, IntegData(subset(fes, rl1), gr), CSys(3, 3, updatecs!), material))
    region2 = FDataDict("femm"=>FEMMDeforLinearMSH8(MR, IntegData(subset(fes, rl2), gr), CSys(3, 3, updatecs!), material))
    
    # File =  "NAFEMS-R0031-1-plate-r1.vtk"
    # vtkexportmesh(File, region1["femm"].integdata.fes.conn, fens.xyz, FinEtools.MeshExportModule.H8)
    # # @async run(`"paraview.exe" $File`)
    # File =  "NAFEMS-R0031-1-plate-r2.vtk"
    # vtkexportmesh(File, region2["femm"].integdata.fes.conn, fens.xyz, FinEtools.MeshExportModule.H8)
    # @async run(`"paraview.exe" $File`)
    
    # The essential boundary conditions are applied on the symmetry planes.
    # First the plane X=0;...
    lx0 = selectnode(fens, box=[0.0 0.0 -Inf Inf -Inf Inf], inflate=tolerance)
    ex0 = FDataDict( "displacement"=>  0.0, "component"=> 1, "node_list"=>lx0 )
    # ... and then the plane Y=0.
    ly0 = selectnode(fens, box=[-Inf Inf 0.0 0.0 -Inf Inf], inflate=tolerance)
    ey0 = FDataDict( "displacement"=>  0.0, "component"=> 2, "node_list"=>ly0 )
    # The transverse displacement is fixed along the line  passing through
    # point B. The nodes are fixed in the box along this line in the Z
    # direction.
    lz0 = selectnode(fens, box=[AB/2 AB/2 -Inf Inf -Inf Inf], inflate=tolerance)
    ez0 = FDataDict( "displacement"=>  0.0, "component"=> 3, "node_list"=>lz0 )
    
    # The traction boundary condition is applied  along  the edge of the
    # mesh passing through point C at the top surface of the strip.   First
    # we extract the boundary of the hexahedral mesh.
    bfes = meshboundary(fes)
    # From  the entire boundary we select those quadrilaterals that lie on the plane
    # X = 0
    xl = selectelem(fens, bfes, box = [0.0 0.0 -Inf Inf -Inf Inf], inflate=tolerance)
    # Now we extract the boundary  of these selected quadrilaterals
    bbfes = meshboundary(subset(bfes, xl))
    # …  And from these  we extract the ones at the top
    zl = selectelem(fens, bbfes, box = [0.0 0.0 -Inf Inf TH TH], inflate=tolerance)
    # Note that  we have to apply only half of the line load given that
    # were modeling  just one quarter of the geometry and we are splitting
    # the line load  with the symmetry plane X=0. Also note that the
    # quadrature rule is one-dimensional  since we are integrating along
    # a curve.
    Trac = FDataDict("traction_vector"=>vec([0.0; 0.0; -q0/2]), "femm"=>FEMMBase(IntegData(subset(bbfes, zl), GaussRule(1, 3))))
    
    modeldata = FDataDict("fens"=>fens, "regions"=>[region1, region2], "essential_bcs"=>[ex0, ey0, ez0], "traction_bcs"=> [Trac])
    modeldata = AlgoDeforLinearModule.linearstatics(modeldata)
    
    modeldata["postprocessing"] = FDataDict("file"=>"NAFEMS-R0031-1-plate")
    modeldata = AlgoDeforLinearModule.exportdeformation(modeldata)
    
    u = modeldata["u"]
    geom = modeldata["geom"]
    
    # The results of the displacement and stresses will be reported at
    # nodes located at the appropriate points.
    nE = selectnode(fens, box=[0.0 0.0 0.0 0.0 0.0 0.0], inflate=tolerance)
    nC = selectnode(fens, box=[0.0 0.0 0.0 0.0 TH TH], inflate=tolerance)
    nD = selectnode(fens, box=[CD CD 0.0 0.0 ts[1] ts[1]], inflate=tolerance)
    n0z = selectnode(fens, box=[0.0 0.0 0.0 0.0 0.0 TH], inflate=tolerance)
    ix = sortperm(geom.values[n0z, 3])
    # println("ix = $(ix)")
    
    cdis = mean(u.values[nE, 3])
    println("")
    println("Normalized Center deflection: $(cdis/wEref)")
    
    extrap = :extraptrend
    # # extrap = :extrapmean
    inspectormeth = :averaging
    # extrap = :default
    # inspectormeth = :invdistance
    
    modeldata["postprocessing"] = FDataDict("file"=>"NAFEMS-R0031-1-plate-sx",
    "quantity"=>:Cauchy, "component"=>1, "outputcsys"=>CSys(3),
    "nodevalmethod"=>inspectormeth, "reportat"=>extrap)
    modeldata = AlgoDeforLinearModule.exportstress(modeldata)
    s = modeldata["postprocessing"]["exported"][1]["field"]
    println("sx@E = $(s.values[nE]/phun("MPa")) [MPa]")
    
    modeldata["postprocessing"] = FDataDict("file"=>"NAFEMS-R0031-1-plate-sxz",
    "quantity"=>:Cauchy, "component"=>5, "outputcsys"=>CSys(3),
    "nodevalmethod"=>inspectormeth, "reportat"=>extrap)
    modeldata = AlgoDeforLinearModule.exportstress(modeldata)
    s = modeldata["postprocessing"]["exported"][1]["field"]
    println("sxz@D_1 = $(s.values[nD]/phun("MPa")) [MPa]")
    s = modeldata["postprocessing"]["exported"][2]["field"]
    println("sxz@D_2 = $(s.values[nD]/phun("MPa")) [MPa]")
    
    
    #
    # s = fieldfromintegpoints(region1["femm"], geom, u, :Cauchy, 1;
    #     outputcsys = CSys(3), nodevalmethod = inspectormeth, reportat = extrap)
    # println("sx@E = $(s.values[nE]/phun("MPa")) [MPa]")
    # sx_z = s.values[n0z]/phun("MPa")
    # println("sx(z)_1 = $(sx_z)")
    #
    # s = fieldfromintegpoints(region1["femm"], geom, u, :Cauchy, 5;
    #     outputcsys = CSys(3), nodevalmethod = inspectormeth, reportat = extrap)
    # println("sxz@D_1 = $(s.values[nD]/phun("MPa")) [MPa]")
    # sxz_z_1 = s.values[n0z]/phun("MPa")
    # println("sxz(z)_1 = $(sxz_z_1)")
    # s = fieldfromintegpoints(region2["femm"], geom, u, :Cauchy, 5;
    #     outputcsys = CSys(3), nodevalmethod = inspectormeth, reportat = extrap)
    # println("sxz@D_2 = $(s.values[nD]/phun("MPa")) [MPa]")
    # sxz_z_2 = s.values[n0z]/phun("MPa")
    # println("sxz(z)_2 = $(sxz_z_2)")
    
    # function _inspector(idat, elnum, conn, xe,  out,  xq)
    #     # xe = coordinates of the nodes of the element
    #     # xq = coordinate of the quadrature point
    #     println("@$(xq): $(out/1.0e6)")
    #     return idat
    # end
    #
    # felist = selectelem(fens, region1["femm"].integdata.fes,
    #     box=[0.0 0.0 0.0 0.0 0.0 0.0], inflate=tolerance, allin = false)
    #
    # inspectintegpoints(region1["femm"], geom, u, felist,
    #     _inspector, 0, quantity=:Cauchy, outputcsys = CSys(3))
    #
    # femm = deepcopy(region1["femm"])
    # femm.integdata.fes = subset(femm.integdata.fes, felist)
    # associategeometry!(femm, geom)
    # s = fieldfromintegpoints(femm, geom, u, :Cauchy, 5;
    #     outputcsys = CSys(3), nodevalmethod = inspectormeth, reportat = extrap)
    # println("sxz@D_1 = $(s.values[nD]/phun("MPa")) [MPa]")
    
    # felist = selectelem(fens, region2["femm"].integdata.fes,
    #     box=[0.0 0.0 0.0 0.0 0.0 TH], inflate=tolerance, allin = false)
    #
    # inspectintegpoints(region2["femm"], geom, u, felist,
    #     _inspector, 0, quantity=:Cauchy, outputcsys = CSys(3))
    
    
    println("Done")
    true
    
end # NAFEMS_R0031_1


function NAFEMS_R0031_1_H20()
    println("""
    Laminated Strip Under Three-Point Bending
    """)
    
    # Determine the central transverse displacement in a simply-supported seven
    # layer symmetric strip with a central line load. A 0/90/0/90/0/90/0
    # material lay-up is specified with the center ply being four times as
    # thick as the others.
    # Reference: NAFEMS Report R0031, Test No.1, 17-Dec-1998.
    
    # Because of the symmetries of the geometry and load, only the
    # first-quadrant   (in XY) quarter of the plate is modeled.
    
    # The coordinate system is centered at point E (at the difference with
    # respect to the original benchmark definition).  The  load is applied
    # along a curve passing through point C. The simple support is applied
    # along the curve passing through point B.
    
    # We realize the simple supports along the lines  A, B and the line load at
    # point C  are illegal from the point of view of convergence.  No
    # convergence can be hoped for as the stress underneath the load and above
    # the simple supports  is infinite in the limit (these locations are stress
    # singularities).   However, for relatively coarse meshes the results away
    # from the singularities are still meaningful.
    
    # The target quantities are displacement at the bottom surface at point E,
    # the tensile axial stress at the same point,  and of the transverse shear
    # stress at point D  in between the bottommost two layers (See figure 1).
    
    t0 = time()
    # Orthotropic material parameters of the material of the layers
    E1s = 100.0*phun("GPa")
    E2s = E3s = 5.0*phun("GPa")
    nu12s = 0.4
    nu13s = 0.3
    nu23s = 0.3
    G12s = 3.0*phun("GPa")
    G13s = G23s = 2.0*phun("GPa")
    CTE1 = 3.0e-6
    CTE2 = 2.0e-5
    CTE3 = 2.0e-5
    
    AB = 30.0*phun("mm") # span between simple supports
    CD = 4.0*phun("mm") # distance between the point D and the center
    OH = 10.0*phun("mm") # overhang
    W = 10.0*phun("mm") # width of the strip
    
    # Here we define the layout and the thicknesses of the layers.
    angles = vec([0 90 0 90 0 90 0]);
    ts = vec([0.1  0.1  0.1  0.4  0.1  0.1  0.1])*phun("mm"); # layer thicknesses
    TH = sum(ts); # total thickness of the plate
    
    tolerance = 0.0001*TH
    
    # The line load is in the negative Z direction.
    q0 = 10*phun("N/mm"); #    line load
    
    # Reference deflection under the load is
    wEref = -1.06*phun("mm");
    
    # The reference tensile stress at the bottom of the lowest layer is
    sigma11Eref = 684*phun("MPa");
    
    # Because we model the first-quadrant quarter of the plate using
    # coordinate axes centered  at the point E  the shear at the point D is
    # positive instead of negative as in the benchmark where the coordinate
    # system is located at the outer corner of the strip.
    sigma13Dref=4.1*phun("MPa");
    
    Refinement = 8
    # We select 8 elements spanwise and 2 elements widthwise.  The overhang
    # of the plate is given one element.
    nL = Refinement * 4; nO = Refinement * 1; nW = Refinement * 1;
    
    # Each layer is modeled with a single element.
    nts= Refinement * ones(Int, length(angles));# number of elements per layer
    
    xs = unique(vcat(collect(linspace(0,AB/2,nL+1)), [CD],
    collect(linspace(AB/2,AB/2+OH,nO+1))))
    xs = xs[sortperm(xs)]
    ys = collect(linspace(0,W/2,nW+1));
    
    fens,fes = H8layeredplatex(xs, ys, ts, nts)
    fens,fes = H8toH20(fens,fes)
    
    # This is the material  model
    MR = DeforModelRed3D
    material = MatDeforElastOrtho(MR,
    0.0, E1s, E2s, E3s,
    nu12s, nu13s, nu23s,
    G12s, G13s, G23s,
    CTE1, CTE2, CTE3)
    
    # The material coordinate system function is defined as:
    function updatecs!(csmatout::FFltMat, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)
        rotmat3!(csmatout, angles[fe_label]/180.0*pi* [0.0; 0.0; 1.0]);
    end
    
    # The vvolume integrals are evaluated using this rule
    gr = GaussRule(3, 3)
    
    # We will create two regions, one for the layers with 0°  orientation,
    # and one for the layers with 90° orientation.
    rl1 = vcat(selectelem(fens, fes, label = 1), selectelem(fens, fes, label = 3),
    selectelem(fens, fes, label = 5), selectelem(fens, fes, label = 7))
    rl2 = vcat(selectelem(fens, fes, label = 2), selectelem(fens, fes, label = 4),
    selectelem(fens, fes, label = 6))
    region1 = FDataDict("femm"=>FEMMDeforLinear(MR, IntegData(subset(fes, rl1), gr), CSys(3, 3, updatecs!), material))
    region2 = FDataDict("femm"=>FEMMDeforLinear(MR, IntegData(subset(fes, rl2), gr), CSys(3, 3, updatecs!), material))
    
    # File =  "NAFEMS-R0031-1-plate-r1.vtk"
    # vtkexportmesh(File, region1["femm"].integdata.fes.conn, fens.xyz, FinEtools.MeshExportModule.H8)
    # # @async run(`"paraview.exe" $File`)
    # File =  "NAFEMS-R0031-1-plate-r2.vtk"
    # vtkexportmesh(File, region2["femm"].integdata.fes.conn, fens.xyz, FinEtools.MeshExportModule.H8)
    # @async run(`"paraview.exe" $File`)
    
    # The essential boundary conditions are applied on the symmetry planes.
    # First the plane X=0;...
    lx0 = selectnode(fens, box=[0.0 0.0 -Inf Inf -Inf Inf], inflate=tolerance)
    ex0 = FDataDict( "displacement"=>  0.0, "component"=> 1, "node_list"=>lx0 )
    # ... and then the plane Y=0.
    ly0 = selectnode(fens, box=[-Inf Inf 0.0 0.0 -Inf Inf], inflate=tolerance)
    ey0 = FDataDict( "displacement"=>  0.0, "component"=> 2, "node_list"=>ly0 )
    # The transverse displacement is fixed along the line  passing through
    # point B. The nodes are fixed in the box along this line in the Z
    # direction.
    lz0 = selectnode(fens, box=[AB/2 AB/2 -Inf Inf 0.0 0.0], inflate=tolerance)
    ez0 = FDataDict( "displacement"=>  0.0, "component"=> 3, "node_list"=>lz0 )
    
    # The traction boundary condition is applied  along  the edge of the
    # mesh passing through point C at the top surface of the strip.   First
    # we extract the boundary of the hexahedral mesh.
    bfes = meshboundary(fes)
    # From  the entire boundary we select those quadrilaterals that lie on the plane
    # X = 0
    xl = selectelem(fens, bfes, box = [0.0 0.0 -Inf Inf -Inf Inf], inflate=tolerance)
    # Now we extract the boundary  of these selected quadrilaterals
    bbfes = meshboundary(subset(bfes, xl))
    # …  And from these  we extract the ones at the top
    zl = selectelem(fens, bbfes, box = [0.0 0.0 -Inf Inf TH TH], inflate=tolerance)
    # Note that  we have to apply only half of the line load given that
    # were modeling  just one quarter of the geometry and we are splitting
    # the line load  with the symmetry plane X=0. Also note that the
    # quadrature rule is one-dimensional  since we are integrating along
    # a curve.
    Trac = FDataDict("traction_vector"=>vec([0.0; 0.0; -q0/2]),    "femm"=>FEMMBase(IntegData(subset(bbfes, zl), GaussRule(1, 3))))
    
    modeldata = FDataDict("fens"=>fens, "regions"=>[region1, region2], "essential_bcs"=>[ex0, ey0, ez0], "traction_bcs"=> [Trac])
    modeldata = AlgoDeforLinearModule.linearstatics(modeldata)
    
    modeldata["postprocessing"] = FDataDict("file"=>"NAFEMS-R0031-1-plate")
    modeldata = AlgoDeforLinearModule.exportdeformation(modeldata)
    
    u = modeldata["u"]
    geom = modeldata["geom"]
    
    # The results of the displacement and stresses will be reported at
    # nodes located at the appropriate points.
    nE = selectnode(fens, box=[0.0 0.0 0.0 0.0 0.0 0.0], inflate=tolerance)
    nC = selectnode(fens, box=[0.0 0.0 0.0 0.0 TH TH], inflate=tolerance)
    nD = selectnode(fens, box=[CD CD 0.0 0.0 ts[1] ts[1]], inflate=tolerance)
    
    cdis = mean(u.values[nE, 3])
    println("")
    println("Normalized Center deflection: $(cdis/wEref)")
    
    modeldata["postprocessing"] = FDataDict("file"=>"NAFEMS-R0031-1-plate-sx",
    "quantity"=>:Cauchy, "component"=>1, "outputcsys"=>CSys(3))
    modeldata = AlgoDeforLinearModule.exportstress(modeldata)
    s = modeldata["postprocessing"]["exported"][1]["field"]
    println("sx@E = $(s.values[nE]/phun("MPa")) [MPa]")
    
    modeldata["postprocessing"] = FDataDict("file"=>"NAFEMS-R0031-1-plate-sxz",
    "quantity"=>:Cauchy, "component"=>5, "outputcsys"=>CSys(3))
    modeldata = AlgoDeforLinearModule.exportstress(modeldata)
    s = modeldata["postprocessing"]["exported"][1]["field"]
    println("sxz@D = $(s.values[nD]/phun("MPa")) [MPa]")
    s = modeldata["postprocessing"]["exported"][2]["field"]
    println("sxz@D = $(s.values[nD]/phun("MPa")) [MPa]")
    
    println("Done")
    true
    
end # NAFEMS_R0031_1_H20

function NAFEMS_R0031_2_both()
    println("""
    R0031(2): Wrapped thick cylinder under pressure and thermal loading.
    """)
    
    t0 = time()
    pu = ustring -> phun(ustring; system_of_units = :SIMM)
    # Orthotropic material parameters of the outer cylinder
    E1s = 130.0*pu("GPa")
    E2s = 5.0*pu("GPa")
    E3s = E2s
    nu12s = nu13s = 0.25
    nu23s = 0.0
    G12s = 10.0*pu("GPa")
    G13s = G12s
    G23s = 5.0*pu("GPa")
    CTE1 = 3.0e-6
    CTE2 = 2.0e-5
    CTE3 = 2.0e-5
    # Isotropic material parameters of the inner cylinder
    E = 2.1e5*pu("MPa")
    nu = 0.3
    CTE = 2.0e-5
    
    L = 200.0*pu("mm") # length of the cylinder
    ri = 23.0*pu("mm") # inner radius of the cylinder
    ti = 2.0*pu("mm") # thickness of the inner cylinder
    te = 2.0*pu("mm") # thickness of the outer cylinder
    q0 = 200.0*pu("MPa") # inner pressure
    dT = 130*pu("K") # temperature rise
    
    tolerance = 0.0001*ti
    
    # Generate mesh
    nL = 16 # number of elements lengthwise
    nc = 8 # number of elements circumferentially
    xs = collect(linspace(0.0, L/2, nL+1))
    ys = collect(linspace(0.0, pi/2, nc+1))
    ts = [ti; te];# layer thicknesses
    nts= 6*ones(Int, length(ts));# number of elements per layer
    fens,fes = H8layeredplatex(xs, ys, ts, nts)
    fens,fes = H8toH20(fens,fes)
    bfes = meshboundary(fes)
    # inner surface  for the pressure loading
    intl = selectelem(fens, bfes; facing=true, direction = [0.0 0.0 -1.0])
    # Shape into a cylinder
    for i = 1:count(fens)
        z = fens.xyz[i,1]; a = fens.xyz[i,2]; t = fens.xyz[i,3];
        fens.xyz[i,:] = [(ri+t)*cos(pi/2-a) (ri+t)*sin(pi/2-a) z];
    end
    
    MR = DeforModelRed3D
    outermaterial = MatDeforElastOrtho(MR,
    0.0, E1s, E2s, E3s,
    nu12s, nu13s, nu23s,
    G12s, G13s, G23s,
    CTE1, CTE2, CTE3)
    innermaterial = MatDeforElastIso(MR,
    0.0, E, nu, CTE)
    
    function cylcs!(csmatout::FFltMat, XYZ::FFltMat)
        csmatout[:, 2] = [0.0 0.0 1.0]
        csmatout[:, 3] = XYZ
        csmatout[3, 3] = 0.0
        csmatout[:, 3] = csmatout[:, 3]/norm(csmatout[:, 3])
        csmatout[:, 1] = cross(csmatout[:, 2], csmatout[:, 3])
    end
    
    function updatecs!(csmatout::FFltMat, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)
        cylcs!(csmatout, XYZ)
    end
    
    gr = GaussRule(3, 3)
    
    rli = selectelem(fens, fes, label=1)
    innerregion = FDataDict("femm"=>FEMMDeforLinear(MR, IntegData(subset(fes, rli), gr), innermaterial))
    rle = selectelem(fens, fes, label=2)
    outerregion = FDataDict("femm"=>FEMMDeforLinear(MR, IntegData(subset(fes, rle), gr), CSys(3, 3, updatecs!), outermaterial))
    
    lx0 = selectnode(fens, box=[0.0 0.0 -Inf Inf -Inf Inf], inflate=tolerance)
    ly0 = selectnode(fens, box=[-Inf Inf 0.0 0.0 -Inf Inf], inflate=tolerance)
    lz0 = selectnode(fens, box=[-Inf Inf -Inf Inf 0.0 0.0], inflate=tolerance)
    
    ex0 = FDataDict( "displacement"=>  0.0, "component"=> 1, "node_list"=>lx0 )
    ey0 = FDataDict( "displacement"=>  0.0, "component"=> 2, "node_list"=>ly0 )
    ez0 = FDataDict( "displacement"=>  0.0, "component"=> 3, "node_list"=>lz0 )
    
    function getpr!(forceout::FFltVec, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)
        csmatout = zeros(3, 3)
        cylcs!(csmatout, XYZ)
        copy!(forceout, q0*csmatout[:, 3])
    end
    
    Trac = FDataDict("traction_vector"=>getpr!, "femm"=>FEMMBase(IntegData(subset(bfes, intl), GaussRule(2, 3))))
    
    modeldata = FDataDict("fens"=>fens, "regions"=>[innerregion, outerregion], "essential_bcs"=>[ex0, ey0, ez0], "traction_bcs"=>[Trac], "temperature_change"=>FDataDict("temperature"=>dT))
    modeldata = AlgoDeforLinearModule.linearstatics(modeldata)
    
    u = modeldata["u"]
    geom = modeldata["geom"]
    # File =  "NAFEMS-R0031-2-plate.vtk"
    # vtkexportmesh(File, fes.conn, geom.values, FinEtools.MeshExportModule.H20;
    #     scalars = [("Layer", fes.label)], vectors = [("displacement", u.values)])
    # @async run(`"paraview.exe" $File`)
    
    modeldata["postprocessing"] = FDataDict("file"=>"NAFEMS-R0031-2", "outputcsys"=>CSys(3, 3, updatecs!), "quantity"=>:Cauchy, "component"=>6)
    modeldata = AlgoDeforLinearModule.exportstress(modeldata)
    modeldata["postprocessing"] = FDataDict("file"=>"NAFEMS-R0031-2-elem", "outputcsys"=>CSys(3, 3, updatecs!), "quantity"=>:Cauchy, "component"=>6)
    modeldata = AlgoDeforLinearModule.exportstresselementwise(modeldata)
    
    println("Done")
    true
    
end # NAFEMS_R0031_2_both


function NAFEMS_R0031_2_pressure()
    println("""
    R0031(2): Wrapped thick cylinder under pressure and thermal loading.
    """)
    
    t0 = time()
    # Orthotropic material parameters of the outer cylinder
    E1s = 130.0*phun("GPa")
    E2s = 5.0*phun("GPa")
    E3s = E2s
    nu12s = nu13s = 0.25
    nu23s = 0.0
    G12s = 10.0*phun("GPa")
    G13s = G12s
    G23s = 5.0*phun("GPa")
    CTE1 = 3.0e-6
    CTE2 = 2.0e-5
    CTE3 = 2.0e-5
    # Isotropic material parameters of the inner cylinder
    E = 2.1e5*phun("MPa")
    nu = 0.3
    CTE = 2.0e-5
    
    L = 200.0*phun("mm") # length of the cylinder
    ri = 23.0*phun("mm") # inner radius of the cylinder
    ti = 2.0*phun("mm") # thickness of the inner cylinder
    te = 2.0*phun("mm") # thickness of the outer cylinder
    q0 = 200.0*phun("MPa") # inner pressure
    dT = 0*phun("K") # NO temperature rise
    
    tolerance = 0.0001*ti
    
    # Generate mesh
    nL = 18 # number of elements lengthwise
    nc = 18 # number of elements circumferentially
    xs = collect(linspace(0.0, L/2, nL+1))
    ys = collect(linspace(0.0, pi/2, nc+1))
    ts = [ti; te];# layer thicknesses
    nts= 6*ones(Int, length(ts));# number of elements per layer
    fens,fes = H8layeredplatex(xs, ys, ts, nts)
    fens,fes = H8toH20(fens,fes)
    bfes = meshboundary(fes)
    # inner surface  for the pressure loading
    intl = selectelem(fens, bfes; facing=true, direction = [0.0 0.0 -1.0])
    # Shape into a cylinder
    for i = 1:count(fens)
        z = fens.xyz[i,1]; a = fens.xyz[i,2]; t = fens.xyz[i,3];
        fens.xyz[i,:] = [(ri+t)*cos(pi/2-a) (ri+t)*sin(pi/2-a) z];
    end
    
    MR = DeforModelRed3D
    outermaterial = MatDeforElastOrtho(MR,
    0.0, E1s, E2s, E3s,
    nu12s, nu13s, nu23s,
    G12s, G13s, G23s,
    CTE1, CTE2, CTE3)
    innermaterial = MatDeforElastIso(MR,
    0.0, E, nu, CTE)
    
    function cylcs!(csmatout::FFltMat, XYZ::FFltMat)
        csmatout[:, 2] = [0.0 0.0 1.0]
        radial = XYZ; radial[3] = 0.0
        csmatout[:, 3] = radial/norm(radial)
        csmatout[:, 1] = cross(csmatout[:, 2], csmatout[:, 3])
    end
    
    function updatecs!(csmatout::FFltMat, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)
        cylcs!(csmatout, XYZ)
    end
    
    gr = GaussRule(3, 3)
    
    rli = selectelem(fens, fes, label=1)
    innerregion = FDataDict("femm"=>FEMMDeforLinear(MR,
    IntegData(subset(fes, rli), gr), innermaterial))
    rle = selectelem(fens, fes, label=2)
    outerregion = FDataDict("femm"=>FEMMDeforLinear(MR,
    IntegData(subset(fes, rle), gr), CSys(3, 3, updatecs!), outermaterial))
    
    lx0 = selectnode(fens, box=[0.0 0.0 -Inf Inf -Inf Inf], inflate=tolerance)
    ly0 = selectnode(fens, box=[-Inf Inf 0.0 0.0 -Inf Inf], inflate=tolerance)
    lz0 = selectnode(fens, box=[-Inf Inf -Inf Inf 0.0 0.0], inflate=tolerance)
    
    ex0 = FDataDict( "displacement"=>  0.0, "component"=> 1, "node_list"=>lx0 )
    ey0 = FDataDict( "displacement"=>  0.0, "component"=> 2, "node_list"=>ly0 )
    ez0 = FDataDict( "displacement"=>  0.0, "component"=> 3, "node_list"=>lz0 )
    
    function getpr!(forceout::FFltVec, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)
        csmatout = zeros(3, 3)
        cylcs!(csmatout, XYZ)
        copy!(forceout, q0*csmatout[:, 3])
    end
    
    Trac = FDataDict("traction_vector"=>getpr!, "femm"=>FEMMBase(IntegData(subset(bfes, intl), GaussRule(2, 3))))
    
    modeldata = FDataDict("fens"=>fens,
    "regions"=>[innerregion, outerregion],
    "essential_bcs"=>[ex0, ey0, ez0],
    "traction_bcs"=>[Trac],
    "temperature_change"=>FDataDict("temperature"=>dT)
    )
    modeldata = AlgoDeforLinearModule.linearstatics(modeldata)
    
    u = modeldata["u"]
    geom = modeldata["geom"]
    # lcenter = selectnode(fens, box=[a/2 a/2  b/2 b/2 -Inf Inf], inflate=tolerance)
    # cdis = abs(mean(u.values[lcenter, 3]))
    # println("")
    # println("Normalized Center deflection: $(cdis/wc_analytical)")
    
    File =  "NAFEMS-R0031-2-plate.vtk"
    vtkexportmesh(File, fes.conn, geom.values, FinEtools.MeshExportModule.H20; scalars = [("Layer", fes.label)], vectors = [("displacement", u.values)])
    @async run(`"paraview.exe" $File`)
    
    println("Done")
    true
    
end # NAFEMS_R0031_2_pressure


function NAFEMS_R0031_3()
    println("""
    NAFEMS publication R0031/3 Composite plate test.
    Simply supported on all four edges.  Uniform transverse  loading.
    The modeled part is one quarter of the full plate here.
    """)
    
    # This is a test recommended by the National Agency for Finite Element Methods
    # and Standards (U.K.): Test R0031/3 from NAFEMS publication R0031, “Composites
    # Benchmarks,” February 1995.
    t0 = time()
    # Skin (face) material parameters
    E1s = 1.0e7*phun("psi")
    E2s = 0.4e7*phun("psi")
    E3s = 0.4e7*phun("psi")
    nu12s = 0.3
    nu13s = 0.3
    nu23s = 0.3
    G12s = 0.1875e7*phun("psi")
    G13s = 0.1875e7*phun("psi")
    G23s = 0.1875e7*phun("psi")
    # Core material parameters
    E1c = 10.*phun("psi")
    E2c = 10.*phun("psi")
    E3c = 10e4.*phun("psi")
    nu12c = 0.
    nu13c = 0.
    nu23c = 0.
    G12c = 10.*phun("psi")
    G13c = 3.0e4*phun("psi")
    G23c = 1.2e4*phun("psi")
    L = 10.0*phun("in") # side of the square plate
    nL = 8 # number of elements along the side of the plate
    tolerance = 0.0001*phun("in")
    xs = collect(linspace(0.0, L/2, nL+1))
    ys = collect(linspace(0.0, L/2, nL+1))
    ts = [0.028; 0.75; 0.028]*phun("in")
    nts = [2; 3;  2; ] # number of elements through the thickness
    tmag = 100*phun("psi")
    
    # Generate mesh
    fens,fes = H8layeredplatex(xs, ys, ts, nts)
    fens,fes = H8toH20(fens,fes)
    
    MR = DeforModelRed3D
    skinmaterial = MatDeforElastOrtho(MR,
    0.0, E1s, E2s, E3s,
    nu12s, nu13s, nu23s,
    G12s, G13s, G23s,
    0.0, 0.0, 0.0)
    corematerial = MatDeforElastOrtho(MR,
    0.0, E1c, E2c, E3c,
    nu12c, nu13c, nu23c,
    G12c, G13c, G23c,
    0.0, 0.0, 0.0)
    
    gr = GaussRule(3, 3)
    
    rl1 = selectelem(fens, fes, label=1)
    skinbot = FDataDict("femm"=>FEMMDeforLinear(MR, IntegData(subset(fes, rl1), gr), skinmaterial))
    
    rl3 = selectelem(fens, fes, label=3)
    skintop = FDataDict("femm"=>FEMMDeforLinear(MR, IntegData(subset(fes, rl3), gr), skinmaterial))
    
    rl2 = selectelem(fens, fes, label=2)
    core = FDataDict("femm"=>FEMMDeforLinear(MR, IntegData(subset(fes, rl2), gr), corematerial))
    
    lx0 = selectnode(fens, box=[0.0 0.0 -Inf Inf -Inf Inf], inflate=tolerance)
    lxL2 = selectnode(fens, box=[L/2 L/2 -Inf Inf -Inf Inf], inflate=tolerance)
    ly0 = selectnode(fens, box=[-Inf Inf 0.0 0.0 -Inf Inf], inflate=tolerance)
    lyL2 = selectnode(fens, box=[-Inf Inf L/2 L/2 -Inf Inf], inflate=tolerance)
    
    ex0 = FDataDict( "displacement"=>  0.0, "component"=> 3, "node_list"=>lx0 )
    exL2 = FDataDict( "displacement"=>  0.0, "component"=> 1, "node_list"=>lxL2 )
    ey0 = FDataDict( "displacement"=>  0.0, "component"=> 3, "node_list"=>ly0 )
    eyL2 = FDataDict( "displacement"=>  0.0, "component"=> 2, "node_list"=>lyL2 )
    
    bfes = meshboundary(fes)
    ttopl = selectelem(fens, bfes; facing=true, direction = [0.0 0.0 1.0])
    Trac = FDataDict("traction_vector"=>[0.0; 0.0; -tmag], "femm"=>FEMMBase(IntegData(subset(bfes, ttopl), GaussRule(2, 3))))
    
    modeldata = FDataDict("fens"=>fens, "regions"=>[skinbot, core, skintop], "essential_bcs"=>[ex0, exL2, ey0, eyL2], "traction_bcs"=> [Trac])
    modeldata = AlgoDeforLinearModule.linearstatics(modeldata)
    
    u = modeldata["u"]
    geom = modeldata["geom"]
    lcenter = selectnode(fens, box=[L/2 L/2  L/2 L/2 -Inf Inf], inflate=tolerance)
    cdis = mean(u.values[lcenter, 3])/phun("in")
    println("Center node displacements $(cdis) [in]; NAFEMS-R0031-3 lists –0.123	[in]")
    println("")
    
    File =  "NAFEMS-R0031-3-plate.vtk"
    vtkexportmesh(File, fes.conn, geom.values, FinEtools.MeshExportModule.H20; scalars = [("Layer", fes.label)], vectors = [("displacement", u.values)])
    @async run(`"paraview.exe" $File`)
    
    println("Done")
    true
    
end # NAFEMS_R0031_3

function allrun()
    println("#####################################################") 
    println("# NAFEMS_R0031_1 ")
    NAFEMS_R0031_1()
    println("#####################################################") 
    println("# NAFEMS_R0031_1_H20 ")
    NAFEMS_R0031_1_H20()
    println("#####################################################") 
    println("# NAFEMS_R0031_2_pressure ")
    NAFEMS_R0031_2_pressure()
    println("#####################################################") 
    println("# NAFEMS_R0031_2_both ")
    NAFEMS_R0031_2_both()
    println("#####################################################") 
    println("# NAFEMS_R0031_3 ")
    NAFEMS_R0031_3()
    return true
end # function allrun

end # module NAFEMS_R0031_examples
