module LE11NAFEMS_examples
using FinEtools
using FinEtools.AlgoDeforLinearModule

function LE11NAFEMS_Q8_algo()
    # NAFEMS LE11 benchmark with Q8 elements.
    # # This is a test recommended by the National Agency for Finite Element
    # # Methods and Standards (U.K.): Test LE11 from NAFEMS Publication TNSB,
    # # Rev. 3, “The Standard NAFEMS Benchmarks,” October 1990.
    # #
    # # Target solution: Direct stress,   =  –105 MPa at point A.
    #function  LE11NAFEMS()
    # Parameters:
    Ea =  210000*phun("MEGA*Pa")
    nua =  0.3;
    alphaa = 2.3e-4;              # thermal expansion coefficient
    sigmaA = -105*phun("MEGA*Pa")
    nref =  1;                        # how many times should we refine the mesh?
    X = [1.     0.;#A
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
    tolerance  = 1.e-6*phun("M")
    ##
    # Note that the material object needs to be created with the proper
    # model-dimension reduction in mind.  In this case that is the axial symmetry
    # assumption.
    MR  =  DeforModelRed2DAxisymm
    
    fens = FENodeSet(X);
    fes = FESetQ4([1 2 4 3; 3 4 6 5; 5 6 8 7; 7 8 10 9; 9 10 12 11; 11 12 14 13; 13 14 16 15]);
    for ref = 1:nref
        fens,fes = Q4refine(fens,fes);
        list = selectnode(fens,distance = 1.0+0.1/2^nref, from = [0. 0.], inflate = tolerance);
        fens.xyz[list,:] =  FinEtools.MeshUtilModule.ontosphere(fens.xyz[list,:],1.0);
    end
    fens,fes = Q4toQ8(fens,fes)
    list = selectnode(fens,distance = 1.0+0.1/2^nref, from = [0. 0.], inflate = tolerance);
    fens.xyz[list,:] =  FinEtools.MeshUtilModule.ontosphere(fens.xyz[list,:],1.0);
    
    
    # EBC's
    l1 = selectnode(fens,box = [-Inf Inf 0 0],inflate = tolerance)
    e1 = FDataDict("node_list"=>l1, "component"=>2, "displacement"=>x -> 0.0)
    l1 = selectnode(fens,box = [-Inf Inf 1.79  1.79],inflate = tolerance)
    e2 = FDataDict("node_list"=>l1, "component"=>2, "displacement"=>x -> 0.0)
    
    # Temperature field
    dtemp = FDataDict("temperature"=>x -> x[1] + x[2])
    
    # Property and material
    material = MatDeforElastIso(MR, 0.0, Ea, nua, alphaa)
    
    femm  =  FEMMDeforLinear(MR, IntegData(fes, GaussRule(2, 3), true), material)
    
    # Make region 1
    region = FDataDict("femm"=>femm)
    # Make model data
    modeldata =  FDataDict(
    "fens"=> fens, "regions"=>  [region],
    "essential_bcs"=>[e1, e2], "temperature_change"=>dtemp)
    
    # Call the solver
    modeldata = AlgoDeforLinearModule.linearstatics(modeldata)
    geom = modeldata["geom"]
    u = modeldata["u"]
    dT = modeldata["temp"]
    
    nA  = selectnode(fens,box = FFlt[1.0  1.0 0.0 0.0], inflate = tolerance);
    
    fld =  fieldfromintegpoints(femm, geom, u, dT, :Cauchy, 2)
    
    
    File  =   "LE11NAFEMS_Q8_sigmay.vtk"
    vtkexportmesh(File, fens, fes; scalars = [("sigmay", fld.values)],
    vectors = [("u", u.values)])
    println("range of  sigmay = $((minimum(fld.values), maximum(fld.values)))")
    @async run(`"paraview.exe" $File`)
    
    
    sA  =  fld.values[nA]/phun("MEGA*Pa")
    sAn  =  fld.values[nA]/sigmaA
    println("Stress at point A: $(sA) i. e.  $( sAn*100  )% of reference value")
    
    fen2fe  = FENodeToFEMap(fes.conn, nnodes(geom))
    function inspector(idat, elnum, conn, xe,  out,  xq)
        println("loc = $(  xq  ) : $(  transpose(out)/phun("MEGA*Pa")  )")
        return idat
    end
    
    inspectintegpoints(femm, geom, u, dT,  fen2fe.map[nA[1]],
    inspector, []; quantity = :Cauchy)
    
end # LE11NAFEMS_Q8_algo


function LE11NAFEMS_Q8_algo2()
    # NAFEMS LE11 benchmark with Q8 elements.
    # # This is a test recommended by the National Agency for Finite Element
    # # Methods and Standards (U.K.): Test LE11 from NAFEMS Publication TNSB,
    # # Rev. 3, “The Standard NAFEMS Benchmarks,” October 1990.
    # #
    # # Target solution: Direct stress,   =  –105 MPa at point A.
    #function  LE11NAFEMS()
    # Parameters:
    Ea =  210000*phun("MEGA*Pa")
    nua =  0.3;
    alphaa = 2.3e-4;              # thermal expansion coefficient
    sigmaA = -105*phun("MEGA*Pa")
    nref =  1;                        # how many times should we refine the mesh?
    X = [1.     0.;#A
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
    tolerance  = 1.e-6*phun("M")
    ##
    # Note that the material object needs to be created with the proper
    # model-dimension reduction in mind.  In this case that is the axial symmetry
    # assumption.
    MR  =  DeforModelRed2DAxisymm
    
    fens = FENodeSet(X);
    fes = FESetQ4([1 2 4 3; 3 4 6 5; 5 6 8 7; 7 8 10 9; 9 10 12 11; 11 12 14 13; 13 14 16 15]);
    for ref = 1:nref
        fens,fes = Q4refine(fens,fes);
        list = selectnode(fens,distance = 1.0+0.1/2^nref, from = [0. 0.], inflate = tolerance);
        fens.xyz[list,:] =  FinEtools.MeshUtilModule.ontosphere(fens.xyz[list,:],1.0);
    end
    fens,fes = Q4toQ8(fens,fes)
    list = selectnode(fens,distance = 1.0+0.1/2^nref, from = [0. 0.], inflate = tolerance);
    fens.xyz[list,:] =  FinEtools.MeshUtilModule.ontosphere(fens.xyz[list,:],1.0);
    
    
    # EBC's
    l1 = selectnode(fens,box = [-Inf Inf 0 0],inflate = tolerance)
    e1 = FDataDict("node_list"=>l1, "component"=>2, "displacement"=>x -> 0.0)
    l1 = selectnode(fens,box = [-Inf Inf 1.79  1.79],inflate = tolerance)
    e2 = FDataDict("node_list"=>l1, "component"=>2, "displacement"=>x -> 0.0)
    
    # Temperature field
    dtemp = FDataDict("temperature"=>x -> x[1] + x[2])
    
    # Property and material
    material = MatDeforElastIso(MR, 0.0, Ea, nua, alphaa)
    
    femm  =  FEMMDeforLinear(MR, IntegData(fes, GaussRule(2, 3), true), material)
    
    # Make region 1
    region = FDataDict("femm"=>femm)
    # Make model data
    modeldata =  FDataDict(
    "fens"=> fens, "regions"=>  [region],
    "essential_bcs"=>[e1, e2], "temperature_change"=>dtemp)
    
    # Call the solver
    modeldata = AlgoDeforLinearModule.linearstatics(modeldata)
    geom = modeldata["geom"]
    u = modeldata["u"]
    dT = modeldata["temp"]
    
    modeldata["postprocessing"] = FDataDict("boundary_only"=> true,
    "file"=>"LE11NAFEMS_Q8_deformation.vtk")
    modeldata = AlgoDeforLinearModule.exportdeformation(modeldata)
    @async run(`"paraview.exe" $(modeldata["postprocessing"]["exported"][1]["file"])`)
    
    nA  = selectnode(fens,box = FFlt[1.0  1.0 0.0 0.0], inflate = tolerance);
    
    modeldata["postprocessing"] = FDataDict("boundary_only"=> true,
    "file"=>"LE11NAFEMS_Q8_sigmay.vtk", "quantity"=>:Cauchy,
    "component"=>2)
    modeldata = AlgoDeforLinearModule.exportstress(modeldata)
    modeldata["postprocessing"] = FDataDict("boundary_only"=> false,
    "file"=>"LE11NAFEMS_Q8_sigmay.vtk", "quantity"=>:Cauchy,
    "component"=>2)
    modeldata = AlgoDeforLinearModule.exportstress(modeldata)
    @async run(`"paraview.exe" $(modeldata["postprocessing"]["exported"][1]["file"])`)
    fld =  modeldata["postprocessing"]["exported"][1]["field"]
    
    sA  =  fld.values[nA]/phun("MEGA*Pa")
    sAn  =  fld.values[nA]/sigmaA
    println("Stress at point A: $(sA) i. e.  $( sAn*100  )% of reference value")
    
    # Loop over only those elements that share the node nA
    fen2fe  = FENodeToFEMap(fes.conn, nnodes(geom))
    function inspector(idat, elnum, conn, xe,  out,  xq)
        println("loc = $(  xq  ) : $(  transpose(out)/phun("MEGA*Pa")  )")
        return idat
    end
    
    inspectintegpoints(femm, geom, u, dT,  fen2fe.map[nA[1]],
    inspector, []; quantity = :Cauchy)
    
    modeldata["postprocessing"] = FDataDict("boundary_only"=> false,
    "file"=>"LE11NAFEMS_Q8_sigmay_ew.vtk", "quantity"=>:Cauchy,
    "component"=>2)
    modeldata = AlgoDeforLinearModule.exportstresselementwise(modeldata)
    @async run(`"paraview.exe" $(modeldata["postprocessing"]["exported"][1]["file"])`)
    
end # LE11NAFEMS_Q8_algo2


function LE11NAFEMS_Q8_export_stress()
    # NAFEMS LE11 benchmark with Q8 elements.
    # # This is a test recommended by the National Agency for Finite Element
    # # Methods and Standards (U.K.): Test LE11 from NAFEMS Publication TNSB,
    # # Rev. 3, “The Standard NAFEMS Benchmarks,” October 1990.
    # #
    # # Target solution: Direct stress,   =  –105 MPa at point A.
    #function  LE11NAFEMS()
    # Parameters:
    Ea =  210000*phun("MEGA*Pa")
    nua =  0.3;
    alphaa = 2.3e-4;              # thermal expansion coefficient
    sigmaA = -105*phun("MEGA*Pa")
    nref = 2;                        # how many times should we refine the mesh?
    X = [1.     0.;#A
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
    tolerance  = 1.e-6*phun("M")
    ##
    # Note that the material object needs to be created with the proper
    # model-dimension reduction in mind.  In this case that is the axial symmetry
    # assumption.
    MR  =  DeforModelRed2DAxisymm
    
    fens = FENodeSet(X);
    fes = FESetQ4([1 2 4 3; 3 4 6 5; 5 6 8 7; 7 8 10 9; 9 10 12 11; 11 12 14 13; 13 14 16 15]);
    for ref = 1:nref
        fens,fes = Q4refine(fens,fes);
        list = selectnode(fens,distance = 1.0+0.1/2^nref, from = [0. 0.], inflate = tolerance);
        fens.xyz[list,:] =  FinEtools.MeshUtilModule.ontosphere(fens.xyz[list,:],1.0);
    end
    fens,fes = Q4toQ8(fens,fes)
    list = selectnode(fens,distance = 1.0+0.1/2^nref, from = [0. 0.], inflate = tolerance);
    fens.xyz[list,:] =  FinEtools.MeshUtilModule.ontosphere(fens.xyz[list,:],1.0);
    
    #     File  =   "mesh.vtk"
    # vtkexportmesh(File, fens, fes)
    
    # now we create the geometry and displacement fields
    geom  =  NodalField(fens.xyz)
    u  =  NodalField(zeros(size(fens.xyz,1),2)) # displacement field
    
    # Apply EBC's
    lbottom = selectnode(fens,box = [-Inf Inf 0 0],inflate = tolerance)
    setebc!(u, lbottom, true, 2, 00.0)
    ltop = selectnode(fens,box = [-Inf Inf 1.79  1.79],inflate = tolerance)
    setebc!(u, ltop, true, 2, 00.0)
    applyebc!(u)
    numberdofs!(u)
    
    # Temperature field
    dT  = NodalField(reshape(fens.xyz[:,1]+fens.xyz[:,2],size(fens.xyz,1),1));
    
    # Property and material
    material = MatDeforElastIso(MR, 0.0, Ea, nua, alphaa)
    
    femm  =  FEMMDeforLinear(MR, IntegData(fes, GaussRule(2, 3), true), material)
    
    K  = stiffness(femm, geom, u)
    F  =  thermalstrainloads(femm, geom, u, dT)
    #K = cholesky(K)
    U =   K\F
    scattersysvec!(u, U[:])
    
    nA  = selectnode(fens,box = FFlt[1.0  1.0 0.0 0.0], inflate = tolerance);
    
    fld =  fieldfromintegpoints(femm, geom, u, dT, :Cauchy, 2)
    
    
    File  =   "LE11NAFEMS_Q8_sigmay.vtk"
    vtkexportmesh(File, fens, fes; scalars = [("sigmay", fld.values)],
    vectors = [("u", u.values)])
    println("range of  sigmay = $((minimum(fld.values), maximum(fld.values)))")
    @async run(`"paraview.exe" $File`)
    
    
    sA  =  fld.values[nA]/phun("MEGA*Pa")
    sAn  =  fld.values[nA]/sigmaA
    println("Stress at point A: $(sA) i. e.  $( sAn*100  )% of reference value")
    
    fen2fe  = FENodeToFEMap(fes.conn, nnodes(geom))
    function inspector(idat, elnum, conn, xe,  out,  xq)
        # println("loc = $(  xq  ) : $(  transpose(out)/phun("MEGA*Pa")  )")
        return idat
    end
    
    inspectintegpoints(femm, geom, u, dT,  fen2fe.map[nA[1]],
    inspector, []; quantity = :Cauchy)
    
    fld =  fieldfromintegpoints(femm, geom, u, dT, :Pressure, 1)
    File  =   "LE11NAFEMS_Q8_pressure.vtk"
    vtkexportmesh(File, fens, fes; scalars = [("pressure", fld.values)],
    vectors = [("u", u.values)])
    println("range of  pressure = $((minimum(fld.values), maximum(fld.values)))")
    @async run(`"paraview.exe" $File`)
    
    
    fld =  fieldfromintegpoints(femm, geom, u, dT, :vm, 1)
    File  =   "LE11NAFEMS_Q8_vm.vtk"
    vtkexportmesh(File, fens, fes; scalars = [("pressure", fld.values)],
    vectors = [("u", u.values)])
    println("range of von Mises = $((minimum(fld.values), maximum(fld.values)))")
    @async run(`"paraview.exe" $File`)
    
    AE = AbaqusExporter("LE11NAFEMS_Q8_export_stress");
    HEADING(AE, "NAFEMS LE11 benchmark with Q8 elements.");
    COMMENT(AE, "sigmaA = -105 MPa ");
    PART(AE, "part1");
    END_PART(AE);
    ASSEMBLY(AE, "ASSEM1");
    INSTANCE(AE, "INSTNC1", "PART1");
    NODE(AE, fens.xyz);
    ELEMENT(AE, "cax8", "AllElements", 1, connasarray(fes))
    NSET_NSET(AE, "ltop", ltop)
    NSET_NSET(AE, "lbottom", lbottom)
    ORIENTATION(AE, "GlobalOrientation", vec([1. 0 0]), vec([0 1. 0]));
    SOLID_SECTION(AE, "elasticity", "GlobalOrientation", "AllElements");
    SURFACE_SECTION(AE, "TractionElements")
    END_INSTANCE(AE);
    END_ASSEMBLY(AE);
    MATERIAL(AE, "elasticity")
    ELASTIC(AE, Ea, nua)
    EXPANSION(AE, alphaa)
    STEP_PERTURBATION_STATIC(AE)
    BOUNDARY(AE, "ASSEM1.INSTNC1.ltop", 2)
    BOUNDARY(AE, "ASSEM1.INSTNC1.lbottom", 2)
    TEMPERATURE(AE, "ASSEM1.INSTNC1.", collect(1:count(fens)), vec(dT.values))
    END_STEP(AE)
    close(AE)
end # LE11NAFEMS_Q8_export_stress


function LE11NAFEMS_Q8_stress()
    # NAFEMS LE11 benchmark with Q8 elements.
    # # This is a test recommended by the National Agency for Finite Element
    # # Methods and Standards (U.K.): Test LE11 from NAFEMS Publication TNSB,
    # # Rev. 3, “The Standard NAFEMS Benchmarks,” October 1990.
    # #
    # # Target solution: Direct stress,   =  –105 MPa at point A.
    #function  LE11NAFEMS()
    # Parameters:
    Ea =  210000*phun("MEGA*Pa")
    nua =  0.3;
    alphaa = 2.3e-4;              # thermal expansion coefficient
    sigmaA = -105*phun("MEGA*Pa")
    nref =  1;                        # how many times should we refine the mesh?
    X = [1.     0.;#A
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
    tolerance  = 1.e-6*phun("M")
    ##
    # Note that the material object needs to be created with the proper
    # model-dimension reduction in mind.  In this case that is the axial symmetry
    # assumption.
    MR  =  DeforModelRed2DAxisymm
    
    fens = FENodeSet(X);
    fes = FESetQ4([1 2 4 3; 3 4 6 5; 5 6 8 7; 7 8 10 9; 9 10 12 11; 11 12 14 13; 13 14 16 15]);
    for ref = 1:nref
        fens,fes = Q4refine(fens,fes);
        list = selectnode(fens,distance = 1.0+0.1/2^nref, from = [0. 0.], inflate = tolerance);
        fens.xyz[list,:] =  FinEtools.MeshUtilModule.ontosphere(fens.xyz[list,:],1.0);
    end
    fens,fes = Q4toQ8(fens,fes)
    list = selectnode(fens,distance = 1.0+0.1/2^nref, from = [0. 0.], inflate = tolerance);
    fens.xyz[list,:] =  FinEtools.MeshUtilModule.ontosphere(fens.xyz[list,:],1.0);
    
    #     File  =   "mesh.vtk"
    # vtkexportmesh(File, fens, fes)
    
    # now we create the geometry and displacement fields
    geom  =  NodalField(fens.xyz)
    u  =  NodalField(zeros(size(fens.xyz,1),2)) # displacement field
    
    # Apply EBC's
    l1 = selectnode(fens,box = [-Inf Inf 0 0],inflate = tolerance)
    setebc!(u, l1, true, 2, 00.0)
    l1 = selectnode(fens,box = [-Inf Inf 1.79  1.79],inflate = tolerance)
    setebc!(u, l1, true, 2, 00.0)
    applyebc!(u)
    numberdofs!(u)
    
    # Temperature field
    dT  = NodalField(reshape(fens.xyz[:,1]+fens.xyz[:,2],size(fens.xyz,1),1));
    
    
    # Property and material
    material = MatDeforElastIso(MR, 0.0, Ea, nua, alphaa)
    
    femm  =  FEMMDeforLinear(MR, IntegData(fes, GaussRule(2, 3), true), material)
    
    K  = stiffness(femm, geom, u)
    F  =  thermalstrainloads(femm, geom, u, dT)
    #K = cholesky(K)
    U =   K\F
    scattersysvec!(u, U[:])
    
    nA  = selectnode(fens,box = FFlt[1.0  1.0 0.0 0.0], inflate = tolerance);
    
    fld =  fieldfromintegpoints(femm, geom, u, dT, :Cauchy, 2)
    
    
    File  =   "LE11NAFEMS_Q8_sigmay.vtk"
    vtkexportmesh(File, fens, fes; scalars = [("sigmay", fld.values)],
    vectors = [("u", u.values)])
    println("range of  sigmay = $((minimum(fld.values), maximum(fld.values)))")
    @async run(`"paraview.exe" $File`)
    
    
    sA  =  fld.values[nA]/phun("MEGA*Pa")
    sAn  =  fld.values[nA]/sigmaA
    println("Stress at point A: $(sA) i. e.  $( sAn*100  )% of reference value")
    
    fen2fe  = FENodeToFEMap(fes.conn, nnodes(geom))
    function inspector(idat, elnum, conn, xe,  out,  xq)
        println("loc = $(  xq  ) : $(  transpose(out)/phun("MEGA*Pa")  )")
        return idat
    end
    
    inspectintegpoints(femm, geom, u, dT,  fen2fe.map[nA[1]],
    inspector, []; quantity = :Cauchy)
    
end # LE11NAFEMS_Q8_stress

function allrun()
    println("#####################################################") 
    println("# LE11NAFEMS_Q8_algo ")
    LE11NAFEMS_Q8_algo()
    println("#####################################################") 
    println("# LE11NAFEMS_Q8_algo2 ")
    LE11NAFEMS_Q8_algo2()
    println("#####################################################") 
    println("# LE11NAFEMS_Q8_export_stress ")
    LE11NAFEMS_Q8_export_stress()
    println("#####################################################") 
    println("# LE11NAFEMS_Q8_stress ")
    LE11NAFEMS_Q8_stress()
    return true
end # function allrun

end # module LE11NAFEMS_examples
