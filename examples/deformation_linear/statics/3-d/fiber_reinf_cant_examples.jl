module fiber_reinf_cant_examples
using FinEtools
using FinEtools.AlgoDeforLinearModule
# using IterativeSolvers
using Statistics: mean
using LinearAlgebra: Symmetric, cholesky

function fiber_reinf_cant_iso()
    println("""
    Cantilever example.  Isotropic material.
    @article{
    author = {Krysl, P.},
    title = {Mean-strain 8-node hexahedron with optimized energy-sampling stabilization},
    journal = {Finite Elements in Analysis and Design},
    volume = {108}, pages = {41-53}, DOI = {10.1016/j.finel.2015.09.008}, year = {2016}
    }
    """)
    
    t0 = time()
    # # Orthotropic material parameters of the external cylinder
    # E1s = 130.0*phun("GPa")
    # E2s = 5.0*phun("GPa")
    # E3s = E2s
    # nu12s = nu13s = 0.25
    # nu23s = 0.0
    # G12s = 10.0*phun("GPa")
    # G13s = G12s
    # G23s = 5.0*phun("GPa")
    # CTE1 = 3.0e-6
    # CTE2 = 2.0e-5
    # CTE3 = 2.0e-5
    # Isotropic material
    E = 1.0e9*phun("Pa")
    nu = 0.25
    CTE = 0.0
    
    # Reference value for  the vertical deflection of the tip
    uz_ref =  -7.516310912734678e-06;
    
    a = 90.0*phun("mm") # length of the cantilever
    b = 10.0*phun("mm") # width of the cross-section
    t = 20.0*phun("mm") # height of the cross-section
    q0 = -1000.0*phun("Pa") # shear traction
    dT = 0*phun("K") # temperature rise
    
    tolerance = 0.00001*t
    
    # Generate mesh
    n = 16
    na = n # number of elements lengthwise
    nb = n # number of elements through the wwith
    nt = n # number of elements through the thickness
    xs = collect(linearspace(0.0, a, na+1))
    ys = collect(linearspace(0.0, b, nb+1))
    ts = collect(linearspace(0.0, t, nt+1))
    fens,fes = H8blockx(xs, ys, ts)
    fens,fes = H8toH20(fens,fes)
    bfes = meshboundary(fes)
    # end cross-section surface  for the shear loading
    sshearl = selectelem(fens, bfes; facing=true, direction = [+1.0 0.0 0.0])
    
    MR = DeforModelRed3D
    # externalmaterial = MatDeforElastOrtho(MR,
    #   0.0, E1s, E2s, E3s,
    #   nu12s, nu13s, nu23s,
    #   G12s, G13s, G23s,
    #   CTE1, CTE2, CTE3)
    material = MatDeforElastIso(MR,
    0.0, E, nu, CTE)
    
    # Material orientation matrix
    csmat = [i==j ? one(FFlt) : zero(FFlt) for i=1:3, j=1:3]
    
    function updatecs!(csmatout::FFltMat, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)
        copyto!(csmatout, csmat)
    end
    
    gr = GaussRule(3, 2)
    
    region = FDataDict("femm"=>FEMMDeforLinear(MR, IntegData(fes, gr), material))
    
    lx0 = selectnode(fens, box=[0.0 0.0 -Inf Inf -Inf Inf], inflate=tolerance)
    
    ex01 = FDataDict( "displacement"=>  0.0, "component"=> 1, "node_list"=>lx0 )
    ex02 = FDataDict( "displacement"=>  0.0, "component"=> 2, "node_list"=>lx0 )
    ex03 = FDataDict( "displacement"=>  0.0, "component"=> 3, "node_list"=>lx0 )
    
    function getshr!(forceout::FFltVec, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)
        copyto!(forceout, q0*[0.0; 0.0; 1.0])
    end
    
    Trac = FDataDict("traction_vector"=>getshr!, "femm"=>FEMMBase(IntegData(subset(bfes, sshearl), GaussRule(2, 3))))
    
    modeldata = FDataDict("fens"=>fens,
    "regions"=>[region],
    "essential_bcs"=>[ex01, ex02, ex03],
    "traction_bcs"=>[Trac],
    "temperature_change"=>FDataDict("temperature"=>dT)
    )
    modeldata = AlgoDeforLinearModule.linearstatics(modeldata)
    
    u = modeldata["u"]
    geom = modeldata["geom"]
    
    Tipl = selectnode(fens, box=[a a b b 0. 0.], inflate=tolerance)
    utip = mean(u.values[Tipl, 3])
    println(" Normalized deflection: $(utip/uz_ref)")
    println("Solution: $(  time()-t0 )")
    
    # File =  "NAFEMS-R0031-2-plate.vtk"
    # vtkexportmesh(File, fes.conn, geom.values, FinEtools.MeshExportModule.H20;
    #     scalars = [("Layer", fes.label)], vectors = [("displacement", u.values)])
    # @async run(`"paraview.exe" $File`)
    
    modeldata["postprocessing"] = FDataDict("file"=>"fiber_reinf_cant_iso",
    "outputcsys"=>CSys(3, 3, updatecs!), "quantity"=>:Cauchy, "component"=>5)
    # modeldata = AlgoDeforLinearModule.exportstress(modeldata)
    modeldata = AlgoDeforLinearModule.exportdeformation(modeldata)
    
    println("Done: $(  time()-t0 )")
    true
    
end # fiber_reinf_cant_iso


function fiber_reinf_cant_iso_stresses()
    println("""
    Cantilever example.  Isotropic material.
    @article{
    author = {Krysl, P.},
    title = {Mean-strain 8-node hexahedron with optimized energy-sampling stabilization},
    journal = {Finite Elements in Analysis and Design},
    volume = {108}, pages = {41-53}, DOI = {10.1016/j.finel.2015.09.008}, year = {2016}
    }
    """)
    
    # Isotropic material
    E = 1.0e9*phun("Pa")
    nu = 0.25
    CTE = 0.0
    
    # Reference value for  the vertical deflection of the tip
    uz_ref =  -7.516310912734678e-06;
    
    a = 90.0*phun("mm") # length of the cantilever
    b = 10.0*phun("mm") # width of the cross-section
    t = 20.0*phun("mm") # height of the cross-section
    q0 = -1000.0*phun("Pa") # shear traction
    dT = 0*phun("K") # temperature rise
    
    tolerance = 0.00001*t
    
    # Generate mesh
    for n = [2 4 8 16 32]
        
        na = n # number of elements lengthwise
        nb = n # number of elements through the wwith
        nt = n # number of elements through the thickness
        xs = collect(linearspace(0.0, a, na+1))
        ys = collect(linearspace(0.0, b, nb+1))
        ts = collect(linearspace(0.0, t, nt+1))
        fens,fes = H8blockx(xs, ys, ts)
        # fens,fes = H8toH20(fens,fes)
        bfes = meshboundary(fes)
        # end cross-section surface  for the shear loading
        sshearl = selectelem(fens, bfes; facing=true, direction = [+1.0 0.0 0.0])
        
        MR = DeforModelRed3D
        # externalmaterial = MatDeforElastOrtho(MR,
        #   0.0, E1s, E2s, E3s,
        #   nu12s, nu13s, nu23s,
        #   G12s, G13s, G23s,
        #   CTE1, CTE2, CTE3)
        material = MatDeforElastIso(MR,
        0.0, E, nu, CTE)
        
        # Material orientation matrix
        csmat = [i==j ? one(FFlt) : zero(FFlt) for i=1:3, j=1:3]
        
        function updatecs!(csmatout::FFltMat, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)
            copyto!(csmatout, csmat)
        end
        
        gr = GaussRule(3, 2)
        
        region = FDataDict("femm"=>FEMMDeforLinearMSH8(MR,
        IntegData(fes, gr), material))
        
        lx0 = selectnode(fens, box=[0.0 0.0 -Inf Inf -Inf Inf], inflate=tolerance)
        
        ex01 = FDataDict( "displacement"=>  0.0, "component"=> 1, "node_list"=>lx0 )
        ex02 = FDataDict( "displacement"=>  0.0, "component"=> 2, "node_list"=>lx0 )
        ex03 = FDataDict( "displacement"=>  0.0, "component"=> 3, "node_list"=>lx0 )
        
        function getshr!(forceout::FFltVec, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)
            copyto!(forceout, q0*[0.0; 0.0; 1.0])
        end
        
        Trac = FDataDict("traction_vector"=>getshr!,
        "femm"=>FEMMBase(IntegData(subset(bfes, sshearl), GaussRule(2, 3))))
        
        modeldata = FDataDict("fens"=>fens,
        "regions"=>[region],
        "essential_bcs"=>[ex01, ex02, ex03],
        "traction_bcs"=>[Trac],
        "temperature_change"=>FDataDict("temperature"=>dT)
        )
        modeldata = AlgoDeforLinearModule.linearstatics(modeldata)
        
        u = modeldata["u"]
        geom = modeldata["geom"]
        
        Tipl = selectnode(fens, box=[a a b b 0. 0.], inflate=tolerance)
        utip = mean(u.values[Tipl, 3])
        println(" Normalized deflection: $(utip/uz_ref)")
        
        stressfields = NodalField[]
        for c = 1:6
            modeldata["postprocessing"] = FDataDict("file"=>"fiber_reinf_cant_iso",
            "outputcsys"=>CSys(3, 3, updatecs!), "quantity"=>:Cauchy, "component"=>c)
            modeldata = AlgoDeforLinearModule.exportstress(modeldata)
            push!(stressfields, modeldata["postprocessing"]["exported"][1]["field"])
        end
        
        # println("@__FILE__ = $(@__FILE__)")
        # jldopen("fiber_reinf_cant_iso_stresses" * "n=$(n)" * ".jld", "w") do file
        #     write(file, "n", n)
        #     write(file, "fens", fens)
        #     write(file, "fes", fes)
        #     write(file, "geom", geom)
        #     write(file, "u", u)
        #     # write(file, "femm", region["femm"])
        #     write(file, "integrationrule", region["femm"].integdata.integration_rule)
        #     write(file, "stressfields", stressfields)
        #     write(file, "tolerance", tolerance)
        # end
        
    end
    true
    
end # fiber_reinf_cant_iso_stresses


function fiber_reinf_cant_iso_stresses_MST10()
    elementtag = "MST10"
    println("""
    Cantilever example.  Isotropic material.
    ############### $(elementtag) ###############
    @article{
    author = {Krysl, P.},
    title = {Mean-strain 8-node hexahedron with optimized energy-sampling stabilization},
    journal = {Finite Elements in Analysis and Design},
    volume = {108}, pages = {41-53}, DOI = {10.1016/j.finel.2015.09.008}, year = {2016}
    }
    """)
    
    # Isotropic material
    E = 1.0e9*phun("Pa")
    nu = 0.25
    CTE = 0.0
    
    # Reference value for  the vertical deflection of the tip
    uz_ref =  -7.516310912734678e-06;
    
    a = 90.0*phun("mm") # length of the cantilever
    b = 10.0*phun("mm") # width of the cross-section
    t = 20.0*phun("mm") # height of the cross-section
    q0 = -1000.0*phun("Pa") # shear traction
    dT = 0*phun("K") # temperature rise
    
    tolerance = 0.00001*t
    
    convergencestudy = FDataDict[]
    for n = [1 2 4 8 16]
        
        na = n # number of elements lengthwise
        nb = n # number of elements through the wwith
        nt = n # number of elements through the thickness
        xs = collect(linearspace(0.0, a, na+1))
        ys = collect(linearspace(0.0, b, nb+1))
        ts = collect(linearspace(0.0, t, nt+1))
        fens,fes = T10blockx(xs, ys, ts)
        bfes = meshboundary(fes)
        # end cross-section surface  for the shear loading
        sshearl = selectelem(fens, bfes; facing=true, direction = [+1.0 0.0 0.0])
        
        MR = DeforModelRed3D
        # externalmaterial = MatDeforElastOrtho(MR,
        #   0.0, E1s, E2s, E3s,
        #   nu12s, nu13s, nu23s,
        #   G12s, G13s, G23s,
        #   CTE1, CTE2, CTE3)
        material = MatDeforElastIso(MR,
        0.0, E, nu, CTE)
        
        # Material orientation matrix
        csmat = [i==j ? one(FFlt) : zero(FFlt) for i=1:3, j=1:3]
        
        function updatecs!(csmatout::FFltMat, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)
            copyto!(csmatout, csmat)
        end
        
        gr = SimplexRule(3, 4)
        
        region = FDataDict("femm"=>FEMMDeforLinearMST10(MR,
        IntegData(fes, gr), material))
        
        lx0 = selectnode(fens, box=[0.0 0.0 -Inf Inf -Inf Inf], inflate=tolerance)
        
        ex01 = FDataDict( "displacement"=>  0.0, "component"=> 1, "node_list"=>lx0 )
        ex02 = FDataDict( "displacement"=>  0.0, "component"=> 2, "node_list"=>lx0 )
        ex03 = FDataDict( "displacement"=>  0.0, "component"=> 3, "node_list"=>lx0 )
        
        function getshr!(forceout::FFltVec, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)
            copyto!(forceout, q0*[0.0; 0.0; 1.0])
        end
        
        Trac = FDataDict("traction_vector"=>getshr!,
        "femm"=>FEMMBase(IntegData(subset(bfes, sshearl), SimplexRule(2, 3))))
        
        modeldata = FDataDict("fens"=>fens,
        "regions"=>[region],
        "essential_bcs"=>[ex01, ex02, ex03],
        "traction_bcs"=>[Trac],
        "temperature_change"=>FDataDict("temperature"=>dT)
        )
        modeldata = AlgoDeforLinearModule.linearstatics(modeldata)
        
        u = modeldata["u"]
        geom = modeldata["geom"]
        
        Tipl = selectnode(fens, box=[a a b b 0. 0.], inflate=tolerance)
        utip = mean(u.values[Tipl, 3])
        println(" Normalized deflection: $(utip/uz_ref)")
        
        modeldata["postprocessing"] = FDataDict("file"=>"fiber_reinf_cant_iso_stresses_$(elementtag)",
        "outputcsys"=>CSys(3, 3, updatecs!), "quantity"=>:Cauchy,
        "component"=>[5])
        modeldata = AlgoDeforLinearModule.exportstresselementwise(modeldata)
        
        modeldata["postprocessing"] = FDataDict("file"=>"fiber_reinf_cant_iso_stresses_$(elementtag)",
        "outputcsys"=>CSys(3, 3, updatecs!), "quantity"=>:Cauchy,
        "component"=>collect(1:6))
        modeldata = AlgoDeforLinearModule.exportstresselementwise(modeldata)
        stressfields = ElementalField[modeldata["postprocessing"]["exported"][1]["field"]]
        
        push!(convergencestudy, FDataDict(
        "elementsize"=> 1.0 / n,
        "fens"=>fens,
        "fes"=>fes,
        "geom"=>geom,
        "u"=>u,
        "femm"=>region["femm"],
        "integrationrule"=>region["femm"].integdata.integration_rule,
        "stressfields"=>stressfields,
        "tolerance"=>tolerance)
        )
    end
    
    # File = "fiber_reinf_cant_iso_stresses_$(elementtag)"
    # open(File * ".jls", "w") do file
    #     serialize(file, convergencestudy)
    # end
    
    true
    
end # fiber_reinf_cant_iso_stresses_MST10


function fiber_reinf_cant_iso_stresses_T10()
    elementtag = "T10"
    println("""
    Cantilever example.  Isotropic material.
    ############### $(elementtag) ###############
    @article{
    author = {Krysl, P.},
    title = {Mean-strain 8-node hexahedron with optimized energy-sampling stabilization},
    journal = {Finite Elements in Analysis and Design},
    volume = {108}, pages = {41-53}, DOI = {10.1016/j.finel.2015.09.008}, year = {2016}
    }
    """)
    
    # Isotropic material
    E = 1.0e9*phun("Pa")
    nu = 0.25
    CTE = 0.0
    
    # Reference value for  the vertical deflection of the tip
    uz_ref =  -7.516310912734678e-06;
    
    a = 90.0*phun("mm") # length of the cantilever
    b = 10.0*phun("mm") # width of the cross-section
    t = 20.0*phun("mm") # height of the cross-section
    q0 = -1000.0*phun("Pa") # shear traction
    dT = 0*phun("K") # temperature rise
    
    tolerance = 0.00001*t
    
    convergencestudy = FDataDict[]
    for n = [1 2 4 8 16]
        
        na = n # number of elements lengthwise
        nb = n # number of elements through the wwith
        nt = n # number of elements through the thickness
        xs = collect(linearspace(0.0, a, na+1))
        ys = collect(linearspace(0.0, b, nb+1))
        ts = collect(linearspace(0.0, t, nt+1))
        fens,fes = T10blockx(xs, ys, ts)
        bfes = meshboundary(fes)
        # end cross-section surface  for the shear loading
        sshearl = selectelem(fens, bfes; facing=true, direction = [+1.0 0.0 0.0])
        
        MR = DeforModelRed3D
        # externalmaterial = MatDeforElastOrtho(MR,
        #   0.0, E1s, E2s, E3s,
        #   nu12s, nu13s, nu23s,
        #   G12s, G13s, G23s,
        #   CTE1, CTE2, CTE3)
        material = MatDeforElastIso(MR,
        0.0, E, nu, CTE)
        
        # Material orientation matrix
        csmat = [i==j ? one(FFlt) : zero(FFlt) for i=1:3, j=1:3]
        
        function updatecs!(csmatout::FFltMat, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)
            copyto!(csmatout, csmat)
        end
        
        gr = SimplexRule(3, 4)
        
        region = FDataDict("femm"=>FEMMDeforLinear(MR,
        IntegData(fes, gr), material))
        
        lx0 = selectnode(fens, box=[0.0 0.0 -Inf Inf -Inf Inf], inflate=tolerance)
        
        ex01 = FDataDict( "displacement"=>  0.0, "component"=> 1, "node_list"=>lx0 )
        ex02 = FDataDict( "displacement"=>  0.0, "component"=> 2, "node_list"=>lx0 )
        ex03 = FDataDict( "displacement"=>  0.0, "component"=> 3, "node_list"=>lx0 )
        
        function getshr!(forceout::FFltVec, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)
            copyto!(forceout, q0*[0.0; 0.0; 1.0])
        end
        
        Trac = FDataDict("traction_vector"=>getshr!,
        "femm"=>FEMMBase(IntegData(subset(bfes, sshearl), SimplexRule(2, 3))))
        
        modeldata = FDataDict("fens"=>fens,
        "regions"=>[region],
        "essential_bcs"=>[ex01, ex02, ex03],
        "traction_bcs"=>[Trac],
        "temperature_change"=>FDataDict("temperature"=>dT)
        )
        modeldata = AlgoDeforLinearModule.linearstatics(modeldata)
        
        u = modeldata["u"]
        geom = modeldata["geom"]
        
        Tipl = selectnode(fens, box=[a a b b 0. 0.], inflate=tolerance)
        utip = mean(u.values[Tipl, 3])
        println(" Normalized deflection: $(utip/uz_ref)")
        
        modeldata["postprocessing"] = FDataDict("file"=>"fiber_reinf_cant_iso_stresses_$(elementtag)",
        "outputcsys"=>CSys(3, 3, updatecs!), "quantity"=>:Cauchy,
        "component"=>[5])
        modeldata = AlgoDeforLinearModule.exportstresselementwise(modeldata)
        
        modeldata["postprocessing"] = FDataDict("file"=>"fiber_reinf_cant_iso_stresses_$(elementtag)",
        "outputcsys"=>CSys(3, 3, updatecs!), "quantity"=>:Cauchy,
        "component"=>collect(1:6))
        modeldata = AlgoDeforLinearModule.exportstresselementwise(modeldata)
        stressfields = ElementalField[modeldata["postprocessing"]["exported"][1]["field"]]
        
        push!(convergencestudy, FDataDict(
        "elementsize"=> 1.0 / n,
        "fens"=>fens,
        "fes"=>fes,
        "geom"=>geom,
        "u"=>u,
        "femm"=>region["femm"],
        "integrationrule"=>region["femm"].integdata.integration_rule,
        "stressfields"=>stressfields,
        "tolerance"=>tolerance)
        )
    end
    
    # File = "fiber_reinf_cant_iso_stresses_$(elementtag)"
    # open(File * ".jls", "w") do file
    #     serialize(file, convergencestudy)
    # end
    
    true
    
end # fiber_reinf_cant_iso_stresses_T10


function fiber_reinf_cant_yn_strong()
    println("""
    Cantilever example.  Strongly orthotropic material. Orientation "y".
    @article{
    author = {Krysl, P.},
    title = {Mean-strain 8-node hexahedron with optimized energy-sampling stabilization},
    journal = {Finite Elements in Analysis and Design},
    volume = {108}, pages = {41-53}, DOI = {10.1016/j.finel.2015.09.008}, year = {2016}
    }
    """)
    
    t0 = time()
    # # Orthotropic material
    E1s = 100000.0*phun("GPa")
    E2s = 1.0*phun("GPa")
    E3s = E2s
    nu23s = nu12s = nu13s = 0.25
    G12s = 0.2*phun("GPa")
    G23s = G13s = G12s
    CTE1 = 0.0
    CTE2 = 0.0
    CTE3 = 0.0
    # # Isotropic material
    # E = 1.0e9*phun("Pa")
    # nu = 0.25
    # CTE = 0.0
    
    # Reference value for  the vertical deflection of the tip
    uz_ref = -1.027498445054843e-05;
    
    a = 90.0*phun("mm") # length of the cantilever
    b = 10.0*phun("mm") # width of the cross-section
    t = 20.0*phun("mm") # height of the cross-section
    q0 = -1000.0*phun("Pa") # shear traction
    dT = 0*phun("K") # temperature rise
    
    tolerance = 0.00001*t
    
    # Generate mesh
    n = 4
    na = 8*n # number of elements lengthwise
    nb = n # number of elements through the wwith
    nt = n # number of elements through the thickness
    xs = collect(linearspace(0.0, a, na+1))
    ys = collect(linearspace(0.0, b, nb+1))
    ts = collect(linearspace(0.0, t, nt+1))
    fens,fes = H8blockx(xs, ys, ts)
    fens,fes = H8toH20(fens,fes)
    bfes = meshboundary(fes)
    # end cross-section surface  for the shear loading
    sshearl = selectelem(fens, bfes; facing=true, direction = [+1.0 0.0 0.0])
    
    MR = DeforModelRed3D
    material = MatDeforElastOrtho(MR,
    0.0, E1s, E2s, E3s,
    nu12s, nu13s, nu23s,
    G12s, G13s, G23s,
    CTE1, CTE2, CTE3)
    # material = MatDeforElastIso(MR,
    #   0.0, E, nu, CTE)
    
    # Material orientation matrix
    csmat = zeros(3, 3)
    rotmat3!(csmat, -45.0/180.0*pi*[0,1,0])
    
    function updatecs!(csmatout::FFltMat, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)
        copyto!(csmatout, csmat)
    end
    
    gr = GaussRule(3, 2)
    
    region = FDataDict("femm"=>FEMMDeforLinear(MR, IntegData(fes, gr), CSys(3, 3, updatecs!), material))
    
    lx0 = selectnode(fens, box=[0.0 0.0 -Inf Inf -Inf Inf], inflate=tolerance)
    
    ex01 = FDataDict( "displacement"=>  0.0, "component"=> 1, "node_list"=>lx0 )
    ex02 = FDataDict( "displacement"=>  0.0, "component"=> 2, "node_list"=>lx0 )
    ex03 = FDataDict( "displacement"=>  0.0, "component"=> 3, "node_list"=>lx0 )
    
    function getshr!(forceout::FFltVec, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)
        copyto!(forceout, q0*[0.0; 0.0; 1.0])
    end
    
    Trac = FDataDict("traction_vector"=>getshr!, "femm"=>FEMMBase(IntegData(subset(bfes, sshearl), GaussRule(2, 3))))
    
    modeldata = FDataDict("fens"=>fens,
    "regions"=>[region],
    "essential_bcs"=>[ex01, ex02, ex03],
    "traction_bcs"=>[Trac],
    "temperature_change"=>FDataDict("temperature"=>dT)
    )
    modeldata = AlgoDeforLinearModule.linearstatics(modeldata)
    
    u = modeldata["u"]
    geom = modeldata["geom"]
    
    Tipl = selectnode(fens, box=[a a b b 0. 0.], inflate=tolerance)
    utip = mean(u.values[Tipl, 3])
    println("Deflection $utip, normalized: $(utip/uz_ref)")
    println("Solution: $(  time()-t0 )")
    
    # File =  "NAFEMS-R0031-2-plate.vtk"
    # vtkexportmesh(File, fes.conn, geom.values, FinEtools.MeshExportModule.H20;
    #     scalars = [("Layer", fes.label)], vectors = [("displacement", u.values)])
    # @async run(`"paraview.exe" $File`)
    
    modeldata["postprocessing"] = FDataDict("file"=>"fiber_reinf_cant_yn_strong", "outputcsys"=>CSys(3, 3, updatecs!), "quantity"=>:Cauchy, "component"=>5)
    modeldata = AlgoDeforLinearModule.exportstress(modeldata)
    File = modeldata["postprocessing"]["exported"][1]["file"]
    @async run(`"paraview.exe" $File`)
    
    # modeldata = AlgoDeforLinearModule.exportdeformation(modeldata)
    # File = modeldata["postprocessing"]["exported"][1]["file"]
    # @async run(`"paraview.exe" $File`)
    
    println("Done: $(  time()-t0 )")
    true
    
end # fiber_reinf_cant_yn_strong


function fiber_reinf_cant_yn_strong_no_algo()
    println("""
    Cantilever example.  Strongly orthotropic material. Orientation "y".
    @article{
    author = {Krysl, P.},
    title = {Mean-strain 8-node hexahedron with optimized energy-sampling stabilization},
    journal = {Finite Elements in Analysis and Design},
    volume = {108}, pages = {41-53}, DOI = {10.1016/j.finel.2015.09.008}, year = {2016}
    }
    """)
    
    t0 = time()
    pu = ustring -> phun(ustring; system_of_units = :SIMM)
    
    # # Orthotropic material
    E1s = 100000.0*pu("GPa")
    E2s = 1.0*pu("GPa")
    E3s = E2s
    nu23s = nu12s = nu13s = 0.25
    G12s = 0.2*pu("GPa")
    G23s = G13s = G12s
    CTE1 = 0.0
    CTE2 = 0.0
    CTE3 = 0.0
    # # Isotropic material
    # E = 1.0e9*pu("Pa")
    # nu = 0.25
    # CTE = 0.0
    
    # Reference value for  the vertical deflection of the tip
    uz_ref = -1.027498445054843e-05*pu("m");
    
    a = 90.0*pu("mm") # length of the cantilever
    b = 10.0*pu("mm") # width of the cross-section
    t = 20.0*pu("mm") # height of the cross-section
    q0 = -1000.0*pu("Pa") # shear traction
    dT = 0*pu("K") # temperature rise
    
    tolerance = 0.00001*t
    
    # Generate mesh
    n = 10
    na = n # number of elements lengthwise
    nb = n # number of elements through the wwith
    nt = n # number of elements through the thickness
    xs = collect(linearspace(0.0, a, na+1))
    ys = collect(linearspace(0.0, b, nb+1))
    ts = collect(linearspace(0.0, t, nt+1))
    println("fens,fes = H8blockx(xs, ys, ts)")
    @time fens,fes = H8blockx(xs, ys, ts)
    println("fens,fes = H8toH20(fens,fes)")
    @time fens,fes = H8toH20(fens,fes)
    println("bfes = meshboundary(fes)")
    @time bfes = meshboundary(fes)
    # end cross-section surface  for the shear loading
    sshearl = selectelem(fens, bfes; facing=true, direction = [+1.0 0.0 0.0])
    
    MR = DeforModelRed3D
    material = MatDeforElastOrtho(MR,
    0.0, E1s, E2s, E3s,
    nu12s, nu13s, nu23s,
    G12s, G13s, G23s,
    CTE1, CTE2, CTE3)
    # material = MatDeforElastIso(MR,
    #   0.0, E, nu, CTE)
    
    # Material orientation matrix
    csmat = zeros(3, 3)
    rotmat3!(csmat, -45.0/180.0*pi*[0,1,0])
    
    function updatecs!(csmatout::FFltMat, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)
        copyto!(csmatout, csmat)
    end
    
    gr = GaussRule(3, 2)
    
    femm = FEMMDeforLinear(MR, IntegData(fes, gr), CSys(3, 3, updatecs!), material)

    lx0 = selectnode(fens, box=[0.0 0.0 -Inf Inf -Inf Inf], inflate=tolerance)
    
    geom = NodalField(fens.xyz)
    u = NodalField(zeros(size(fens.xyz,1),3)) # displacement field
    nnodes(geom)
    
    setebc!(u, lx0, true, 1, zeros(size(lx0)))
    setebc!(u, lx0, true, 2, zeros(size(lx0)))
    setebc!(u, lx0, true, 3, zeros(size(lx0)))
    applyebc!(u)
    
    S = connectionmatrix(femm, nnodes(geom))
    
    numberdofs!(u)
    
    function getshr!(forceout::FFltVec, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)
        copyto!(forceout, q0*[0.0; 0.0; 1.0])
    end
    
    Tracfemm = FEMMBase(IntegData(subset(bfes, sshearl), GaussRule(2, 3)))
    
    println("K = stiffness(femm, geom, u)")
    @time K = stiffness(femm, geom, u)
    fi = ForceIntensity(FFlt, 3, getshr!);
    println("F =  distribloads(Tracfemm, geom, u, fi, 2);")
    @time F =  distribloads(Tracfemm, geom, u, fi, 2);
    
    println("K = cholesky(K)")
    K = (K + K')/2;
    @time K = cholesky(Symmetric(K))
    println("U = K\\F")
    @time U = K\F
    # println("U = cg(K, F; tol=1e-3, maxiter=2000)")
    # @time U = cg(K, F; tol=1e-3, maxiter=2000)
    scattersysvec!(u, U[:])
    
    Tipl = selectnode(fens, box=[a a b b 0. 0.], inflate=tolerance)
    utip = mean(u.values[Tipl, 3])
    println("Deflection $utip, normalized: $(utip/uz_ref)")
    println("Solution: $(  time()-t0 )")
    
    println("Done: $(  time()-t0 )")
    true
    
end # fiber_reinf_cant_yn_strong_no_algo

function allrun()
    println("#####################################################") 
    println("# fiber_reinf_cant_iso ")
    fiber_reinf_cant_iso()
    println("#####################################################") 
    println("# fiber_reinf_cant_iso_stresses ")
    fiber_reinf_cant_iso_stresses()
    println("#####################################################") 
    println("# fiber_reinf_cant_iso_stresses_MST10 ")
    fiber_reinf_cant_iso_stresses_MST10()
    println("#####################################################") 
    println("# fiber_reinf_cant_iso_stresses_T10 ")
    fiber_reinf_cant_iso_stresses_T10()
    println("#####################################################") 
    println("# fiber_reinf_cant_yn_strong ")
    fiber_reinf_cant_yn_strong()
    println("#####################################################") 
    println("# fiber_reinf_cant_yn_strong_no_algo ")
    fiber_reinf_cant_yn_strong_no_algo()
    return true
end # function allrun

end # module fiber_reinf_cant_examples
