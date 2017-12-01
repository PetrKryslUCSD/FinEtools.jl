module twisted_beam_examples
using FinEtools
using FinEtools.AlgoDeforLinearModule
using FinEtools.MeshExportModule

function twisted_beam_algo()
    println("""    The initially twisted cantilever beam is one of the standard test problems for verifying the finite-element accuracy [1]. The beam is clamped at one end and loaded either with unit in-plane or
    unit out-of-plane force at the other. The centroidal axis of the beam is
    straight at the undeformed  configuration, while its cross-sections are
    twisted about the centroidal axis from 0 at the clamped end to pi/2 at
    the free end.
    
    Reference:
    Zupan D, Saje M (2004) On "A proposed standard set of problems to test
    finite element accuracy": the twisted beam. Finite Elements in Analysis
    and Design 40: 1445-1451.
    """)
    E = 0.29e8;
    nu = 0.22;
    W = 1.1;
    L = 12.;
    t =  0.32;
    nl = 2; nt = 1; nw = 1; ref = 7;
    p =   1/W/t;
    #  Loading in the Z direction
    loadv = [0;0;p]; dir = 3; uex = 0.005424534868469; # Harder: 5.424e-3;
    #   Loading in the Y direction
    #loadv = [0;p;0]; dir = 2; uex = 0.001753248285256; # Harder: 1.754e-3;
    tolerance  = t/1000;
    
    fens,fes  = H20block(L,W,t, nl*ref,nw*ref,nt*ref)
    
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
    el1femm  = FEMMBase(IntegData(subset(boundaryfes,Toplist), GaussRule(2, 2)))
    flux1 = FDataDict("femm"=>el1femm, "traction_vector"=>loadv)
    
    
    # Make the region
    MR = DeforModelRed3D
    material = MatDeforElastIso(MR, 00.0, E, nu, 0.0)
    region1 = FDataDict("femm"=>FEMMDeforLinear(MR, IntegData(fes, GaussRule(3,2)),
    material))
    
    # Make model data
    modeldata =  FDataDict(
    "fens"=> fens, "regions"=>  [region1],
    "essential_bcs"=>[e1, e2, e3], "traction_bcs"=>  [flux1])
    
    # Call the solver
    modeldata = AlgoDeforLinearModule.linearstatics(modeldata)
    geom = modeldata["geom"]
    u = modeldata["u"]
    
    # Extract the solution
    nl = selectnode(fens, box = [L L -100*W 100*W -100*W 100*W],inflate = tolerance);
    theutip = mean(u.values[nl,:],1)
    println("displacement  = $(theutip[dir]) as compared to converged $uex")
    println("normalized displacement  = $(theutip[dir]/uex*100) %")
    
    # Write out mesh with displacements
    modeldata["postprocessing"] = FDataDict("file"=>"twisted_beam")
    modeldata = AlgoDeforLinearModule.exportdeformation(modeldata)
    
    # Write out mesh with stresses
    modeldata["postprocessing"] = FDataDict("file"=>"twisted_beam",
    "quantity"=> :Cauchy, "component"=> :xy)
    modeldata = AlgoDeforLinearModule.exportstress(modeldata)
    
    # Write out mesh with stresses
    modeldata["postprocessing"] = FDataDict("file"=>"twisted_beam",
    "quantity"=> :Cauchy, "component"=> :xz)
    modeldata = AlgoDeforLinearModule.exportstress(modeldata)
    
    # Write out mesh with von Mises stresses
    modeldata["postprocessing"] = FDataDict("file"=>"twisted_beam",
    "quantity"=> :vm)
    modeldata = AlgoDeforLinearModule.exportstress(modeldata)
    
    # Write out mesh with von Mises stresses, elementwise
    modeldata["postprocessing"] = FDataDict("file"=>"twisted_beam-ew",
    "quantity"=> :vm)
    modeldata = AlgoDeforLinearModule.exportstresselementwise(modeldata)
    vm  = modeldata["postprocessing"]["exported"][1]["field"]
    println("extremes of von Mises: $([minimum(vm.values),   maximum(vm.values)])")
    
    # Write out mesh with von Mises stresses, elementwise
    modeldata["postprocessing"] = FDataDict("file"=>"twisted_beam-ew",
    "quantity"=> :Cauchy, "component"=> :xz)
    modeldata = AlgoDeforLinearModule.exportstresselementwise(modeldata)
    
    println("Done")
    true
end # twisted_beam_algo


function twisted_beam_algo_stress()
    println("""The initially twisted cantilever beam is one of the standard test problems for verifying the finite-element accuracy [1]. The beam is clamped at one end and loaded either with unit in-plane or unit ut-of-plane force at the other. The centroidal axis of the beam is
    straight at the undeformed  configuration, while its cross-sections are
    twisted about the centroidal axis from 0 at the clamped end to pi/2 at
    the free end.
    
    Reference:
    Zupan D, Saje M (2004) On "A proposed standard set of problems to test
    finite element accuracy": the twisted beam. Finite Elements in Analysis
    and Design 40: 1445-1451.
    """)
    E = 0.29e8;
    nu = 0.22;
    W = 1.1;
    L = 12.;
    t =  0.32;
    nl = 2; nt = 1; nw = 1; ref = 4;
    p =   1/W/t;
    #  Loading in the Z direction
    loadv = [0;0;p]; dir = 3; uex = 0.005424534868469; # Harder: 5.424e-3;
    #   Loading in the Y direction
    #loadv = [0;p;0]; dir = 2; uex = 0.001753248285256; # Harder: 1.754e-3;
    tolerance  = t/1000;
    
    fens,fes  = H20block(L,W,t, nl*ref,nw*ref,nt*ref)
    
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
    el1femm  = FEMMBase(IntegData(subset(boundaryfes,Toplist), GaussRule(2, 2)))
    flux1 = FDataDict("femm"=>el1femm, "traction_vector"=>loadv)
    
    
    # Make the region
    MR = DeforModelRed3D
    material = MatDeforElastIso(MR, 00.0, E, nu, 0.0)
    region1 = FDataDict("femm"=>FEMMDeforLinear(MR, IntegData(fes, GaussRule(3,2)),
    material))
    
    # Make model data
    modeldata =  FDataDict(
    "fens"=> fens, "regions"=>  [region1],
    "essential_bcs"=>[e1, e2, e3], "traction_bcs"=>  [flux1])
    
    # Call the solver
    modeldata = AlgoDeforLinearModule.linearstatics(modeldata)
    geom = modeldata["geom"]
    u = modeldata["u"]
    
    # Extract the solution
    nl = selectnode(fens, box = [L L -100*W 100*W -100*W 100*W],inflate = tolerance);
    theutip = mean(u.values[nl,:],1)
    println("displacement  = $(theutip[dir]) as compared to converged $uex")
    println("normalized displacement  = $(theutip[dir]/uex*100) %")
    
    # Write out mesh with displacements
    modeldata["postprocessing"] = FDataDict("file"=>"twisted_beam")
    modeldata = AlgoDeforLinearModule.exportdeformation(modeldata)
    
    # Write out mesh with stresses
    modeldata["postprocessing"] = FDataDict("file"=>"twisted_beam",
    "quantity"=> :Cauchy, "component"=> :xy)
    modeldata = AlgoDeforLinearModule.exportstress(modeldata)
    
    # Write out mesh with stresses
    modeldata["postprocessing"] = FDataDict("file"=>"twisted_beam",
    "quantity"=> :Cauchy, "component"=> :xz)
    modeldata = AlgoDeforLinearModule.exportstress(modeldata)
    
    # Write out mesh with von Mises stresses
    modeldata["postprocessing"] = FDataDict("file"=>"twisted_beam",
    "quantity"=> :vm)
    modeldata = AlgoDeforLinearModule.exportstress(modeldata)
    
    # Write out mesh with von Mises stresses, elementwise
    modeldata["postprocessing"] = FDataDict("file"=>"twisted_beam-ew",
    "quantity"=> :vm)
    modeldata = AlgoDeforLinearModule.exportstresselementwise(modeldata)
    vm  = modeldata["postprocessing"]["exported"][1]["field"]
    println("extremes of von Mises: $([minimum(vm.values),   maximum(vm.values)])")
    
    # Write out mesh with von Mises stresses, elementwise
    modeldata["postprocessing"] = FDataDict("file"=>"twisted_beam-ew",
    "quantity"=> :Cauchy, "component"=> :xz)
    modeldata = AlgoDeforLinearModule.exportstresselementwise(modeldata)
    
    # Write out mesh with principal stresses, elementwise
    modeldata["postprocessing"] = FDataDict("file"=>"twisted_beam-principal-1-ew",
    "quantity"=> :princCauchy, "component"=> 1)
    modeldata = AlgoDeforLinearModule.exportstresselementwise(modeldata)
    ps  = modeldata["postprocessing"]["exported"][1]["field"]
    println("extremes of first principal stress: $([minimum(ps.values),   maximum(ps.values)])")
    
    # Write out mesh with principal stresses, elementwise
    modeldata["postprocessing"] = FDataDict("file"=>"twisted_beam-principal-3-ew",
    "quantity"=> :princCauchy, "component"=> 3)
    modeldata = AlgoDeforLinearModule.exportstresselementwise(modeldata)
    ps  = modeldata["postprocessing"]["exported"][1]["field"]
    println("extremes of third principal stress: $([minimum(ps.values),   maximum(ps.values)])")
    
    # Write out mesh with principal stresses, elementwise
    modeldata["postprocessing"] = FDataDict("file"=>"twisted_beam-press-ew",
    "quantity"=> :pressure, "component"=> 1)
    modeldata = AlgoDeforLinearModule.exportstresselementwise(modeldata)
    ps  = modeldata["postprocessing"]["exported"][1]["field"]
    println("extremes of pressure: $([minimum(ps.values),   maximum(ps.values)])")
    
    println("Done")
    true
end # twisted_beam_algo_stress


function twisted_beam_export()
    println("""
    Refer to twisted_beam.jl.
    
    This example EXPORTS the model to Abaqus. Import the .inp file
    into Abaqus and run the job.
    """)
    E = 0.29e8;
    nu = 0.22;
    W = 1.1;
    L = 12.;
    t =  0.32;
    nl = 2; nt = 1; nw = 1; ref = 5;
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
    el1femm  = FEMMBase(IntegData(subset(boundaryfes,Toplist), GaussRule(2, 2)))
    flux1 = FDataDict("femm"=>el1femm, "traction_vector"=>loadv)
    
    
    # Make the region
    MR = DeforModelRed3D
    material = MatDeforElastIso(MR, 00.0, E, nu, 0.0)
    region1 = FDataDict("femm"=>FEMMDeforLinearMSH8(MR, IntegData(fes, GaussRule(3,2)),
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
    ELEMENT(AE, "c3d8rh", "AllElements", 1, region1["femm"].integdata.fes.conn)
    ELEMENT(AE, "SFM3D4", "TractionElements",
    1+count(region1["femm"].integdata.fes), flux1["femm"].integdata.fes.conn)
    NSET_NSET(AE, "l1", l1)
    ORIENTATION(AE, "GlobalOrientation", vec([1. 0 0]), vec([0 1. 0]));
    SOLID_SECTION(AE, "elasticity", "GlobalOrientation", "AllElements", "Hourglassctl");
    SURFACE_SECTION(AE, "TractionElements")
    END_INSTANCE(AE);
    END_ASSEMBLY(AE);
    MATERIAL(AE, "elasticity")
    ELASTIC(AE, E, nu)
    SECTION_CONTROLS(AE, "Hourglassctl", "HOURGLASS=ENHANCED")
    STEP_PERTURBATION_STATIC(AE)
    BOUNDARY(AE, "ASSEM1.INSTNC1.l1", 1)
    BOUNDARY(AE, "ASSEM1.INSTNC1.l1", 2)
    BOUNDARY(AE, "ASSEM1.INSTNC1.l1", 3)
    DLOAD(AE, "ASSEM1.INSTNC1.TractionElements", vec(flux1["traction_vector"]))
    END_STEP(AE)
    close(AE)
    
    true
end # twisted_beam_export


function twisted_beam_export_nb()
    E = 0.29e8;
    nu = 0.22;
    W = 1.1;
    L = 12.;
    t =  0.32;
    nl = 2; nt = 1; nw = 1; ref = 4;
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
    el1femm  = FEMMBase(IntegData(subset(boundaryfes,Toplist), GaussRule(2, 2)))
    flux1 = FDataDict("femm"=>el1femm, "traction_vector"=>loadv)
    
    
    # Make the region
    MR = DeforModelRed3D
    material = MatDeforElastIso(MR, 00.0, E, nu, 0.0)
    region1 = FDataDict("femm"=>FEMMDeforLinearMSH8(MR, IntegData(fes, GaussRule(3,2)),
    material))
    
    # Make model data
    modeldata =  FDataDict(
    "fens"=> fens, "regions"=>  [region1],
    "essential_bcs"=>[e1, e2, e3], "traction_bcs"=>  [flux1])
    
    
    AE = AbaqusExporter("twisted_beam");
    # AE.ios = STDOUT;
    HEADING(AE, "Twisted beam example");
    PART(AE, "part1");
    END_PART(AE);
    ASSEMBLY(AE, "ASSEM1");
    INSTANCE(AE, "INSTNC1", "PART1");
    NODE(AE, fens.xyz);
    ELEMENT(AE, "c3d8rh", "AllElements", 1, region1["femm"].integdata.fes.conn)
    ELEMENT(AE, "SFM3D4", "TractionElements",
    1+count(region1["femm"].integdata.fes), flux1["femm"].integdata.fes.conn)
    NSET_NSET(AE, "l1", l1)
    ORIENTATION(AE, "GlobalOrientation", vec([1. 0 0]), vec([0 1. 0]));
    SOLID_SECTION(AE, "elasticity", "GlobalOrientation", "AllElements", "Hourglassctl");
    SURFACE_SECTION(AE, "TractionElements")
    END_INSTANCE(AE);
    END_ASSEMBLY(AE);
    MATERIAL(AE, "elasticity")
    ELASTIC(AE, E, nu)
    SECTION_CONTROLS(AE, "Hourglassctl", "HOURGLASS=ENHANCED")
    STEP_PERTURBATION_STATIC(AE)
    BOUNDARY(AE, "ASSEM1.INSTNC1.l1", 1)
    BOUNDARY(AE, "ASSEM1.INSTNC1.l1", 2)
    BOUNDARY(AE, "ASSEM1.INSTNC1.l1", 3)
    DLOAD(AE, "ASSEM1.INSTNC1.TractionElements", vec(flux1["traction_vector"]))
    END_STEP(AE)
    close(AE)
    
end # twisted_beam_export_nb


function twisted_beam_msh8()
    println("""
    The initially twisted cantilever beam is one of the standard test    problems for verifying the finite-element accuracy [1]. The beam is        clamped at one end and loaded either with unit in-plane or        unit out-of-plane force at the other. The centroidal axis of the beam is        straight at the undeformed  configuration, while its cross-sections are        twisted about the centroidal axis from 0 at the clamped end to pi/2 at        the free end.
    
    Reference:
    Zupan D, Saje M (2004) On "A proposed standard set of problems to test
    finite element accuracy": the twisted beam. Finite Elements in Analysis
    and Design 40: 1445-1451.
    """)
    E = 0.29e8;
    nu = 0.22;
    W = 1.1;
    L = 12.;
    t =  0.32;
    nl = 2; nt = 1; nw = 1; ref = 7;
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
    el1femm  = FEMMBase(IntegData(subset(boundaryfes,Toplist), GaussRule(2, 2)))
    flux1 = FDataDict("femm"=>el1femm, "traction_vector"=>loadv)
    
    
    # Make the region
    MR = DeforModelRed3D
    material = MatDeforElastIso(MR, 00.0, E, nu, 0.0)
    region1 = FDataDict("femm"=>FEMMDeforLinearMSH8(MR, IntegData(fes, GaussRule(3,2)),
    material))
    
    # Make model data
    modeldata =  FDataDict(
    "fens"=> fens, "regions"=>  [region1],
    "essential_bcs"=>[e1, e2, e3], "traction_bcs"=>  [flux1])
    
    # Call the solver
    modeldata = AlgoDeforLinearModule.linearstatics(modeldata)
    geom = modeldata["geom"]
    u = modeldata["u"]
    
    # Extract the solution
    nl = selectnode(fens, box = [L L -100*W 100*W -100*W 100*W],inflate = tolerance);
    theutip = mean(u.values[nl,:],1)
    println("displacement  = $(theutip[dir]) as compared to converged $uex")
    
    # Write out mesh with displacements
    modeldata["postprocessing"] = FDataDict("file"=>"twisted_beam_msh8")
    modeldata = AlgoDeforLinearModule.exportdeformation(modeldata)
    
    # Write out mesh with stresses
    modeldata["postprocessing"] = FDataDict("file"=>"twisted_beam_msh8",
    "quantity"=> :Cauchy, "component"=> :xy)
    modeldata = AlgoDeforLinearModule.exportstress(modeldata)
    
    # Write out mesh with stresses
    modeldata["postprocessing"] = FDataDict("file"=>"twisted_beam_msh8",
    "quantity"=> :Cauchy, "component"=> :xz)
    modeldata = AlgoDeforLinearModule.exportstress(modeldata)
    
    # Write out mesh with von Mises stresses
    modeldata["postprocessing"] = FDataDict("file"=>"twisted_beam_msh8",
    "quantity"=> :vm)
    modeldata = AlgoDeforLinearModule.exportstress(modeldata)
    vm  = modeldata["postprocessing"]["exported"][1]["field"]
    println("extremes of vm, nodal: $([minimum(vm.values),   maximum(vm.values)])")
    
    # Write out mesh with von Mises stresses, elementwise
    modeldata["postprocessing"] = FDataDict("file"=>"twisted_beam_msh8-ew",
    "quantity"=> :vm)
    modeldata = AlgoDeforLinearModule.exportstresselementwise(modeldata)
    vm  = modeldata["postprocessing"]["exported"][1]["field"]
    println("extremes of vm, elemental: $([minimum(vm.values),   maximum(vm.values)])")
    
    # Write out mesh with von Mises stresses, elementwise
    modeldata["postprocessing"] = FDataDict("file"=>"twisted_beam_msh8-ew",
    "quantity"=> :Cauchy, "component"=> :xz)
    modeldata = AlgoDeforLinearModule.exportstresselementwise(modeldata)
    
    println("Done")
    true
end # twisted_beam_msh8


function twisted_beam_msh8_algo_stress()
    println("""
    The initially twisted cantilever beam is one of the standard test    problems for verifying the finite-element accuracy [1]. The beam is        clamped at one end and loaded either with unit in-plane or        unit out-of-plane force at the other. The centroidal axis of the beam s
    straight at the undeformed  configuration, while its cross-sections are
    twisted about the centroidal axis from 0 at the clamped end to pi/2 at
    the free end.
    
    Reference:
    Zupan D, Saje M (2004) On "A proposed standard set of problems to test
    finite element accuracy": the twisted beam. Finite Elements in Analysis
    and Design 40: 1445-1451.
    """)
    E = 0.29e8;
    nu = 0.22;
    W = 1.1;
    L = 12.;
    t =  0.32;
    nl = 2; nt = 1; nw = 1; ref = 4;
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
    el1femm  = FEMMBase(IntegData(subset(boundaryfes,Toplist), GaussRule(2, 2)))
    flux1 = FDataDict("femm"=>el1femm, "traction_vector"=>loadv)
    
    
    # Make the region
    MR = DeforModelRed3D
    material = MatDeforElastIso(MR, 00.0, E, nu, 0.0)
    region1 = FDataDict("femm"=>FEMMDeforLinearMSH8(MR, IntegData(fes, GaussRule(3,2)),    material))
    
    # Make model data
    modeldata =  FDataDict(
    "fens"=> fens, "regions"=>  [region1],
    "essential_bcs"=>[e1, e2, e3], "traction_bcs"=>  [flux1])
    
    # Call the solver
    modeldata = AlgoDeforLinearModule.linearstatics(modeldata)
    geom = modeldata["geom"]
    u = modeldata["u"]
    
    # Extract the solution
    nl = selectnode(fens, box = [L L -100*W 100*W -100*W 100*W],inflate = tolerance);
    theutip = mean(u.values[nl,:],1)
    println("displacement  = $(theutip[dir]) as compared to converged $uex")
    println("normalized displacement  = $(theutip[dir]/uex*100) %")
    
    # Write out mesh with displacements
    modeldata["postprocessing"] = FDataDict("file"=>"twisted_beam")
    modeldata = AlgoDeforLinearModule.exportdeformation(modeldata)
    
    # Write out mesh with stresses
    modeldata["postprocessing"] = FDataDict("file"=>"twisted_beam",
    "quantity"=> :Cauchy, "component"=> :xy)
    modeldata = AlgoDeforLinearModule.exportstress(modeldata)
    
    # Write out mesh with stresses
    modeldata["postprocessing"] = FDataDict("file"=>"twisted_beam",
    "quantity"=> :Cauchy, "component"=> :xz)
    modeldata = AlgoDeforLinearModule.exportstress(modeldata)
    
    # Write out mesh with von Mises stresses
    modeldata["postprocessing"] = FDataDict("file"=>"twisted_beam",
    "quantity"=> :vm)
    modeldata = AlgoDeforLinearModule.exportstress(modeldata)
    
    # Write out mesh with von Mises stresses, elementwise
    modeldata["postprocessing"] = FDataDict("file"=>"twisted_beam-ew",
    "quantity"=> :vm)
    modeldata = AlgoDeforLinearModule.exportstresselementwise(modeldata)
    vm  = modeldata["postprocessing"]["exported"][1]["field"]
    println("extremes of von Mises: $([minimum(vm.values),   maximum(vm.values)])")
    
    # Write out mesh with von Mises stresses, elementwise
    modeldata["postprocessing"] = FDataDict("file"=>"twisted_beam-ew",
    "quantity"=> :Cauchy, "component"=> :xz)
    modeldata = AlgoDeforLinearModule.exportstresselementwise(modeldata)
    
    # Write out mesh with principal stresses, elementwise
    modeldata["postprocessing"] = FDataDict("file"=>"twisted_beam-principal-1-ew",
    "quantity"=> :princCauchy, "component"=> 1)
    modeldata = AlgoDeforLinearModule.exportstresselementwise(modeldata)
    ps  = modeldata["postprocessing"]["exported"][1]["field"]
    println("extremes of first principal stress: $([minimum(ps.values),   maximum(ps.values)])")
    
    # Write out mesh with principal stresses, elementwise
    modeldata["postprocessing"] = FDataDict("file"=>"twisted_beam-principal-3-ew",
    "quantity"=> :princCauchy, "component"=> 3)
    modeldata = AlgoDeforLinearModule.exportstresselementwise(modeldata)
    ps  = modeldata["postprocessing"]["exported"][1]["field"]
    println("extremes of third principal stress: $([minimum(ps.values),   maximum(ps.values)])")
    
    # Write out mesh with principal stresses, elementwise
    modeldata["postprocessing"] = FDataDict("file"=>"twisted_beam-press-ew",
    "quantity"=> :pressure, "component"=> 1)
    modeldata = AlgoDeforLinearModule.exportstresselementwise(modeldata)
    ps  = modeldata["postprocessing"]["exported"][1]["field"]
    println("extremes of pressure: $([minimum(ps.values),   maximum(ps.values)])")
    
    println("Done")
    true
end # twisted_beam_msh8_algo_stress

function allrun()
    println("#####################################################") 
    println("# twisted_beam_algo ")
    twisted_beam_algo()
    println("#####################################################") 
    println("# twisted_beam_algo_stress ")
    twisted_beam_algo_stress()
    println("#####################################################") 
    println("# twisted_beam_export ")
    twisted_beam_export()
    println("#####################################################") 
    println("# twisted_beam_export_nb ")
    twisted_beam_export_nb()
    println("#####################################################") 
    println("# twisted_beam_msh8 ")
    twisted_beam_msh8()
    println("#####################################################") 
    println("# twisted_beam_msh8_algo_stress ")
    twisted_beam_msh8_algo_stress()
    return true
end # function allrun

end # module twisted_beam_examples
