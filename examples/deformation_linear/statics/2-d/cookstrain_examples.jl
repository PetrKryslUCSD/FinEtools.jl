module cookstrain_examples
using FinEtools
using FinEtools.AlgoDeforLinearModule
using FinEtools.MeshExportModule

function cookstrain_algo_export()
    E = 1.0;
    nu = 1.0/3;
    width = 48.0; height = 44.0; thickness  = 1.0;
    free_height  = 16.0;
    Mid_edge  = [48.0, 52.0];# Location of tracked  deflection
    magn = 1.0/(free_height*thickness);# Density of applied load
    convutip = 23.97;
    n = 30;# number of elements per side
    tolerance = minimum([width, height])/n/1000.;#Geometrical tolerance
    
    fens,fes = T3block(width, height, n, n)
    
    # Reshape into a trapezoidal panel
    for i = 1:count(fens)
        fens.xyz[i,2] = fens.xyz[i,2]+(fens.xyz[i,1]/width)*(height -fens.xyz[i,2]/height*(height-free_height));
    end
    
    # Clamped edge of the membrane
    l1 = selectnode(fens; box=[0.,0.,-Inf, Inf], inflate = tolerance)
    ess1 = FDataDict("displacement"=>  0.0, "component"=> 1, "node_list"=>l1)
    ess2 = FDataDict("displacement"=>  0.0, "component"=> 2, "node_list"=>l1)
    
    # Traction on the opposite edge
    boundaryfes =  meshboundary(fes);
    Toplist  = selectelem(fens, boundaryfes, box= [width, width, -Inf, Inf ], inflate=  tolerance);
    el1femm = FEMMBase(IntegData(subset(boundaryfes, Toplist), GaussRule(1, 2), thickness))
    flux1 = FDataDict("traction_vector"=>[0.0,+magn],
    "femm"=>el1femm
    )
    
    # Make the region
    MR = DeforModelRed2DStrain
    material = MatDeforElastIso(MR,  0.0, E, nu, 0.0)
    region1 = FDataDict("femm"=>FEMMDeforLinear(MR,
    IntegData(fes, TriRule(1), thickness), material))
    
    modeldata = FDataDict("fens"=>fens,
    "regions"=>[region1],
    "essential_bcs"=>[ess1, ess2],
    "traction_bcs"=>[flux1]
    )
    
    # Call the solver
    modeldata = AlgoDeforLinearModule.linearstatics(modeldata)
    
    u = modeldata["u"]
    geom = modeldata["geom"]
    
    # Extract the solution
    nl = selectnode(fens, box=[Mid_edge[1],Mid_edge[1],Mid_edge[2],Mid_edge[2]],
    inflate=tolerance);
    theutip = u.values[nl,:]
    println("displacement =$(theutip[2]) as compared to converged $convutip")
    
    modeldata["postprocessing"] = FDataDict("file"=>"cookstress-ew",
    "quantity"=>:Cauchy, "component"=>:xy)
    modeldata = AlgoDeforLinearModule.exportstresselementwise(modeldata)
    fld = modeldata["postprocessing"]["exported"][1]["field"]
    println("range of Cauchy_xy = $((minimum(fld.values), maximum(fld.values)))")
    File = modeldata["postprocessing"]["exported"][1]["file"]
    @async run(`"paraview.exe" $File`)
    
    modeldata["postprocessing"] = FDataDict("file"=>"cookstress-ew-vm",
    "quantity"=>:vm, "component"=>1)
    modeldata = AlgoDeforLinearModule.exportstresselementwise(modeldata)
    fld = modeldata["postprocessing"]["exported"][1]["field"]
    println("range of vm = $((minimum(fld.values), maximum(fld.values)))")
    File = modeldata["postprocessing"]["exported"][1]["file"]
    @async run(`"paraview.exe" $File`)
    
    modeldata["postprocessing"] = FDataDict("file"=>"cookstress-ew-pressure",
    "quantity"=>:pressure, "component"=>1)
    modeldata = AlgoDeforLinearModule.exportstresselementwise(modeldata)
    fld = modeldata["postprocessing"]["exported"][1]["field"]
    println("range of pressure = $((minimum(fld.values), maximum(fld.values)))")
    File = modeldata["postprocessing"]["exported"][1]["file"]
    @async run(`"paraview.exe" $File`)
    
    modeldata["postprocessing"] = FDataDict("file"=>"cookstress-ew-princ1",
    "quantity"=>:princCauchy, "component"=>1)
    modeldata = AlgoDeforLinearModule.exportstresselementwise(modeldata)
    fld = modeldata["postprocessing"]["exported"][1]["field"]
    println("range of princCauchy Max = $((minimum(fld.values), maximum(fld.values)))")
    File = modeldata["postprocessing"]["exported"][1]["file"]
    @async run(`"paraview.exe" $File`)
    
    modeldata["postprocessing"] = FDataDict("file"=>"cookstress-ew-princ3",
    "quantity"=>:princCauchy, "component"=>3)
    modeldata = AlgoDeforLinearModule.exportstresselementwise(modeldata)
    fld = modeldata["postprocessing"]["exported"][1]["field"]
    println("range of princCauchy Min = $((minimum(fld.values), maximum(fld.values)))")
    File = modeldata["postprocessing"]["exported"][1]["file"]
    @async run(`"paraview.exe" $File`)
    
    AE = AbaqusExporter("Cookstress_algo_stress");
    HEADING(AE, "Cook trapezoidal panel, plane stress");
    COMMENT(AE, "Converged free mid-edge displacement = 23.97");
    PART(AE, "part1");
    END_PART(AE);
    ASSEMBLY(AE, "ASSEM1");
    INSTANCE(AE, "INSTNC1", "PART1");
    NODE(AE, fens.xyz);
    COMMENT(AE, "We are assuming three node triangles in plane-stress");
    COMMENT(AE, "CPE3 are pretty poor-accuracy elements, but here we don't care about it.");
    @assert nodesperelem(modeldata["regions"][1]["femm"].integdata.fes) == 3
    ELEMENT(AE, "CPE3", "AllElements", connasarray(modeldata["regions"][1]["femm"].integdata.fes))
    NSET_NSET(AE, "clamped", modeldata["essential_bcs"][1]["node_list"])
    ORIENTATION(AE, "GlobalOrientation", vec([1. 0 0]), vec([0 1. 0]));
    SOLID_SECTION(AE, "elasticity", "GlobalOrientation", "AllElements", thickness);
    END_INSTANCE(AE);
    END_ASSEMBLY(AE);
    MATERIAL(AE, "elasticity")
    ELASTIC(AE, E, nu)
    STEP_PERTURBATION_STATIC(AE)
    BOUNDARY(AE, "ASSEM1.INSTNC1.clamped", 1)
    BOUNDARY(AE, "ASSEM1.INSTNC1.clamped", 2)
    bfes = modeldata["traction_bcs"][1]["femm"].integdata.fes
    COMMENT(AE, "Concentrated loads: we are assuming that the elements on the boundary");
    COMMENT(AE, "have two nodes each and also that they are the same length.");
    COMMENT(AE, "Then the concentrated loads below will be correctly lumped.");
    nl = connectednodes(bfes)
    F = zeros(count(modeldata["fens"]))
    for ix = 1:count(bfes)
        for jx = 1:2
            F[bfes.conn[ix][jx]] += 1.0/n/2/thickness
        end
    end
    for ixxxx = 1:length(F)
        if F[ixxxx] != 0.0
            CLOAD(AE, "ASSEM1.INSTNC1.$(ixxxx)", 2, F[ixxxx])
        end
    end
    END_STEP(AE)
    close(AE)
    
    true
    
end # cookstrain_algo_export


function cookstrain_algo_export_ortho()
    E = 1.0;
    nu = 1.0/3;
    width = 48.0; height = 44.0; thickness  = 1.0;
    free_height  = 16.0;
    Mid_edge  = [48.0, 52.0];# Location of tracked  deflection
    magn = 1.0/(free_height*thickness);# Density of applied load
    convutip = 23.97;
    n = 30;# number of elements per side
    tolerance = minimum([width, height])/n/1000.;#Geometrical tolerance
    
    fens,fes = T3block(width, height, n, n)
    
    # Reshape into a trapezoidal panel
    for i = 1:count(fens)
        fens.xyz[i,2] = fens.xyz[i,2]+(fens.xyz[i,1]/width)*(height -fens.xyz[i,2]/height*(height-free_height));
    end
    
    # Clamped edge of the membrane
    l1 = selectnode(fens; box=[0.,0.,-Inf, Inf], inflate = tolerance)
    ess1 = FDataDict("displacement"=>  0.0, "component"=> 1, "node_list"=>l1)
    ess2 = FDataDict("displacement"=>  0.0, "component"=> 2, "node_list"=>l1)
    
    # Traction on the opposite edge
    boundaryfes =  meshboundary(fes);
    Toplist  = selectelem(fens, boundaryfes, box= [width, width, -Inf, Inf ], inflate=  tolerance);
    el1femm = FEMMBase(IntegData(subset(boundaryfes, Toplist), GaussRule(1, 2), thickness))
    flux1 = FDataDict("traction_vector"=>[0.0,+magn],
    "femm"=>el1femm
    )
    
    # Make the region
    MR = DeforModelRed2DStrain
    material = MatDeforElastOrtho(MR,  0.0, E, nu, 0.0)
    region1 = FDataDict("femm"=>FEMMDeforLinear(MR,
    IntegData(fes, TriRule(1), thickness), material))
    
    modeldata = FDataDict("fens"=>fens,
    "regions"=>[region1],
    "essential_bcs"=>[ess1, ess2],
    "traction_bcs"=>[flux1]
    )
    
    # Call the solver
    modeldata = AlgoDeforLinearModule.linearstatics(modeldata)
    
    u = modeldata["u"]
    geom = modeldata["geom"]
    
    # Extract the solution
    nl = selectnode(fens, box=[Mid_edge[1],Mid_edge[1],Mid_edge[2],Mid_edge[2]],
    inflate=tolerance);
    theutip = u.values[nl,:]
    println("displacement =$(theutip[2]) as compared to converged $convutip")
    
    modeldata["postprocessing"] = FDataDict("file"=>"cookstress-ew",
    "quantity"=>:Cauchy, "component"=>:xy)
    modeldata = AlgoDeforLinearModule.exportstresselementwise(modeldata)
    fld = modeldata["postprocessing"]["exported"][1]["field"]
    println("range of Cauchy_xy = $((minimum(fld.values), maximum(fld.values)))")
    File = modeldata["postprocessing"]["exported"][1]["file"]
    @async run(`"paraview.exe" $File`)
    
    modeldata["postprocessing"] = FDataDict("file"=>"cookstress-ew-vm",
    "quantity"=>:vm, "component"=>1)
    modeldata = AlgoDeforLinearModule.exportstresselementwise(modeldata)
    fld = modeldata["postprocessing"]["exported"][1]["field"]
    println("range of vm = $((minimum(fld.values), maximum(fld.values)))")
    File = modeldata["postprocessing"]["exported"][1]["file"]
    @async run(`"paraview.exe" $File`)
    
    modeldata["postprocessing"] = FDataDict("file"=>"cookstress-ew-pressure",
    "quantity"=>:pressure, "component"=>1)
    modeldata = AlgoDeforLinearModule.exportstresselementwise(modeldata)
    fld = modeldata["postprocessing"]["exported"][1]["field"]
    println("range of pressure = $((minimum(fld.values), maximum(fld.values)))")
    File = modeldata["postprocessing"]["exported"][1]["file"]
    @async run(`"paraview.exe" $File`)
    
    modeldata["postprocessing"] = FDataDict("file"=>"cookstress-ew-princ1",
    "quantity"=>:princCauchy, "component"=>1)
    modeldata = AlgoDeforLinearModule.exportstresselementwise(modeldata)
    fld = modeldata["postprocessing"]["exported"][1]["field"]
    println("range of princCauchy Max = $((minimum(fld.values), maximum(fld.values)))")
    File = modeldata["postprocessing"]["exported"][1]["file"]
    @async run(`"paraview.exe" $File`)
    
    modeldata["postprocessing"] = FDataDict("file"=>"cookstress-ew-princ3",
    "quantity"=>:princCauchy, "component"=>3)
    modeldata = AlgoDeforLinearModule.exportstresselementwise(modeldata)
    fld = modeldata["postprocessing"]["exported"][1]["field"]
    println("range of princCauchy Min = $((minimum(fld.values), maximum(fld.values)))")
    File = modeldata["postprocessing"]["exported"][1]["file"]
    @async run(`"paraview.exe" $File`)
    
    AE = AbaqusExporter("Cookstress_algo_stress");
    HEADING(AE, "Cook trapezoidal panel, plane stress");
    COMMENT(AE, "Converged free mid-edge displacement = 23.97");
    PART(AE, "part1");
    END_PART(AE);
    ASSEMBLY(AE, "ASSEM1");
    INSTANCE(AE, "INSTNC1", "PART1");
    NODE(AE, fens.xyz);
    COMMENT(AE, "We are assuming three node triangles in plane-stress");
    COMMENT(AE, "CPE3 are pretty poor-accuracy elements, but here we don't care about it.");
    @assert nodesperelem(modeldata["regions"][1]["femm"].integdata.fes) == 3
    ELEMENT(AE, "CPE3", "AllElements", connasarray(modeldata["regions"][1]["femm"].integdata.fes))
    NSET_NSET(AE, "clamped", modeldata["essential_bcs"][1]["node_list"])
    ORIENTATION(AE, "GlobalOrientation", vec([1. 0 0]), vec([0 1. 0]));
    SOLID_SECTION(AE, "elasticity", "GlobalOrientation", "AllElements", thickness);
    END_INSTANCE(AE);
    END_ASSEMBLY(AE);
    MATERIAL(AE, "elasticity")
    ELASTIC(AE, E, nu)
    STEP_PERTURBATION_STATIC(AE)
    BOUNDARY(AE, "ASSEM1.INSTNC1.clamped", 1)
    BOUNDARY(AE, "ASSEM1.INSTNC1.clamped", 2)
    bfes = modeldata["traction_bcs"][1]["femm"].integdata.fes
    COMMENT(AE, "Concentrated loads: we are assuming that the elements on the boundary");
    COMMENT(AE, "have two nodes each and also that they are the same length.");
    COMMENT(AE, "Then the concentrated loads below will be correctly lumped.");
    nl = connectednodes(bfes)
    F = zeros(count(modeldata["fens"]))
    for ix = 1:count(bfes)
        for jx = 1:2
            F[bfes.conn[ix][jx]] += 1.0/n/2/thickness
        end
    end
    for ixxxx = 1:length(F)
        if F[ixxxx] != 0.0
            CLOAD(AE, "ASSEM1.INSTNC1.$(ixxxx)", 2, F[ixxxx])
        end
    end
    END_STEP(AE)
    close(AE)
    
    true
    
end # cookstrain_algo_export_ortho

function allrun()
    println("#####################################################") 
    println("# cookstrain_algo_export ")
    cookstrain_algo_export()
    println("#####################################################") 
    println("# cookstrain_algo_export_ortho ")
    cookstrain_algo_export_ortho()
    return true
end # function allrun

end # module cookstrain_examples
