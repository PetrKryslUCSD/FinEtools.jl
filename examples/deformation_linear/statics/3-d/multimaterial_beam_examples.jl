module multimaterial_beam_examples
using FinEtools
using FinEtools.AlgoDeforLinearModule
using FinEtools.MeshExportModule


function multimaterial_beam_algo()
    println("""
    Multi-material beam. Rubber-like and metal-like halves,
    clamped, with shear traction at free end.
    """)
    E1 = 0.29e3;
    nu1 = 0.49;
    E2 = 0.4e4;
    nu2 = 0.3;
    W = 4.1;
    L = 12.;
    t = 6.5;
    nl = 2; nt = 1; nw = 1; ref = 9;
    p =   200.0/W/t;
    #  Loading in the Z direction
    loadv = [0;0;p]; dir = 3;
    tolerance  = t/1000;
    
    fens,fes  = H20block(L,W,t, nl*ref,nw*ref,nt*ref)
    
    # Clamped end of the beam
    l1  = selectnode(fens; box = [0 0 -Inf Inf -Inf Inf], inflate  =  tolerance)
    e1 = FDataDict("node_list"=>l1, "component"=>1, "displacement"=>0.0)
    e2 = FDataDict("node_list"=>l1, "component"=>2, "displacement"=>0.0)
    e3 = FDataDict("node_list"=>l1, "component"=>3, "displacement"=>0.0)
    
    # Traction on the opposite edge
    boundaryfes  =   meshboundary(fes);
    Toplist   = selectelem(fens,boundaryfes, box =  [L L -Inf Inf -Inf Inf], inflate =   tolerance);
    el1femm  =   FEMMBase(IntegData(subset(boundaryfes,Toplist), GaussRule(2, 2)))
    flux1 = FDataDict("femm"=>el1femm, "traction_vector"=>loadv)
    
    r1list   = selectelem(fens,fes, box =  [0 L/2. -Inf Inf -Inf Inf], inflate =   tolerance);
    r2list   = selectelem(fens,fes, box =  [L/2. L -Inf Inf -Inf Inf], inflate =   tolerance);
    
    # Model reduction type
    MR = DeforModelRed3D
    
    # Make region 1
    region1 = FDataDict("femm"=>FEMMDeforLinear(MR,
    IntegData(subset(fes,r1list), GaussRule(3,2)),
    MatDeforElastIso(MR, 0.0, E1, nu1, 0.0)))
    
    # Make region 2
    region2 = FDataDict("femm"=>FEMMDeforLinear(MR,
    IntegData(subset(fes,r2list), GaussRule(3,2)),
    MatDeforElastIso(MR, 0.0, E2, nu2, 0.0)))
    
    # Make model data
    modeldata =  FDataDict(
    "fens"=> fens, "regions"=>  [region1, region2],
    "essential_bcs"=>[e1, e2, e3], "traction_bcs"=>  [flux1])
    
    # Call the solver
    modeldata = AlgoDeforLinearModule.linearstatics(modeldata)
    geom = modeldata["geom"]
    u = modeldata["u"]
    
    # Write out mesh with displacements
    modeldata["postprocessing"] = FDataDict("file"=>"multimaterial_beam")
    modeldata = AlgoDeforLinearModule.exportdeformation(modeldata)
    
    # Write out mesh with stresses
    modeldata["postprocessing"] = FDataDict("file"=>"multimaterial_beam_xy",
    "quantity"=> :Cauchy, "component"=> :xy)
    modeldata = AlgoDeforLinearModule.exportstress(modeldata)
    
    # Write out mesh with stresses
    modeldata["postprocessing"] = FDataDict("file"=>"multimaterial_beam_xz",
    "quantity"=> :Cauchy, "component"=> :xz)
    modeldata = AlgoDeforLinearModule.exportstress(modeldata)
    
    # Write out mesh with von Mises stresses
    modeldata["postprocessing"] = FDataDict("file"=>"multimaterial_beam_vm",
    "quantity"=> :vm)
    modeldata = AlgoDeforLinearModule.exportstress(modeldata)
    
    # Write out mesh with von Mises stresses, elementwise
    modeldata["postprocessing"] = FDataDict("file"=>"multimaterial_beam_vm-ew",
    "quantity"=> :vm)
    modeldata = AlgoDeforLinearModule.exportstresselementwise(modeldata)
    
end # multimaterial_beam_algo

function allrun()
    println("#####################################################") 
    println("# multimaterial_beam_algo ")
    multimaterial_beam_algo()
    return true
end # function allrun

end # module multimaterial_beam_examples
