module rltb_examples
using FinEtools
using FinEtools.AlgoBaseModule: evalconvergencestudy
using FinEtools.AlgoDeforLinearModule: linearstatics, exportstresselementwise, exportstress
using Statistics: mean
using LinearAlgebra: Symmetric, cholesky

# Isotropic material
E=1000.0;
nu=0.4999; #Taylor data
W=2.5;
H=5.0;
L= 50.0;
htol = minimum([L,H,W])/1000;
uzex =-12.6;
magn = 0.2*uzex/4;
Force =magn*W*H*2;
CTE = 0.0
mult = 8
n = 3 #

function getfrcL!(forceout::FFltVec, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)
    copyto!(forceout, [0.0; 0.0; magn])
end

function rltb_H8_by_hand()
    elementtag = "H8"
    println("""
    Taylor Cantilever example. Element: $(elementtag)
    """)

    fens,fes = H8block(W, L, H, n, mult*n, 2*n)
    bfes = meshboundary(fes)
    # end cross-section surface  for the shear loading
    sectionL = selectelem(fens, bfes; facing=true, direction = [0.0 +1.0 0.0])
    # 0 cross-section surface  for the reactions
    section0 = selectelem(fens, bfes; facing=true, direction = [0.0 -1.0 0.0])
    # 0 cross-section surface  for the reactions
    sectionlateral = selectelem(fens, bfes; facing=true, direction = [1.0 0.0 0.0])

    MR = DeforModelRed3D
    material = MatDeforElastIso(MR, 0.0, E, nu, CTE)
    
    # Material orientation matrix
    csmat = [i==j ? one(FFlt) : zero(FFlt) for i=1:3, j=1:3]
    
    function updatecs!(csmatout::FFltMat, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)
        copyto!(csmatout, csmat)
    end
        
    femm = FEMMDeforLinear(MR, IntegData(fes, GaussRule(3, 2)), material)
    
    geom = NodalField(fens.xyz)
    u = NodalField(zeros(size(fens.xyz,1), 3)) # displacement field
    
    lx0 = connectednodes(subset(bfes, section0))
    setebc!(u,lx0,true,1,0.0)
    setebc!(u,lx0,true,2,0.0)
    setebc!(u,lx0,true,3,0.0)
    lx1 = connectednodes(subset(bfes, sectionlateral))
    setebc!(u,lx1,true,1,0.0)
    applyebc!(u)
    numberdofs!(u)

    fi = ForceIntensity(FFlt, 3, getfrcL!);
    el2femm = FEMMBase(IntegData(subset(bfes, sectionL), GaussRule(2, 2)))
    F2 = distribloads(el2femm, geom, u, fi, 2);
    @show sum(F2)
    associategeometry!(femm, geom)
    K = stiffness(femm, geom, u)
    K=cholesky(K)
    U = K\(F2)
    scattersysvec!(u,U[:])
    @show length(U)
    Tipl = selectnode(fens, box=[0 W L L 0 H], inflate=htol)
    utip = mean(u.values[Tipl, 3], dims=1)
    println("Deflection: $(utip), compared to $(uzex)")

    File =  "rltb_H8_by_hand.vtk"
    vtkexportmesh(File, fens, fes;  vectors=[("u", u.values)])
    @async run(`"paraview.exe" $File`)

    # modeldata["postprocessing"] = FDataDict("file"=>"hughes_cantilever_stresses_$(elementtag)", "outputcsys"=>CSys(3, 3, updatecs!), "quantity"=>:Cauchy, "component"=>[5])
    # modeldata = exportstresselementwise(modeldata)
    
    # modeldata["postprocessing"] = FDataDict("file"=>"hughes_cantilever_stresses_$(elementtag)",
    # "outputcsys"=>CSys(3, 3, updatecs!), "quantity"=>:Cauchy,
    # "component"=>collect(1:6))
    # modeldata = exportstresselementwise(modeldata)
    # stressfields = ElementalField[modeldata["postprocessing"]["exported"][1]["field"]]
    
    true
    
end # hughes_cantilever_stresses_MST10

function rltb_H20_by_hand()
    elementtag = "H20"
    println("""
    Taylor Cantilever example. Element: $(elementtag)
    """)

    fens,fes = H20block(W, L, H, n, mult*n, 2*n)
    bfes = meshboundary(fes)
    # end cross-section surface  for the shear loading
    sectionL = selectelem(fens, bfes; facing=true, direction = [0.0 +1.0 0.0])
    # 0 cross-section surface  for the reactions
    section0 = selectelem(fens, bfes; facing=true, direction = [0.0 -1.0 0.0])
    # 0 cross-section surface  for the reactions
    sectionlateral = selectelem(fens, bfes; facing=true, direction = [1.0 0.0 0.0])

    MR = DeforModelRed3D
    material = MatDeforElastIso(MR, 0.0, E, nu, CTE)
    
    # Material orientation matrix
    csmat = [i==j ? one(FFlt) : zero(FFlt) for i=1:3, j=1:3]
    
    function updatecs!(csmatout::FFltMat, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)
        copyto!(csmatout, csmat)
    end
        
    femm = FEMMDeforLinear(MR, IntegData(fes, GaussRule(3, 2)), material)
    
    geom = NodalField(fens.xyz)
    u = NodalField(zeros(size(fens.xyz,1), 3)) # displacement field
    
    lx0 = connectednodes(subset(bfes, section0))
    setebc!(u,lx0,true,1,0.0)
    setebc!(u,lx0,true,2,0.0)
    setebc!(u,lx0,true,3,0.0)
    lx1 = connectednodes(subset(bfes, sectionlateral))
    setebc!(u,lx1,true,1,0.0)
    applyebc!(u)
    numberdofs!(u)

    fi = ForceIntensity(FFlt, 3, getfrcL!);
    el2femm = FEMMBase(IntegData(subset(bfes, sectionL), GaussRule(2, 2)))
    F2 = distribloads(el2femm, geom, u, fi, 2);
    @show sum(F2)
    associategeometry!(femm, geom)
    K = stiffness(femm, geom, u)
    K=cholesky(K)
    U = K\(F2)
    scattersysvec!(u,U[:])
    @show length(U)
    Tipl = selectnode(fens, box=[0 W L L 0 H], inflate=htol)
    utip = mean(u.values[Tipl, 3], dims=1)
    println("Deflection: $(utip), compared to $(uzex)")

    File =  "rltb_H20_by_hand.vtk"
    vtkexportmesh(File, fens, fes;  vectors=[("u", u.values)])
    @async run(`"paraview.exe" $File`)

    # modeldata["postprocessing"] = FDataDict("file"=>"hughes_cantilever_stresses_$(elementtag)", "outputcsys"=>CSys(3, 3, updatecs!), "quantity"=>:Cauchy, "component"=>[5])
    # modeldata = exportstresselementwise(modeldata)
    
    # modeldata["postprocessing"] = FDataDict("file"=>"hughes_cantilever_stresses_$(elementtag)",
    # "outputcsys"=>CSys(3, 3, updatecs!), "quantity"=>:Cauchy,
    # "component"=>collect(1:6))
    # modeldata = exportstresselementwise(modeldata)
    # stressfields = ElementalField[modeldata["postprocessing"]["exported"][1]["field"]]
    
    true
    
end # hughes_cantilever_stresses_MST10

function allrun()
    println("#####################################################") 
    println("# hughes_cantilever_stresses_MST10 ")
    hughes_cantilever_stresses_MST10()
    println("#####################################################") 
    println("# hughes_cantilever_stresses_MST10_incompressible ")
    hughes_cantilever_stresses_MST10_incompressible()
    println("#####################################################") 
    println("# hughes_cantilever_stresses_nodal_MST10 ")
    hughes_cantilever_stresses_nodal_MST10()
    println("#####################################################") 
    println("# hughes_cantilever_stresses_nodal_T10 ")
    hughes_cantilever_stresses_nodal_T10()
    println("#####################################################") 
    println("# hughes_cantilever_stresses_T10 ")
    hughes_cantilever_stresses_T10()
    println("#####################################################") 
    println("# hughes_cantilever_stresses_T10_incompressible ")
    hughes_cantilever_stresses_T10_incompressible()
    return true
end # function allrun

end # module hughes_cantilever_examples
