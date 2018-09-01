module hughes_cantilever_examples
using FinEtools
using FinEtools.AlgoBaseModule: evalconvergencestudy
using FinEtools.AlgoDeforLinearModule: linearstatics, exportstresselementwise, exportstress
using Statistics: mean
using LinearAlgebra: Symmetric, cholesky

# Example from TJR Hughes, The Finite Element Method: Linear Static and Dynamic Finite Element Analysis, 1987.
# Two-dimensional plane-elasticity Solution of beam deflection. Shear force P at one end, and corresponding reactions (tractions) at the other end. The Analytical solution is plane-strain.
# Find out: where does the exact solution come from?

# Isotropic material
E=1.0;
P=1.0;
L=16.0;
c=2.0;
h=1.0;
W=2/3*h*c^3;
CTE = 0.0
tolerance = 0.00001*c

function getfrcL!(forceout::FFltVec, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)
    copyto!(forceout, [0.0; 0.0; P/2/W*(c^2-XYZ[3]^2)])
end

function getfrc0!(forceout::FFltVec, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)
    copyto!(forceout, [P*L/W*XYZ[3]; 0.0; -P/2/W*(c^2-XYZ[3]^2)])
end

function evaluateerrors(filebase, modeldatasequence)
    println("")
    println("Stress RMS error")
    for md = modeldatasequence
        md["targetfields"] = [e["field"] for e in md["postprocessing"]["exported"]]
    end
    elementsizes, errornorms, p = evalconvergencestudy(modeldatasequence)
    
    println("Normalized Approximate Error = $(errornorms)")
    
    f = log.(vec(errornorms))
    A = hcat(log.(vec(elementsizes[1:end-1])), ones(size(f)))
    p = A \ f
    println("Linear log-log fit: p = $(p)")
    
    csvFile = filebase * "_Stress" * ".CSV"
    savecsv(csvFile,
    elementsizes=vec(elementsizes[1:end-1]),
    elementsizes2=vec(elementsizes[1:end-1].^2),
    elementsizes3=vec(elementsizes[1:end-1].^3),
    errornorms=vec(errornorms)
    )
    println("Wrote $csvFile")
    
    println("")
    println("Displacement RMS error")
    for md = modeldatasequence
        md["targetfields"] = [md["u"] for r in md["regions"]]
    end
    elementsizes, errornorms, p = evalconvergencestudy(modeldatasequence)
    
    println("Normalized Approximate Error = $(errornorms)")
    
    f = log.(vec(errornorms))
    A = hcat(log.(vec(elementsizes[1:end-1])), ones(size(f)))
    p = A \ f
    println("Linear log-log fit: p = $(p)")
    
    csvFile = filebase * "_Displ" * ".CSV"
    savecsv(csvFile,
    elementsizes=vec(elementsizes[1:end-1]),
    elementsizes2=vec(elementsizes[1:end-1].^2),
    elementsizes3=vec(elementsizes[1:end-1].^3),
    errornorms=vec(errornorms)
    )
    println("Wrote $csvFile")
end

function hughes_cantilever_stresses_H8_by_hand()
    elementtag = "H8"
    println("""
    Cantilever example.  Hughes 1987. Element: $(elementtag)
    """)
    nu=0.3; # COMPRESSIBLE
    E1=E/(1-nu^2);  # Plane strain
    nu1=nu/(1-nu); I=h*(2*c)^3/12;
    exactux(x,y) = (P/(6*E1*I)*(-y)*(3*(L^2-(L-x)^2)+(2+nu1)*(y^2-c^2)));
    exactuy(x,y) = (P/(6*E1*I)*(((L-x)^3-L^3)-((4+5*nu1)*c^2+3*L^2)*(L-x-L)+3*nu1*(L-x)*y^2));

    n = 2 #
        
    nL = 3*n # number of elements lengthwise
    nc = 2*n # number of elements through the depth
    nh = 1 # number of elements through the thickness
    xs = collect(linearspace(0.0, L, nL+1))
    ys = collect(linearspace(0.0, h, nh+1))
    zs = collect(linearspace(-c, +c, nc+1))
    fens,fes = H8blockx(xs, ys, zs)
    # fens,fes = H8toH20(fens,fes)
    bfes = meshboundary(fes)
    # end cross-section surface  for the shear loading
    sshearL = selectelem(fens, bfes; facing=true, direction = [+1.0 0.0 0.0])
    # 0 cross-section surface  for the reactions
    sshear0 = selectelem(fens, bfes; facing=true, direction = [-1.0 0.0 0.0])
    
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
    
    lx0 = selectnode(fens, box=[0.0 0.0 0.0 0.0 0.0 0.0], inflate=tolerance)
    setebc!(u,lx0,true,1,0.0)
    setebc!(u,lx0,true,2,0.0)
    setebc!(u,lx0,true,3,0.0)
    lx1 = selectnode(fens, box=[0.0 0.0 0.0 0.0 c c], inflate=tolerance)
    lx2 = selectnode(fens, box=[0.0 0.0 0.0 0.0 -c -c], inflate=tolerance)
    setebc!(u,vcat(lx1, lx2),true,1,0.0)
    ly1 = selectnode(fens, box=[-Inf Inf 0.0 0.0 -Inf Inf], inflate=tolerance)
    ly2 = selectnode(fens, box=[-Inf Inf h h -Inf Inf], inflate=tolerance)
    setebc!(u,vcat(ly1, ly2),true,2,0.0)
    applyebc!(u)
    numberdofs!(u)

    fi = ForceIntensity(FFlt, 3, getfrc0!);
    el1femm = FEMMBase(IntegData(subset(bfes, sshear0), GaussRule(2, 2)))
    F1 = distribloads(el1femm, geom, u, fi, 2);
    fi = ForceIntensity(FFlt, 3, getfrcL!);
    el2femm = FEMMBase(IntegData(subset(bfes, sshearL), GaussRule(2, 2)))
    F2 = distribloads(el2femm, geom, u, fi, 2);

    associategeometry!(femm, geom)
    K = stiffness(femm, geom, u)

    K = cholesky(K)
    U = K\(F1 + F2)
    scattersysvec!(u,U[:])
    
    Tipl = selectnode(fens, box=[L L 0.0 0.0 0. 0.], inflate=tolerance)
    utip = mean(u.values[Tipl, 3])
    println("Deflection: $(utip), compared to $(exactuy(L,0.0))")

    File =  "hughes_cantilever_stresses_H8_by_hand.vtk"
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
    
end # hughes_cantilever_stresses_H8_by_hand

function hughes_cantilever_stresses_H20_by_hand()
    elementtag = "H20"
    println("""
    Cantilever example.  Hughes 1987. Element: $(elementtag)
    """)
    nu=0.3; # COMPRESSIBLE
    E1=E/(1-nu^2);  # Plane strain
    nu1=nu/(1-nu); I=h*(2*c)^3/12;
    exactux(x,y) = (P/(6*E1*I)*(-y)*(3*(L^2-(L-x)^2)+(2+nu1)*(y^2-c^2)));
    exactuy(x,y) = (P/(6*E1*I)*(((L-x)^3-L^3)-((4+5*nu1)*c^2+3*L^2)*(L-x-L)+3*nu1*(L-x)*y^2));

    n = 8 #
        
    nL = 3*n # number of elements lengthwise
    nc = 2*n # number of elements through the depth
    nh = 1 # number of elements through the thickness
    xs = collect(linearspace(0.0, L, nL+1))
    ys = collect(linearspace(0.0, h, nh+1))
    zs = collect(linearspace(-c, +c, nc+1))
    fens,fes = H8blockx(xs, ys, zs)
    fens,fes = H8toH20(fens,fes)
    bfes = meshboundary(fes)
    # end cross-section surface  for the shear loading
    sshearL = selectelem(fens, bfes; facing=true, direction = [+1.0 0.0 0.0])
    # 0 cross-section surface  for the reactions
    sshear0 = selectelem(fens, bfes; facing=true, direction = [-1.0 0.0 0.0])
    
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
    
    lx0 = selectnode(fens, box=[0.0 0.0 0.0 0.0 0.0 0.0], inflate=tolerance)
    setebc!(u,lx0,true,1,0.0)
    setebc!(u,lx0,true,2,0.0)
    setebc!(u,lx0,true,3,0.0)
    lx1 = selectnode(fens, box=[0.0 0.0 0.0 0.0 c c], inflate=tolerance)
    lx2 = selectnode(fens, box=[0.0 0.0 0.0 0.0 -c -c], inflate=tolerance)
    setebc!(u,vcat(lx1, lx2),true,1,0.0)
    ly1 = selectnode(fens, box=[-Inf Inf 0.0 0.0 -Inf Inf], inflate=tolerance)
    ly2 = selectnode(fens, box=[-Inf Inf h h -Inf Inf], inflate=tolerance)
    setebc!(u,vcat(ly1, ly2),true,2,0.0)
    applyebc!(u)
    numberdofs!(u)
    u

    fi = ForceIntensity(FFlt, 3, getfrc0!);
    el1femm = FEMMBase(IntegData(subset(bfes, sshear0), GaussRule(2, 2)))
    F1 = distribloads(el1femm, geom, u, fi, 2);
    fi = ForceIntensity(FFlt, 3, getfrcL!);
    el2femm = FEMMBase(IntegData(subset(bfes, sshearL), GaussRule(2, 2)))
    F2 = distribloads(el2femm, geom, u, fi, 2);

    associategeometry!(femm, geom)
    K = stiffness(femm, geom, u)

    K = cholesky(K)
    U = K\(F1 + F2)
    scattersysvec!(u,U[:])
    
    Tipl = selectnode(fens, box=[L L 0.0 0.0 0. 0.], inflate=tolerance)
    utip = mean(u.values[Tipl, 3])
    println("Deflection: $(utip), compared to $(exactuy(L,0.0))")

    File =  "hughes_cantilever_stresses_H8_by_hand.vtk"
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

function hughes_cantilever_stresses_T10_by_hand()
    elementtag = "T10"
    println("""
    Cantilever example.  Hughes 1987. Element: $(elementtag)
    """)
    nu=0.3; # COMPRESSIBLE
    E1=E/(1-nu^2);  # Plane strain
    nu1=nu/(1-nu); I=h*(2*c)^3/12;
    exactux(x,y) = (P/(6*E1*I)*(-y)*(3*(L^2-(L-x)^2)+(2+nu1)*(y^2-c^2)));
    exactuy(x,y) = (P/(6*E1*I)*(((L-x)^3-L^3)-((4+5*nu1)*c^2+3*L^2)*(L-x-L)+3*nu1*(L-x)*y^2));

    n = 2 #
        
    nL = 3*n # number of elements lengthwise
    nc = 2*n # number of elements through the depth
    nh = n # number of elements through the thickness
    xs = collect(linearspace(0.0, L, nL+1))
    ys = collect(linearspace(0.0, h, nh+1))
    zs = collect(linearspace(-c, +c, nc+1))
    fens,fes = T10blockx(xs, ys, zs)
    bfes = meshboundary(fes)
    # end cross-section surface  for the shear loading
    sshearL = selectelem(fens, bfes; facing=true, direction = [+1.0 0.0 0.0])
    # 0 cross-section surface  for the reactions
    sshear0 = selectelem(fens, bfes; facing=true, direction = [-1.0 0.0 0.0])
    
    MR = DeforModelRed3D
    material = MatDeforElastIso(MR, 0.0, E, nu, CTE)
    
    # Material orientation matrix
    csmat = [i==j ? one(FFlt) : zero(FFlt) for i=1:3, j=1:3]
    
    function updatecs!(csmatout::FFltMat, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)
        copyto!(csmatout, csmat)
    end
        
    gr = SimplexRule(3, 4)
    
    femm = FEMMDeforLinear(MR, IntegData(fes, gr), material)
    
    geom = NodalField(fens.xyz)
    u = NodalField(zeros(size(fens.xyz,1), 3)) # displacement field
    
    lx0 = selectnode(fens, box=[0.0 0.0 0.0 0.0 0.0 0.0], inflate=tolerance)
    setebc!(u,lx0,true,1,0.0)
    setebc!(u,lx0,true,2,0.0)
    setebc!(u,lx0,true,3,0.0)
    lx1 = selectnode(fens, box=[0.0 0.0 0.0 0.0 c c], inflate=tolerance)
    lx2 = selectnode(fens, box=[0.0 0.0 0.0 0.0 -c -c], inflate=tolerance)
    setebc!(u,vcat(lx1, lx2),true,1,0.0)
    ly1 = selectnode(fens, box=[-Inf Inf 0.0 0.0 -Inf Inf], inflate=tolerance)
    ly2 = selectnode(fens, box=[-Inf Inf h h -Inf Inf], inflate=tolerance)
    setebc!(u,vcat(ly1, ly2),true,2,0.0)
    applyebc!(u)
    numberdofs!(u)

    
    fi = ForceIntensity(FFlt, 3, getfrc0!);
    el1femm = FEMMBase(IntegData(subset(bfes, sshear0), SimplexRule(2, 3)))
    F1 = distribloads(el1femm, geom, u, fi, 2);
    fi = ForceIntensity(FFlt, 3, getfrcL!);
    el2femm = FEMMBase(IntegData(subset(bfes, sshearL), SimplexRule(2, 3)))
    F2 = distribloads(el2femm, geom, u, fi, 2);

    associategeometry!(femm, geom)
    K = stiffness(femm, geom, u)
    K=cholesky(K)
    U = K\(F1 + F2)
    scattersysvec!(u,U[:])
    
    Tipl = selectnode(fens, box=[L L 0.0 0.0 0. 0.], inflate=tolerance)
    utip = mean(u.values[Tipl, 3])
    println("Deflection: $(utip), compared to $(exactuy(L,0.0))")

    File =  "hughes_cantilever_stresses_T10_by_hand.vtk"
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

function hughes_cantilever_stresses_MST10_by_hand()
    elementtag = "MST10"
    println("""
    Cantilever example.  Hughes 1987. Element: $(elementtag)
    """)
    nu=0.3; # COMPRESSIBLE
    E1=E/(1-nu^2);  # Plane strain
    nu1=nu/(1-nu); I=h*(2*c)^3/12;
    exactux(x,y) = (P/(6*E1*I)*(-y)*(3*(L^2-(L-x)^2)+(2+nu1)*(y^2-c^2)));
    exactuy(x,y) = (P/(6*E1*I)*(((L-x)^3-L^3)-((4+5*nu1)*c^2+3*L^2)*(L-x-L)+3*nu1*(L-x)*y^2));

    n = 2 #
        
    nL = 3*n # number of elements lengthwise
    nc = 2*n # number of elements through the depth
    nh = n # number of elements through the thickness
    xs = collect(linearspace(0.0, L, nL+1))
    ys = collect(linearspace(0.0, h, nh+1))
    zs = collect(linearspace(-c, +c, nc+1))
    fens,fes = T10blockx(xs, ys, zs)
    bfes = meshboundary(fes)
    # end cross-section surface  for the shear loading
    sshearL = selectelem(fens, bfes; facing=true, direction = [+1.0 0.0 0.0])
    # 0 cross-section surface  for the reactions
    sshear0 = selectelem(fens, bfes; facing=true, direction = [-1.0 0.0 0.0])
    
    MR = DeforModelRed3D
    material = MatDeforElastIso(MR, 0.0, E, nu, CTE)
    
    # Material orientation matrix
    csmat = [i==j ? one(FFlt) : zero(FFlt) for i=1:3, j=1:3]
    
    function updatecs!(csmatout::FFltMat, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)
        copyto!(csmatout, csmat)
    end
        
    gr = SimplexRule(3, 4)
    
    femm = FEMMDeforLinearMST10(MR, IntegData(fes, gr), material)
    
    geom = NodalField(fens.xyz)
    u = NodalField(zeros(size(fens.xyz,1), 3)) # displacement field
    
    lx0 = selectnode(fens, box=[0.0 0.0 0.0 0.0 0.0 0.0], inflate=tolerance)
    setebc!(u,lx0,true,1,0.0)
    setebc!(u,lx0,true,2,0.0)
    setebc!(u,lx0,true,3,0.0)
    lx1 = selectnode(fens, box=[0.0 0.0 0.0 0.0 c c], inflate=tolerance)
    lx2 = selectnode(fens, box=[0.0 0.0 0.0 0.0 -c -c], inflate=tolerance)
    setebc!(u,vcat(lx1, lx2),true,1,0.0)
    ly1 = selectnode(fens, box=[-Inf Inf 0.0 0.0 -Inf Inf], inflate=tolerance)
    ly2 = selectnode(fens, box=[-Inf Inf h h -Inf Inf], inflate=tolerance)
    setebc!(u,vcat(ly1, ly2),true,2,0.0)
    applyebc!(u)
    numberdofs!(u)

    
    fi = ForceIntensity(FFlt, 3, getfrc0!);
    el1femm = FEMMBase(IntegData(subset(bfes, sshear0), SimplexRule(2, 3)))
    F1 = distribloads(el1femm, geom, u, fi, 2);
    fi = ForceIntensity(FFlt, 3, getfrcL!);
    el2femm = FEMMBase(IntegData(subset(bfes, sshearL), SimplexRule(2, 3)))
    F2 = distribloads(el2femm, geom, u, fi, 2);

    associategeometry!(femm, geom)
    K =stiffness(femm, geom, u)
    K=cholesky(K)
    U=  K\(F1 + F2)
    scattersysvec!(u,U[:])
    
    Tipl = selectnode(fens, box=[L L 0.0 0.0 0. 0.], inflate=tolerance)
    utip = mean(u.values[Tipl, 3])
    println("Deflection: $(utip), compared to $(exactuy(L,0.0))")
    
    # modeldata["postprocessing"] = FDataDict("file"=>"hughes_cantilever_stresses_$(elementtag)", "outputcsys"=>CSys(3, 3, updatecs!), "quantity"=>:Cauchy, "component"=>[5])
    # modeldata = exportstresselementwise(modeldata)
    
    # modeldata["postprocessing"] = FDataDict("file"=>"hughes_cantilever_stresses_$(elementtag)",
    # "outputcsys"=>CSys(3, 3, updatecs!), "quantity"=>:Cauchy,
    # "component"=>collect(1:6))
    # modeldata = exportstresselementwise(modeldata)
    # stressfields = ElementalField[modeldata["postprocessing"]["exported"][1]["field"]]
    
    true
    
end # hughes_cantilever_stresses_MST10

function hughes_cantilever_stresses_MST10()
    elementtag = "MST10"
    println("""
    Cantilever example.  Hughes 1987. Element: $(elementtag)
    """)
    nu=0.3; # COMPRESSIBLE
    E1=E/(1-nu^2);  # Plane strain
    nu1=nu/(1-nu); I=h*(2*c)^3/12;
    exactux(x,y) = (P/(6*E1*I)*(-y)*(3*(L^2-(L-x)^2)+(2+nu1)*(y^2-c^2)));
    exactuy(x,y) = (P/(6*E1*I)*(((L-x)^3-L^3)-((4+5*nu1)*c^2+3*L^2)*(L-x-L)+3*nu1*(L-x)*y^2));

    modeldatasequence = FDataDict[]
    for n = [1 2 4 8] #
        
        nL = 3*n # number of elements lengthwise
        nc = 2*n # number of elements through the depth
        nh = n # number of elements through the thickness
        xs = collect(linearspace(0.0, L, nL+1))
        ys = collect(linearspace(0.0, h, nh+1))
        zs = collect(linearspace(-c, +c, nc+1))
        fens,fes = T10blockx(xs, ys, zs)
        bfes = meshboundary(fes)
        # end cross-section surface  for the shear loading
        sshearL = selectelem(fens, bfes; facing=true, direction = [+1.0 0.0 0.0])
        # 0 cross-section surface  for the reactions
        sshear0 = selectelem(fens, bfes; facing=true, direction = [-1.0 0.0 0.0])
        
        MR = DeforModelRed3D
        material = MatDeforElastIso(MR, 0.0, E, nu, CTE)
        
        # Material orientation matrix
        csmat = [i==j ? one(FFlt) : zero(FFlt) for i=1:3, j=1:3]
        
        function updatecs!(csmatout::FFltMat, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)
            copyto!(csmatout, csmat)
        end
        
        gr = SimplexRule(3, 4)
        
        region = FDataDict("femm"=>FEMMDeforLinearMST10(MR, IntegData(fes, gr), material))
        
        lx0 = selectnode(fens, box=[0.0 0.0 0.0 0.0 0.0 0.0], inflate=tolerance)
        # println("lx0 = $(lx0)")
        ex01 = FDataDict( "displacement"=>  0.0, "component"=> 1, "node_list"=>lx0 )
        ex02 = FDataDict( "displacement"=>  0.0, "component"=> 2, "node_list"=>lx0 )
        ex03 = FDataDict( "displacement"=>  0.0, "component"=> 3, "node_list"=>lx0 )
        lx1 = selectnode(fens, box=[0.0 0.0 0.0 0.0 c c], inflate=tolerance)
        lx2 = selectnode(fens, box=[0.0 0.0 0.0 0.0 -c -c], inflate=tolerance)
        # println("vcat(lx1, lx2) = $(vcat(lx1, lx2))")
        ex04 = FDataDict( "displacement"=>  0.0, "component"=> 1, "node_list"=>vcat(lx1, lx2) )
        ly1 = selectnode(fens, box=[-Inf Inf 0.0 0.0 -Inf Inf], inflate=tolerance)
        ly2 = selectnode(fens, box=[-Inf Inf h h -Inf Inf], inflate=tolerance)
        # println("vcat(ly1, ly2) = $(vcat(ly1, ly2))")
        ey01 = FDataDict( "displacement"=>  0.0, "component"=> 2, "node_list"=>vcat(ly1, ly2) )
        
        Trac0 = FDataDict("traction_vector"=>getfrc0!, "femm"=>FEMMBase(IntegData(subset(bfes, sshear0), SimplexRule(2, 3))))
        TracL = FDataDict("traction_vector"=>getfrcL!, "femm"=>FEMMBase(IntegData(subset(bfes, sshearL), SimplexRule(2, 3))))
        
        modeldata = FDataDict("fens"=>fens,
        "regions"=>[region],
        "essential_bcs"=>[ex01, ex02, ex03, ex04, ey01],
        "traction_bcs"=>[Trac0, TracL],
        "temperature_change"=>FDataDict("temperature"=>0.0)
        )
        modeldata = linearstatics(modeldata)
        
        u = modeldata["u"]
        geom = modeldata["geom"]
        
        Tipl = selectnode(fens, box=[L L 0.0 0.0 0. 0.], inflate=tolerance)
        utip = mean(u.values[Tipl, 3])
        println("Deflection: $(utip), compared to $(exactuy(L,0.0))")
        
        modeldata["postprocessing"] = FDataDict("file"=>"hughes_cantilever_stresses_$(elementtag)", "outputcsys"=>CSys(3, 3, updatecs!), "quantity"=>:Cauchy, "component"=>[5])
        modeldata = exportstresselementwise(modeldata)
        
        modeldata["postprocessing"] = FDataDict("file"=>"hughes_cantilever_stresses_$(elementtag)",
        "outputcsys"=>CSys(3, 3, updatecs!), "quantity"=>:Cauchy,
        "component"=>collect(1:6))
        modeldata = exportstresselementwise(modeldata)
        stressfields = ElementalField[modeldata["postprocessing"]["exported"][1]["field"]]
        
        modeldata["geometricaltolerance"] = tolerance
        modeldata["elementsize"] = 1.0 / n
        push!(modeldatasequence, modeldata)
    end
    
    filebase = "hughes_cantilever_stresses_$(elementtag)"
    evaluateerrors(filebase, modeldatasequence)
    
    true
    
end # hughes_cantilever_stresses_MST10


function hughes_cantilever_stresses_MST10_incompressible()
    elementtag = "MST10"
    println("""
    Cantilever example.  Hughes 1987. Element: $(elementtag)
    """)
    nu=0.499999999; # INCOMPRESSIBLE
    E1=E/(1-nu^2);  # Plane strain
    nu1=nu/(1-nu); I=h*(2*c)^3/12;
    exactux(x,y) = (P/(6*E1*I)*(-y)*(3*(L^2-(L-x)^2)+(2+nu1)*(y^2-c^2)));
    exactuy(x,y) = (P/(6*E1*I)*(((L-x)^3-L^3)-((4+5*nu1)*c^2+3*L^2)*(L-x-L)+3*nu1*(L-x)*y^2));

    modeldatasequence = FDataDict[]
    for n = [1 2 4 8] #
        
        nL = 3*n # number of elements lengthwise
        nc = 2*n # number of elements through the wwith
        nh = n # number of elements through the thickness
        xs = collect(linearspace(0.0, L, nL+1))
        ys = collect(linearspace(0.0, h, nh+1))
        zs = collect(linearspace(-c, +c, nc+1))
        fens,fes = T10blockx(xs, ys, zs)
        bfes = meshboundary(fes)
        # end cross-section surface  for the shear loading
        sshearL = selectelem(fens, bfes; facing=true, direction = [+1.0 0.0 0.0])
        # 0 cross-section surface  for the reactions
        sshear0 = selectelem(fens, bfes; facing=true, direction = [-1.0 0.0 0.0])
        
        MR = DeforModelRed3D
        material = MatDeforElastIso(MR,  0.0, E, nu, CTE)
        
        # Material orientation matrix
        csmat = [i==j ? one(FFlt) : zero(FFlt) for i=1:3, j=1:3]
        
        function updatecs!(csmatout::FFltMat, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)
            copyto!(csmatout, csmat)
        end
        
        gr = SimplexRule(3, 4)
        
        region = FDataDict("femm"=>FEMMDeforLinearMST10(MR,
        IntegData(fes, gr), material))
        
        lx0 = selectnode(fens, box=[0.0 0.0 0.0 0.0 0.0 0.0], inflate=tolerance)
        # println("lx0 = $(lx0)")
        ex01 = FDataDict( "displacement"=>  0.0, "component"=> 1, "node_list"=>lx0 )
        ex02 = FDataDict( "displacement"=>  0.0, "component"=> 2, "node_list"=>lx0 )
        ex03 = FDataDict( "displacement"=>  0.0, "component"=> 3, "node_list"=>lx0 )
        lx1 = selectnode(fens, box=[0.0 0.0 0.0 0.0 c c], inflate=tolerance)
        lx2 = selectnode(fens, box=[0.0 0.0 0.0 0.0 -c -c], inflate=tolerance)
        # println("vcat(lx1, lx2) = $(vcat(lx1, lx2))")
        ex04 = FDataDict( "displacement"=>  0.0, "component"=> 1, "node_list"=>vcat(lx1, lx2) )
        ly1 = selectnode(fens, box=[-Inf Inf 0.0 0.0 -Inf Inf], inflate=tolerance)
        ly2 = selectnode(fens, box=[-Inf Inf h h -Inf Inf], inflate=tolerance)
        # println("vcat(ly1, ly2) = $(vcat(ly1, ly2))")
        ey01 = FDataDict( "displacement"=>  0.0, "component"=> 2, "node_list"=>vcat(ly1, ly2) )
        
        Trac0 = FDataDict("traction_vector"=>getfrc0!,
        "femm"=>FEMMBase(IntegData(subset(bfes, sshear0), SimplexRule(2, 3))))
        TracL = FDataDict("traction_vector"=>getfrcL!,
        "femm"=>FEMMBase(IntegData(subset(bfes, sshearL), SimplexRule(2, 3))))
        
        modeldata = FDataDict("fens"=>fens,
        "regions"=>[region],
        "essential_bcs"=>[ex01, ex02, ex03, ex04, ey01],
        "traction_bcs"=>[Trac0, TracL],
        "temperature_change"=>FDataDict("temperature"=>0.0)
        )
        modeldata = linearstatics(modeldata)
        
        u = modeldata["u"]
        geom = modeldata["geom"]
        
        Tipl = selectnode(fens, box=[L L 0.0 0.0 0. 0.], inflate=tolerance)
        utip = mean(u.values[Tipl, 3])
        println("Deflection: $(utip), compared to $(exactuy(L,0.0))")
        
        modeldata["postprocessing"] = FDataDict("file"=>"hughes_cantilever_stresses_$(elementtag)",
        "outputcsys"=>CSys(3, 3, updatecs!), "quantity"=>:Cauchy,
        "component"=>[5])
        modeldata = exportstresselementwise(modeldata)
        
        modeldata["postprocessing"] = FDataDict("file"=>"hughes_cantilever_stresses_$(elementtag)",
        "outputcsys"=>CSys(3, 3, updatecs!), "quantity"=>:Cauchy,
        "component"=>collect(1:6))
        modeldata = exportstresselementwise(modeldata)
        stressfields = ElementalField[modeldata["postprocessing"]["exported"][1]["field"]]
        
        modeldata["geometricaltolerance"] = tolerance
        modeldata["elementsize"] = 1.0 / n
        push!(modeldatasequence, modeldata)
    end
    
    filebase = "hughes_cantilever_stresses_incompressible_$(elementtag)"
    evaluateerrors(filebase, modeldatasequence)
    
    true
    
end # hughes_cantilever_stresses_MST10_incompressible


function hughes_cantilever_stresses_nodal_MST10()
    elementtag = "MST10"
    println("""
    Cantilever example.  Hughes 1987. Element: $(elementtag)
    """)
    nu=0.3; # COMPRESSIBLE
    E1=E/(1-nu^2);  # Plane strain
    nu1=nu/(1-nu); I=h*(2*c)^3/12;
    exactux(x,y) = (P/(6*E1*I)*(-y)*(3*(L^2-(L-x)^2)+(2+nu1)*(y^2-c^2)));
    exactuy(x,y) = (P/(6*E1*I)*(((L-x)^3-L^3)-((4+5*nu1)*c^2+3*L^2)*(L-x-L)+3*nu1*(L-x)*y^2));

    modeldatasequence = FDataDict[]
    for n = [1 2 4 8] #
        
        nL = 3*n # number of elements lengthwise
        nc = 2*n # number of elements through the wwith
        nh = n # number of elements through the thickness
        xs = collect(linearspace(0.0, L, nL+1))
        ys = collect(linearspace(0.0, h, nh+1))
        zs = collect(linearspace(-c, +c, nc+1))
        fens,fes = T10blockx(xs, ys, zs)
        bfes = meshboundary(fes)
        # end cross-section surface  for the shear loading
        sshearL = selectelem(fens, bfes; facing=true, direction = [+1.0 0.0 0.0])
        # 0 cross-section surface  for the reactions
        sshear0 = selectelem(fens, bfes; facing=true, direction = [-1.0 0.0 0.0])
        
        MR = DeforModelRed3D
        material = MatDeforElastIso(MR,   0.0, E, nu, CTE)
        
        # Material orientation matrix
        csmat = [i==j ? one(FFlt) : zero(FFlt) for i=1:3, j=1:3]
        
        function updatecs!(csmatout::FFltMat, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)
            copyto!(csmatout, csmat)
        end
        
        gr = SimplexRule(3, 4) # rule for tetrahedral meshes
        
        region = FDataDict("femm"=>FEMMDeforLinearMST10(MR,
        IntegData(fes, gr), material))
        
        lx0 = selectnode(fens, box=[0.0 0.0 0.0 0.0 0.0 0.0], inflate=tolerance)
        # println("lx0 = $(lx0)")
        ex01 = FDataDict( "displacement"=>  0.0, "component"=> 1, "node_list"=>lx0 )
        ex02 = FDataDict( "displacement"=>  0.0, "component"=> 2, "node_list"=>lx0 )
        ex03 = FDataDict( "displacement"=>  0.0, "component"=> 3, "node_list"=>lx0 )
        lx1 = selectnode(fens, box=[0.0 0.0 0.0 0.0 c c], inflate=tolerance)
        lx2 = selectnode(fens, box=[0.0 0.0 0.0 0.0 -c -c], inflate=tolerance)
        # println("vcat(lx1, lx2) = $(vcat(lx1, lx2))")
        ex04 = FDataDict( "displacement"=>  0.0, "component"=> 1, "node_list"=>vcat(lx1, lx2) )
        ly1 = selectnode(fens, box=[-Inf Inf 0.0 0.0 -Inf Inf], inflate=tolerance)
        ly2 = selectnode(fens, box=[-Inf Inf h h -Inf Inf], inflate=tolerance)
        # println("vcat(ly1, ly2) = $(vcat(ly1, ly2))")
        ey01 = FDataDict( "displacement"=>  0.0, "component"=> 2, "node_list"=>vcat(ly1, ly2) )
        
        
        Trac0 = FDataDict("traction_vector"=>getfrc0!,
        "femm"=>FEMMBase(IntegData(subset(bfes, sshear0), SimplexRule(2, 3))))
        TracL = FDataDict("traction_vector"=>getfrcL!,
        "femm"=>FEMMBase(IntegData(subset(bfes, sshearL), SimplexRule(2, 3))))
        
        modeldata = FDataDict("fens"=>fens,
        "regions"=>[region],
        "essential_bcs"=>[ex01, ex02, ex03, ex04, ey01],
        "traction_bcs"=>[Trac0, TracL],
        "temperature_change"=>FDataDict("temperature"=>0.0)
        )
        modeldata = linearstatics(modeldata)
        
        u = modeldata["u"]
        geom = modeldata["geom"]
        
        Tipl = selectnode(fens, box=[L L 0.0 0.0 0. 0.], inflate=tolerance)
        utip = mean(u.values[Tipl, 3])
        println("Deflection: $(utip), compared to $(exactuy(L,0.0))")
        
        modeldata["postprocessing"] = FDataDict("file"=>"hughes_cantilever_stresses_nodal_$(elementtag)",
        "outputcsys"=>CSys(3, 3, updatecs!), "quantity"=>:Cauchy,
        "component"=>collect(1:6), "nodevalmethod"=>:averaging, "reportat"=>:extraptrend)
        modeldata = exportstress(modeldata)
        stressfields = NodalField[modeldata["postprocessing"]["exported"][1]["field"]]
        
        modeldata["geometricaltolerance"] = tolerance
        modeldata["elementsize"] = 1.0 / n
        push!(modeldatasequence, modeldata)
    end
    
    filebase = "hughes_cantilever_stresses_nodal_$(elementtag)"
    evaluateerrors(filebase, modeldatasequence)
    
    true
    
end # hughes_cantilever_stresses_nodal_MST10


function hughes_cantilever_stresses_nodal_T10()
    elementtag = "T10"
    println("""
    Cantilever example.  Hughes 1987. Element: $(elementtag)
    """)
    nu=0.3; # COMPRESSIBLE
    E1=E/(1-nu^2);  # Plane strain
    nu1=nu/(1-nu); I=h*(2*c)^3/12;
    exactux(x,y) = (P/(6*E1*I)*(-y)*(3*(L^2-(L-x)^2)+(2+nu1)*(y^2-c^2)));
    exactuy(x,y) = (P/(6*E1*I)*(((L-x)^3-L^3)-((4+5*nu1)*c^2+3*L^2)*(L-x-L)+3*nu1*(L-x)*y^2));

    modeldatasequence = FDataDict[]
    for n = [1 2 4 8] #
        
        nL = 3*n # number of elements lengthwise
        nc = 2*n # number of elements through the wwith
        nh = n # number of elements through the thickness
        xs = collect(linearspace(0.0, L, nL+1))
        ys = collect(linearspace(0.0, h, nh+1))
        zs = collect(linearspace(-c, +c, nc+1))
        fens,fes = T10blockx(xs, ys, zs)
        bfes = meshboundary(fes)
        # end cross-section surface  for the shear loading
        sshearL = selectelem(fens, bfes; facing=true, direction = [+1.0 0.0 0.0])
        # 0 cross-section surface  for the reactions
        sshear0 = selectelem(fens, bfes; facing=true, direction = [-1.0 0.0 0.0])
        
        MR = DeforModelRed3D
        material = MatDeforElastIso(MR,
        0.0, E, nu, CTE)
        
        # Material orientation matrix
        csmat = [i==j ? 1.0 : 0.0 for i=1:3, j=1:3]
        
        function updatecs!(csmatout::FFltMat, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)
            copyto!(csmatout, csmat)
        end
        
        gr = SimplexRule(3, 4)
        
        region = FDataDict("femm"=>FEMMDeforLinear(MR,
        IntegData(fes, gr), material))
        
        lx0 = selectnode(fens, box=[0.0 0.0 0.0 0.0 0.0 0.0], inflate=tolerance)
        # println("lx0 = $(lx0)")
        ex01 = FDataDict( "displacement"=>  0.0, "component"=> 1, "node_list"=>lx0 )
        ex02 = FDataDict( "displacement"=>  0.0, "component"=> 2, "node_list"=>lx0 )
        ex03 = FDataDict( "displacement"=>  0.0, "component"=> 3, "node_list"=>lx0 )
        lx1 = selectnode(fens, box=[0.0 0.0 0.0 0.0 c c], inflate=tolerance)
        lx2 = selectnode(fens, box=[0.0 0.0 0.0 0.0 -c -c], inflate=tolerance)
        # println("vcat(lx1, lx2) = $(vcat(lx1, lx2))")
        ex04 = FDataDict( "displacement"=>  0.0, "component"=> 1, "node_list"=>vcat(lx1, lx2) )
        ly1 = selectnode(fens, box=[-Inf Inf 0.0 0.0 -Inf Inf], inflate=tolerance)
        ly2 = selectnode(fens, box=[-Inf Inf h h -Inf Inf], inflate=tolerance)
        # println("vcat(ly1, ly2) = $(vcat(ly1, ly2))")
        ey01 = FDataDict( "displacement"=>  0.0, "component"=> 2, "node_list"=>vcat(ly1, ly2) )
        
        Trac0 = FDataDict("traction_vector"=>getfrc0!,
        "femm"=>FEMMBase(IntegData(subset(bfes, sshear0), SimplexRule(2, 3))))
        TracL = FDataDict("traction_vector"=>getfrcL!,
        "femm"=>FEMMBase(IntegData(subset(bfes, sshearL), SimplexRule(2, 3))))
        
        modeldata = FDataDict("fens"=>fens,
        "regions"=>[region],
        "essential_bcs"=>[ex01, ex02, ex03, ex04, ey01],
        "traction_bcs"=>[Trac0, TracL],
        "temperature_change"=>FDataDict("temperature"=>0.0)
        )
        modeldata = linearstatics(modeldata)
        
        u = modeldata["u"]
        geom = modeldata["geom"]
        
        Tipl = selectnode(fens, box=[L L 0.0 0.0 0. 0.], inflate=tolerance)
        utip = mean(u.values[Tipl, 3])
        println("Deflection: $(utip), compared to $(exactuy(L,0.0))")
        
        modeldata["postprocessing"] = FDataDict("file"=>"hughes_cantilever_stresses_nodal_$(elementtag)",
        "outputcsys"=>CSys(3, 3, updatecs!), "quantity"=>:Cauchy,
        "component"=>collect(1:6))
        modeldata = exportstress(modeldata)
        stressfields = NodalField[modeldata["postprocessing"]["exported"][1]["field"]]
        
        modeldata["geometricaltolerance"] = tolerance
        modeldata["elementsize"] = 1.0 / n
        push!(modeldatasequence, modeldata)
    end
    
    filebase = "hughes_cantilever_stresses_nodal_$(elementtag)"
    evaluateerrors(filebase, modeldatasequence)
    
    true
    
end # hughes_cantilever_stresses_nodal_T10


function hughes_cantilever_stresses_T10()
    elementtag = "T10"
    println("""
    Cantilever example.  Hughes 1987. Element: $(elementtag)
    """)
    nu=0.3; # COMPRESSIBLE
    E1=E/(1-nu^2);  # Plane strain
    nu1=nu/(1-nu); I=h*(2*c)^3/12;
    exactux(x,y) = (P/(6*E1*I)*(-y)*(3*(L^2-(L-x)^2)+(2+nu1)*(y^2-c^2)));
    exactuy(x,y) = (P/(6*E1*I)*(((L-x)^3-L^3)-((4+5*nu1)*c^2+3*L^2)*(L-x-L)+3*nu1*(L-x)*y^2));

    modeldatasequence = FDataDict[]
    for n = [1 2 4 8] #
        
        nL = 3*n # number of elements lengthwise
        nc = 2*n # number of elements through the wwith
        nh = n # number of elements through the thickness
        xs = collect(linearspace(0.0, L, nL+1))
        ys = collect(linearspace(0.0, h, nh+1))
        zs = collect(linearspace(-c, +c, nc+1))
        fens,fes = T10blockx(xs, ys, zs)
        bfes = meshboundary(fes)
        # end cross-section surface  for the shear loading
        sshearL = selectelem(fens, bfes; facing=true, direction = [+1.0 0.0 0.0])
        # 0 cross-section surface  for the reactions
        sshear0 = selectelem(fens, bfes; facing=true, direction = [-1.0 0.0 0.0])
        
        MR = DeforModelRed3D
        material = MatDeforElastIso(MR,
        0.0, E, nu, CTE)
        
        # Material orientation matrix
        csmat = [i==j ? one(FFlt) : zero(FFlt) for i=1:3, j=1:3]
        
        function updatecs!(csmatout::FFltMat, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)
            copyto!(csmatout, csmat)
        end
        
        gr = SimplexRule(3, 4)
        
        region = FDataDict("femm"=>FEMMDeforLinear(MR, IntegData(fes, gr), material))
        
        lx0 = selectnode(fens, box=[0.0 0.0 0.0 0.0 0.0 0.0], inflate=tolerance)
        # println("lx0 = $(lx0)")
        ex01 = FDataDict( "displacement"=>  0.0, "component"=> 1, "node_list"=>lx0 )
        ex02 = FDataDict( "displacement"=>  0.0, "component"=> 2, "node_list"=>lx0 )
        ex03 = FDataDict( "displacement"=>  0.0, "component"=> 3, "node_list"=>lx0 )
        lx1 = selectnode(fens, box=[0.0 0.0 0.0 0.0 c c], inflate=tolerance)
        lx2 = selectnode(fens, box=[0.0 0.0 0.0 0.0 -c -c], inflate=tolerance)
        # println("vcat(lx1, lx2) = $(vcat(lx1, lx2))")
        ex04 = FDataDict( "displacement"=>  0.0, "component"=> 1, "node_list"=>vcat(lx1, lx2) )
        ly1 = selectnode(fens, box=[-Inf Inf 0.0 0.0 -Inf Inf], inflate=tolerance)
        ly2 = selectnode(fens, box=[-Inf Inf h h -Inf Inf], inflate=tolerance)
        # println("vcat(ly1, ly2) = $(vcat(ly1, ly2))")
        ey01 = FDataDict( "displacement"=>  0.0, "component"=> 2, "node_list"=>vcat(ly1, ly2) )
        
        Trac0 = FDataDict("traction_vector"=>getfrc0!,
        "femm"=>FEMMBase(IntegData(subset(bfes, sshear0), SimplexRule(2, 3))))
        TracL = FDataDict("traction_vector"=>getfrcL!,
        "femm"=>FEMMBase(IntegData(subset(bfes, sshearL), SimplexRule(2, 3))))
        
        modeldata = FDataDict("fens"=>fens,
        "regions"=>[region],
        "essential_bcs"=>[ex01, ex02, ex03, ex04, ey01],
        "traction_bcs"=>[Trac0, TracL],
        "temperature_change"=>FDataDict("temperature"=>0.0)
        )
        modeldata = linearstatics(modeldata)
        
        u = modeldata["u"]
        geom = modeldata["geom"]
        
        Tipl = selectnode(fens, box=[L L 0.0 0.0 0. 0.], inflate=tolerance)
        utip = mean(u.values[Tipl, 3])
        println("Deflection: $(utip), compared to $(exactuy(L,0.0))")
        
        modeldata["postprocessing"] = FDataDict("file"=>"hughes_cantilever_stresses_$(elementtag)",
        "outputcsys"=>CSys(3, 3, updatecs!), "quantity"=>:Cauchy,
        "component"=>[5])
        modeldata = exportstresselementwise(modeldata)
        
        modeldata["postprocessing"] = FDataDict("file"=>"hughes_cantilever_stresses_$(elementtag)",
        "outputcsys"=>CSys(3, 3, updatecs!), "quantity"=>:Cauchy,
        "component"=>collect(1:6))
        modeldata = exportstresselementwise(modeldata)
        stressfields = ElementalField[modeldata["postprocessing"]["exported"][1]["field"]]
        
        modeldata["geometricaltolerance"] = tolerance
        modeldata["elementsize"] = 1.0 / n
        push!(modeldatasequence, modeldata)
    end
    
    filebase = "hughes_cantilever_stresses_$(elementtag)"
    evaluateerrors(filebase, modeldatasequence)
    
    true
    
end # hughes_cantilever_stresses_T10


function hughes_cantilever_stresses_T10_incompressible()
    elementtag = "T10"
    println("""
    Cantilever example.  Hughes 1987. Element: $(elementtag)
    """)
        nu=0.499999999; # INCOMPRESSIBLE
        E1=E/(1-nu^2);  # Plane strain
        nu1=nu/(1-nu); I=h*(2*c)^3/12;
        exactux(x,y) = (P/(6*E1*I)*(-y)*(3*(L^2-(L-x)^2)+(2+nu1)*(y^2-c^2)));
        exactuy(x,y) = (P/(6*E1*I)*(((L-x)^3-L^3)-((4+5*nu1)*c^2+3*L^2)*(L-x-L)+3*nu1*(L-x)*y^2));

    modeldatasequence = FDataDict[]
    for n = [1 2 4 8] #
        
        nL = 3*n # number of elements lengthwise
        nc = 2*n # number of elements through the wwith
        nh = n # number of elements through the thickness
        xs = collect(linearspace(0.0, L, nL+1))
        ys = collect(linearspace(0.0, h, nh+1))
        zs = collect(linearspace(-c, +c, nc+1))
        fens,fes = T10blockx(xs, ys, zs)
        bfes = meshboundary(fes)
        # end cross-section surface  for the shear loading
        sshearL = selectelem(fens, bfes; facing=true, direction = [+1.0 0.0 0.0])
        # 0 cross-section surface  for the reactions
        sshear0 = selectelem(fens, bfes; facing=true, direction = [-1.0 0.0 0.0])
        
        MR = DeforModelRed3D
        material = MatDeforElastIso(MR,  0.0, E, nu, CTE)
        
        # Material orientation matrix
        csmat = [i==j ? one(FFlt) : zero(FFlt) for i=1:3, j=1:3]
        
        function updatecs!(csmatout::FFltMat, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)
            copyto!(csmatout, csmat)
        end
        
        gr = SimplexRule(3, 4)
        
        region = FDataDict("femm"=>FEMMDeforLinear(MR,
        IntegData(fes, gr), material))
        
        lx0 = selectnode(fens, box=[0.0 0.0 0.0 0.0 0.0 0.0], inflate=tolerance)
        # println("lx0 = $(lx0)")
        ex01 = FDataDict( "displacement"=>  0.0, "component"=> 1, "node_list"=>lx0 )
        ex02 = FDataDict( "displacement"=>  0.0, "component"=> 2, "node_list"=>lx0 )
        ex03 = FDataDict( "displacement"=>  0.0, "component"=> 3, "node_list"=>lx0 )
        lx1 = selectnode(fens, box=[0.0 0.0 0.0 0.0 c c], inflate=tolerance)
        lx2 = selectnode(fens, box=[0.0 0.0 0.0 0.0 -c -c], inflate=tolerance)
        # println("vcat(lx1, lx2) = $(vcat(lx1, lx2))")
        ex04 = FDataDict( "displacement"=>  0.0, "component"=> 1, "node_list"=>vcat(lx1, lx2) )
        ly1 = selectnode(fens, box=[-Inf Inf 0.0 0.0 -Inf Inf], inflate=tolerance)
        ly2 = selectnode(fens, box=[-Inf Inf h h -Inf Inf], inflate=tolerance)
        # println("vcat(ly1, ly2) = $(vcat(ly1, ly2))")
        ey01 = FDataDict( "displacement"=>  0.0, "component"=> 2, "node_list"=>vcat(ly1, ly2) )
        
        Trac0 = FDataDict("traction_vector"=>getfrc0!,
        "femm"=>FEMMBase(IntegData(subset(bfes, sshear0), SimplexRule(2, 3))))
        TracL = FDataDict("traction_vector"=>getfrcL!,
        "femm"=>FEMMBase(IntegData(subset(bfes, sshearL), SimplexRule(2, 3))))
        
        modeldata = FDataDict("fens"=>fens,
        "regions"=>[region],
        "essential_bcs"=>[ex01, ex02, ex03, ex04, ey01],
        "traction_bcs"=>[Trac0, TracL],
        "temperature_change"=>FDataDict("temperature"=>0.0)
        )
        modeldata = linearstatics(modeldata)
        
        u = modeldata["u"]
        geom = modeldata["geom"]
        
        Tipl = selectnode(fens, box=[L L 0.0 0.0 0. 0.], inflate=tolerance)
        utip = mean(u.values[Tipl, 3])
        println("Deflection: $(utip), compared to $(exactuy(L,0.0))")
        
        modeldata["postprocessing"] = FDataDict("file"=>"hughes_cantilever_stresses_$(elementtag)",
        "outputcsys"=>CSys(3, 3, updatecs!), "quantity"=>:Cauchy,
        "component"=>[5])
        modeldata = exportstresselementwise(modeldata)
        
        modeldata["postprocessing"] = FDataDict("file"=>"hughes_cantilever_stresses_$(elementtag)",
        "outputcsys"=>CSys(3, 3, updatecs!), "quantity"=>:Cauchy,
        "component"=>collect(1:6))
        modeldata = exportstresselementwise(modeldata)
        stressfields = ElementalField[modeldata["postprocessing"]["exported"][1]["field"]]
        
        modeldata["geometricaltolerance"] = tolerance
        modeldata["elementsize"] = 1.0 / n
        push!(modeldatasequence, modeldata)
    end
    
    filebase = "hughes_cantilever_stresses_incompressible_$(elementtag)"
    evaluateerrors(filebase, modeldatasequence)
    
    true
    
end # hughes_cantilever_stresses_T10_incompressible

function allrun()
    println("#####################################################") 
    println("# hughes_cantilever_stresses_H8_by_hand ")
    hughes_cantilever_stresses_H8_by_hand()
    println("#####################################################") 
    println("# hughes_cantilever_stresses_H20_by_hand ")
    hughes_cantilever_stresses_H20_by_hand()
    println("#####################################################") 
    println("# hughes_cantilever_stresses_T10_by_hand ")
    hughes_cantilever_stresses_T10_by_hand()
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
