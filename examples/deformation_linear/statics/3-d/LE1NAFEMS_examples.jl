module LE1NAFEMS_examples
using FinEtools
using FinEtools.MeshExportModule
using FinEtools.MeshImportModule: import_ABAQUS


# function LE1NAFEMS_compare_meshes()
#     Thick0 = 0.1*phun("m")/2.0 # to account for the symmetry reduction
    
#     ref  = 0
#     Thickness = Thick0
#     tolerance = Thickness/2^ref/300.; # Geometrical tolerance
    
    
#     fens,fes = T10block(1.0, pi/2, Thickness, 2^ref*5, 2^ref*6, 1; orientation = :b)
#     for i=1:count(fens)
#         t=fens.xyz[i,1]; a=fens.xyz[i,2]; z=fens.xyz[i,3]
#         fens.xyz[i,:]=[(t*3.25+(1-t)*2)*cos(a), (t*2.75+(1-t)*1)*sin(a), z];
#     end
#     println("$((count(fens), count(fes)))")
    
#     output = import_ABAQUS("LE1AbaqusExport-C3D10HS-5-6-1.inp")
#     fens1, fes1 = output["fens"], output["fesets"][1]
#     println("$((count(fens1), count(fes1[1])))")
    
#     fens, newfes1, fes2 = mergemeshes(fens,fes, fens1,fes1[1], tolerance)
#     # fes = cat(fes2, newfes1)
#     # println("$((count(fens), count(fes)))")
    
#     File =  "a1.vtk"
#     vtkexportmesh(File, newfes1.conn, fens.xyz,
#     FinEtools.MeshExportModule.T10)
#     @async run(`"paraview.exe" $File`)
#     File =  "a2.vtk"
#     vtkexportmesh(File, fes2.conn, fens.xyz,
#     FinEtools.MeshExportModule.T10)
#     @async run(`"paraview.exe" $File`)
    
    
#     #
#     # fens,fes = H8block(1.0, pi/2, Thickness, 2^ref*5, 2^ref*6, 1)
#     # for i=1:count(fens)
#     #     t=fens.xyz[i,1]; a=fens.xyz[i,2]; z=fens.xyz[i,3]
#     #     fens.xyz[i,:]=[(t*3.25+(1-t)*2)*cos(a), (t*2.75+(1-t)*1)*sin(a), z];
#     # end
#     # println("$((count(fens), count(fes)))")
#     #
#     # fens1,fes1 = import_ABAQUS("LE1AbaqusExport-C3D8S-80-96-1.inp")
#     # println("$((count(fens1), count(fes1[1])))")
#     #
#     #  fens, newfes1, fes2 = mergemeshes(fens,fes, fens1,fes1[1], tolerance)
#     #  fes = cat(fes2, newfes1)
#     #  println("$((count(fens), count(fes)))")
#     #
#     # File =  "a.vtk"
#     # vtkexportmesh(File, fes.conn, fens.xyz,
#     #                FinEtools.MeshExportModule.H8)
#     # @async run(`"paraview.exe" $File`)
    
    
#     true
    
# end # LE1NAFEMS_compare_meshes


function LE1NAFEMS_MSH8()
    println("LE1NAFEMS, 3D version."        )
    t0 = time()
    
    E = 210e3*phun("MEGA*PA");# 210e3 MPa
    nu = 0.3;
    p = 10*phun("MEGA*PA");# 10 MPA Outward pressure on the outside ellipse
    sigma_yD = 92.7*phun("MEGA*PA");# tensile stress at [2.0, 0.0] meters
    Radius = 1.0*phun("m")
    Thickness = 0.1*phun("m")
    n = 2; # number of elements per side
    tolerance = 1.0/n/1000.; # Geometrical tolerance
    
    fens,fes = Q4block(1.0, pi/2, n, n*2)
    fens,fes  = H8extrudeQ4(fens, fes,
    1, (xyz, layer)->[xyz[1], xyz[2], (layer)*Thickness]);
    
    bdryfes = meshboundary(fes);
    icl = selectelem(fens, bdryfes, box=[1.0, 1.0, 0.0, pi/2, 0.0, Thickness], inflate=tolerance);
    for i=1:count(fens)
        t=fens.xyz[i,1]; a=fens.xyz[i,2]; z=fens.xyz[i,3]
        fens.xyz[i,:]=[(t*3.25+(1-t)*2)*cos(a), (t*2.75+(1-t)*1)*sin(a), z];
    end
    
    
    geom = NodalField(fens.xyz)
    u = NodalField(zeros(size(fens.xyz,1),3)) # displacement field
    
    l1 =selectnode(fens; box=[0.0, Inf, 0.0, 0.0, 0.0, Thickness], inflate = tolerance)
    setebc!(u,l1,true, 2, 0.0)
    l1 =selectnode(fens; box=[0.0, 0.0, 0.0, Inf, 0.0, Thickness], inflate = tolerance)
    setebc!(u,l1,true, 1, 0.0)
    l1 =selectnode(fens; box=[0.0, Inf, 0.0, Inf, 0.0, 0.0], inflate = tolerance)
    setebc!(u,l1,true, 3, 0.0)
    
    applyebc!(u)
    numberdofs!(u)
    
    
    el1femm =  FEMMBase(IntegData(subset(bdryfes,icl), GaussRule(2, 2)))
    function pfun(forceout::FVec{T}, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt) where {T}
        pt= [2.75/3.25*XYZ[1], 3.25/2.75*XYZ[2], 0.0]
        forceout .=    vec(p*pt/norm(pt));
        return forceout
    end
    fi = ForceIntensity(FFlt, 3, pfun);
    F2 = distribloads(el1femm, geom, u, fi, 2);
    
    # Note that the material object needs to be created with the proper
    # model-dimension reduction in mind.  In this case that is the fully three-dimensional solid.
    MR = DeforModelRed3D
    
    material = MatDeforElastIso(MR, E, nu)
    
    femm = FEMMDeforLinearMSH8(MR, IntegData(fes, GaussRule(3, 2)), material)
    
    # The geometry field now needs to be associated with the FEMM
    femm = associategeometry!(femm, geom)
    
    K = stiffness(femm, geom, u)
    K = cholesky(K)
    U = K\(F2)
    scattersysvec!(u, U[:])
    
    nl = selectnode(fens, box=[2.0, 2.0, 0.0, 0.0, 0.0, 0.0],inflate=tolerance);
    thecorneru = zeros(FFlt,1,3)
    gathervalues_asmat!(u, thecorneru, nl);
    thecorneru = thecorneru/phun("mm")
    println("$(time()-t0) [s];  displacement =$(thecorneru) [MM] as compared to reference [-0.10215,0] [MM]")
    
    
    fld= fieldfromintegpoints(femm, geom, u, :Cauchy, 2;
    nodevalmethod = :averaging, reportat = :meanonly)
    println("Sigma_y =$(fld.values[nl,1][1]/phun("MPa")) as compared to reference sigma_yD = $(sigma_yD/phun("MPa")) [MPa]")
    fld= fieldfromintegpoints(femm, geom, u, :Cauchy, 2;
    nodevalmethod = :averaging, reportat = :extrapmean)
    println("Sigma_y =$(fld.values[nl,1][1]/phun("MPa")) as compared to reference sigma_yD = $(sigma_yD/phun("MPa")) [MPa]")
    fld= fieldfromintegpoints(femm, geom, u, :Cauchy, 2;
    nodevalmethod = :averaging, reportat = :extraptrend)
    println("Sigma_y =$(fld.values[nl,1][1]/phun("MPa")) as compared to reference sigma_yD = $(sigma_yD/phun("MPa")) [MPa]")
    
    fld= fieldfromintegpoints(femm, geom, u, :Cauchy, 2; nodevalmethod = :invdistance, reportat = :meanonly)
    println("Sigma_y =$(fld.values[nl,1][1]/phun("MPa")) as compared to reference sigma_yD = $(sigma_yD/phun("MPa")) [MPa]")
    fld= fieldfromintegpoints(femm, geom, u, :Cauchy, 2; nodevalmethod = :averaging, reportat = :extrapmean)
    println("Sigma_y =$(fld.values[nl,1][1]/phun("MPa")) as compared to reference sigma_yD = $(sigma_yD/phun("MPa")) [MPa]")
    fld= fieldfromintegpoints(femm, geom, u, :Cauchy, 2; nodevalmethod = :averaging, reportat = :extraptrend)
    println("Sigma_y =$(fld.values[nl,1][1]/phun("MPa")) as compared to reference sigma_yD = $(sigma_yD/phun("MPa")) [MPa]")
    
    println("$(n), $(fld.values[nl,1][1]/phun("MPa"))")
    
    File =  "a.vtk"
    vtkexportmesh(File, fes.conn, geom.values,
    FinEtools.MeshExportModule.H8; vectors=[("u", u.values)],
    scalars=[("sigmay", fld.values)])
    @async run(`"paraview.exe" $File`)
    true
    
end # LE1NAFEMS_MSH8


function LE1NAFEMS_MSH8_convergence()
    # Example from the "Improved Stress Recovery for Mean-strain Finite Elements" paper by Sivapuram and Krysl, 2017
    println("LE1NAFEMS, 3D version. MSH8"        )
    t0 = time()
    
    E = 210e3*phun("MEGA*PA");# 210e3 MPa
    nu = 0.3;
    p = 10*phun("MEGA*PA");# 10 MPA Outward pressure on the outside ellipse
    sigma_yD = 92.7*phun("MEGA*PA");# tensile stress at [2.0, 0.0] meters
    Thick0 = 0.1*phun("m")/2.0 # to account for the symmetry reduction
    
    sigyderrs = Dict{Symbol, FFltVec}()
    
    nnodes = []
    for extrapolation in [:extraptrend :extrapmean]
        sigyderrs[extrapolation] = FFltVec[]
        nnodes = []
        for ref in 0:1:4
            # Thickness = Thick0
            Thickness = Thick0/2^ref
            tolerance = Thickness/2^ref/1000.; # Geometrical tolerance
            
            fens,fes = H8block(1.0, pi/2, Thickness, 2^ref*5, 2^ref*6, 1)
            
            bdryfes = meshboundary(fes);
            icl = selectelem(fens, bdryfes, box=[1.0, 1.0, 0.0, pi/2, 0.0, Thickness], inflate=tolerance);
            for i=1:count(fens)
                t=fens.xyz[i,1]; a=fens.xyz[i,2]; z=fens.xyz[i,3]
                fens.xyz[i,:]=[(t*3.25+(1-t)*2)*cos(a), (t*2.75+(1-t)*1)*sin(a), z];
            end
            
            
            geom = NodalField(fens.xyz)
            u = NodalField(zeros(size(fens.xyz,1),3)) # displacement field
            
            l1 =selectnode(fens; box=[0.0, Inf, 0.0, 0.0, 0.0, Thickness], inflate = tolerance)
            setebc!(u,l1,true, 2, 0.0)
            l1 =selectnode(fens; box=[0.0, 0.0, 0.0, Inf, 0.0, Thickness], inflate = tolerance)
            setebc!(u,l1,true, 1, 0.0)
            l1 =selectnode(fens; box=[0.0, Inf, 0.0, Inf, 0.0, 0.0], inflate = tolerance)
            setebc!(u,l1,true, 3, 0.0)
            
            applyebc!(u)
            numberdofs!(u)
            
            
            el1femm =  FEMMBase(IntegData(subset(bdryfes,icl), GaussRule(2, 2)))
            function pfun(forceout::FVec{T}, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt) where {T}
                pt= [2.75/3.25*XYZ[1], 3.25/2.75*XYZ[2], 0.0]
                forceout .=    vec(p*pt/norm(pt));
                return forceout
            end
            fi = ForceIntensity(FFlt, 3, pfun);
            F2 = distribloads(el1femm, geom, u, fi, 2);
            
            # Note that the material object needs to be created with the proper
            # model-dimension reduction in mind.  In this case that is the fully three-dimensional solid.
            MR = DeforModelRed3D
            
            material = MatDeforElastIso(MR, E, nu)
            
            femm = FEMMDeforLinearMSH8(MR, IntegData(fes, GaussRule(3, 2)), material)
            
            # The geometry field now needs to be associated with the FEMM
            femm = associategeometry!(femm, geom)
            
            K = stiffness(femm, geom, u)
            K = cholesky(K)
            U = K\(F2)
            scattersysvec!(u, U[:])
            
            nl = selectnode(fens, box=[2.0, 2.0, 0.0, 0.0, 0.0, Thickness], inflate=tolerance);
            thecorneru = zeros(FFlt,length(nl),3)
            gathervalues_asmat!(u, thecorneru, nl);
            thecorneru = mean(thecorneru, 1)[1]/phun("mm")
            println("displacement =$(thecorneru) vs -0.10215 [MM]")
            
            fld = fieldfromintegpoints(femm, geom, u, :Cauchy, 2; nodevalmethod = :averaging, reportat = extrapolation)
            sigyd = mean(fld.values[nl,1], 1)[1]
            println("Sigma_y =$(sigyd/phun("MPa")) vs $(sigma_yD/phun("MPa")) [MPa]")
            
            println("$extrapolation, $(count(fes)), $(sigyd/phun("MPa"))")
            push!(nnodes, count(fes))
            push!(sigyderrs[extrapolation], (sigyd/sigma_yD - 1.0))
            # File =  "a.vtk"
            # vtkexportmesh(File, fes.conn, geom.values,
            # FinEtools.MeshExportModule.H8; vectors=[("u", u.values)],
            # scalars=[("sigmay", fld.values)])
            # @async run(`"paraview.exe" $File`)
        end
    end
    
    File = "LE1NAFEMS_MSH8_convergence.CSV"
    savecsv(File, nnodes=vec(nnodes), sigyderrtrend=vec(sigyderrs[:extraptrend]), sigyderrmean=vec(sigyderrs[:extrapmean]))
    
end # LE1NAFEMS_MSH8_convergence


function LE1NAFEMS_MSH8_export()
    println("LE1NAFEMS, 3D version."        )
    t0 = time()
    
    E = 210e3*phun("MEGA*PA");# 210e3 MPa
    nu = 0.3;
    p = 10*phun("MEGA*PA");# 10 MPA Outward pressure on the outside ellipse
    sigma_yD = 92.7*phun("MEGA*PA");# tensile stress at [2.0, 0.0] meters
    Radius = 1.0*phun("m")
    Thickness = 0.1*phun("m")
    n = 64; # number of elements per side
    tolerance = 1.0/n/1000.; # Geometrical tolerance
    
    fens,fes = Q4block(1.0, pi/2, n, n*2)
    fens,fes  = H8extrudeQ4(fens, fes,
    1, (xyz, layer)->[xyz[1], xyz[2], (layer)*Thickness]);
    
    bdryfes = meshboundary(fes);
    icl = selectelem(fens, bdryfes, box=[1.0, 1.0, 0.0, pi/2, 0.0, Thickness], inflate=tolerance);
    for i=1:count(fens)
        t=fens.xyz[i,1]; a=fens.xyz[i,2]; z=fens.xyz[i,3]
        fens.xyz[i,:]=[(t*3.25+(1-t)*2)*cos(a), (t*2.75+(1-t)*1)*sin(a), z];
    end
    
    
    geom = NodalField(fens.xyz)
    u = NodalField(zeros(size(fens.xyz,1),3)) # displacement field
    
    l1 =selectnode(fens; box=[0.0, Inf, 0.0, 0.0, 0.0, Thickness], inflate = tolerance)
    setebc!(u,l1,true, 2, 0.0)
    l1 =selectnode(fens; box=[0.0, 0.0, 0.0, Inf, 0.0, Thickness], inflate = tolerance)
    setebc!(u,l1,true, 1, 0.0)
    l1 =selectnode(fens; box=[0.0, Inf, 0.0, Inf, 0.0, 0.0], inflate = tolerance)
    setebc!(u,l1,true, 3, 0.0)
    
    applyebc!(u)
    numberdofs!(u)
    
    
    el1femm =  FEMMBase(IntegData(subset(bdryfes,icl), GaussRule(2, 2)))
    function pfun(forceout::FVec{T}, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt) where {T}
        pt= [2.75/3.25*XYZ[1], 3.25/2.75*XYZ[2], 0.0]
        forceout .=    vec(p*pt/norm(pt));
        return forceout
    end
    fi = ForceIntensity(FFlt, 3, pfun);
    F2 = distribloads(el1femm, geom, u, fi, 2);
    
    # Note that the material object needs to be created with the proper
    # model-dimension reduction in mind.  In this case that is the fully three-dimensional solid.
    MR = DeforModelRed3D
    
    material = MatDeforElastIso(MR, E, nu)
    
    femm = FEMMDeforLinearMSH8(MR, IntegData(fes, GaussRule(3, 2)), material)
    
    # The geometry field now needs to be associated with the FEMM
    femm = associategeometry!(femm, geom)
    
    K = stiffness(femm, geom, u)
    K = cholesky(K)
    U = K\(F2)
    scattersysvec!(u, U[:])
    
    nl = selectnode(fens, box=[2.0, 2.0, 0.0, 0.0, 0.0, 0.0],inflate=tolerance);
    thecorneru = zeros(FFlt,1,3)
    gathervalues_asmat!(u, thecorneru, nl);
    thecorneru = thecorneru/phun("mm")
    println("$(time()-t0) [s];  displacement =$(thecorneru) [MM] as compared to reference [-0.10215,0] [MM]")
    
    
    fld= fieldfromintegpoints(femm, geom, u, :Cauchy, 2; nodevalmethod = :averaging, reportat = :extraptrend)
    println("Sigma_y =$(fld.values[nl,1][1]/phun("MPa")) as compared to reference sigma_yD = $(sigma_yD/phun("MPa")) [MPa]")
    
    println("$(n), $(fld.values[nl,1][1]/phun("MPa"))")
    
    File =  "a.vtk"
    vtkexportmesh(File, fes.conn, geom.values,
    FinEtools.MeshExportModule.H8; vectors=[("u", u.values)],
    scalars=[("sigmay", fld.values)])
    @async run(`"paraview.exe" $File`)
    true
    
end # LE1NAFEMS_MSH8_export


function LE1NAFEMS_MST10_convergence()
    # Example from the "Improved Stress Recovery for Mean-strain Finite Elements" paper by Sivapuram and Krysl, 2017
    println("LE1NAFEMS, 3D version. MST10"        )
    t0 = time()
    
    E = 210e3*phun("MEGA*PA");# 210e3 MPa
    nu = 0.3;
    p = 10*phun("MEGA*PA");# 10 MPA Outward pressure on the outside ellipse
    sigma_yD = 92.7*phun("MEGA*PA");# tensile stress at [2.0, 0.0] meters
    Thick0 = 0.1*phun("m")/2.0 # to account for the symmetry reduction
    
    sigyderrs = Dict{Symbol, FFltVec}()
    
    nnodes = []
    for extrapolation in [:extraptrend :extrapmean]
        sigyderrs[extrapolation] = FFltVec[]
        nnodes = []
        for ref in 0:1:4
            # Thickness = Thick0
            Thickness = Thick0/2^ref
            tolerance = Thickness/2^ref/1000.; # Geometrical tolerance
            
            fens,fes = T10block(1.0, pi/2, Thickness, 2^ref*5, 2^ref*6, 1; orientation = :b)
            
            bdryfes = meshboundary(fes);
            icl = selectelem(fens, bdryfes, box=[1.0, 1.0, 0.0, pi/2, 0.0, Thickness], inflate=tolerance);
            for i=1:count(fens)
                t=fens.xyz[i,1]; a=fens.xyz[i,2]; z=fens.xyz[i,3]
                fens.xyz[i,:]=[(t*3.25+(1-t)*2)*cos(a), (t*2.75+(1-t)*1)*sin(a), z];
            end
            
            
            geom = NodalField(fens.xyz)
            u = NodalField(zeros(size(fens.xyz,1),3)) # displacement field
            
            l1 =selectnode(fens; box=[0.0, Inf, 0.0, 0.0, 0.0, Thickness], inflate = tolerance)
            setebc!(u,l1,true, 2, 0.0)
            l1 =selectnode(fens; box=[0.0, 0.0, 0.0, Inf, 0.0, Thickness], inflate = tolerance)
            setebc!(u,l1,true, 1, 0.0)
            l1 =selectnode(fens; box=[0.0, Inf, 0.0, Inf, 0.0, 0.0], inflate = tolerance)
            setebc!(u,l1,true, 3, 0.0)
            
            applyebc!(u)
            numberdofs!(u)
            
            
            el1femm =  FEMMBase(IntegData(subset(bdryfes,icl), TriRule(3)))
            function pfun(forceout::FVec{T}, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt) where {T}
                pt= [2.75/3.25*XYZ[1], 3.25/2.75*XYZ[2], 0.0]
                forceout .=    vec(p*pt/norm(pt));
                return forceout
            end
            fi = ForceIntensity(FFlt, 3, pfun);
            F2 = distribloads(el1femm, geom, u, fi, 2);
            
            # Note that the material object needs to be created with the proper
            # model-dimension reduction in mind.  In this case that is the fully three-dimensional solid.
            MR = DeforModelRed3D
            
            material = MatDeforElastIso(MR, E, nu)
            
            femm = FEMMDeforLinearMST10(MR, IntegData(fes, TetRule(4)), material)
            
            # The geometry field now needs to be associated with the FEMM
            femm = associategeometry!(femm, geom)
            
            K = stiffness(femm, geom, u)
            K = cholesky(K)
            U = K\(F2)
            scattersysvec!(u, U[:])
            
            nl = selectnode(fens, box=[2.0, 2.0, 0.0, 0.0, 0.0, Thickness], inflate=tolerance);
            thecorneru = zeros(FFlt,length(nl),3)
            gathervalues_asmat!(u, thecorneru, nl);
            thecorneru = mean(thecorneru, 1)[1]/phun("mm")
            println("displacement =$(thecorneru) vs -0.10215 [MM]")
            
            fld = fieldfromintegpoints(femm, geom, u, :Cauchy, 2; nodevalmethod = :averaging, reportat = extrapolation)
            sigyd = mean(fld.values[nl,1], 1)[1]
            println("Sigma_y =$(sigyd/phun("MPa")) vs $(sigma_yD/phun("MPa")) [MPa]")
            
            println("$extrapolation, $(count(fes)), $(sigyd/phun("MPa"))")
            push!(nnodes, count(fens))
            push!(sigyderrs[extrapolation], abs(sigyd/sigma_yD - 1.0))
            # File =  "a.vtk"
            # vtkexportmesh(File, fes.conn, geom.values,
            # FinEtools.MeshExportModule.H8; vectors=[("u", u.values)],
            # scalars=[("sigmay", fld.values)])
            # @async run(`"paraview.exe" $File`)
        end
    end
    
    File = "LE1NAFEMS_MST10_convergence.CSV"
    savecsv(File, nnodes=vec(nnodes), sigyderrtrend=vec(sigyderrs[:extraptrend]), sigyderrdefault=vec(sigyderrs[:extrapmean]))
    
end # LE1NAFEMS_MST10_convergence


function LE1NAFEMS_MST10_one()
    println("LE1NAFEMS, 3D version."        )
    t0 = time()
    
    E = 210e3*phun("MEGA*PA");# 210e3 MPa
    nu = 0.3;
    p = 10*phun("MEGA*PA");# 10 MPA Outward pressure on the outside ellipse
    sigma_yD = 92.7*phun("MEGA*PA");# tensile stress at [2.0, 0.0] meters
    Thick0 = 0.1*phun("m")/2.0 # to account for the symmetry reduction
    
    sigyderrs = Dict{Symbol, FFltVec}()
    
    for extrapolation in [:extraptrend]
        sigyderrs[extrapolation] = FFltVec[]
        nnodes = []
        for ref in 1:1
            Thickness = Thick0
            # Thickness = Thick0/2^ref
            tolerance = Thickness/2^ref/1000.; # Geometrical tolerance
            
            fens,fes = T10block(1.0, pi/2, Thickness, 2^ref*5, 2^ref*6, 1; orientation = :b)
            
            bdryfes = meshboundary(fes);
            icl = selectelem(fens, bdryfes, box=[1.0, 1.0, 0.0, pi/2, 0.0, Thickness], inflate=tolerance);
            for i=1:count(fens)
                t=fens.xyz[i,1]; a=fens.xyz[i,2]; z=fens.xyz[i,3]
                fens.xyz[i,:]=[(t*3.25+(1-t)*2)*cos(a), (t*2.75+(1-t)*1)*sin(a), z];
            end
            
            
            geom = NodalField(fens.xyz)
            u = NodalField(zeros(size(fens.xyz,1),3)) # displacement field
            
            l1 =selectnode(fens; box=[0.0, Inf, 0.0, 0.0, 0.0, Thickness], inflate = tolerance)
            setebc!(u,l1,true, 2, 0.0)
            l1 =selectnode(fens; box=[0.0, 0.0, 0.0, Inf, 0.0, Thickness], inflate = tolerance)
            setebc!(u,l1,true, 1, 0.0)
            l1 =selectnode(fens; box=[0.0, Inf, 0.0, Inf, 0.0, 0.0], inflate = tolerance)
            setebc!(u,l1,true, 3, 0.0)
            
            applyebc!(u)
            numberdofs!(u)
            
            
            el1femm =  FEMMBase(IntegData(subset(bdryfes,icl), TriRule(3)))
            function pfun(forceout::FVec{T}, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt) where {T}
                pt= [2.75/3.25*XYZ[1], 3.25/2.75*XYZ[2], 0.0]
                forceout .=    vec(p*pt/norm(pt));
                return forceout
            end
            fi = ForceIntensity(FFlt, 3, pfun);
            F2 = distribloads(el1femm, geom, u, fi, 2);
            
            # Note that the material object needs to be created with the proper
            # model-dimension reduction in mind.  In this case that is the fully three-dimensional solid.
            MR = DeforModelRed3D
            
            material = MatDeforElastIso(MR, E, nu)
            
            femm = FEMMDeforLinearMST10(MR, IntegData(fes, TetRule(4)), material)
            
            # The geometry field now needs to be associated with the FEMM
            femm = associategeometry!(femm, geom)
            
            K = stiffness(femm, geom, u)
            K = cholesky(K)
            U = K\(F2)
            scattersysvec!(u, U[:])
            
            nl = selectnode(fens, box=[2.0, 2.0, 0.0, 0.0, 0.0, Thickness], inflate=tolerance);
            thecorneru = zeros(FFlt,length(nl),3)
            gathervalues_asmat!(u, thecorneru, nl);
            thecorneru = mean(thecorneru, 1)[1]/phun("mm")
            println("displacement =$(thecorneru) vs -0.10215 [MM]")
            
            fld = fieldfromintegpoints(femm, geom, u, :Cauchy, 2; nodevalmethod = :averaging, reportat = extrapolation)
            sigyd = mean(fld.values[nl,1], 1)[1]
            println("Sigma_y =$(sigyd/phun("MPa")) vs $(sigma_yD/phun("MPa")) [MPa]")
            
            println("$extrapolation, $(count(fes)), $(sigyd/phun("MPa"))")
            push!(nnodes, count(fens))
            push!(sigyderrs[extrapolation], abs(sigyd/sigma_yD - 1.0))
            File =  "LE1NAFEMS_MST10_a.vtk"
            vtkexportmesh(File, fes.conn, geom.values,
            FinEtools.MeshExportModule.T10; vectors=[("u", u.values)],
            scalars=[("sigmay", fld.values)])
            @async run(`"paraview.exe" $File`)
        end
    end
    
    #
    #
    # using DataFrames
    # using CSV
    #
    # df = DataFrame(nnodes=vec(nnodes),
    #     sigyderrtrend=vec(sigyderrs[:extraptrend]),
    #     sigyderrdefault=vec(sigyderrs[:extrapmean]))
    # File = "LE1NAFEMS_MST10_convergence.CSV"
    # CSV.write(File, df)
    # @async run(`"paraview.exe" $File`)
    
end # LE1NAFEMS_MST10_one

# using SymPy
# function LE1NAFEMS_MST10_single()
    
    
#     @vars ua ub uc ud
#     @vars va vb vc vd
#     @vars wa wb wc wd
#     @vars x y z
#     u(x, y, z) = ua*x +  ub*y +  uc*z +  ud
#     v(x, y, z) = va*x +  vb*y +  vc*z +  vd
#     w(x, y, z) = wa*x +  wb*y +  wc*z +  wd
    
#     ex = diff(u(x, y, z), x)
#     ey = diff(v(x, y, z), y)
#     ez = diff(w(x, y, z), z)
#     gxy = diff(v(x, y, z), x) + diff(u(x, y, z), y)
#     gxz = diff(w(x, y, z), x) + diff(u(x, y, z), z)
#     gyz = diff(w(x, y, z), y) + diff(v(x, y, z), z)
    
#     display(ex)
#     display(ey)
#     display(ez)
#     display(gxy)
#     display(gxz)
#     display(gyz)
    
# end # LE1NAFEMS_MST10_single


function LE1NAFEMS_MST10_stresses_nodal()
    elementtag = "MST10"
    println("LE1NAFEMS, 3D version. Element: $(elementtag)")
    
    E = 210e3*phun("MEGA*PA");# 210e3 MPa
    nu = 0.3;
    p = 10*phun("MEGA*PA");# 10 MPA Outward pressure on the outside ellipse
    sigma_yD = 92.7*phun("MEGA*PA");# tensile stress at [2.0, 0.0] meters
    Thick0 = 0.1*phun("m")/2.0 # to account for the symmetry reduction
    
    for extrapolation in [:extraptrend :extrapmean]
        convergencestudy = FDataDict[]
        for ref in  [0, 1, 2, 3, 4]
            Thickness = Thick0
            tolerance = Thickness/2^ref/1000.; # Geometrical tolerance
            
            fens,fes = T10block(1.0, pi/2, Thickness, 2^ref*2, 2^ref*3, 1; orientation = :b)
            
            bdryfes = meshboundary(fes);
            icl = selectelem(fens, bdryfes, box=[1.0, 1.0, 0.0, pi/2, 0.0, Thickness], inflate=tolerance);
            for i=1:count(fens)
                t=fens.xyz[i,1]; a=fens.xyz[i,2]; z=fens.xyz[i,3]
                fens.xyz[i,:]=[(t*3.25+(1-t)*2)*cos(a), (t*2.75+(1-t)*1)*sin(a), z];
            end
            
            geom = NodalField(fens.xyz)
            u = NodalField(zeros(size(fens.xyz,1),3)) # displacement field
            
            l1 =selectnode(fens; box=[0.0, Inf, 0.0, 0.0, 0.0, Thickness], inflate = tolerance)
            setebc!(u,l1,true, 2, 0.0)
            l1 =selectnode(fens; box=[0.0, 0.0, 0.0, Inf, 0.0, Thickness], inflate = tolerance)
            setebc!(u,l1,true, 1, 0.0)
            l1 =selectnode(fens; box=[0.0, Inf, 0.0, Inf, 0.0, 0.0], inflate = tolerance)
            setebc!(u,l1,true, 3, 0.0)
            
            applyebc!(u)
            numberdofs!(u)
            
            el1femm =  FEMMBase(IntegData(subset(bdryfes,icl), TriRule(3)))
            function pfun(forceout::FVec{T}, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt) where {T}
                pt= [2.75/3.25*XYZ[1], 3.25/2.75*XYZ[2], 0.0]
                forceout .=    vec(p*pt/norm(pt));
                return forceout
            end
            fi = ForceIntensity(FFlt, 3, pfun);
            F2 = distribloads(el1femm, geom, u, fi, 2);
            
            # Note that the material object needs to be created with the proper
            # model-dimension reduction in mind.  In this case that is the fully three-dimensional solid.
            MR = DeforModelRed3D
            
            material = MatDeforElastIso(MR, E, nu)
            
            femm = FEMMDeforLinearMST10(MR, IntegData(fes, TetRule(4)), material)
            
            # The geometry field now needs to be associated with the FEMM
            femm = associategeometry!(femm, geom)
            
            K = stiffness(femm, geom, u)
            K = cholesky(K)
            U = K\(F2)
            scattersysvec!(u, U[:])
            
            nl = selectnode(fens, box=[2.0, 2.0, 0.0, 0.0, 0.0, Thickness], inflate=tolerance);
            thecorneru = zeros(FFlt,length(nl),3)
            gathervalues_asmat!(u, thecorneru, nl);
            thecorneru = mean(thecorneru, 1)[1]/phun("mm")
            println("displacement =$(thecorneru) vs -0.10215 [MM]")
            
            fld = fieldfromintegpoints(femm, geom, u, :Cauchy, 2;
            nodevalmethod = :averaging, reportat = extrapolation)
            sigyd = mean(fld.values[nl,1], 1)[1]
            println("Sigma_y =$(sigyd/phun("MPa")) vs $(sigma_yD/phun("MPa")) [MPa]")
            
            stressfield = fieldfromintegpoints(femm, geom, u, :Cauchy, collect(1:6);
            nodevalmethod = :averaging, reportat = extrapolation)
            
            # File =  "LE1NAFEMS_MST10_sigma.vtk"
            # vtkexportmesh(File, fes.conn, geom.values,
            #     FinEtools.MeshExportModule.T10; vectors=[("u", u.values)],
            #     scalars=[("sig", stressfield.values)])
            # @async run(`"paraview.exe" $File`)
            
            push!(convergencestudy, FDataDict(
            "elementsize"=> 1.0 / (2^ref),
            "fens"=>fens,
            "fes"=>fes,
            "geom"=>geom,
            "u"=>u,
            "femm"=>femm,
            "integrationrule"=>femm.integdata.integration_rule,
            "stressfields"=>[stressfield],
            "tolerance"=>tolerance)
            )
        end # for ref
        
        # File = "LE1NAFEMS_MST10_stresses_nodal_convergence_$(elementtag)_$(extrapolation)"
        # open(File * ".jls", "w") do file
        #     serialize(file, convergencestudy)
        # end
        
    end # for extrapolation
    
end # LE1NAFEMS_MST10_stresses_nodal


function LE1NAFEMS_MST10_S_convergence()
    println("LE1NAFEMS, 3D version."        )
    t0 = time()
    
    E = 210e3*phun("MEGA*PA");# 210e3 MPa
    nu = 0.3;
    p = 10*phun("MEGA*PA");# 10 MPA Outward pressure on the outside ellipse
    sigma_yD = 92.7*phun("MEGA*PA");# tensile stress at [2.0, 0.0] meters
    Thick0 = 0.1*phun("m")/2.0 # to account for the symmetry reduction
    
    sigyderrs = Dict{Symbol, FFltVec}()
    
    nnodes = []
    for extrapolation in [:extraptrend :extrapmean]
        sigyderrs[extrapolation] = FFltVec[]
        nnodes = []
        for ref in 0:1:4
            Thickness = Thick0
            # Thickness = Thick0/2^ref
            tolerance = Thickness/2^ref/1000.; # Geometrical tolerance
            
            fens,fes = T4block(1.0, pi/2, Thickness, 2^ref*5, 2^ref*6, 1)
            
            bdryfes = meshboundary(fes);
            icl = selectelem(fens, bdryfes, box=[1.0, 1.0, 0.0, pi/2, 0.0, Thickness], inflate=tolerance);
            for i=1:count(fens)
                t=fens.xyz[i,1]; a=fens.xyz[i,2]; z=fens.xyz[i,3]
                fens.xyz[i,:]=[(t*3.25+(1-t)*2)*cos(a), (t*2.75+(1-t)*1)*sin(a), z];
            end
            fens,fes = T4toT10(fens,fes)
            
            geom = NodalField(fens.xyz)
            u = NodalField(zeros(size(fens.xyz,1),3)) # displacement field
            
            l1 =selectnode(fens; box=[0.0, Inf, 0.0, 0.0, 0.0, Thickness], inflate = tolerance)
            setebc!(u,l1,true, 2, 0.0)
            l1 =selectnode(fens; box=[0.0, 0.0, 0.0, Inf, 0.0, Thickness], inflate = tolerance)
            setebc!(u,l1,true, 1, 0.0)
            l1 =selectnode(fens; box=[0.0, Inf, 0.0, Inf, 0.0, 0.0], inflate = tolerance)
            setebc!(u,l1,true, 3, 0.0)
            
            applyebc!(u)
            numberdofs!(u)
            
            
            el1femm =  FEMMBase(IntegData(subset(bdryfes,icl), TriRule(3)))
            function pfun(forceout::FVec{T}, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt) where {T}
                pt= [2.75/3.25*XYZ[1], 3.25/2.75*XYZ[2], 0.0]
                forceout .=    vec(p*pt/norm(pt));
                return forceout
            end
            fi = ForceIntensity(FFlt, 3, pfun);
            F2 = distribloads(el1femm, geom, u, fi, 2);
            
            # Note that the material object needs to be created with the proper
            # model-dimension reduction in mind.  In this case that is the fully three-dimensional solid.
            MR = DeforModelRed3D
            
            material = MatDeforElastIso(MR, E, nu)
            
            femm = FEMMDeforLinearMST10(MR, IntegData(fes, TetRule(4)), material)
            
            # The geometry field now needs to be associated with the FEMM
            femm = associategeometry!(femm, geom)
            
            K = stiffness(femm, geom, u)
            K = cholesky(K)
            U = K\(F2)
            scattersysvec!(u, U[:])
            
            nl = selectnode(fens, box=[2.0, 2.0, 0.0, 0.0, 0.0, Thickness], inflate=tolerance);
            thecorneru = zeros(FFlt,length(nl),3)
            gathervalues_asmat!(u, thecorneru, nl);
            thecorneru = mean(thecorneru, 1)[1]/phun("mm")
            println("displacement =$(thecorneru) vs -0.10215 [MM]")
            
            fld = fieldfromintegpoints(femm, geom, u, :Cauchy, 2; nodevalmethod = :averaging, reportat = extrapolation)
            sigyd = mean(fld.values[nl,1], 1)[1]
            println("Sigma_y =$(sigyd/phun("MPa")) vs $(sigma_yD/phun("MPa")) [MPa]")
            
            println("$extrapolation, $(count(fes)), $(sigyd/phun("MPa"))")
            push!(nnodes, count(fens))
            push!(sigyderrs[extrapolation], abs(sigyd/sigma_yD - 1.0))
            # File =  "a.vtk"
            # vtkexportmesh(File, fes.conn, geom.values,
            # FinEtools.MeshExportModule.H8; vectors=[("u", u.values)],
            # scalars=[("sigmay", fld.values)])
            # @async run(`"paraview.exe" $File`)
        end
    end
    
    File = "LE1NAFEMS_MST10_S_convergence.CSV"
    savecsv(File, nnodes=vec(nnodes), sigyderrtrend=vec(sigyderrs[:extraptrend]), sigyderrdefault=vec(sigyderrs[:extrapmean]))
    
end # LE1NAFEMS_MST10_S_convergence


function LE1NAFEMS_T10_stresses_nodal()
    elementtag = "T10"
    println("LE1NAFEMS, 3D version. Element: $(elementtag)")
    
    E = 210e3*phun("MEGA*PA");# 210e3 MPa
    nu = 0.3;
    p = 10*phun("MEGA*PA");# 10 MPA Outward pressure on the outside ellipse
    sigma_yD = 92.7*phun("MEGA*PA");# tensile stress at [2.0, 0.0] meters
    Thick0 = 0.1*phun("m")/2.0 # to account for the symmetry reduction
    
    convergencestudy = FDataDict[]
    for ref in  [0, 1, 2, 3, 4]
        Thickness = Thick0
        tolerance = Thickness/2^ref/1000.; # Geometrical tolerance
        
        fens,fes = T10block(1.0, pi/2, Thickness, 2^ref*2, 2^ref*3, 1; orientation = :b)
        
        bdryfes = meshboundary(fes);
        icl = selectelem(fens, bdryfes, box=[1.0, 1.0, 0.0, pi/2, 0.0, Thickness], inflate=tolerance);
        for i=1:count(fens)
            t=fens.xyz[i,1]; a=fens.xyz[i,2]; z=fens.xyz[i,3]
            fens.xyz[i,:]=[(t*3.25+(1-t)*2)*cos(a), (t*2.75+(1-t)*1)*sin(a), z];
        end
        
        
        geom = NodalField(fens.xyz)
        u = NodalField(zeros(size(fens.xyz,1),3)) # displacement field
        
        l1 =selectnode(fens; box=[0.0, Inf, 0.0, 0.0, 0.0, Thickness], inflate = tolerance)
        setebc!(u,l1,true, 2, 0.0)
        l1 =selectnode(fens; box=[0.0, 0.0, 0.0, Inf, 0.0, Thickness], inflate = tolerance)
        setebc!(u,l1,true, 1, 0.0)
        l1 =selectnode(fens; box=[0.0, Inf, 0.0, Inf, 0.0, 0.0], inflate = tolerance)
        setebc!(u,l1,true, 3, 0.0)
        
        applyebc!(u)
        numberdofs!(u)
        
        
        el1femm =  FEMMBase(IntegData(subset(bdryfes,icl), TriRule(3)))
        function pfun(forceout::FVec{T}, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt) where {T}
            pt= [2.75/3.25*XYZ[1], 3.25/2.75*XYZ[2], 0.0]
            forceout .=    vec(p*pt/norm(pt));
            return forceout
        end
        fi = ForceIntensity(FFlt, 3, pfun);
        F2 = distribloads(el1femm, geom, u, fi, 2);
        
        # Note that the material object needs to be created with the proper
        # model-dimension reduction in mind.  In this case that is the fully three-dimensional solid.
        MR = DeforModelRed3D
        
        material = MatDeforElastIso(MR, E, nu)
        
        femm = FEMMDeforLinear(MR, IntegData(fes, TetRule(4)), material)
        
        # The geometry field now needs to be associated with the FEMM
        femm = associategeometry!(femm, geom)
        
        K = stiffness(femm, geom, u)
        K = cholesky(K)
        U = K\(F2)
        scattersysvec!(u, U[:])
        
        nl = selectnode(fens, box=[2.0, 2.0, 0.0, 0.0, 0.0, Thickness], inflate=tolerance);
        thecorneru = zeros(FFlt,length(nl),3)
        gathervalues_asmat!(u, thecorneru, nl);
        thecorneru = mean(thecorneru, 1)[1]/phun("mm")
        println("displacement =$(thecorneru) vs -0.10215 [MM]")
        
        fld = fieldfromintegpoints(femm, geom, u, :Cauchy, 2)
        sigyd = mean(fld.values[nl,1], 1)[1]
        println("Sigma_y =$(sigyd/phun("MPa")) vs $(sigma_yD/phun("MPa")) [MPa]")
        
        stressfield = fieldfromintegpoints(femm, geom, u, :Cauchy, collect(1:6))
        
        # File =  "LE1NAFEMS_T10_sigma.vtk"
        # vtkexportmesh(File, fes.conn, geom.values,
        #     FinEtools.MeshExportModule.T10; vectors=[("u", u.values)],
        #     scalars=[("sig", stressfield.values)])
        # @async run(`"paraview.exe" $File`)
        
        push!(convergencestudy, FDataDict(
        "elementsize"=> 1.0 / (2^ref),
        "fens"=>fens,
        "fes"=>fes,
        "geom"=>geom,
        "u"=>u,
        "femm"=>femm,
        "integrationrule"=>femm.integdata.integration_rule,
        "stressfields"=>[stressfield],
        "tolerance"=>tolerance)
        )
    end # for ref 
    
    # File = "LE1NAFEMS_T10_stresses_nodal_convergence_$(elementtag)"
    # open(File * ".jls", "w") do file
    #     serialize(file, convergencestudy)
    # end
    
end # LE1NAFEMS_T10_stresses_nodal

function allrun()
    println("#####################################################") 
    println("# LE1NAFEMS_MSH8 ")
    LE1NAFEMS_MSH8()
    println("#####################################################") 
    println("# LE1NAFEMS_MSH8_convergence ")
    LE1NAFEMS_MSH8_convergence()
    println("#####################################################") 
    println("# LE1NAFEMS_MSH8_export ")
    LE1NAFEMS_MSH8_export()
    println("#####################################################") 
    println("# LE1NAFEMS_MST10_convergence ")
    LE1NAFEMS_MST10_convergence()
    println("#####################################################") 
    println("# LE1NAFEMS_MST10_one ")
    LE1NAFEMS_MST10_one()
    # println("#####################################################") 
    # println("# LE1NAFEMS_MST10_single ")
    # LE1NAFEMS_MST10_single()
    println("#####################################################") 
    println("# LE1NAFEMS_MST10_stresses_nodal ")
    LE1NAFEMS_MST10_stresses_nodal()
    println("#####################################################") 
    println("# LE1NAFEMS_MST10_S_convergence ")
    LE1NAFEMS_MST10_S_convergence()
    println("#####################################################") 
    println("# LE1NAFEMS_T10_stresses_nodal ")
    LE1NAFEMS_T10_stresses_nodal()
    return true
end # function allrun

end # module LE1NAFEMS_examples
