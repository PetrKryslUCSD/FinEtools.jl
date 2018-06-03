module plate_w_hole_examples
using FinEtools
using FinEtools.MeshExportModule

function plate_w_hole_H20_stress()
    E = 2.4*phun("MEGA*PA");# 210e3 MPa
    nu = 0.3;
    Re = 0.3*phun("M"); # outer radius
    Ri= 0.1*phun("M"); # hole radius
    H = 0.1*phun("M") # thickness of the plate
    nRadial, nCircumferential, nThickness = 6, 8, 1;
    sigma0=1*phun("MEGA*PA");
    
    function sigmaxx(x)
        local r = norm(x[1:2]);
        local th = atan2(x[2],x[1]);
        return sigma0*(1-Ri^2/r^2*(3/2*cos(2*th)+cos(4*th))+3/2*Ri^4/r^4*cos(4*th));
    end
    function sigmayy(x)
        local r=norm(x[1:2]);
        local th = atan2(x[2],x[1]);
        return -sigma0*(Ri^2/r^2*(1/2*cos(2*th)-cos(4*th))+3/2*Ri^4/r^4*cos(4*th));
    end
    function sigmaxy(x)
        local r=norm(x[1:2]);
        local th = atan2(x[2],x[1]);
        return -sigma0*(Ri^2/r^2*(1/2*sin(2*th)+sin(4*th))-3/2*Ri^4/r^4*sin(4*th));
    end
    function sigmarr(x)
        local r = norm(x[1:2]);
        local th = atan2(x[2],x[1]);
        return sigma0/2*(1-Ri^2/r^2) + sigma0/2*(1-4*Ri^2/r^2+3*Ri^4/r^4)*cos(2*th)
    end
    function sigmatt(x)
        local r = norm(x[1:2]);
        local th = atan2(x[2],x[1]);
        return sigma0/2*(1+Ri^2/r^2) - sigma0/2*(1+3*Ri^4/r^4)*cos(2*th)
    end
    function sigmart(x)
        local r = norm(x[1:2]);
        local th = atan2(x[2],x[1]);
        return -sigma0/2*(1+2*Ri^2/r^2-3*Ri^4/r^4)*sin(2*th)
    end
    
    convergencestudy = FDataDict[]
    for ref in 0:1:2
        println("ref = $(ref)")
        # Thickness = H
        Thickness = H/2^ref
        tolerance = Thickness/2^ref/1000.; # Geometrical tolerance
        
        fens,fes = H20block(1.0, pi/2, Thickness, 2^ref*nRadial, 2^ref*nCircumferential, 2^ref*nThickness)
        
        bdryfes = meshboundary(fes);
        icl = selectelem(fens, bdryfes, box=[1.0, 1.0, 0.0, pi/2, 0.0, Thickness], inflate=tolerance);
        
        for i=1:count(fens)
            t=fens.xyz[i,1]; a=fens.xyz[i,2]; z=fens.xyz[i,3]
            fens.xyz[i,:] = [(t*Re+(1-t)*Ri)*cos(a), (t*Re+(1-t)*Ri)*sin(a), z];
        end
        
        geom = NodalField(fens.xyz)
        u = NodalField(zeros(size(fens.xyz,1),3)) # displacement field
        
        l1 =selectnode(fens; box=[0.0, Inf, 0.0, 0.0, 0.0, Thickness], inflate = tolerance)
        setebc!(u,l1,true, 2, 0.0)
        l1 =selectnode(fens; box=[0.0, 0.0, 0.0, Inf, 0.0, Thickness], inflate = tolerance)
        setebc!(u,l1,true, 1, 0.0)
        # Plane-stress constraint: assume the plane z=0 is the plane of symmetry of the plate
        l1 =selectnode(fens; box=[0.0, Inf, 0.0, Inf, 0.0, 0.0], inflate = tolerance)
        setebc!(u,l1,true, 3, 0.0)
        # If this was enabled, the plane-strain  constraint would be enforced.
        # l1 =selectnode(fens; box=[0.0, Inf, 0.0, Inf, Thickness, Thickness], inflate = tolerance)
        # setebc!(u,l1,true, 3, 0.0)
        
        applyebc!(u)
        numberdofs!(u)
        el1femm =  FEMMBase(IntegData(subset(bdryfes,icl), GaussRule(2, 3)))
        function pfun(forceout::FVec{T}, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt) where {T}
            local r = sqrt(XYZ[1]^2 + XYZ[2]^2)
            nx = XYZ[1]/r; ny = XYZ[2]/r
            forceout[1] = sigmarr(XYZ) * nx - sigmart(XYZ) * ny
            forceout[2] = sigmarr(XYZ) * ny + sigmart(XYZ) * nx
            forceout[3] = 0.0
            return forceout
        end
        fi = ForceIntensity(FFlt, 3, pfun);
        F2 = distribloads(el1femm, geom, u, fi, 2);
        
        MR = DeforModelRed3D
        
        material = MatDeforElastIso(MR, E, nu)
        
        femm = FEMMDeforLinear(MR, IntegData(fes, GaussRule(3, 2)), material)
        
        # The geometry field now needs to be associated with the FEMM
        femm = associategeometry!(femm, geom)
        
        K = stiffness(femm, geom, u)
        K = cholesky(K)
        U = K\(F2)
        scattersysvec!(u, U[:])
        
        stressfields = elemfieldfromintegpoints(femm, geom, u, :Cauchy, collect(1:6))
        
        push!(convergencestudy, FDataDict(
        "elementsize"=> 1.0 / 2^ref,
        "fens"=>fens,
        "fes"=>fes,
        "geom"=>geom,
        "u"=>u,
        "femm"=>femm,
        "stressfields"=>[stressfields],
        "tolerance"=>tolerance)
        )
    end # for ref in 0:1:5
    
    # File = "mplate_w_hole_H20m_stress"
    # open(File * ".jls", "w") do file
    #     serialize(file, convergencestudy)
    # end
    
end # plate_w_hole_H20_stress


function plate_w_hole_MSH8_convergence()
    E = 2.4*phun("MEGA*PA");# 210e3 MPa
    nu = 0.49995;
    Re = 0.3*phun("M"); # outer radius
    Ri= 0.1*phun("M"); # hole radius
    H = 0.1*phun("M") # thickness of the plate
    nRadial, nCircumferential=3, 5;
    sigma0=1*phun("MEGA*PA");
    
    function sigmaxx(x)
        local r = norm(x[1:2]);
        local th = atan2(x[2],x[1]);
        return sigma0*(1-Ri^2/r^2*(3/2*cos(2*th)+cos(4*th))+3/2*Ri^4/r^4*cos(4*th));
    end
    function sigmayy(x)
        local r=norm(x[1:2]);
        local th = atan2(x[2],x[1]);
        return -sigma0*(Ri^2/r^2*(1/2*cos(2*th)-cos(4*th))+3/2*Ri^4/r^4*cos(4*th));
    end
    function sigmaxy(x)
        local r=norm(x[1:2]);
        local th = atan2(x[2],x[1]);
        return -sigma0*(Ri^2/r^2*(1/2*sin(2*th)+sin(4*th))-3/2*Ri^4/r^4*sin(4*th));
    end
    function sigmarr(x)
        local r = norm(x[1:2]);
        local th = atan2(x[2],x[1]);
        return sigma0/2*(1-Ri^2/r^2) + sigma0/2*(1-4*Ri^2/r^2+3*Ri^4/r^4)*cos(2*th)
    end
    function sigmatt(x)
        local r = norm(x[1:2]);
        local th = atan2(x[2],x[1]);
        return sigma0/2*(1+Ri^2/r^2) - sigma0/2*(1+3*Ri^4/r^4)*cos(2*th)
    end
    function sigmart(x)
        local r = norm(x[1:2]);
        local th = atan2(x[2],x[1]);
        return -sigma0/2*(1+2*Ri^2/r^2-3*Ri^4/r^4)*sin(2*th)
    end
    
    sigxderrs = Dict{Symbol, FFltVec}()
    sigyderrs = Dict{Symbol, FFltVec}()
    numelements = []
    numnodes = []
    for extrap in [:extraptrend :extrapmean]
        sigxderrs[extrap] = FFltVec[]
        sigyderrs[extrap] = FFltVec[]
        numelements = []
        numnodes = []
        for ref in 0:1:5
            # Thickness = H
            Thickness = H/2^ref
            tolerance = Thickness/2^ref/1000.; # Geometrical tolerance
            
            fens,fes = H8block(1.0, pi/2, Thickness, 2^ref*nRadial, 2^ref*nCircumferential, 1)
            
            bdryfes = meshboundary(fes);
            icl = selectelem(fens, bdryfes, box=[1.0, 1.0, 0.0, pi/2, 0.0, Thickness], inflate=tolerance);
            
            for i=1:count(fens)
                t=fens.xyz[i,1]; a=fens.xyz[i,2]; z=fens.xyz[i,3]
                fens.xyz[i,:] = [(t*Re+(1-t)*Ri)*cos(a), (t*Re+(1-t)*Ri)*sin(a), z];
            end
            
            # File =  "a.vtk"
            # vtkexportmesh(File, fes.conn, fens.xyz,
            #     FinEtools.MeshExportModule.H8)
            # @async run(`"paraview.exe" $File`)
            
            geom = NodalField(fens.xyz)
            u = NodalField(zeros(size(fens.xyz,1),3)) # displacement field
            
            l1 =selectnode(fens; box=[0.0, Inf, 0.0, 0.0, 0.0, Thickness], inflate = tolerance)
            setebc!(u,l1,true, 2, 0.0)
            l1 =selectnode(fens; box=[0.0, 0.0, 0.0, Inf, 0.0, Thickness], inflate = tolerance)
            setebc!(u,l1,true, 1, 0.0)
            # Plane-stress constraint: assume the plane z=0 is the plane of symmetry of the plate
            l1 =selectnode(fens; box=[0.0, Inf, 0.0, Inf, 0.0, 0.0], inflate = tolerance)
            setebc!(u,l1,true, 3, 0.0)
            # If this was enabled, the plane-strain  constraint would be enforced.
            # l1 =selectnode(fens; box=[0.0, Inf, 0.0, Inf, Thickness, Thickness], inflate = tolerance)
            # setebc!(u,l1,true, 3, 0.0)
            
            applyebc!(u)
            numberdofs!(u)
            el1femm =  FEMMBase(IntegData(subset(bdryfes,icl), GaussRule(2, 2)))
            function pfun(forceout::FVec{T}, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt) where {T}
                local r = sqrt(XYZ[1]^2 + XYZ[2]^2)
                nx = XYZ[1]/r; ny = XYZ[2]/r
                # local sx, sy, txy
                # sx, sy, txy = sigmaxx(XYZ), sigmayy(XYZ), sigmaxy(XYZ)
                # sn = sx * nx^2 + sy * ny^2 + 2 * nx * ny * txy
                # tn = -(sx - sy) * nx * ny + (nx^2 - ny^2) * txy
                # forceout[1] = sn * nx - tn * ny
                # forceout[2] = sn * ny + tn * nx
                # forceout[3] = 0.0
                forceout[1] = sigmarr(XYZ) * nx - sigmart(XYZ) * ny
                forceout[2] = sigmarr(XYZ) * ny + sigmart(XYZ) * nx
                forceout[3] = 0.0
                return forceout
            end
            fi = ForceIntensity(FFlt, 3, pfun);
            F2 = distribloads(el1femm, geom, u, fi, 2);
            
            MR = DeforModelRed3D
            
            material = MatDeforElastIso(MR, E, nu)
            
            femm = FEMMDeforLinearMSH8(MR, IntegData(fes, GaussRule(3, 2)), material)
            
            # The geometry field now needs to be associated with the FEMM
            femm = associategeometry!(femm, geom)
            
            K = stiffness(femm, geom, u)
            K = cholesky(K)
            U = K\(F2)
            scattersysvec!(u, U[:])
            
            nlA = selectnode(fens, box=[Ri, Ri, 0.0, 0.0, 0.0, Thickness], inflate=tolerance);
            nlB = selectnode(fens, box=[0.0, 0.0, Ri, Ri, 0.0, Thickness], inflate=tolerance);
            # thecorneru = zeros(FFlt,length(nlA),3)
            # gathervalues_asmat!(u, thecorneru, nl);
            # thecorneru = mean(thecorneru, 1)[1]/phun("mm")
            # println("displacement = $(thecorneru) vs -0.10215 [MM]")
            
            println("Extrapolation: $( extrap )--------------- ")
            sigx = fieldfromintegpoints(femm, geom, u, :Cauchy, 1;
            nodevalmethod = :averaging, reportat = extrap)
            sigy = fieldfromintegpoints(femm, geom, u, :Cauchy, 2;
            nodevalmethod = :averaging, reportat = extrap)
            sigyA = mean(sigy.values[nlA,1], 1)[1]
            sigyAtrue = sigmatt([Ri, 0.0, 0.0])
            println("sig_y@A =$(sigyA/phun("MPa")) vs $(sigyAtrue/phun("MPa")) [MPa]")
            sigxB = mean(sigx.values[nlB,1], 1)[1]
            sigxBtrue = sigmatt([0.0, Ri, 0.0])
            println("sig_x@B =$(sigxB/phun("MPa")) vs $(sigxBtrue/phun("MPa")) [MPa]")
            push!(numnodes, count(fens))
            push!(numelements, count(fes))
            push!(sigxderrs[extrap], abs(sigxB/sigxBtrue - 1.0))
            push!(sigyderrs[extrap], abs(sigyA/sigyAtrue - 1.0))
            # File =  "a.vtk"
            # vtkexportmesh(File, fes.conn, geom.values,
            # FinEtools.MeshExportModule.H8; vectors=[("u", u.values)],
            # scalars=[("sigmax", sigx.values/phun("MEGA*PA"))])
            # @async run(`"paraview.exe" $File`)
        end
    end
    
    File = "plate_w_hole_MSH8_convergence.CSV"
    savecsv(File, numelements=vec(numelements), numnodes=vec(numnodes), sigxderrtrend=vec(sigxderrs[:extraptrend]), sigxderrdefault=vec(sigxderrs[:extrapmean]), sigyderrtrend=vec(sigyderrs[:extraptrend]), sigyderrdefault=vec(sigyderrs[:extrapmean]))
    
end # plate_w_hole_MSH8_convergence


function plate_w_hole_MSH8_PE_convergence()
    E = 2.4*phun("MEGA*PA");# 210e3 MPa
    nu = 0.49995;
    Re = 0.3*phun("M"); # outer radius
    Ri= 0.1*phun("M"); # hole radius
    H = 0.1*phun("M") # thickness of the plate
    nRadial, nCircumferential=3, 5;
    sigma0=1*phun("MEGA*PA");
    
    function sigmaxx(x)
        local r = norm(x[1:2]);
        local th = atan2(x[2],x[1]);
        return sigma0*(1-Ri^2/r^2*(3/2*cos(2*th)+cos(4*th))+3/2*Ri^4/r^4*cos(4*th));
    end
    function sigmayy(x)
        local r=norm(x[1:2]);
        local th = atan2(x[2],x[1]);
        return -sigma0*(Ri^2/r^2*(1/2*cos(2*th)-cos(4*th))+3/2*Ri^4/r^4*cos(4*th));
    end
    function sigmaxy(x)
        local r=norm(x[1:2]);
        local th = atan2(x[2],x[1]);
        return -sigma0*(Ri^2/r^2*(1/2*sin(2*th)+sin(4*th))-3/2*Ri^4/r^4*sin(4*th));
    end
    function sigmarr(x)
        local r = norm(x[1:2]);
        local th = atan2(x[2],x[1]);
        return sigma0/2*(1-Ri^2/r^2) + sigma0/2*(1-4*Ri^2/r^2+3*Ri^4/r^4)*cos(2*th)
    end
    function sigmatt(x)
        local r = norm(x[1:2]);
        local th = atan2(x[2],x[1]);
        return sigma0/2*(1+Ri^2/r^2) - sigma0/2*(1+3*Ri^4/r^4)*cos(2*th)
    end
    function sigmart(x)
        local r = norm(x[1:2]);
        local th = atan2(x[2],x[1]);
        return -sigma0/2*(1+2*Ri^2/r^2-3*Ri^4/r^4)*sin(2*th)
    end
    
    sigxderrs = Dict{Symbol, FFltVec}()
    sigyderrs = Dict{Symbol, FFltVec}()
    numelements = []
    numnodes = []
    for extrapolation in [:extraptrend :extrapmean]
        sigxderrs[extrapolation] = FFltVec[]
        sigyderrs[extrapolation] = FFltVec[]
        numelements = []
        numnodes = []
        for ref in 0:1:5
            # Thickness = H
            Thickness = H/2^ref
            tolerance = Thickness/2^ref/1000.; # Geometrical tolerance
            
            fens,fes = H8block(1.0, pi/2, Thickness, 2^ref*nRadial, 2^ref*nCircumferential, 1)
            
            bdryfes = meshboundary(fes);
            icl = selectelem(fens, bdryfes, box=[1.0, 1.0, 0.0, pi/2, 0.0, Thickness], inflate=tolerance);
            
            for i=1:count(fens)
                t=fens.xyz[i,1]; a=fens.xyz[i,2]; z=fens.xyz[i,3]
                fens.xyz[i,:] = [(t*Re+(1-t)*Ri)*cos(a), (t*Re+(1-t)*Ri)*sin(a), z];
            end
            
            # File =  "a.vtk"
            # vtkexportmesh(File, fes.conn, fens.xyz,
            #     FinEtools.MeshExportModule.H8)
            # @async run(`"paraview.exe" $File`)
            
            geom = NodalField(fens.xyz)
            u = NodalField(zeros(size(fens.xyz,1),3)) # displacement field
            
            l1 =selectnode(fens; box=[0.0, Inf, 0.0, 0.0, 0.0, Thickness], inflate = tolerance)
            setebc!(u,l1,true, 2, 0.0)
            l1 =selectnode(fens; box=[0.0, 0.0, 0.0, Inf, 0.0, Thickness], inflate = tolerance)
            setebc!(u,l1,true, 1, 0.0)
            # Plane-stress constraint: assume the plane z=0 is the plane of symmetry of the plate
            l1 =selectnode(fens; box=[0.0, Inf, 0.0, Inf, 0.0, 0.0], inflate = tolerance)
            setebc!(u,l1,true, 3, 0.0)
            # If this was enabled, the plane-strain  constraint would be enforced.
            l1 =selectnode(fens; box=[0.0, Inf, 0.0, Inf, Thickness, Thickness], inflate = tolerance)
            setebc!(u,l1,true, 3, 0.0)
            
            applyebc!(u)
            numberdofs!(u)
            el1femm =  FEMMBase(IntegData(subset(bdryfes,icl), GaussRule(2, 2)))
            function pfun(forceout::FVec{T}, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt) where {T}
                local r = sqrt(XYZ[1]^2 + XYZ[2]^2)
                nx = XYZ[1]/r; ny = XYZ[2]/r
                # local sx, sy, txy
                # sx, sy, txy = sigmaxx(XYZ), sigmayy(XYZ), sigmaxy(XYZ)
                # sn = sx * nx^2 + sy * ny^2 + 2 * nx * ny * txy
                # tn = -(sx - sy) * nx * ny + (nx^2 - ny^2) * txy
                # forceout[1] = sn * nx - tn * ny
                # forceout[2] = sn * ny + tn * nx
                # forceout[3] = 0.0
                forceout[1] = sigmarr(XYZ) * nx - sigmart(XYZ) * ny
                forceout[2] = sigmarr(XYZ) * ny + sigmart(XYZ) * nx
                forceout[3] = 0.0
                return forceout
            end
            fi = ForceIntensity(FFlt, 3, pfun);
            F2 = distribloads(el1femm, geom, u, fi, 2);
            
            MR = DeforModelRed3D
            
            material = MatDeforElastIso(MR, E, nu)
            
            femm = FEMMDeforLinearMSH8(MR, IntegData(fes, GaussRule(3, 2)), material)
            
            # The geometry field now needs to be associated with the FEMM
            femm = associategeometry!(femm, geom)
            
            K = stiffness(femm, geom, u)
            K = cholesky(K)
            U = K\(F2)
            scattersysvec!(u, U[:])
            
            nlA = selectnode(fens, box=[Ri, Ri, 0.0, 0.0, 0.0, Thickness], inflate=tolerance);
            nlB = selectnode(fens, box=[0.0, 0.0, Ri, Ri, 0.0, Thickness], inflate=tolerance);
            # thecorneru = zeros(FFlt,length(nlA),3)
            # gathervalues_asmat!(u, thecorneru, nl);
            # thecorneru = mean(thecorneru, 1)[1]/phun("mm")
            # println("displacement = $(thecorneru) vs -0.10215 [MM]")
            
            println("Extrapolation: $( extrapolation )")
            sigx = fieldfromintegpoints(femm, geom, u, :Cauchy, 1;
            nodevalmethod = :averaging, reportat = extrapolation)
            sigy = fieldfromintegpoints(femm, geom, u, :Cauchy, 2;
            nodevalmethod = :averaging, reportat = extrapolation)
            sigyA = mean(sigy.values[nlA,1], 1)[1]
            sigyAtrue = sigmatt([Ri, 0.0, 0.0])
            println("sig_y@A =$(sigyA/phun("MPa")) vs $(sigyAtrue/phun("MPa")) [MPa]")
            sigxB = mean(sigx.values[nlB,1], 1)[1]
            sigxBtrue = sigmatt([0.0, Ri, 0.0])
            println("sig_x@B =$(sigxB/phun("MPa")) vs $(sigxBtrue/phun("MPa")) [MPa]")
            push!(numnodes, count(fens))
            push!(numelements, count(fes))
            push!(sigxderrs[extrapolation], abs(sigxB/sigxBtrue - 1.0))
            push!(sigyderrs[extrapolation], abs(sigyA/sigyAtrue - 1.0))
            # File =  "a.vtk"
            # vtkexportmesh(File, fes.conn, geom.values,
            # FinEtools.MeshExportModule.H8; vectors=[("u", u.values)],
            # scalars=[("sigmax", sigx.values/phun("MEGA*PA"))])
            # @async run(`"paraview.exe" $File`)
        end
    end
    
    File = "plate_w_hole_PE_MSH8_convergence.CSV"
    savecsv(File, numelements=vec(numelements), numnodes=vec(numnodes), sigxderrtrend=vec(sigxderrs[:extraptrend]), sigxderrdefault=vec(sigxderrs[:extrapmean]), sigyderrtrend=vec(sigyderrs[:extraptrend]), sigyderrdefault=vec(sigyderrs[:extrapmean]))
    
end # plate_w_hole_MSH8_PE_convergence


function plate_w_hole_MST10_convergence()
    E = 2.4*phun("MEGA*PA");# 210e3 MPa
    nu = 0.49995;
    Re = 0.3*phun("M"); # outer radius
    Ri= 0.1*phun("M"); # hole radius
    H = 0.1*phun("M") # thickness of the plate
    nRadial, nCircumferential=3, 5;
    sigma0=1*phun("MEGA*PA");
    
    function sigmaxx(x)
        local r = norm(x[1:2]);
        local th = atan2(x[2],x[1]);
        return sigma0*(1-Ri^2/r^2*(3/2*cos(2*th)+cos(4*th))+3/2*Ri^4/r^4*cos(4*th));
    end
    function sigmayy(x)
        local r=norm(x[1:2]);
        local th = atan2(x[2],x[1]);
        return -sigma0*(Ri^2/r^2*(1/2*cos(2*th)-cos(4*th))+3/2*Ri^4/r^4*cos(4*th));
    end
    function sigmaxy(x)
        local r=norm(x[1:2]);
        local th = atan2(x[2],x[1]);
        return -sigma0*(Ri^2/r^2*(1/2*sin(2*th)+sin(4*th))-3/2*Ri^4/r^4*sin(4*th));
    end
    function sigmarr(x)
        local r = norm(x[1:2]);
        local th = atan2(x[2],x[1]);
        return sigma0/2*(1-Ri^2/r^2) + sigma0/2*(1-4*Ri^2/r^2+3*Ri^4/r^4)*cos(2*th)
    end
    function sigmatt(x)
        local r = norm(x[1:2]);
        local th = atan2(x[2],x[1]);
        return sigma0/2*(1+Ri^2/r^2) - sigma0/2*(1+3*Ri^4/r^4)*cos(2*th)
    end
    function sigmart(x)
        local r = norm(x[1:2]);
        local th = atan2(x[2],x[1]);
        return -sigma0/2*(1+2*Ri^2/r^2-3*Ri^4/r^4)*sin(2*th)
    end
    
    sigxderrs = Dict{Symbol, FFltVec}()
    sigyderrs = Dict{Symbol, FFltVec}()
    numelements = []
    numnodes = []
    for extrapolation in [:extraptrend :extrapmean]
        sigxderrs[extrapolation] = FFltVec[]
        sigyderrs[extrapolation] = FFltVec[]
        numelements = []
        numnodes = []
        for ref in 0:1:5
            # Thickness = H
            Thickness = H/2^ref
            tolerance = Thickness/2^ref/1000.; # Geometrical tolerance
            
            fens,fes = T10block(1.0, pi/2, Thickness, 2^ref*nRadial, 2^ref*nCircumferential, 1)
            
            bdryfes = meshboundary(fes);
            icl = selectelem(fens, bdryfes, box=[1.0, 1.0, 0.0, pi/2, 0.0, Thickness], inflate=tolerance);
            
            for i=1:count(fens)
                t=fens.xyz[i,1]; a=fens.xyz[i,2]; z=fens.xyz[i,3]
                fens.xyz[i,:] = [(t*Re+(1-t)*Ri)*cos(a), (t*Re+(1-t)*Ri)*sin(a), z];
            end
            
            # File =  "a.vtk"
            # vtkexportmesh(File, fes.conn, fens.xyz,
            #     FinEtools.MeshExportModule.H8)
            # @async run(`"paraview.exe" $File`)
            
            geom = NodalField(fens.xyz)
            u = NodalField(zeros(size(fens.xyz,1),3)) # displacement field
            
            l1 =selectnode(fens; box=[0.0, Inf, 0.0, 0.0, 0.0, Thickness], inflate = tolerance)
            setebc!(u,l1,true, 2, 0.0)
            l1 =selectnode(fens; box=[0.0, 0.0, 0.0, Inf, 0.0, Thickness], inflate = tolerance)
            setebc!(u,l1,true, 1, 0.0)
            # Plane-stress constraint: assume the plane z=0 is the plane of symmetry of the plate
            l1 =selectnode(fens; box=[0.0, Inf, 0.0, Inf, 0.0, 0.0], inflate = tolerance)
            setebc!(u,l1,true, 3, 0.0)
            # If this was enabled, the plane-strain  constraint would be enforced.
            # l1 =selectnode(fens; box=[0.0, Inf, 0.0, Inf, Thickness, Thickness], inflate = tolerance)
            # setebc!(u,l1,true, 3, 0.0)
            
            applyebc!(u)
            numberdofs!(u)
            el1femm =  FEMMBase(IntegData(subset(bdryfes,icl), TriRule(3)))
            function pfun(forceout::FVec{T}, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt) where {T}
                local r = sqrt(XYZ[1]^2 + XYZ[2]^2)
                nx = XYZ[1]/r; ny = XYZ[2]/r
                # local sx, sy, txy
                # sx, sy, txy = sigmaxx(XYZ), sigmayy(XYZ), sigmaxy(XYZ)
                # sn = sx * nx^2 + sy * ny^2 + 2 * nx * ny * txy
                # tn = -(sx - sy) * nx * ny + (nx^2 - ny^2) * txy
                # forceout[1] = sn * nx - tn * ny
                # forceout[2] = sn * ny + tn * nx
                # forceout[3] = 0.0
                forceout[1] = sigmarr(XYZ) * nx - sigmart(XYZ) * ny
                forceout[2] = sigmarr(XYZ) * ny + sigmart(XYZ) * nx
                forceout[3] = 0.0
                return forceout
            end
            fi = ForceIntensity(FFlt, 3, pfun);
            F2 = distribloads(el1femm, geom, u, fi, 2);
            
            MR = DeforModelRed3D
            
            material = MatDeforElastIso(MR, E, nu)
            
            femm = FEMMDeforLinearMST10(MR, IntegData(fes, TetRule(4)), material)
            
            # The geometry field now needs to be associated with the FEMM
            femm = associategeometry!(femm, geom)
            
            K = stiffness(femm, geom, u)
            K = cholesky(K)
            U = K\(F2)
            scattersysvec!(u, U[:])
            
            nlA = selectnode(fens, box=[Ri, Ri, 0.0, 0.0, 0.0, Thickness], inflate=tolerance);
            nlB = selectnode(fens, box=[0.0, 0.0, Ri, Ri, 0.0, Thickness], inflate=tolerance);
            # thecorneru = zeros(FFlt,length(nlA),3)
            # gathervalues_asmat!(u, thecorneru, nl);
            # thecorneru = mean(thecorneru, 1)[1]/phun("mm")
            # println("displacement = $(thecorneru) vs -0.10215 [MM]")
            
            println("Extrapolation: $( extrapolation )")
            sigx = fieldfromintegpoints(femm, geom, u, :Cauchy, 1;            nodevalmethod = :averaging, reportat = extrapolation)
            sigy = fieldfromintegpoints(femm, geom, u, :Cauchy, 2;            nodevalmethod = :averaging, reportat = extrapolation)
            sigyA = mean(sigy.values[nlA,1], 1)[1]
            sigyAtrue = sigmatt([Ri, 0.0, 0.0])
            println("sig_y@A =$(sigyA/phun("MPa")) vs $(sigyAtrue/phun("MPa")) [MPa]")
            sigxB = mean(sigx.values[nlB,1], 1)[1]
            sigxBtrue = sigmatt([0.0, Ri, 0.0])
            println("sig_x@B =$(sigxB/phun("MPa")) vs $(sigxBtrue/phun("MPa")) [MPa]")
            push!(numnodes, count(fens))
            push!(numelements, count(fes))
            push!(sigxderrs[extrapolation], abs(sigxB/sigxBtrue - 1.0))
            push!(sigyderrs[extrapolation], abs(sigyA/sigyAtrue - 1.0))
            # File =  "a.vtk"
            # vtkexportmesh(File, fes.conn, geom.values,
            # FinEtools.MeshExportModule.H8; vectors=[("u", u.values)],
            # scalars=[("sigmax", sigx.values/phun("MEGA*PA"))])
            # @async run(`"paraview.exe" $File`)
        end
    end
    
    File = "plate_w_hole_MST10_convergence.CSV"
    savecsv(File, numelements=vec(numelements), numnodes=vec(numnodes), sigxderrtrend=vec(sigxderrs[:extraptrend]), sigxderrdefault=vec(sigxderrs[:extrapmean]), sigyderrtrend=vec(sigyderrs[:extraptrend]), sigyderrdefault=vec(sigyderrs[:extrapmean]))
    
end # plate_w_hole_MST10_convergence


function plate_w_hole_MST10_PE_convergence()
    E = 2.4*phun("MEGA*PA");# 210e3 MPa
    nu = 0.49995;
    Re = 0.3*phun("M"); # outer radius
    Ri= 0.1*phun("M"); # hole radius
    H = 0.1*phun("M") # thickness of the plate
    nRadial, nCircumferential=3, 5;
    sigma0=1*phun("MEGA*PA");
    
    function sigmaxx(x)
        local r = norm(x[1:2]);
        local th = atan2(x[2],x[1]);
        return sigma0*(1-Ri^2/r^2*(3/2*cos(2*th)+cos(4*th))+3/2*Ri^4/r^4*cos(4*th));
    end
    function sigmayy(x)
        local r=norm(x[1:2]);
        local th = atan2(x[2],x[1]);
        return -sigma0*(Ri^2/r^2*(1/2*cos(2*th)-cos(4*th))+3/2*Ri^4/r^4*cos(4*th));
    end
    function sigmaxy(x)
        local r=norm(x[1:2]);
        local th = atan2(x[2],x[1]);
        return -sigma0*(Ri^2/r^2*(1/2*sin(2*th)+sin(4*th))-3/2*Ri^4/r^4*sin(4*th));
    end
    function sigmarr(x)
        local r = norm(x[1:2]);
        local th = atan2(x[2],x[1]);
        return sigma0/2*(1-Ri^2/r^2) + sigma0/2*(1-4*Ri^2/r^2+3*Ri^4/r^4)*cos(2*th)
    end
    function sigmatt(x)
        local r = norm(x[1:2]);
        local th = atan2(x[2],x[1]);
        return sigma0/2*(1+Ri^2/r^2) - sigma0/2*(1+3*Ri^4/r^4)*cos(2*th)
    end
    function sigmart(x)
        local r = norm(x[1:2]);
        local th = atan2(x[2],x[1]);
        return -sigma0/2*(1+2*Ri^2/r^2-3*Ri^4/r^4)*sin(2*th)
    end
    
    sigxderrs = Dict{Symbol, FFltVec}()
    sigyderrs = Dict{Symbol, FFltVec}()
    numelements = []
    numnodes = []
    for extrapolation in [:extraptrend :extrapmean]
        sigxderrs[extrapolation] = FFltVec[]
        sigyderrs[extrapolation] = FFltVec[]
        numelements = []
        numnodes = []
        for ref in 0:1:5
            # Thickness = H
            Thickness = H/2^ref
            tolerance = Thickness/2^ref/1000.; # Geometrical tolerance
            
            fens,fes = T10block(1.0, pi/2, Thickness, 2^ref*nRadial, 2^ref*nCircumferential, 1)
            
            bdryfes = meshboundary(fes);
            icl = selectelem(fens, bdryfes, box=[1.0, 1.0, 0.0, pi/2, 0.0, Thickness], inflate=tolerance);
            
            for i=1:count(fens)
                t=fens.xyz[i,1]; a=fens.xyz[i,2]; z=fens.xyz[i,3]
                fens.xyz[i,:] = [(t*Re+(1-t)*Ri)*cos(a), (t*Re+(1-t)*Ri)*sin(a), z];
            end
            
            # File =  "a.vtk"
            # vtkexportmesh(File, fes.conn, fens.xyz,
            #     FinEtools.MeshExportModule.H8)
            # @async run(`"paraview.exe" $File`)
            
            geom = NodalField(fens.xyz)
            u = NodalField(zeros(size(fens.xyz,1),3)) # displacement field
            
            l1 =selectnode(fens; box=[0.0, Inf, 0.0, 0.0, 0.0, Thickness], inflate = tolerance)
            setebc!(u,l1,true, 2, 0.0)
            l1 =selectnode(fens; box=[0.0, 0.0, 0.0, Inf, 0.0, Thickness], inflate = tolerance)
            setebc!(u,l1,true, 1, 0.0)
            # Plane-stress constraint: assume the plane z=0 is the plane of symmetry of the plate
            l1 =selectnode(fens; box=[0.0, Inf, 0.0, Inf, 0.0, 0.0], inflate = tolerance)
            setebc!(u,l1,true, 3, 0.0)
            # If this was enabled, the plane-strain  constraint would be enforced.
            l1 =selectnode(fens; box=[0.0, Inf, 0.0, Inf, Thickness, Thickness], inflate = tolerance)
            setebc!(u,l1,true, 3, 0.0)
            
            applyebc!(u)
            numberdofs!(u)
            el1femm =  FEMMBase(IntegData(subset(bdryfes,icl), TriRule(3)))
            function pfun(forceout::FVec{T}, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt) where {T}
                local r = sqrt(XYZ[1]^2 + XYZ[2]^2)
                nx = XYZ[1]/r; ny = XYZ[2]/r
                # local sx, sy, txy
                # sx, sy, txy = sigmaxx(XYZ), sigmayy(XYZ), sigmaxy(XYZ)
                # sn = sx * nx^2 + sy * ny^2 + 2 * nx * ny * txy
                # tn = -(sx - sy) * nx * ny + (nx^2 - ny^2) * txy
                # forceout[1] = sn * nx - tn * ny
                # forceout[2] = sn * ny + tn * nx
                # forceout[3] = 0.0
                forceout[1] = sigmarr(XYZ) * nx - sigmart(XYZ) * ny
                forceout[2] = sigmarr(XYZ) * ny + sigmart(XYZ) * nx
                forceout[3] = 0.0
                return forceout
            end
            fi = ForceIntensity(FFlt, 3, pfun);
            F2 = distribloads(el1femm, geom, u, fi, 2);
            
            MR = DeforModelRed3D
            
            material = MatDeforElastIso(MR, E, nu)
            
            femm = FEMMDeforLinearMST10(MR, IntegData(fes, TetRule(4)), material)
            
            # The geometry field now needs to be associated with the FEMM
            femm = associategeometry!(femm, geom)
            
            K = stiffness(femm, geom, u)
            K = cholesky(K)
            U = K\(F2)
            scattersysvec!(u, U[:])
            
            nlA = selectnode(fens, box=[Ri, Ri, 0.0, 0.0, 0.0, Thickness], inflate=tolerance);
            nlB = selectnode(fens, box=[0.0, 0.0, Ri, Ri, 0.0, Thickness], inflate=tolerance);
            # thecorneru = zeros(FFlt,length(nlA),3)
            # gathervalues_asmat!(u, thecorneru, nl);
            # thecorneru = mean(thecorneru, 1)[1]/phun("mm")
            # println("displacement = $(thecorneru) vs -0.10215 [MM]")
            
            println("Extrapolation: $( extrapolation )")
            sigx = fieldfromintegpoints(femm, geom, u, :Cauchy, 1;
            nodevalmethod = :averaging, reportat = extrapolation)
            sigy = fieldfromintegpoints(femm, geom, u, :Cauchy, 2;
            nodevalmethod = :averaging, reportat = extrapolation)
            sigyA = mean(sigy.values[nlA,1], 1)[1]
            sigyAtrue = sigmatt([Ri, 0.0, 0.0])
            println("sig_y@A =$(sigyA/phun("MPa")) vs $(sigyAtrue/phun("MPa")) [MPa]")
            sigxB = mean(sigx.values[nlB,1], 1)[1]
            sigxBtrue = sigmatt([0.0, Ri, 0.0])
            println("sig_x@B =$(sigxB/phun("MPa")) vs $(sigxBtrue/phun("MPa")) [MPa]")
            push!(numnodes, count(fens))
            push!(numelements, count(fes))
            push!(sigxderrs[extrapolation], abs(sigxB/sigxBtrue - 1.0))
            push!(sigyderrs[extrapolation], abs(sigyA/sigyAtrue - 1.0))
            # File =  "a.vtk"
            # vtkexportmesh(File, fes.conn, geom.values,
            # FinEtools.MeshExportModule.H8; vectors=[("u", u.values)],
            # scalars=[("sigmax", sigx.values/phun("MEGA*PA"))])
            # @async run(`"paraview.exe" $File`)
        end
    end
    
    File = "plate_w_hole_PE_MST10_convergence.CSV"
    savecsv(File, numelements=vec(numelements), numnodes=vec(numnodes), sigxderrtrend=vec(sigxderrs[:extraptrend]), sigxderrdefault=vec(sigxderrs[:extrapmean]), sigyderrtrend=vec(sigyderrs[:extraptrend]), sigyderrdefault=vec(sigyderrs[:extrapmean]))
    
end # plate_w_hole_MST10_PE_convergence


function plate_w_hole_MST10_stress()
    E = 2.4*phun("MEGA*PA");# 210e3 MPa
    nu = 0.3;
    Re = 0.3*phun("M"); # outer radius
    Ri= 0.1*phun("M"); # hole radius
    H = 0.1*phun("M") # thickness of the plate
    nRadial, nCircumferential, nThickness = 6, 8, 1;
    sigma0=1*phun("MEGA*PA");
    
    function sigmaxx(x)
        local r = norm(x[1:2]);
        local th = atan2(x[2],x[1]);
        return sigma0*(1-Ri^2/r^2*(3/2*cos(2*th)+cos(4*th))+3/2*Ri^4/r^4*cos(4*th));
    end
    function sigmayy(x)
        local r=norm(x[1:2]);
        local th = atan2(x[2],x[1]);
        return -sigma0*(Ri^2/r^2*(1/2*cos(2*th)-cos(4*th))+3/2*Ri^4/r^4*cos(4*th));
    end
    function sigmaxy(x)
        local r=norm(x[1:2]);
        local th = atan2(x[2],x[1]);
        return -sigma0*(Ri^2/r^2*(1/2*sin(2*th)+sin(4*th))-3/2*Ri^4/r^4*sin(4*th));
    end
    function sigmarr(x)
        local r = norm(x[1:2]);
        local th = atan2(x[2],x[1]);
        return sigma0/2*(1-Ri^2/r^2) + sigma0/2*(1-4*Ri^2/r^2+3*Ri^4/r^4)*cos(2*th)
    end
    function sigmatt(x)
        local r = norm(x[1:2]);
        local th = atan2(x[2],x[1]);
        return sigma0/2*(1+Ri^2/r^2) - sigma0/2*(1+3*Ri^4/r^4)*cos(2*th)
    end
    function sigmart(x)
        local r = norm(x[1:2]);
        local th = atan2(x[2],x[1]);
        return -sigma0/2*(1+2*Ri^2/r^2-3*Ri^4/r^4)*sin(2*th)
    end
    
    convergencestudy = FDataDict[]
    for ref in 0:1:2
        println("ref = $(ref)")
        # Thickness = H
        Thickness = H/2^ref
        tolerance = Thickness/2^ref/1000.; # Geometrical tolerance
        
        fens,fes = T10block(1.0, pi/2, Thickness, 2^ref*nRadial, 2^ref*nCircumferential, 2^ref*nThickness)
        
        bdryfes = meshboundary(fes);
        icl = selectelem(fens, bdryfes, box=[1.0, 1.0, 0.0, pi/2, 0.0, Thickness], inflate=tolerance);
        
        for i=1:count(fens)
            t=fens.xyz[i,1]; a=fens.xyz[i,2]; z=fens.xyz[i,3]
            fens.xyz[i,:] = [(t*Re+(1-t)*Ri)*cos(a), (t*Re+(1-t)*Ri)*sin(a), z];
        end
        
        geom = NodalField(fens.xyz)
        u = NodalField(zeros(size(fens.xyz,1),3)) # displacement field
        
        l1 =selectnode(fens; box=[0.0, Inf, 0.0, 0.0, 0.0, Thickness], inflate = tolerance)
        setebc!(u,l1,true, 2, 0.0)
        l1 =selectnode(fens; box=[0.0, 0.0, 0.0, Inf, 0.0, Thickness], inflate = tolerance)
        setebc!(u,l1,true, 1, 0.0)
        # Plane-stress constraint: assume the plane z=0 is the plane of symmetry of the plate
        l1 =selectnode(fens; box=[0.0, Inf, 0.0, Inf, 0.0, 0.0], inflate = tolerance)
        setebc!(u,l1,true, 3, 0.0)
        # If this was enabled, the plane-strain  constraint would be enforced.
        # l1 =selectnode(fens; box=[0.0, Inf, 0.0, Inf, Thickness, Thickness], inflate = tolerance)
        # setebc!(u,l1,true, 3, 0.0)
        
        applyebc!(u)
        numberdofs!(u)
        el1femm =  FEMMBase(IntegData(subset(bdryfes,icl), TriRule(3)))
        function pfun(forceout::FVec{T}, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt) where {T}
            local r = sqrt(XYZ[1]^2 + XYZ[2]^2)
            nx = XYZ[1]/r; ny = XYZ[2]/r
            forceout[1] = sigmarr(XYZ) * nx - sigmart(XYZ) * ny
            forceout[2] = sigmarr(XYZ) * ny + sigmart(XYZ) * nx
            forceout[3] = 0.0
            return forceout
        end
        fi = ForceIntensity(FFlt, 3, pfun);
        F2 = distribloads(el1femm, geom, u, fi, 2);
        
        MR = DeforModelRed3D
        
        material = MatDeforElastIso(MR, E, nu)
        
        femm = FEMMDeforLinearMST10(MR, IntegData(fes, TetRule(4)), material)
        
        # The geometry field now needs to be associated with the FEMM
        femm = associategeometry!(femm, geom)
        
        K = stiffness(femm, geom, u)
        K = cholesky(K)
        U = K\(F2)
        scattersysvec!(u, U[:])
        
        stressfields = elemfieldfromintegpoints(femm, geom, u, :Cauchy, collect(1:6))
        
        push!(convergencestudy, FDataDict(
        "elementsize"=> 1.0 / 2^ref,
        "fens"=>fens,
        "fes"=>fes,
        "geom"=>geom,
        "u"=>u,
        "femm"=>femm,
        "stressfields"=>[stressfields],
        "tolerance"=>tolerance)
        )
    end # for ref in 0:1:5
    
    # File = "mplate_w_hole_MST10m_stress"
    # open(File * ".jls", "w") do file
    #     serialize(file, convergencestudy)
    # end
    
end # plate_w_hole_MST10_stress


function plate_w_hole_RECT_H20_convergence()
    E = 210000*phun("MEGA*PA");# 210e3 MPa
    nu = 0.3;
    Ri= 0.15*phun("M"); # hole radius
    Re = 2*Ri; # outer radius
    H = 0.01*phun("M") # thickness of the plate
    nRadial, nCircumferential=6, 3;
    sigma0=1*phun("MEGA*PA");
    
    function sigmaxx(x)
        local r = norm(vec(x[1:2]));
        local th = atan2(x[2],x[1]);
        return sigma0*(1-Ri^2/r^2*(3/2*cos(2*th)+cos(4*th))+3/2*Ri^4/r^4*cos(4*th));
    end
    function sigmayy(x)
        local r = norm(vec(x[1:2]));
        local th = atan2(x[2],x[1]);
        return -sigma0*(Ri^2/r^2*(1/2*cos(2*th)-cos(4*th))+3/2*Ri^4/r^4*cos(4*th));
    end
    function sigmaxy(x)
        local r = norm(vec(x[1:2]));
        local th = atan2(x[2],x[1]);
        return -sigma0*(Ri^2/r^2*(1/2*sin(2*th)+sin(4*th))-3/2*Ri^4/r^4*sin(4*th));
    end
    function sigmarr(x)
        local r = norm(vec(x[1:2]));
        local th = atan2(x[2],x[1]);
        return sigma0/2*(1-Ri^2/r^2) + sigma0/2*(1-4*Ri^2/r^2+3*Ri^4/r^4)*cos(2*th)
    end
    function sigmatt(x)
        local r = norm(vec(x[1:2]));
        local th = atan2(x[2],x[1]);
        return sigma0/2*(1+Ri^2/r^2) - sigma0/2*(1+3*Ri^4/r^4)*cos(2*th)
    end
    function sigmart(x)
        local r = norm(vec(x[1:2]));
        local th = atan2(x[2],x[1]);
        return -sigma0/2*(1+2*Ri^2/r^2-3*Ri^4/r^4)*sin(2*th)
    end
    
    sigyderrs = Dict{Symbol, FFltVec}()
    
    nelems = []
    for extrapolation in [:meanonly]
        sigyderrs[extrapolation] = FFltVec[]
        nelems = []
        for ref in [1]
            Thickness = H
            # Thickness = H/2^ref
            tolerance = Thickness/2^ref/1000.; # Geometrical tolerance
            
            fens,fes = H8elliphole(Ri, Ri, Re, Re, Thickness,
            2^ref*nCircumferential, 2^ref*nCircumferential, 2^ref*nRadial, 1)
            fens,fes = H8toH20(fens,fes)
            # File =  "a.vtk"
            # vtkexportmesh(File, fes.conn, fens.xyz,
            #     FinEtools.MeshExportModule.H20)
            # @async run(`"paraview.cexe" $File`)
            
            println("My mesh=>$((count(fens), count(fes)))")
            #
            # output = import_ABAQUS("plane_w_hole_m_debug.inp")
            # fens1,fes1 = output["fens"], output["fesets"][1]
            # println("Matlab mesh=>$((count(fens1), count(fes1[1])))")
            #
            #  fens3, newfes1, fes2 = mergemeshes(fens,fes, fens1,fes1[1], tolerance)
            #  fes3 = cat(2, newfes1)
            #  println("Merged mesh=>$((count(fens3), count(fes3)))")
            
            geom = NodalField(fens.xyz)
            u = NodalField(zeros(size(fens.xyz,1),3)) # displacement field
            
            l1 =selectnode(fens; box=[0.0, Inf, 0.0, 0.0, 0.0, Thickness], inflate = tolerance)
            setebc!(u,l1,true, 2, 0.0)
            l1 =selectnode(fens; box=[0.0, 0.0, 0.0, Inf, 0.0, Thickness], inflate = tolerance)
            setebc!(u,l1,true, 1, 0.0)
            l1 =selectnode(fens; box=[0.0, Inf, 0.0, Inf, 0.0, 0.0], inflate = tolerance)
            setebc!(u,l1,true, 3, 0.0)
            # l1 =selectnode(fens; box=[0.0, Inf, 0.0, Inf, Thickness, Thickness], inflate = tolerance)
            # setebc!(u,l1,true, 3, 0.0)
            
            applyebc!(u)
            numberdofs!(u)
            
            
            bdryfes = meshboundary(fes);
            # ixl = selectelem(fens, bdryfes, plane=[1.0, 0.0, 0.0, Re], thickness=tolerance);
            ixl = selectelem(fens, bdryfes, box=[Re, Re, -Inf, +Inf, -Inf, +Inf], inflate = tolerance);
            elxfemm =  FEMMBase(IntegData(subset(bdryfes,ixl), GaussRule(2, 2)))
            function pfunx(forceout::FVec{T}, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt) where {T}
                forceout[1] = sigmaxx(XYZ)
                forceout[2] = sigmaxy(XYZ)
                forceout[3] = 0.0
                return forceout
            end
            fi = ForceIntensity(FFlt, 3, pfunx);
            Fx = distribloads(elxfemm, geom, u, fi, 2);
            # iyl = selectelem(fens, bdryfes, plane=[0.0, 1.0, 0.0, Re], thickness=tolerance);
            iyl = selectelem(fens, bdryfes, box=[-Inf, +Inf, Re, Re, -Inf, +Inf], inflate = tolerance);
            elyfemm =  FEMMBase(IntegData(subset(bdryfes,iyl), GaussRule(2, 2)))
            function pfuny(forceout::FVec{T}, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt) where {T}
                forceout[1] = sigmaxy(XYZ)
                forceout[2] = sigmayy(XYZ)
                forceout[3] = 0.0
                return forceout
            end
            fi = ForceIntensity(FFlt, 3, pfuny);
            Fy = distribloads(elyfemm, geom, u, fi, 2);
            
            MR = DeforModelRed3D
            
            material = MatDeforElastIso(MR, E, nu)
            
            femm = FEMMDeforLinear(MR, IntegData(fes, GaussRule(3, 2)), material)
            
            # The geometry field now needs to be associated with the FEMM
            femm = associategeometry!(femm, geom)
            
            K = stiffness(femm, geom, u)
            K = cholesky(K)
            U = K\(Fx + Fy)
            scattersysvec!(u, U[:])
            println("oof load = $(norm(Fx + Fy, 2))")
            
            nlA = selectnode(fens, box=[Ri, Ri, 0.0, 0.0, 0.0, 00.0], inflate=tolerance);
            pointu = zeros(FFlt,length(nlA),3)
            gathervalues_asmat!(u, pointu, nlA);
            println("disp@A = $(pointu/phun("mm")) [MM]")
            nlB = selectnode(fens, box=[0.0, 0.0, Ri, Ri, 0.0, 0.0], inflate=tolerance);
            pointu = zeros(FFlt,length(nlB),3)
            gathervalues_asmat!(u, pointu, nlB);
            println("disp@B = $(pointu/phun("mm")) [MM]")
            nlC = selectnode(fens, box=[Re, Re, Re, Re, Thickness, Thickness], inflate=tolerance);
            pointu = zeros(FFlt,length(nlC),3)
            gathervalues_asmat!(u, pointu, nlC);
            println("disp@C = $(pointu/phun("mm")) [MM]")
            
            nlAallz = selectnode(fens, box=[Ri, Ri, 0.0, 0.0, 0.0, Thickness], inflate=tolerance);
            nlBallz = selectnode(fens, box=[0.0, 0.0, Ri, Ri, 0.0, Thickness], inflate=tolerance);
            sigx = fieldfromintegpoints(femm, geom, u, :Cauchy, 1; nodevalmethod = :invdistance, reportat = extrapolation)
            sigy = fieldfromintegpoints(femm, geom, u, :Cauchy, 2; nodevalmethod = :invdistance, reportat = extrapolation)
            sigyA = mean(sigy.values[nlAallz,1], 1)[1]
            sigyAtrue = sigmayy([Ri, 0.0, 0.0])
            println("sig_y@A =$(sigyA/phun("MPa")) vs $(sigyAtrue/phun("MPa")) [MPa]")
            sigxB = mean(sigx.values[nlBallz,1], 1)[1]
            sigxBtrue = sigmaxx([0.0, Ri, 0.0])
            println("sig_x@B =$(sigxB/phun("MPa")) vs $(sigxBtrue/phun("MPa")) [MPa]")
            # println("$extrapolation, $(count(fes)), $(sigyd/phun("MPa"))")
            # push!(nelems, count(fes))
            # push!(sigyderrs[extrapolation], abs(sigyd/sigma_yD - 1.0))
            File =  "a.vtk"
            vtkexportmesh(File, fes.conn, geom.values,
            FinEtools.MeshExportModule.H20; vectors=[("u", u.values)],
            scalars=[("sigmax", sigx.values/phun("MEGA*PA")),
            ("sigmay", sigy.values/phun("MEGA*PA"))])
            @async run(`"paraview.exe" $File`)
        end
    end
    
    # df = DataFrame(nelems=vec(nelems),
    #     sigyderrtrend=vec(sigyderrs[:extraptrend]),
    #     sigyderrdefault=vec(sigyderrs[:extrapmean]))
    # File = "LE1NAFEMS_MSH8_convergence.CSV"
    # CSV.write(File, df)
    # @async run(`"paraview.exe" $File`)
    
end # plate_w_hole_RECT_H20_convergence


function plate_w_hole_RECT_MSH8_convergence()
    E = 210000*phun("MEGA*PA");# 210e3 MPa
    nu = 0.3;
    Ri= 0.1*phun("M"); # hole radius
    Re = 2*Ri; # outer radius
    H = 0.01*phun("M") # thickness of the plate
    nRadial, nCircumferential=6, 3;
    sigma0=1*phun("MEGA*PA");
    
    function sigmaxx(x)
        local r = norm(x[1:2]);
        local th = atan2(x[2],x[1]);
        return sigma0*(1-Ri^2/r^2*(3/2*cos(2*th)+cos(4*th))+3/2*Ri^4/r^4*cos(4*th));
    end
    function sigmayy(x)
        local r=norm(x[1:2]);
        local th = atan2(x[2],x[1]);
        return -sigma0*(Ri^2/r^2*(1/2*cos(2*th)-cos(4*th))+3/2*Ri^4/r^4*cos(4*th));
    end
    function sigmaxy(x)
        local r=norm(x[1:2]);
        local th = atan2(x[2],x[1]);
        return -sigma0*(Ri^2/r^2*(1/2*sin(2*th)+sin(4*th))-3/2*Ri^4/r^4*sin(4*th));
    end
    function sigmarr(x)
        local r = norm(x[1:2]);
        local th = atan2(x[2],x[1]);
        return sigma0/2*(1-Ri^2/r^2) + sigma0/2*(1-4*Ri^2/r^2+3*Ri^4/r^4)*cos(2*th)
    end
    function sigmatt(x)
        local r = norm(x[1:2]);
        local th = atan2(x[2],x[1]);
        return sigma0/2*(1+Ri^2/r^2) - sigma0/2*(1+3*Ri^4/r^4)*cos(2*th)
    end
    function sigmart(x)
        local r = norm(x[1:2]);
        local th = atan2(x[2],x[1]);
        return -sigma0/2*(1+2*Ri^2/r^2-3*Ri^4/r^4)*sin(2*th)
    end
    
    sigyderrs = Dict{Symbol, FFltVec}()
    
    nelems = []
    for extrapolation in [:extraptrend :extrapmean]
        sigyderrs[extrapolation] = FFltVec[]
        nelems = []
        for ref in [2]
            # Thickness = H
            Thickness = H/2^ref
            tolerance = Thickness/2^ref/1000.; # Geometrical tolerance
            
            fens,fes = H8elliphole(Ri, Ri, Re, Re, Thickness,
            2^ref*nCircumferential, 2^ref*nCircumferential, 2^ref*nRadial, 1)
            
            # File =  "a.vtk"
            # vtkexportmesh(File, fes.conn, fens.xyz,
            #     FinEtools.MeshExportModule.H8)
            # @async run(`"paraview.exe" $File`)
            
            geom = NodalField(fens.xyz)
            u = NodalField(zeros(size(fens.xyz,1),3)) # displacement field
            
            l1 =selectnode(fens; box=[0.0, Inf, 0.0, 0.0, 0.0, Thickness], inflate = tolerance)
            setebc!(u,l1,true, 2, 0.0)
            l1 =selectnode(fens; box=[0.0, 0.0, 0.0, Inf, 0.0, Thickness], inflate = tolerance)
            setebc!(u,l1,true, 1, 0.0)
            l1 =selectnode(fens; box=[0.0, Inf, 0.0, Inf, 0.0, 0.0], inflate = tolerance)
            setebc!(u,l1,true, 3, 0.0)
            # l1 =selectnode(fens; box=[0.0, Inf, 0.0, Inf, Thickness, Thickness], inflate = tolerance)
            # setebc!(u,l1,true, 3, 0.0)
            
            applyebc!(u)
            numberdofs!(u)
            
            
            bdryfes = meshboundary(fes);
            ixl = selectelem(fens, bdryfes, plane=[1.0, 0.0, 0.0, Re], thickness=tolerance);
            elxfemm =  FEMMBase(IntegData(subset(bdryfes,ixl), GaussRule(2, 2)))
            function pfunx(forceout::FVec{T}, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt) where {T}
                forceout[1] = sigmaxx(XYZ)
                forceout[2] = sigmaxy(XYZ)
                forceout[3] = 0.0
                return forceout
            end
            fi = ForceIntensity(FFlt, 3, pfunx);
            Fx = distribloads(elxfemm, geom, u, fi, 2);
            iyl = selectelem(fens, bdryfes, plane=[0.0, 1.0, 0.0, Re], thickness=tolerance);
            elyfemm =  FEMMBase(IntegData(subset(bdryfes,iyl), GaussRule(2, 2)))
            function pfuny(forceout::FVec{T}, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt) where {T}
                forceout[1] = -sigmaxy(XYZ)
                forceout[2] = sigmayy(XYZ)
                forceout[3] = 0.0
                return forceout
            end
            fi = ForceIntensity(FFlt, 3, pfuny);
            Fy = distribloads(elyfemm, geom, u, fi, 2);
            
            MR = DeforModelRed3D
            
            material = MatDeforElastIso(MR, E, nu)
            
            femm = FEMMDeforLinearMSH8(MR, IntegData(fes, GaussRule(3, 2)), material)
            
            # The geometry field now needs to be associated with the FEMM
            femm = associategeometry!(femm, geom)
            
            K = stiffness(femm, geom, u)
            K = cholesky(K)
            U = K\(Fx + Fy)
            scattersysvec!(u, U[:])
            
            nlA = selectnode(fens, box=[Ri, Ri, 0.0, 0.0, 0.0, Thickness], inflate=tolerance);
            nlB = selectnode(fens, box=[0.0, 0.0, Ri, Ri, 0.0, Thickness], inflate=tolerance);
            # thecorneru = zeros(FFlt,length(nlA),3)
            # gathervalues_asmat!(u, thecorneru, nl);
            # thecorneru = mean(thecorneru, 1)[1]/phun("mm")
            # println("displacement = $(thecorneru) vs -0.10215 [MM]")
            
            println("Extrapolation: $( extrapolation )")
            sigx = fieldfromintegpoints(femm, geom, u, :Cauchy, 1; nodevalmethod = :averaging, reportat = extrapolation)
            sigy = fieldfromintegpoints(femm, geom, u, :Cauchy, 2; nodevalmethod = :averaging, reportat = extrapolation)
            sigyA = mean(sigy.values[nlA,1], 1)[1]
            sigyAtrue = sigmatt([Ri, 0.0, 0.0])
            println("sig_y@A =$(sigyA/phun("MPa")) vs $(sigyAtrue/phun("MPa")) [MPa]")
            sigxB = mean(sigx.values[nlB,1], 1)[1]
            sigxBtrue = sigmatt([0.0, Ri, 0.0])
            println("sig_x@B =$(sigxB/phun("MPa")) vs $(sigxBtrue/phun("MPa")) [MPa]")
            # println("$extrapolation, $(count(fes)), $(sigyd/phun("MPa"))")
            # push!(nelems, count(fes))
            # push!(sigyderrs[extrapolation], abs(sigyd/sigma_yD - 1.0))
            File =  "a.vtk"
            vtkexportmesh(File, fes.conn, geom.values,
            FinEtools.MeshExportModule.H8; vectors=[("u", u.values)],
            scalars=[("sigmax", sigx.values/phun("MEGA*PA")),
            ("sigmay", sigy.values/phun("MEGA*PA"))])
            @async run(`"paraview.exe" $File`)
        end
    end
    
    # df = DataFrame(nelems=vec(nelems),
    #     sigyderrtrend=vec(sigyderrs[:extraptrend]),
    #     sigyderrdefault=vec(sigyderrs[:extrapmean]))
    # File = "LE1NAFEMS_MSH8_convergence.CSV"
    # CSV.write(File, df)
    # @async run(`"paraview.exe" $File`)
    
end # plate_w_hole_RECT_MSH8_convergence


function plate_w_hole_T10_PE_convergence()
    E = 2.4*phun("MEGA*PA");# 210e3 MPa
    nu = 0.49995;
    Re = 0.3*phun("M"); # outer radius
    Ri= 0.1*phun("M"); # hole radius
    H = 0.1*phun("M") # thickness of the plate
    nRadial, nCircumferential=3, 5;
    sigma0=1*phun("MEGA*PA");
    
    function sigmaxx(x)
        local r = norm(x[1:2]);
        local th = atan2(x[2],x[1]);
        return sigma0*(1-Ri^2/r^2*(3/2*cos(2*th)+cos(4*th))+3/2*Ri^4/r^4*cos(4*th));
    end
    function sigmayy(x)
        local r=norm(x[1:2]);
        local th = atan2(x[2],x[1]);
        return -sigma0*(Ri^2/r^2*(1/2*cos(2*th)-cos(4*th))+3/2*Ri^4/r^4*cos(4*th));
    end
    function sigmaxy(x)
        local r=norm(x[1:2]);
        local th = atan2(x[2],x[1]);
        return -sigma0*(Ri^2/r^2*(1/2*sin(2*th)+sin(4*th))-3/2*Ri^4/r^4*sin(4*th));
    end
    function sigmarr(x)
        local r = norm(x[1:2]);
        local th = atan2(x[2],x[1]);
        return sigma0/2*(1-Ri^2/r^2) + sigma0/2*(1-4*Ri^2/r^2+3*Ri^4/r^4)*cos(2*th)
    end
    function sigmatt(x)
        local r = norm(x[1:2]);
        local th = atan2(x[2],x[1]);
        return sigma0/2*(1+Ri^2/r^2) - sigma0/2*(1+3*Ri^4/r^4)*cos(2*th)
    end
    function sigmart(x)
        local r = norm(x[1:2]);
        local th = atan2(x[2],x[1]);
        return -sigma0/2*(1+2*Ri^2/r^2-3*Ri^4/r^4)*sin(2*th)
    end
    
    sigxderrs = Dict{Symbol, FFltVec}()
    sigyderrs = Dict{Symbol, FFltVec}()
    numelements = []
    numnodes = []
    for extrapolation in [:extraptrend :extrapmean]
        sigxderrs[extrapolation] = FFltVec[]
        sigyderrs[extrapolation] = FFltVec[]
        numelements = []
        numnodes = []
        for ref in 0:1:2
            # Thickness = H
            Thickness = H/2^ref
            tolerance = Thickness/2^ref/1000.; # Geometrical tolerance
            
            fens,fes = T10block(1.0, pi/2, Thickness, 2^ref*nRadial, 2^ref*nCircumferential, 1)
            
            bdryfes = meshboundary(fes);
            icl = selectelem(fens, bdryfes, box=[1.0, 1.0, 0.0, pi/2, 0.0, Thickness], inflate=tolerance);
            
            for i=1:count(fens)
                t=fens.xyz[i,1]; a=fens.xyz[i,2]; z=fens.xyz[i,3]
                fens.xyz[i,:] = [(t*Re+(1-t)*Ri)*cos(a), (t*Re+(1-t)*Ri)*sin(a), z];
            end
            
            # File =  "a.vtk"
            # vtkexportmesh(File, fes.conn, fens.xyz,
            #     FinEtools.MeshExportModule.H8)
            # @async run(`"paraview.exe" $File`)
            
            geom = NodalField(fens.xyz)
            u = NodalField(zeros(size(fens.xyz,1),3)) # displacement field
            
            l1 =selectnode(fens; box=[0.0, Inf, 0.0, 0.0, 0.0, Thickness], inflate = tolerance)
            setebc!(u,l1,true, 2, 0.0)
            l1 =selectnode(fens; box=[0.0, 0.0, 0.0, Inf, 0.0, Thickness], inflate = tolerance)
            setebc!(u,l1,true, 1, 0.0)
            # Plane-stress constraint: assume the plane z=0 is the plane of symmetry of the plate
            l1 =selectnode(fens; box=[0.0, Inf, 0.0, Inf, 0.0, 0.0], inflate = tolerance)
            setebc!(u,l1,true, 3, 0.0)
            # If this was enabled, the plane-strain  constraint would be enforced.
            l1 =selectnode(fens; box=[0.0, Inf, 0.0, Inf, Thickness, Thickness], inflate = tolerance)
            setebc!(u,l1,true, 3, 0.0)
            
            applyebc!(u)
            numberdofs!(u)
            el1femm =  FEMMBase(IntegData(subset(bdryfes,icl), TriRule(3)))
            function pfun(forceout::FVec{T}, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt) where {T}
                local r = sqrt(XYZ[1]^2 + XYZ[2]^2)
                nx = XYZ[1]/r; ny = XYZ[2]/r
                # local sx, sy, txy
                # sx, sy, txy = sigmaxx(XYZ), sigmayy(XYZ), sigmaxy(XYZ)
                # sn = sx * nx^2 + sy * ny^2 + 2 * nx * ny * txy
                # tn = -(sx - sy) * nx * ny + (nx^2 - ny^2) * txy
                # forceout[1] = sn * nx - tn * ny
                # forceout[2] = sn * ny + tn * nx
                # forceout[3] = 0.0
                forceout[1] = sigmarr(XYZ) * nx - sigmart(XYZ) * ny
                forceout[2] = sigmarr(XYZ) * ny + sigmart(XYZ) * nx
                forceout[3] = 0.0
                return forceout
            end
            fi = ForceIntensity(FFlt, 3, pfun);
            F2 = distribloads(el1femm, geom, u, fi, 2);
            
            MR = DeforModelRed3D
            
            material = MatDeforElastIso(MR, E, nu)
            
            femm = FEMMDeforLinear(MR, IntegData(fes, TetRule(4)), material)
            
            # The geometry field now needs to be associated with the FEMM
            femm = associategeometry!(femm, geom)
            
            K = stiffness(femm, geom, u)
            K = cholesky(K)
            U = K\(F2)
            scattersysvec!(u, U[:])
            
            nlA = selectnode(fens, box=[Ri, Ri, 0.0, 0.0, 0.0, Thickness], inflate=tolerance);
            nlB = selectnode(fens, box=[0.0, 0.0, Ri, Ri, 0.0, Thickness], inflate=tolerance);
            # thecorneru = zeros(FFlt,length(nlA),3)
            # gathervalues_asmat!(u, thecorneru, nl);
            # thecorneru = mean(thecorneru, 1)[1]/phun("mm")
            # println("displacement = $(thecorneru) vs -0.10215 [MM]")
            
            println("Extrapolation: $( extrapolation )")
            sigx = fieldfromintegpoints(femm, geom, u, :Cauchy, 1; nodevalmethod = :averaging, reportat = extrapolation)
            sigy = fieldfromintegpoints(femm, geom, u, :Cauchy, 2; nodevalmethod = :averaging, reportat = extrapolation)
            sigyA = mean(sigy.values[nlA,1], 1)[1]
            sigyAtrue = sigmatt([Ri, 0.0, 0.0])
            println("sig_y@A =$(sigyA/phun("MPa")) vs $(sigyAtrue/phun("MPa")) [MPa]")
            sigxB = mean(sigx.values[nlB,1], 1)[1]
            sigxBtrue = sigmatt([0.0, Ri, 0.0])
            println("sig_x@B =$(sigxB/phun("MPa")) vs $(sigxBtrue/phun("MPa")) [MPa]")
            push!(numnodes, count(fens))
            push!(numelements, count(fes))
            push!(sigxderrs[extrapolation], abs(sigxB/sigxBtrue - 1.0))
            push!(sigyderrs[extrapolation], abs(sigyA/sigyAtrue - 1.0))
            # File =  "a.vtk"
            # vtkexportmesh(File, fes.conn, geom.values,
            # FinEtools.MeshExportModule.H8; vectors=[("u", u.values)],
            # scalars=[("sigmax", sigx.values/phun("MEGA*PA"))])
            # @async run(`"paraview.exe" $File`)
        end
    end
    
    File = "plate_w_hole_PE_T10_convergence.CSV"
    savecsv(File, numelements=vec(numelements), numnodes=vec(numnodes), sigxderrtrend=vec(sigxderrs[:extraptrend]), sigxderrdefault=vec(sigxderrs[:extrapmean]), sigyderrtrend=vec(sigyderrs[:extraptrend]), sigyderrdefault=vec(sigyderrs[:extrapmean]))
    
end # plate_w_hole_T10_PE_convergence


function plate_w_hole_T10_stress()
    E = 2.4*phun("MEGA*PA");# 210e3 MPa
    nu = 0.3;
    Re = 0.3*phun("M"); # outer radius
    Ri= 0.1*phun("M"); # hole radius
    H = 0.1*phun("M") # thickness of the plate
    nRadial, nCircumferential, nThickness = 6, 8, 1;
    sigma0=1*phun("MEGA*PA");
    
    function sigmaxx(x)
        local r = norm(x[1:2]);
        local th = atan2(x[2],x[1]);
        return sigma0*(1-Ri^2/r^2*(3/2*cos(2*th)+cos(4*th))+3/2*Ri^4/r^4*cos(4*th));
    end
    function sigmayy(x)
        local r=norm(x[1:2]);
        local th = atan2(x[2],x[1]);
        return -sigma0*(Ri^2/r^2*(1/2*cos(2*th)-cos(4*th))+3/2*Ri^4/r^4*cos(4*th));
    end
    function sigmaxy(x)
        local r=norm(x[1:2]);
        local th = atan2(x[2],x[1]);
        return -sigma0*(Ri^2/r^2*(1/2*sin(2*th)+sin(4*th))-3/2*Ri^4/r^4*sin(4*th));
    end
    function sigmarr(x)
        local r = norm(x[1:2]);
        local th = atan2(x[2],x[1]);
        return sigma0/2*(1-Ri^2/r^2) + sigma0/2*(1-4*Ri^2/r^2+3*Ri^4/r^4)*cos(2*th)
    end
    function sigmatt(x)
        local r = norm(x[1:2]);
        local th = atan2(x[2],x[1]);
        return sigma0/2*(1+Ri^2/r^2) - sigma0/2*(1+3*Ri^4/r^4)*cos(2*th)
    end
    function sigmart(x)
        local r = norm(x[1:2]);
        local th = atan2(x[2],x[1]);
        return -sigma0/2*(1+2*Ri^2/r^2-3*Ri^4/r^4)*sin(2*th)
    end
    
    convergencestudy = FDataDict[]
    for ref in 0:1:3
        println("ref = $(ref)")
        # Thickness = H
        Thickness = H/2^ref
        tolerance = Thickness/2^ref/1000.; # Geometrical tolerance
        
        fens,fes = T10block(1.0, pi/2, Thickness, 2^ref*nRadial, 2^ref*nCircumferential, 2^ref*nThickness)
        
        bdryfes = meshboundary(fes);
        icl = selectelem(fens, bdryfes, box=[1.0, 1.0, 0.0, pi/2, 0.0, Thickness], inflate=tolerance);
        
        for i=1:count(fens)
            t=fens.xyz[i,1]; a=fens.xyz[i,2]; z=fens.xyz[i,3]
            fens.xyz[i,:] = [(t*Re+(1-t)*Ri)*cos(a), (t*Re+(1-t)*Ri)*sin(a), z];
        end
        
        geom = NodalField(fens.xyz)
        u = NodalField(zeros(size(fens.xyz,1),3)) # displacement field
        
        l1 =selectnode(fens; box=[0.0, Inf, 0.0, 0.0, 0.0, Thickness], inflate = tolerance)
        setebc!(u,l1,true, 2, 0.0)
        l1 =selectnode(fens; box=[0.0, 0.0, 0.0, Inf, 0.0, Thickness], inflate = tolerance)
        setebc!(u,l1,true, 1, 0.0)
        # Plane-stress constraint: assume the plane z=0 is the plane of symmetry of the plate
        l1 =selectnode(fens; box=[0.0, Inf, 0.0, Inf, 0.0, 0.0], inflate = tolerance)
        setebc!(u,l1,true, 3, 0.0)
        # If this was enabled, the plane-strain  constraint would be enforced.
        # l1 =selectnode(fens; box=[0.0, Inf, 0.0, Inf, Thickness, Thickness], inflate = tolerance)
        # setebc!(u,l1,true, 3, 0.0)
        
        applyebc!(u)
        numberdofs!(u)
        el1femm =  FEMMBase(IntegData(subset(bdryfes,icl), TriRule(3)))
        function pfun(forceout::FVec{T}, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt) where {T}
            local r = sqrt(XYZ[1]^2 + XYZ[2]^2)
            nx = XYZ[1]/r; ny = XYZ[2]/r
            forceout[1] = sigmarr(XYZ) * nx - sigmart(XYZ) * ny
            forceout[2] = sigmarr(XYZ) * ny + sigmart(XYZ) * nx
            forceout[3] = 0.0
            return forceout
        end
        fi = ForceIntensity(FFlt, 3, pfun);
        F2 = distribloads(el1femm, geom, u, fi, 2);
        
        MR = DeforModelRed3D
        
        material = MatDeforElastIso(MR, E, nu)
        
        femm = FEMMDeforLinear(MR, IntegData(fes, TetRule(4)), material)
        
        # The geometry field now needs to be associated with the FEMM
        femm = associategeometry!(femm, geom)
        
        K = stiffness(femm, geom, u)
        K = cholesky(K)
        U = K\(F2)
        scattersysvec!(u, U[:])
        
        stressfields = elemfieldfromintegpoints(femm, geom, u, :Cauchy, collect(1:6))
        
        push!(convergencestudy, FDataDict(
        "elementsize"=> 1.0 / 2^ref,
        "fens"=>fens,
        "fes"=>fes,
        "geom"=>geom,
        "u"=>u,
        "femm"=>femm,
        "stressfields"=>[stressfields],
        "tolerance"=>tolerance)
        )
    end # for ref in 0:1:5
    
    # File = "mplate_w_hole_T10m_stress"
    # open(File * ".jls", "w") do file
    #     serialize(file, convergencestudy)
    # end
    
end # plate_w_hole_T10_stress

function allrun()
    println("#####################################################") 
    println("# plate_w_hole_T10_PE_convergence ")
    plate_w_hole_T10_PE_convergence()
    println("#####################################################") 
    println("# plate_w_hole_T10_stress ")
    plate_w_hole_T10_stress()
    println("#####################################################") 
    println("# plate_w_hole_RECT_MSH8_convergence ")
    plate_w_hole_RECT_MSH8_convergence()
    println("#####################################################") 
    println("# plate_w_hole_H20_stress ")
    plate_w_hole_H20_stress()
    println("#####################################################") 
    println("# plate_w_hole_MSH8_convergence ")
    plate_w_hole_MSH8_convergence()
    println("#####################################################") 
    println("# plate_w_hole_MSH8_PE_convergence ")
    plate_w_hole_MSH8_PE_convergence()
    println("#####################################################") 
    println("# plate_w_hole_MST10_convergence ")
    plate_w_hole_MST10_convergence()
    println("#####################################################") 
    println("# plate_w_hole_MST10_PE_convergence ")
    plate_w_hole_MST10_PE_convergence()
    println("#####################################################") 
    println("# plate_w_hole_MST10_stress ")
    plate_w_hole_MST10_stress()
    println("#####################################################") 
    println("# plate_w_hole_RECT_H20_convergence ")
    plate_w_hole_RECT_H20_convergence()
    return true
end # function allrun

end # module plate_w_hole_examples
