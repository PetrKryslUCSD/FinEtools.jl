module mplate_w_hole_T10_stress
using FinEtools
using FinEtools.MeshExportModule
using ComputeErrorsModule
function test()
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
        el1femm =  FEMMBase(GeoD(subset(bdryfes,icl), TriRule(3)))
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

        femm = FEMMDeforLinear(MR, GeoD(fes, TetRule(4)), material)

        # The geometry field now needs to be associated with the FEMM
        femm = associategeometry!(femm, geom)

        K = stiffness(femm, geom, u)
        K = cholfact(K)
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

    File = "mplate_w_hole_T10m_stress"
    open(File * ".jls", "w") do file
        serialize(file, convergencestudy)
    end

    ComputeErrorsModule.process(File)
end

end # module

using mplate_w_hole_T10_stress
mplate_w_hole_T10_stress.test()
