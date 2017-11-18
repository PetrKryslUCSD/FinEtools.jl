module mplate_w_hole_RECT_MSH8m
using FinEtools
using FinEtools.MeshExportModule
using DataFrames
using CSV
function test()
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
        for ref in [3]
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
            K = cholfact(K)
            U = K\(Fx + Fy)
            scattersysvec!(u, U[:])

            nlA = selectnode(fens, box=[Ri, Ri, 0.0, 0.0, 0.0, Thickness], inflate=tolerance);
            nlB = selectnode(fens, box=[0.0, 0.0, Ri, Ri, 0.0, Thickness], inflate=tolerance);
            # thecorneru = zeros(FFlt,length(nlA),3)
            # gathervalues_asmat!(u, thecorneru, nl);
            # thecorneru = mean(thecorneru, 1)[1]/phun("mm")
            # println("displacement = $(thecorneru) vs -0.10215 [MM]")

            println("Extrapolation: $( extrapolation )")
            sigx = fieldfromintegpoints(femm, geom, u, :Cauchy, 1;
                reportat = extrapolation)
            sigy = fieldfromintegpoints(femm, geom, u, :Cauchy, 2;
                reportat = extrapolation)
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

end
end
using mplate_w_hole_RECT_MSH8m
mplate_w_hole_RECT_MSH8m.test()
