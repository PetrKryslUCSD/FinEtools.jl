module mplate_w_hole_RECT_H20m
using FinEtools
using FinEtools.MeshExportModule
using FinEtools.MeshImportModule: import_ABAQUS
# using DataFrames
# using CSV
function test()
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
    for extrapolation in [:extrapmean]
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
            elxfemm =  FEMMBase(GeoD(subset(bdryfes,ixl), GaussRule(2, 2)))
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
            elyfemm =  FEMMBase(GeoD(subset(bdryfes,iyl), GaussRule(2, 2)))
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

            femm = FEMMDeforLinear(MR, GeoD(fes, GaussRule(3, 2)), material)

            # The geometry field now needs to be associated with the FEMM
            femm = associategeometry!(femm, geom)

            K = stiffness(femm, geom, u)
            K = cholfact(K)
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
            sigx = fieldfromintegpoints(femm, geom, u, :Cauchy, 1;
                reportat = extrapolation)
            sigy = fieldfromintegpoints(femm, geom, u, :Cauchy, 2;
                reportat = extrapolation)
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

end
end
using mplate_w_hole_RECT_H20m
mplate_w_hole_RECT_H20m.test()
