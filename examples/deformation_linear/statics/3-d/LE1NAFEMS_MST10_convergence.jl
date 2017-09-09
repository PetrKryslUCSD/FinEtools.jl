module mLE1NAFEMS_MST10m
using FinEtools
using FinEtools.MeshExportModule
using DataFrames
using CSV
using Base.Test
function test()

    println("LE1NAFEMS, 3D version. MST10"        )
    t0 = time()

    E = 210e3*phun("MEGA*PA");# 210e3 MPa
    nu = 0.3;
    p = 10*phun("MEGA*PA");# 10 MPA Outward pressure on the outside ellipse
    sigma_yD = 92.7*phun("MEGA*PA");# tensile stress at [2.0, 0.0] meters
    Thick0 = 0.1*phun("m")/2.0 # to account for the symmetry reduction

    sigyderrs = Dict{Symbol, FFltVec}()

    nelems = []
    for extrapolation in [:extraptrendpaper :extraptrend :extrapmean]
        sigyderrs[extrapolation] = FFltVec[]
        nelems = []
        for ref in 0:1:4
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


            el1femm =  FEMMBase(GeoD(subset(bdryfes,icl), TriRule(3)))
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

            femm = FEMMDeforLinearMST10(MR, GeoD(fes, TetRule(4)), material)

            # The geometry field now needs to be associated with the FEMM
            femm = associategeometry!(femm, geom)

            K = stiffness(femm, geom, u)
            K = cholfact(K)
            U = K\(F2)
            scattersysvec!(u, U[:])

            nl = selectnode(fens, box=[2.0, 2.0, 0.0, 0.0, 0.0, Thickness], inflate=tolerance);
            thecorneru = zeros(FFlt,length(nl),3)
            gathervalues_asmat!(u, thecorneru, nl);
            thecorneru = mean(thecorneru, 1)[1]/phun("mm")
            println("displacement =$(thecorneru) vs -0.10215 [MM]")

            fld = fieldfromintegpoints(femm, geom, u, :Cauchy, 2;
                tonode = extrapolation)
            sigyd = mean(fld.values[nl,1], 1)[1]
            println("Sigma_y =$(sigyd/phun("MPa")) vs $(sigma_yD/phun("MPa")) [MPa]")

            println("$extrapolation, $(count(fes)), $(sigyd/phun("MPa"))")
            push!(nelems, count(fes))
            push!(sigyderrs[extrapolation], abs(sigyd/sigma_yD - 1.0))
            # File =  "a.vtk"
            # vtkexportmesh(File, fes.conn, geom.values,
            # FinEtools.MeshExportModule.H8; vectors=[("u", u.values)],
            # scalars=[("sigmay", fld.values)])
            # @async run(`"paraview.exe" $File`)
        end
    end



    df = DataFrame(nelems=vec(nelems),
        sigyderrtrendpaper=vec(sigyderrs[:extraptrendpaper]),
        sigyderrtrend=vec(sigyderrs[:extraptrend]),
        sigyderrdefault=vec(sigyderrs[:extrapmean]))
    File = "LE1NAFEMS_MST10_convergence.CSV"
    CSV.write(File, df)
    @async run(`"paraview.exe" $File`)

end
end
using mLE1NAFEMS_MST10m
mLE1NAFEMS_MST10m.test()
