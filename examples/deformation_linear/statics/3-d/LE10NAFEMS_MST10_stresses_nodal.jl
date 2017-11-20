module mLE1NAFEMS_MST10m
using FinEtools
using FinEtools.MeshExportModule
using Main.ComputeErrorsModule
function test()

    elementtag = "MST10"
    println("LE10NAFEMS, 3D version. Element: $(elementtag)")

    E = 210e3*phun("MEGA*PA");# 210e3 MPa
    nu = 0.3;
    qmagn = 1.0*phun("MEGA*PA");# transverse pressure
    sigma_yP = -5.38*phun("MEGA*PA");# tensile stress at [2.0, 0.0] meters
    Ae =3.25*phun("m"); # Major radius of the exterior ellipse
    Be =2.75*phun("m"); # Minor radius of the exterior ellipse
    Ai =2.0*phun("m"); # Major radius of the interior ellipse
    Bi =1.0*phun("m"); # Minor radius of the interior ellipse
    Thickness = 0.6*phun("m")

    for extrapolation in [:extraptrend :extrapmean]
        convergencestudy = FDataDict[]
        for ref in  [0, 1, 2, 3]
            tolerance = Thickness/2^ref/1000.; # Geometrical tolerance

            nr, nc, nt = 2^ref*5, 2^ref*6, 2^ref*2
            @assert nt % 2 == 0 "Number of elements through the thickness must be even"
            fens,fes = T10block(1.0, pi/2, Thickness, nr, nc, nt, orientation=:b)
            
            
            # Select the  boundary faces, on the boundary that is clamped,  and on the part
            # of the boundary that is loaded with the transverse pressure
            bdryfes = meshboundary(fes);
            exteriorbfl = selectelem(fens, bdryfes, box=[1.0, 1.0, 0.0, pi/2, 0.0, Thickness], inflate=tolerance);
            topbfl = selectelem(fens, bdryfes, box=[0.0, 1.0, 0.0, pi/2, Thickness, Thickness], inflate=tolerance);
            
            # Reshape the generated block into the elliptical plate
            for i=1:count(fens)
                r=fens.xyz[i,1]; a=fens.xyz[i,2]; z=fens.xyz[i,3]
                fens.xyz[i,:]=[(r*Ae+(1-r)*Ai)*cos(a) (r*Be+(1-r)*Bi)*sin(a) z];
            end
            
            
            geom = NodalField(fens.xyz)
            u = NodalField(zeros(size(fens.xyz,1),3)) # displacement field
            
            l12 =connectednodes(subset(bdryfes, exteriorbfl)) # external boundary
            setebc!(u, l12, true, 1, 0.0)
            setebc!(u, l12, true, 2, 0.0)
            ll = selectnode(fens; box=[0.0, Inf, 0.0, Inf, Thickness/2.0, Thickness/2.0], inflate = tolerance)
            l3 = intersect(ll, connectednodes(subset(bdryfes, exteriorbfl)))
            setebc!(u, l3, true, 3, 0.0)
            l1 =selectnode(fens; box=[0.0, 0.0, 0.0, Inf, 0.0, Thickness], inflate = tolerance)
            setebc!(u,l1,true, 1, 0.0) # symmetry plane X = 0
            l2 =selectnode(fens; box=[0.0, Inf, 0.0, 0.0, 0.0, Thickness], inflate = tolerance)
            setebc!(u,l2,true, 2, 0.0) # symmetry plane Y = 0
            
            applyebc!(u)
            numberdofs!(u)
            
            el1femm =  FEMMBase(IntegData(subset(bdryfes,topbfl), TriRule(3)))
            function pfun(forceout::FVec{T}, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt) where {T}
                forceout .=  [0.0, 0.0, -qmagn]
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
            K = cholfact(K)
            U = K\(F2)
            scattersysvec!(u, U[:])

            nl = selectnode(fens, box=[Ai,Ai,0,0,Thickness,Thickness],inflate=tolerance);
            thecorneru = zeros(FFlt,1,3)
            gathervalues_asmat!(u, thecorneru, nl);
            thecorneru = thecorneru/phun("mm")
            println("displacement =$(thecorneru) [MM] as compared to reference [-0.030939, 0, -0.10488] [MM]")

            stressfield = fieldfromintegpoints(femm, geom, u, :Cauchy, collect(1:6);
                nodevalmethod = :averaging, reportat = extrapolation)
        
            # File =  "LE10NAFEMS_MST10_sigmay.vtk"
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
        
        File = "LE1NAFEMS_MST10_stresses_nodal_convergence_$(elementtag)_$(extrapolation)"
        open(File * ".jls", "w") do file
            serialize(file, convergencestudy)
        end
        
        ComputeErrorsModule.process(File)
    end # for extrapolation

end
end
using .mLE1NAFEMS_MST10m
mLE1NAFEMS_MST10m.test()
