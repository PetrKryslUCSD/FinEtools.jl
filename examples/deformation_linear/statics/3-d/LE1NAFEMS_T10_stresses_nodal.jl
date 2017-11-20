module mLE1NAFEMS_T10m
using FinEtools
using FinEtools.MeshExportModule
using Main.ComputeErrorsModule
function test()
    
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
        K = cholfact(K)
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
    
    File = "LE1NAFEMS_T10_stresses_nodal_convergence_$(elementtag)"
    open(File * ".jls", "w") do file
        serialize(file, convergencestudy)
    end
    
    ComputeErrorsModule.process(File)

end
end
using .mLE1NAFEMS_T10m
mLE1NAFEMS_T10m.test()
