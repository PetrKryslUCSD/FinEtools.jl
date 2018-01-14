module THICK_plate_2dir_strong_examples
using FinEtools
using FinEtools.AlgoDeforLinearModule
using FinEtools.MeshUtilModule
using FinEtools.AlgoBaseModule

function THICK_plate_2dir_strong_MST10_conv()
    elementtag = "MST10"
    println("""
    Fiber-reinforced cantilever plate: $(elementtag)
    """)
    
    # This example provides three-dimensional finite element model for  the
    # transverse shear stress calculations. The problem consists of a one-, two- or
    # three-layer plate subjected to a sinusoidal distributed load, as
    # described by Pagano (1969). The resulting transverse shear and axial
    # stresses through the thickness of the plate are compared to two existing
    # analytical solutions by Pagano (1969). The first solution is derived from
    # classical laminated plate theory (CPT), while the second is an exact
    # solution from linear elasticity theory.
    
    
    for (extrap, nodevalmeth) = zip([:extrapmean, :extraptrend, :default], [:averaging, :averaging, :invdistance])
        filebase = "THICK_plate_2dir_strong_MST10_conv_$(elementtag)_$(extrap)"
        modeldatasequence = FDataDict[]
        for Refinement = [1, 2, 4, 8]
            
            # Orthotropic material 
            E1 = 1000.0e9*phun("Pa"); E2 = 1000.0e9*phun("Pa"); E3 = 1.0e9*phun("Pa");
            G12 = 0.2e9*phun("Pa");  G13 = G12; G23 = 0.2e9*phun("Pa")
            nu12 = nu13 = nu23 =  0.25;
            # dimensions of the plate
            a = 70.0*phun("mm"); 
            b = 100.0*phun("mm"); 
            t = 50.0*phun("mm"); 
            # Transverse loading
            q0 = 1000.0*phun("Pa")
            # Coefficients of thermal expansion
            CTE1 = CTE2 = CTE3 = 0.0 

            # Here we define the layout and the thicknesses of the layers.
            angles = vec([45.0]);
            nLayers = length(angles)
            ts = t/nLayers * ones(nLayers); # layer thicknesses
            
            tolerance = 0.0001*t
            
            # Select how fine the mesh should be
            nts= Refinement * ones(Int, nLayers);# number of elements per layer
            tnts = sum(nts)
            na, nb = 4 * tnts, 4 * tnts;
            
            xs = collect(linspace(0.0, a, na+1))
            ys = collect(linspace(0.0, b, nb+1))
            fens,fes = T10layeredplatex(xs, ys, ts, nts)
            println("Mesh: na, nb, nts = $na, $nb, $nts")
            println("count(fens) = $(count(fens))")
            
            # This is the material  model
            MR = DeforModelRed3D
            material = MatDeforElastOrtho(MR, 0.0, E1, E2, E3, nu12, nu13, nu23, G12, G13, G23, CTE1, CTE2, CTE3)
            
            # The material coordinate system function is defined as:
            function updatecs!(csmatout::FFltMat, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)
                rotmat3!(csmatout, angles[fe_label]/180.0*pi* [0.0; 0.0; 1.0]);
            end
            
            # The volume integrals are evaluated using this rule
            gr = SimplexRule(3, 4)
            
            # We will create 3 regions, one for each of the layers
            regions = FDataDict[]
            for layer = 1:nLayers
                rls = selectelem(fens, fes, label =  layer)
                push!(regions, FDataDict("femm"=>FEMMDeforLinearMST10(MR, IntegData(subset(fes, rls), gr), CSys(3, 3, updatecs!), material)))
            end
            
            # The essential boundary conditions: clamped face
            lx0 = selectnode(fens, box=[0.0 0.0 -Inf Inf -Inf Inf], inflate=tolerance)
            eclamped1 = FDataDict("displacement"=>  0.0, "component"=> 1, "node_list"=>lx0)
            eclamped2 = FDataDict("displacement"=>  0.0, "component"=> 2, "node_list"=>lx0)
            eclamped3 = FDataDict("displacement"=>  0.0, "component"=> 3, "node_list"=>lx0)
            
            # The traction boundary condition is applied at the free face opposite the clamped face.
            bfes = meshboundary(fes)
            function pfun(forceout::FVec{T}, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt) where {T}
                forceout[1] = 0.0
                forceout[2] = 0.0
                forceout[3] = -q0
                return forceout
            end
            # From  the entire boundary we select those quadrilaterals that lie on the plane
            # Z = thickness
            tl = selectelem(fens, bfes, box = [a a -Inf Inf -Inf Inf], inflate=tolerance)
            Trac = FDataDict("traction_vector"=>pfun, "femm"=>FEMMBase(IntegData(subset(bfes, tl), SimplexRule(2, 3))))
            
            modeldata = FDataDict("fens"=>fens, "regions"=>regions, "essential_bcs"=>[eclamped1, eclamped2, eclamped3],  "traction_bcs"=> [Trac])
            modeldata = AlgoDeforLinearModule.linearstatics(modeldata)
            
            # modeldata["postprocessing"] = FDataDict("file"=>filebase * "-u")
            # modeldata = AlgoDeforLinearModule.exportdeformation(modeldata)
            
            u = modeldata["u"]
            geom = modeldata["geom"]
            
            # The results of the displacement and stresses will be reported at
            # nodes located at the appropriate points.
            nbottomedge = selectnode(fens, box=[a a 0.0 0.0 0.0 0.0], inflate=tolerance)
            
            println("bottom edge deflection: $(mean(u.values[nbottomedge, 3], 1)/phun("mm")) [mm]")
            
            #  Compute  all stresses
            modeldata["postprocessing"] = FDataDict("file"=>filebase * "-s", "quantity"=>:Cauchy, "component"=>collect(1:6), "outputcsys"=>CSys(3), "nodevalmethod"=>nodevalmeth, "reportat"=>extrap)
            modeldata = AlgoDeforLinearModule.exportstress(modeldata)
            
            modeldata["elementsize"] = t/Refinement
            modeldata["geometricaltolerance"] = tolerance
            push!(modeldatasequence, modeldata)
        end # for refinement
        
        println("")
        println("Stress RMS error")
        for md = modeldatasequence
            md["targetfields"] = [e["field"] for e in md["postprocessing"]["exported"]]
        end
        elementsizes, errornorms, p = AlgoBaseModule.evalconvergencestudy(modeldatasequence)
        
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
        elementsizes, errornorms, p = AlgoBaseModule.evalconvergencestudy(modeldatasequence)
        
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
        
        # @async run(`"paraview.exe" $csvFile`)
    end # extrap
    println("Done")
    
end # Pagano_3lay_cyl_bend_MST10_conv

function allrun()
    println("#####################################################") 
    println("# THICK_plate_2dir_strong_MST10_conv ")
    THICK_plate_2dir_strong_MST10_conv()
    return true
end # function allrun

end # module Pagano_examples
