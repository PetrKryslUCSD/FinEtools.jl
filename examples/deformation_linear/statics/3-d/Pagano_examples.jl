module Pagano_examples
using FinEtools
using FinEtools.AlgoDeforLinearModule
using FinEtools.MeshUtilModule
using FinEtools.AlgoBaseModule

function Pagano_3lay_cyl_bend_MST10_conv()
    elementtag = "MST10"
    println("""
    Pagano 3layer cylindrical bending: $(elementtag)
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
        for Span_to_thickness = [4, 50] 
            
            filebase = "Pagano_3lay_cyl_bend_conv_$(elementtag)_$(Span_to_thickness)_$(extrap)"
            modeldatasequence = FDataDict[]
            for Refinement = [1, 2, 4, 8]
                
                # Orthotropic material for the 3 layers
                E1 = 25e6*phun("PSI"); E2 = 1e6*phun("PSI"); E3 = E2;
                G12 = 0.5e6*phun("PSI");  G13 = G12; G23 = 0.2e6*phun("PSI")
                nu12 =  0.25; nu13 =  0.25; nu23 =  0.25;
                T = 2.5*phun("in"); # total thickness of the plate
                L = Span_to_thickness*T;
                h = 1.0/Refinement*phun("in");  # depth of the plate
                q0 = 1*phun("PSI")
                CTE1 =  CTE2 =  CTE3 = 0.0
                
                # Here we define the layout and the thicknesses of the layers.
                angles = vec([0.0 90.0 0.0]);
                nLayers = length(angles)
                ts = T/nLayers * ones(nLayers); # layer thicknesses
                
                tolerance = 0.0001*T
                
                # Select how fine the mesh should be
                nL, nh = Refinement * 2 * 6, 1;
                nts= Refinement * 2 * ones(Int, nLayers);# number of elements per layer
                
                xs = collect(linearspace(0.0, L, nL+1))
                ys = collect(linearspace(0.0, h, nh+1))
                fens,fes = T10layeredplatex(xs, ys, ts, nts)
                display(typeof(fes))
                println("Mesh: nL, nh, nts = $nL, $nh, $nts")
                println("count(fens) = $(count(fens))")
                
                # This is the material  model
                MR = DeforModelRed3D
                skinmaterial = MatDeforElastOrtho(MR,
                0.0, E1, E2, E3,
                nu12, nu13, nu23,
                G12, G13, G23,
                CTE1, CTE2, CTE3)
                
                # The material coordinate system function is defined as:
                function updatecs!(csmatout::FFltMat, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)
                    rotmat3!(csmatout, angles[fe_label]/180.0*pi* [0.0; 0.0; 1.0]);
                end
                
                # The vvolume integrals are evaluated using this rule
                gr = SimplexRule(3, 4)

                # We will create 3 regions, one for each of the layers
                regions = FDataDict[]
                for layer = 1:nLayers
                    rls = selectelem(fens, fes, label =  layer)
                    push!(regions, FDataDict("femm"=>FEMMDeforLinearMST10(MR, IntegData(subset(fes, rls), gr), CSys(3, 3, updatecs!), skinmaterial)))
                end
                
                # The essential boundary conditions are applied to enforce the plane strain constraint.
                ly0 = selectnode(fens, box=[-Inf Inf 0.0 0.0 -Inf Inf], inflate=tolerance)
                lyh = selectnode(fens, box=[-Inf Inf h h -Inf Inf], inflate=tolerance)
                ey = FDataDict("displacement"=>  0.0, "component"=> 2, "node_list"=>vcat(ly0, lyh))
                # The transverse displacement is fixed at the two ends.
                lz0 = selectnode(fens, box=[0.0 0.0 -Inf Inf -Inf Inf], inflate=tolerance)
                lzL = selectnode(fens, box=[L L -Inf Inf -Inf Inf], inflate=tolerance)
                ez = FDataDict("displacement"=>  0.0, "component"=> 3, "node_list"=>vcat(lz0, lzL))
                ex = FDataDict("displacement"=>  0.0, "component"=> 1, "node_list"=>[1])
                
                # The traction boundary condition is applied at the top of the plate.
                bfes = meshboundary(fes)
                function pfun(forceout::FVec{T}, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt) where {T}
                    forceout[1] = 0.0
                    forceout[2] = 0.0
                    forceout[3] = -q0*sin(pi*XYZ[1]/L)
                    return forceout
                end
                # From  the entire boundary we select those quadrilaterals that lie on the plane
                # Z = thickness
                tl = selectelem(fens, bfes, box = [-Inf Inf -Inf Inf T T], inflate=tolerance)
                Trac = FDataDict("traction_vector"=>pfun,
                "femm"=>FEMMBase(IntegData(subset(bfes, tl), SimplexRule(2, 3))))
                
                modeldata = FDataDict("fens"=>fens,
                "regions"=>regions,
                "essential_bcs"=>[ex, ey, ez],
                "traction_bcs"=> [Trac]
                )
                modeldata = AlgoDeforLinearModule.linearstatics(modeldata)
                
                # modeldata["postprocessing"] = FDataDict("file"=>filebase * "-u")
                # modeldata = AlgoDeforLinearModule.exportdeformation(modeldata)
                
                u = modeldata["u"]
                geom = modeldata["geom"]
                
                # The results of the displacement and stresses will be reported at
                # nodes located at the appropriate points.
                ntopcenter = selectnode(fens, box=[L/2 L/2 0.0 h T T], inflate=tolerance)
                ncenterline = selectnode(fens, box=[L/2 L/2 0.0 0.0 0.0 T], inflate=tolerance)
                nx0line = selectnode(fens, box=[0.0 0.0 0.0 0.0 0.0 T], inflate=tolerance)
                
                zclo = sortperm(vec(geom.values[ncenterline, 3]))
                ncenterline = ncenterline[zclo]
                centerz = geom.values[ncenterline, 3]
                
                println("Top Center deflection: $(mean(u.values[ntopcenter, 3], 1)/phun("in")) [in]")
                
                #  Compute  all stresses
                modeldata["postprocessing"] = FDataDict("file"=>filebase * "-s", "quantity"=>:Cauchy, "component"=>collect(1:6), "outputcsys"=>CSys(3), "nodevalmethod"=>nodevalmeth, "reportat"=>extrap)
                modeldata = AlgoDeforLinearModule.exportstress(modeldata)
                
                modeldata["elementsize"] = T/Refinement
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
        end # Span_to_thickness 
    end # extrap
    println("Done")
    
end # Pagano_3lay_cyl_bend_MST10_conv


function Pagano_3lay_cyl_bend_H8_conv()
    elementtag = "H8"
    println("""
    Pagano 3layer cylindrical bending: $(elementtag)
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
        for Span_to_thickness = [4, 50] 
            
            filebase = "Pagano_3lay_cyl_bend_conv_$(elementtag)_$(Span_to_thickness)_$(extrap)"
            modeldatasequence = FDataDict[]
            for Refinement = [1, 2, 4, 8]
                
                # Orthotropic material for the 3 layers
                E1 = 25e6*phun("PSI"); E2 = 1e6*phun("PSI"); E3 = E2;
                G12 = 0.5e6*phun("PSI");  G13 = G12; G23 = 0.2e6*phun("PSI")
                nu12 =  0.25; nu13 =  0.25; nu23 =  0.25;
                T = 2.5*phun("in"); # total thickness of the plate
                L = Span_to_thickness*T;
                h = 1.0/Refinement*phun("in");  # depth of the plate
                q0 = 1*phun("PSI")
                CTE1 =  CTE2 =  CTE3 = 0.0
                
                # Here we define the layout and the thicknesses of the layers.
                angles = vec([0.0 90.0 0.0]);
                nLayers = length(angles)
                ts = T/nLayers * ones(nLayers); # layer thicknesses
                
                tolerance = 0.0001*T
                
                # Select how fine the mesh should be
                nL, nh = Refinement * 2 * 6, 1;
                nts= Refinement * 2 * ones(Int, nLayers);# number of elements per layer
                
                xs = collect(linearspace(0.0, L, nL+1))
                ys = collect(linearspace(0.0, h, nh+1))
                fens,fes = H8layeredplatex(xs, ys, ts, nts)
                
                println("Mesh: nL, nh, nts = $nL, $nh, $nts")
                println("count(fens) = $(count(fens))")
                
                # This is the material  model
                MR = DeforModelRed3D
                skinmaterial = MatDeforElastOrtho(MR,
                0.0, E1, E2, E3,
                nu12, nu13, nu23,
                G12, G13, G23,
                CTE1, CTE2, CTE3)
                
                # The material coordinate system function is defined as:
                function updatecs!(csmatout::FFltMat, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)
                    rotmat3!(csmatout, angles[fe_label]/180.0*pi* [0.0; 0.0; 1.0]);
                end
                
                # The vvolume integrals are evaluated using this rule
                gr = GaussRule(3, 2)
              
                # We will create 3 regions, one for each of the layers
                regions = FDataDict[]
                for layer = 1:nLayers
                    rls = selectelem(fens, fes, label =  layer)
                    push!(regions, FDataDict("femm"=>FEMMDeforLinear(MR, IntegData(subset(fes, rls), gr), CSys(3, 3, updatecs!), skinmaterial)))
                end
                
                # The essential boundary conditions are applied to enforce the plane strain constraint.
                ly0 = selectnode(fens, box=[-Inf Inf 0.0 0.0 -Inf Inf], inflate=tolerance)
                lyh = selectnode(fens, box=[-Inf Inf h h -Inf Inf], inflate=tolerance)
                ey = FDataDict("displacement"=>  0.0, "component"=> 2, "node_list"=>vcat(ly0, lyh))
                # The transverse displacement is fixed at the two ends.
                lz0 = selectnode(fens, box=[0.0 0.0 -Inf Inf -Inf Inf], inflate=tolerance)
                lzL = selectnode(fens, box=[L L -Inf Inf -Inf Inf], inflate=tolerance)
                ez = FDataDict("displacement"=>  0.0, "component"=> 3, "node_list"=>vcat(lz0, lzL))
                ex = FDataDict("displacement"=>  0.0, "component"=> 1, "node_list"=>[1])
                
                # The traction boundary condition is applied at the top of the plate.
                bfes = meshboundary(fes)
                function pfun(forceout::FVec{T}, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt) where {T}
                    forceout[1] = 0.0
                    forceout[2] = 0.0
                    forceout[3] = -q0*sin(pi*XYZ[1]/L)
                    return forceout
                end
                # From  the entire boundary we select those quadrilaterals that lie on the plane
                # Z = thickness
                tl = selectelem(fens, bfes, box = [-Inf Inf -Inf Inf T T], inflate=tolerance)
                Trac = FDataDict("traction_vector"=>pfun,                "femm"=>FEMMBase(IntegData(subset(bfes, tl), GaussRule(2, 2))))
                
                modeldata = FDataDict("fens"=>fens, "regions"=>regions, "essential_bcs"=>[ex, ey, ez], "traction_bcs"=> [Trac])
                modeldata = AlgoDeforLinearModule.linearstatics(modeldata)
                
                # modeldata["postprocessing"] = FDataDict("file"=>filebase * "-u")
                # modeldata = AlgoDeforLinearModule.exportdeformation(modeldata)
                
                u = modeldata["u"]
                geom = modeldata["geom"]
                
                # The results of the displacement and stresses will be reported at
                # nodes located at the appropriate points.
                ntopcenter = selectnode(fens, box=[L/2 L/2 0.0 h T T], inflate=tolerance)
                ncenterline = selectnode(fens, box=[L/2 L/2 0.0 0.0 0.0 T], inflate=tolerance)
                nx0line = selectnode(fens, box=[0.0 0.0 0.0 0.0 0.0 T], inflate=tolerance)
                
                zclo = sortperm(vec(geom.values[ncenterline, 3]))
                ncenterline = ncenterline[zclo]
                centerz = geom.values[ncenterline, 3]
                
                println("Top Center deflection: $(mean(u.values[ntopcenter, 3], 1)/phun("in")) [in]")
                
                #  Compute  all stresses
                modeldata["postprocessing"] = FDataDict("file"=>filebase * "-s",
                "quantity"=>:Cauchy, "component"=>collect(1:6), "outputcsys"=>CSys(3),
                "nodevalmethod"=>nodevalmeth, "reportat"=>extrap)
                modeldata = AlgoDeforLinearModule.exportstress(modeldata)
                
                modeldata["elementsize"] = T/Refinement
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
        end # Span_to_thickness 
    end # extrap
    println("Done")
    
end # Pagano_3lay_cyl_bend_H8_conv


function Pagano_3lay_cyl_bend_MSH8_conv()
    elementtag = "MSH8"
    println("""
    Pagano 3layer cylindrical bending: $(elementtag)
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
        for Span_to_thickness = [4, 50] 
            
            filebase = "Pagano_3lay_cyl_bend_conv_$(elementtag)_$(Span_to_thickness)_$(extrap)"
            modeldatasequence = FDataDict[]
            for Refinement = [1, 2, 4, 8]
                
                # Orthotropic material for the 3 layers
                E1 = 25e6*phun("PSI"); E2 = 1e6*phun("PSI"); E3 = E2;
                G12 = 0.5e6*phun("PSI");  G13 = G12; G23 = 0.2e6*phun("PSI")
                nu12 =  0.25; nu13 =  0.25; nu23 =  0.25;
                T = 2.5*phun("in"); # total thickness of the plate
                L = Span_to_thickness*T;
                h = 1.0/Refinement*phun("in");  # depth of the plate
                q0 = 1*phun("PSI")
                CTE1 =  CTE2 =  CTE3 = 0.0
                
                # Here we define the layout and the thicknesses of the layers.
                angles = vec([0.0 90.0 0.0]);
                nLayers = length(angles)
                ts = T/nLayers * ones(nLayers); # layer thicknesses
                
                tolerance = 0.0001*T
                
                # Select how fine the mesh should be
                nL, nh = Refinement * 2 * 6, 1;
                nts= Refinement * 2 * ones(Int, nLayers);# number of elements per layer
                
                xs = collect(linearspace(0.0, L, nL+1))
                ys = collect(linearspace(0.0, h, nh+1))
                fens,fes = H8layeredplatex(xs, ys, ts, nts)
                
                println("Mesh: nL, nh, nts = $nL, $nh, $nts")
                println("count(fens) = $(count(fens))")
                
                # This is the material  model
                MR = DeforModelRed3D
                skinmaterial = MatDeforElastOrtho(MR,
                0.0, E1, E2, E3,
                nu12, nu13, nu23,
                G12, G13, G23,
                CTE1, CTE2, CTE3)
                
                # The material coordinate system function is defined as:
                function updatecs!(csmatout::FFltMat, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)
                    rotmat3!(csmatout, angles[fe_label]/180.0*pi* [0.0; 0.0; 1.0]);
                end
                
                # The vvolume integrals are evaluated using this rule
                gr = GaussRule(3, 2)
                
                # We will create 3 regions, one for each of the layers
                regions = FDataDict[]
                for layer = 1:nLayers
                    rls = selectelem(fens, fes, label =  layer)
                    push!(regions, FDataDict("femm"=>FEMMDeforLinearMSH8(MR, IntegData(subset(fes, rls), gr), CSys(3, 3, updatecs!), skinmaterial)))
                end
                
                # The essential boundary conditions are applied to enforce the plane strain constraint.
                ly0 = selectnode(fens, box=[-Inf Inf 0.0 0.0 -Inf Inf], inflate=tolerance)
                lyh = selectnode(fens, box=[-Inf Inf h h -Inf Inf], inflate=tolerance)
                ey = FDataDict("displacement"=>  0.0, "component"=> 2, "node_list"=>vcat(ly0, lyh))
                # The transverse displacement is fixed at the two ends.
                lz0 = selectnode(fens, box=[0.0 0.0 -Inf Inf -Inf Inf], inflate=tolerance)
                lzL = selectnode(fens, box=[L L -Inf Inf -Inf Inf], inflate=tolerance)
                ez = FDataDict("displacement"=>  0.0, "component"=> 3, "node_list"=>vcat(lz0, lzL))
                ex = FDataDict("displacement"=>  0.0, "component"=> 1, "node_list"=>[1])
                
                # The traction boundary condition is applied at the top of the plate.
                bfes = meshboundary(fes)
                function pfun(forceout::FVec{T}, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt) where {T}
                    forceout[1] = 0.0
                    forceout[2] = 0.0
                    forceout[3] = -q0*sin(pi*XYZ[1]/L)
                    return forceout
                end
                # From  the entire boundary we select those quadrilaterals that lie on the plane
                # Z = thickness
                tl = selectelem(fens, bfes, box = [-Inf Inf -Inf Inf T T], inflate=tolerance)
                Trac = FDataDict("traction_vector"=>pfun, "femm"=>FEMMBase(IntegData(subset(bfes, tl), GaussRule(2, 2))))
                
                modeldata = FDataDict("fens"=>fens, "regions"=>regions, "essential_bcs"=>[ex, ey, ez], "traction_bcs"=> [Trac])
                modeldata = AlgoDeforLinearModule.linearstatics(modeldata)
                
                # modeldata["postprocessing"] = FDataDict("file"=>filebase * "-u")
                # modeldata = AlgoDeforLinearModule.exportdeformation(modeldata)
                
                u = modeldata["u"]
                geom = modeldata["geom"]
                
                # The results of the displacement and stresses will be reported at
                # nodes located at the appropriate points.
                ntopcenter = selectnode(fens, box=[L/2 L/2 0.0 h T T], inflate=tolerance)
                ncenterline = selectnode(fens, box=[L/2 L/2 0.0 0.0 0.0 T], inflate=tolerance)
                nx0line = selectnode(fens, box=[0.0 0.0 0.0 0.0 0.0 T], inflate=tolerance)
                
                zclo = sortperm(vec(geom.values[ncenterline, 3]))
                ncenterline = ncenterline[zclo]
                centerz = geom.values[ncenterline, 3]
                
                println("Top Center deflection: $(mean(u.values[ntopcenter, 3], 1)/phun("in")) [in]")
                
                #  Compute  all stresses
                modeldata["postprocessing"] = FDataDict("file"=>filebase * "-s", "quantity"=>:Cauchy, "component"=>collect(1:6), "outputcsys"=>CSys(3), "nodevalmethod"=>nodevalmeth, "reportat"=>extrap)
                modeldata = AlgoDeforLinearModule.exportstress(modeldata)
                
                modeldata["elementsize"] = T/Refinement
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
        end # Span_to_thickness 
    end # extrap
    println("Done")
    
end # Pagano_3lay_cyl_bend_MSH8_conv


function Pagano_3lay_cyl_bend_T10_conv()
    elementtag = "T10"
    println("""
    Pagano 3layer cylindrical bending: $(elementtag)
    """)
    
    # This example provides three-dimensional finite element model for  the
    # transverse shear stress calculations. The problem consists of a one-, two- or
    # three-layer plate subjected to a sinusoidal distributed load, as
    # described by Pagano (1969). The resulting transverse shear and axial
    # stresses through the thickness of the plate are compared to two existing
    # analytical solutions by Pagano (1969). The first solution is derived from
    # classical laminated plate theory (CPT), while the second is an exact
    # solution from linear elasticity theory.
    
    
    for (extrap, nodevalmeth) = zip([:default], [:invdistance])
        for Span_to_thickness = [4, 50] 
            
            filebase = "Pagano_3lay_cyl_bend_conv_$(elementtag)_$(Span_to_thickness)_$(extrap)"
            modeldatasequence = FDataDict[]
            for Refinement = [1, 2, 4, 8]
                
                # Orthotropic material for the 3 layers
                E1 = 25e6*phun("PSI"); E2 = 1e6*phun("PSI"); E3 = E2;
                G12 = 0.5e6*phun("PSI");  G13 = G12; G23 = 0.2e6*phun("PSI")
                nu12 =  0.25; nu13 =  0.25; nu23 =  0.25;
                T = 2.5*phun("in"); # total thickness of the plate
                L = Span_to_thickness*T;
                h = 1.0/Refinement*phun("in");  # depth of the plate
                q0 = 1*phun("PSI")
                CTE1 =  CTE2 =  CTE3 = 0.0
                
                # Here we define the layout and the thicknesses of the layers.
                angles = vec([0.0 90.0 0.0]);
                nLayers = length(angles)
                ts = T/nLayers * ones(nLayers); # layer thicknesses
                
                tolerance = 0.0001*T
                
                # Select how fine the mesh should be
                nL, nh = Refinement * 2 * 6, 1;
                nts= Refinement * 2 * ones(Int, nLayers);# number of elements per layer
                
                xs = collect(linearspace(0.0, L, nL+1))
                ys = collect(linearspace(0.0, h, nh+1))
                fens,fes = T10layeredplatex(xs, ys, ts, nts)
                
                println("Mesh: nL, nh, nts = $nL, $nh, $nts")
                println("count(fens) = $(count(fens))")
                
                # This is the material  model
                MR = DeforModelRed3D
                skinmaterial = MatDeforElastOrtho(MR,
                0.0, E1, E2, E3,
                nu12, nu13, nu23,
                G12, G13, G23,
                CTE1, CTE2, CTE3)
                
                # The material coordinate system function is defined as:
                function updatecs!(csmatout::FFltMat, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)
                    rotmat3!(csmatout, angles[fe_label]/180.0*pi* [0.0; 0.0; 1.0]);
                end
                
                # The vvolume integrals are evaluated using this rule
                gr = SimplexRule(3, 4)
                
                # We will create 3 regions, one for each of the layers
                regions = FDataDict[]
                for layer = 1:nLayers
                    rls = selectelem(fens, fes, label =  layer)
                    push!(regions, FDataDict("femm"=>FEMMDeforLinear(MR, IntegData(subset(fes, rls), gr), CSys(3, 3, updatecs!), skinmaterial)))
                end
                
                # The essential boundary conditions are applied to enforce the plane strain constraint.
                ly0 = selectnode(fens, box=[-Inf Inf 0.0 0.0 -Inf Inf], inflate=tolerance)
                lyh = selectnode(fens, box=[-Inf Inf h h -Inf Inf], inflate=tolerance)
                ey = FDataDict("displacement"=>  0.0, "component"=> 2, "node_list"=>vcat(ly0, lyh))
                # The transverse displacement is fixed at the two ends.
                lz0 = selectnode(fens, box=[0.0 0.0 -Inf Inf -Inf Inf], inflate=tolerance)
                lzL = selectnode(fens, box=[L L -Inf Inf -Inf Inf], inflate=tolerance)
                ez = FDataDict("displacement"=>  0.0, "component"=> 3, "node_list"=>vcat(lz0, lzL))
                ex = FDataDict("displacement"=>  0.0, "component"=> 1, "node_list"=>[1])
                
                # The traction boundary condition is applied at the top of the plate.
                bfes = meshboundary(fes)
                function pfun(forceout::FVec{T}, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt) where {T}
                    forceout[1] = 0.0
                    forceout[2] = 0.0
                    forceout[3] = -q0*sin(pi*XYZ[1]/L)
                    return forceout
                end
                # From  the entire boundary we select those quadrilaterals that lie on the plane
                # Z = thickness
                tl = selectelem(fens, bfes, box = [-Inf Inf -Inf Inf T T], inflate=tolerance)
                Trac = FDataDict("traction_vector"=>pfun,
                "femm"=>FEMMBase(IntegData(subset(bfes, tl), SimplexRule(2, 3))))
                
                modeldata = FDataDict("fens"=>fens,
                "regions"=>regions,
                "essential_bcs"=>[ex, ey, ez],
                "traction_bcs"=> [Trac]
                )
                modeldata = AlgoDeforLinearModule.linearstatics(modeldata)
                
                # modeldata["postprocessing"] = FDataDict("file"=>filebase * "-u")
                # modeldata = AlgoDeforLinearModule.exportdeformation(modeldata)
                
                u = modeldata["u"]
                geom = modeldata["geom"]
                
                # The results of the displacement and stresses will be reported at
                # nodes located at the appropriate points.
                ntopcenter = selectnode(fens, box=[L/2 L/2 0.0 h T T], inflate=tolerance)
                ncenterline = selectnode(fens, box=[L/2 L/2 0.0 0.0 0.0 T], inflate=tolerance)
                nx0line = selectnode(fens, box=[0.0 0.0 0.0 0.0 0.0 T], inflate=tolerance)
                
                zclo = sortperm(vec(geom.values[ncenterline, 3]))
                ncenterline = ncenterline[zclo]
                centerz = geom.values[ncenterline, 3]
                
                println("Top Center deflection: $(mean(u.values[ntopcenter, 3], 1)/phun("in")) [in]")
                
                #  Compute  all stresses
                modeldata["postprocessing"] = FDataDict("file"=>filebase * "-s",  "quantity"=>:Cauchy, "component"=>collect(1:6), "outputcsys"=>CSys(3), "nodevalmethod"=>nodevalmeth, "reportat"=>extrap)
                modeldata = AlgoDeforLinearModule.exportstress(modeldata)
                
                modeldata["elementsize"] = T/Refinement
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
        end # Span_to_thickness 
    end # extrap
    println("Done")
    
end # Pagano_3lay_cyl_bend_T10_conv

function allrun()
    println("#####################################################") 
    println("# Pagano_3lay_cyl_bend_MST10_conv ")
    Pagano_3lay_cyl_bend_MST10_conv()
    println("#####################################################") 
    println("# Pagano_3lay_cyl_bend_MSH8_conv ")
    Pagano_3lay_cyl_bend_MSH8_conv()
    println("#####################################################") 
    println("# Pagano_3lay_cyl_bend_H8_conv ")
    Pagano_3lay_cyl_bend_H8_conv()
    println("#####################################################") 
    println("# Pagano_3lay_cyl_bend_MSH8_conv ")
    Pagano_3lay_cyl_bend_MSH8_conv()
    println("#####################################################") 
    println("# Pagano_3lay_cyl_bend_T10_conv ")
    Pagano_3lay_cyl_bend_T10_conv()
    return true
end # function allrun

end # module Pagano_examples
