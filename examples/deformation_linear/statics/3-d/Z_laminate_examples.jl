module Z_laminate_examples
using FinEtools
using FinEtools.AlgoDeforLinearModule


function Z_laminate_u_ss()
    println("""
    % Three-dimensional Elasticity Solution for Uniformly Loaded Cross-ply
    % Laminates and Sandwich Plates
    % Ashraf M. Zenkour, Journal of Sandwich Structures and Materials 2007 9: 213-238
    % DOI: 10.1177/1099636207065675
    """)
    
    t0 = time()
    # Lamina material parameters
    E1s = 25.0e6*phun("psi")
    E2s = 1.0e6*phun("psi")
    E3s = E2s
    nu12s = nu13s = nu23s = 0.25
    G12s = 0.5e6*phun("psi")
    G13s = G12s
    G23s = 0.2e6*phun("psi")
    
    a = 200.0*phun("mm") # side of the square plate
    b = 600.0*phun("mm") # side of the square plate
    q0 = 1.0*phun("psi");
    # The below values come from Table 2
    # h = a/4; wc_analytical = 3.65511/(100*E3s*h^3/a^4/q0);
    # h = a/10; wc_analytical = 1.16899/(100*E3s*h^3/a^4/q0);
    # h = a/50; wc_analytical = 0.66675/(100*E3s*h^3/a^4/q0);
    h = a/100; wc_analytical = 0.65071/(100*E3s*h^3/a^4/q0);
    angles =[0,90,0];
    nLayers = length(angles);
    tolerance = 0.0001*h
    
    
    
    # Generate mesh
    na = 10 # number of elements along the side of the plate
    nb = 30 # number of elements along the side of the plate
    xs = collect(linspace(0.0, a, na+1))
    ys = collect(linspace(0.0, b, nb+1))
    ts = h/nLayers*ones(nLayers);# layer thicknesses
    nts= 3*ones(Int, nLayers);# number of elements per layer
    fens,fes = H8layeredplatex(xs, ys, ts, nts)
    fens,fes = H8toH20(fens,fes)
    
    MR = DeforModelRed3D
    laminamaterial = MatDeforElastOrtho(MR,
    0.0, E1s, E2s, E3s,
    nu12s, nu13s, nu23s,
    G12s, G13s, G23s,
    0.0, 0.0, 0.0)
    
    function updatecs!(csmatout::FFltMat, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)
        rotmat3!(csmatout, angles[fe_label]/180.*pi* [0.;0.;1.]);
    end
    
    gr = GaussRule(3, 3)
    
    region = FDataDict("femm"=>FEMMDeforLinear(MR, IntegData(fes, gr), CSys(3, 3, updatecs!), laminamaterial))
    
    lx0 = selectnode(fens, box=[0.0 0.0 -Inf Inf -Inf Inf], inflate=tolerance)
    lxa = selectnode(fens, box=[a a -Inf Inf -Inf Inf], inflate=tolerance)
    ly0 = selectnode(fens, box=[-Inf Inf 0.0 0.0 -Inf Inf], inflate=tolerance)
    lyb = selectnode(fens, box=[-Inf Inf b b -Inf Inf], inflate=tolerance)
    
    ex02 = FDataDict( "displacement"=>  0.0, "component"=> 2, "node_list"=>lx0 )
    ex03 = FDataDict( "displacement"=>  0.0, "component"=> 3, "node_list"=>lx0 )
    exa2 = FDataDict( "displacement"=>  0.0, "component"=> 2, "node_list"=>lxa )
    exa3 = FDataDict( "displacement"=>  0.0, "component"=> 3, "node_list"=>lxa )
    ey01 = FDataDict( "displacement"=>  0.0, "component"=> 1, "node_list"=>ly0 )
    ey03 = FDataDict( "displacement"=>  0.0, "component"=> 3, "node_list"=>ly0 )
    eyb1 = FDataDict( "displacement"=>  0.0, "component"=> 1, "node_list"=>lyb )
    eyb3 = FDataDict( "displacement"=>  0.0, "component"=> 3, "node_list"=>lyb )
    
    bfes = meshboundary(fes)
    ttopl = selectelem(fens, bfes; facing=true, direction = [0.0 0.0 1.0])
    Trac = FDataDict("traction_vector"=>[0.0; 0.0; -q0],
    "femm"=>FEMMBase(IntegData(subset(bfes, ttopl), GaussRule(2, 3))))
    
    modeldata = FDataDict("fens"=>fens,
    "regions"=>[region],
    "essential_bcs"=>[ex02, ex03, exa2, exa3, ey01, ey03, eyb1, eyb3],
    "traction_bcs"=> [Trac]
    )
    modeldata = AlgoDeforLinearModule.linearstatics(modeldata)
    
    u = modeldata["u"]
    geom = modeldata["geom"]
    lcenter = selectnode(fens, box=[a/2 a/2  b/2 b/2 -Inf Inf], inflate=tolerance)
    cdis = abs(mean(u.values[lcenter, 3]))
    println("")
    println("Normalized Center deflection: $(cdis/wc_analytical)")
    
    File =  "Z_laminate_u_ss.vtk"
    vtkexportmesh(File, fes.conn, geom.values, FinEtools.MeshExportModule.H20;
    scalars = [("Layer", fes.label)], vectors = [("displacement", u.values)])
    @async run(`"paraview.exe" $File`)
    
    println("Done")
    true
    
end # Z_laminate_u_ss

function allrun()
    println("#####################################################") 
    println("# Z_laminate_u_ss ")
    Z_laminate_u_ss()
    return true
end # function allrun

end # module Z_laminate_examples
