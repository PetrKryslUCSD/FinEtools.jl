module unit_cube_mode_examples
using FinEtools
using FinEtools.MeshExportModule
using LinearAlgebra
using Arpack

function unit_cube_modes()
    println("""
    Vibration modes of unit cube  of almost incompressible material.
    %
    Reference: Puso MA, Solberg J (2006) A stabilized nodally integrated
    tetrahedral. International Journal for Numerical Methods in
    Engineering 67: 841-867.
    """)
    t0 = time()
    
    E = 1*phun("PA");
    nu = 0.499;
    rho = 1*phun("KG/M^3");
    a = 1*phun("M"); b = a; h =  a;
    n1 = 10;# How many element edges per side?
    na =  n1; nb =  n1; nh  = n1;
    neigvs = 20                   # how many eigenvalues
    OmegaShift = (0.01*2*pi)^2;
    
    MR = DeforModelRed3D
    fens,fes  = H20block(a,b,h, na,nb,nh)
    
    geom = NodalField(fens.xyz)
    u = NodalField(zeros(size(fens.xyz,1),3)) # displacement field
    
    numberdofs!(u)
    
    material = MatDeforElastIso(MR, rho, E, nu, 0.0)
    
    femm = FEMMDeforLinear(MR, IntegData(fes, GaussRule(3,2)), material)
    
    @time K  = stiffness(femm, geom, u)
    femm = FEMMDeforLinear(MR, IntegData(fes, GaussRule(3,3)), material)
    @time M = mass(femm, geom, u)
    d,v,nev,nconv = eigs(K+OmegaShift*M, M; nev=neigvs, which=:SM)
    d = d .- OmegaShift;
    fs = real(sqrt.(complex(d)))/(2*pi)
    println("Eigenvalues: $fs [Hz]")
    
    mode = 17
    scattersysvec!(u, v[:,mode])
    File =  "unit_cube_modes.vtk"
    vtkexportmesh(File, fens, fes; vectors=[("mode$mode", u.values)])
    
    true
    
end # unit_cube_modes


function unit_cube_modes_algo()
    println("""
    % Vibration modes of unit cube  of almost incompressible material.
    %
    % Reference: Puso MA, Solberg J (2006) A stabilized nodally integrated
    % tetrahedral. International Journal for Numerical Methods in
    % Engineering 67: 841-867.""")
    t0 = time()
    
    
    E = 1*phun("PA");
    nu = 0.499;
    rho= 1*phun("KG/M^3");
    a=1*phun("M"); b=a; h= a;
    n1=2;# How many element edges per side?
    na= n1; nb= n1; nh =n1;
    neigvs=20                   # how many eigenvalues
    omega_shift=(0.1*2*pi)^2;
    
    fens,fes =H20block(a,b,h, na,nb,nh)
    
    # Make the region
    MR = DeforModelRed3D
    material = MatDeforElastIso(MR, rho, E, nu, 0.0)
    region1 = FDataDict("femm"=>FEMMDeforLinear(MR, IntegData(fes, GaussRule(3,2)),
    material), "femm_mass"=>FEMMDeforLinear(MR, IntegData(fes, GaussRule(3,3)),
    material))
    
    # Make model data
    modeldata =  FDataDict(
    "fens"=> fens, "regions"=>  [region1],
    "omega_shift"=>omega_shift, "neigvs"=>neigvs)
    
    # Solve
    modeldata = FinEtools.AlgoDeforLinearModule.modal(modeldata)
    
    fs = modeldata["omega"]/(2*pi)
    println("Eigenvalues: $fs [Hz]")
    
    modeldata["postprocessing"] = FDataDict("file"=>"unit_cube_mode",
    "mode"=>10)
    modeldata=FinEtools.AlgoDeforLinearModule.exportmode(modeldata)
    @async run(`"paraview.exe" $(modeldata["postprocessing"]["file"]*"1.vtk")`)
    
    true
    
end # unit_cube_modes_algo


function unit_cube_modes_export()
    println("""
    Vibration modes of unit cube  of almost incompressible material.
    
    This example EXPORTS the model to Abaqus.
    
    Reference: Puso MA, Solberg J (2006) A stabilized nodally integrated
    tetrahedral. International Journal for Numerical Methods in
    Engineering 67: 841-867.
    """)
    t0 = time()
    
    
    E = 1*phun("PA");
    nu = 0.499;
    rho = 1*phun("KG/M^3");
    a = 1*phun("M"); b = a; h =  a;
    n1 = 5;# How many element edges per side?
    na =  n1; nb =  n1; nh  = n1;
    neigvs = 20                   # how many eigenvalues
    OmegaShift = (0.01*2*pi)^2;
    
    MR = DeforModelRed3D
    fens,fes  = H20block(a,b,h, na,nb,nh)
    
    geom = NodalField(fens.xyz)
    u = NodalField(zeros(size(fens.xyz,1),3)) # displacement field
    
    numberdofs!(u)
    
    material=MatDeforElastIso(MR, rho, E, nu, 0.0)
    
    femm = FEMMDeforLinear(MR, IntegData(fes, GaussRule(3,2)), material)
    
    K =stiffness(femm, geom, u)
    femm = FEMMDeforLinear(MR, IntegData(fes, GaussRule(3,3)), material)
    M =mass(femm, geom, u)
    d,v,nev,nconv = eigs(K+OmegaShift*M, M; nev=neigvs, which=:SM)
    d = d .- OmegaShift;
    fs = real(sqrt.(complex(d)))/(2*pi)
    println("Eigenvalues: $fs [Hz]")
    
    mode = 7
    scattersysvec!(u, v[:,mode])
    File =  "unit_cube_modes.vtk"
    vtkexportmesh(File, fens, fes; vectors=[("mode$mode", u.values)])
    @async run(`"paraview.exe" $File`)
    
    AE = AbaqusExporter("unit_cube_modes_h20");
    # AE.ios = STDOUT;
    HEADING(AE, "Vibration modes of unit cube  of almost incompressible material.");
    COMMENT(AE, "The  first six frequencies are rigid body modes.");
    COMMENT(AE, "The  first nonzero frequency (7) should be around 0.26 Hz");
    PART(AE, "part1");
    END_PART(AE);
    ASSEMBLY(AE, "ASSEM1");
    INSTANCE(AE, "INSTNC1", "PART1");
    NODE(AE, fens.xyz);
    COMMENT(AE, "The hybrid form of the serendipity hexahedron is chosen because");
    COMMENT(AE, "the material is  nearly incompressible.");
    ELEMENT(AE, "c3d20rh", "AllElements", 1, connasarray(fes))
    ORIENTATION(AE, "GlobalOrientation", vec([1. 0 0]), vec([0 1. 0]));
    SOLID_SECTION(AE, "elasticity", "GlobalOrientation", "AllElements");
    END_INSTANCE(AE);
    END_ASSEMBLY(AE);
    MATERIAL(AE, "elasticity")
    ELASTIC(AE, E, nu)
    DENSITY(AE, rho)
    STEP_FREQUENCY(AE, neigvs)
    END_STEP(AE)
    close(AE)
    
    true
    
end # unit_cube_modes_export


function unit_cube_modes_msh8_algo()
    println("""
    % Vibration modes of unit cube  of almost incompressible material.
    % Mean-strain hexahedron.
    % Reference: Puso MA, Solberg J (2006) A stabilized nodally integrated
    % tetrahedral. International Journal for Numerical Methods in
    % Engineering 67: 841-867.""")
    t0 = time()
    
    
    E = 1*phun("PA");
    nu = 0.499;
    rho= 1*phun("KG/M^3");
    a=1*phun("M"); b=a; h= a;
    n1=8 # How many element edges per side?
    na= n1; nb= n1; nh =n1;
    neigvs=20                   # how many eigenvalues
    omega_shift=(0.1*2*pi)^2;
    
    fens,fes = H8block(a,b,h, na,nb,nh)
    
    # Make the region
    MR = DeforModelRed3D
    material = MatDeforElastIso(MR, rho, E, nu, 0.0)
    region1 = FDataDict("femm"=>FEMMDeforLinearMSH8(MR, IntegData(fes, GaussRule(3,2)),
    material), "femm_mass"=>FEMMDeforLinearMSH8(MR, IntegData(fes, GaussRule(3,3)),
    material))
    
    # Make model data
    modeldata =  FDataDict(
    "fens"=> fens, "regions"=>  [region1],
    "omega_shift"=>omega_shift, "neigvs"=>neigvs)
    
    # Solve
    modeldata = FinEtools.AlgoDeforLinearModule.modal(modeldata)
    
    fs = modeldata["omega"]/(2*pi)
    println("Eigenvalues: $fs [Hz]")
    
    modeldata["postprocessing"] = FDataDict("file"=>"unit_cube_mode",
    "mode"=>10)
    modeldata=FinEtools.AlgoDeforLinearModule.exportmode(modeldata)
    @async run(`"paraview.exe" $(modeldata["postprocessing"]["file"]*"1.vtk")`)
    
    true
    
end # unit_cube_modes_msh8_algo

function allrun()
    println("#####################################################") 
    println("# unit_cube_modes ")
    unit_cube_modes()
    println("#####################################################") 
    println("# unit_cube_modes_algo ")
    unit_cube_modes_algo()
    println("#####################################################") 
    println("# unit_cube_modes_export ")
    unit_cube_modes_export()
    println("#####################################################") 
    println("# unit_cube_modes_msh8_algo ")
    unit_cube_modes_msh8_algo()
    return true
end # function allrun

end # module unit_cube_mode_examples
