module alum_cyl_mode_examples
using FinEtools
using FinEtools.MeshExportModule
using Arpack

# Mesh alum_cyl.inp
# Abaqus with the Standard isoparametric C3D4 tetrahedron:
C3D4 = [1 0
2 0
3 6.63586E-005
4 0.000171053
5 0.000211299
6 0.000244378
7 2564.63
8 2568.09
9 2597.26
10 4094.38
11 4714.36
12 4717.19
13 5181.98
14 6865.13
15 6868.17
16 6962.86
17 6965.67
18 7024.97
19 7029.44
20 7108.54]
# Abaqus with the standard quadratic tetrahedron:
C3D10 = [1 0
2 0
3 0
4 0.000139365
5 0.000221551
6 0.000291805
7 2546.81
8 2546.81
9 2560.69
10 4100
11 4693.55
12 4693.56
13 5121.57
14 6841.21
15 6841.24
16 6914.22
17 6914.23
18 6950.64
19 6950.66
20 7000.64]

function alum_cyl_mode_nice_t4()
    
    E = 70000*phun("MPa");
    nu = 0.33;
    rho = 2700*phun("KG/M^3");
    radius = 0.5*phun("ft");
    neigvs = 20                   # how many eigenvalues
    OmegaShift = (10.0*2*pi)^2;
    stabfact = 0.0062
    
    MR = DeforModelRed3D
    output = import_ABAQUS("alum_cyl.inp")
    fens, fes = output["fens"], output["fesets"][1]
    fens.xyz .*= phun("mm") # The input is provided in SI(mm) units
    fens, fes = T10toT4(fens, fes)
     
    geom = NodalField(fens.xyz)
    u = NodalField(zeros(size(fens.xyz,1),3)) # displacement field
    
    numberdofs!(u)
    
    material = MatDeforElastIso(MR, rho, E, nu, 0.0)
    
    femm = FEMMDeforLinearNICET4(MR, IntegData(fes, NodalSimplexRule(3)), material, stabfact)
    associategeometry!(femm,  geom)
    K  = stiffness(femm, geom, u)
    M = mass(femm, geom, u)
    d,v,nev,nconv = eigs(K+OmegaShift*M, M; nev=neigvs, which=:SM)
    d = d .- OmegaShift;
    fs = real(sqrt.(complex(d)))/(2*pi)
    println("Eigenvalues: $fs [Hz]")
    
    true
    
end # alum_cyl_modes

function alum_cyl_mode_esnice_t4()
    
    E = 70000*phun("MPa");
    nu = 0.33;
    rho = 2700*phun("KG/M^3");
    radius = 0.5*phun("ft");
    neigvs = 20                   # how many eigenvalues
    OmegaShift = (10.0*2*pi)^2;
    
    MR = DeforModelRed3D
    output = import_ABAQUS("alum_cyl.inp")
    fens, fes = output["fens"], output["fesets"][1]
    fens.xyz .*= phun("mm") # The input is provided in SI(mm) units
    fens, fes = T10toT4(fens, fes)
     
    geom = NodalField(fens.xyz)
    u = NodalField(zeros(size(fens.xyz,1),3)) # displacement field
    
    numberdofs!(u)
    
    material = MatDeforElastIso(MR, rho, E, nu, 0.0)
    
    femm = FEMMDeforLinearESNICET4(MR, IntegData(fes, NodalSimplexRule(3)), material)
    associategeometry!(femm,  geom)
    K  = stiffness(femm, geom, u)
    M = mass(femm, geom, u)
    d,v,nev,nconv = eigs(K+OmegaShift*M, M; nev=neigvs, which=:SM)
    d = d .- OmegaShift;
    fs = real(sqrt.(complex(d)))/(2*pi)
    println("Eigenvalues: $fs [Hz]")
    
    true
    
end # alum_cyl_modes


function alum_cyl_modes_algo()
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
    
    modeldata["postprocessing"] = FDataDict("file"=>"alum_cyl_mode",
    "mode"=>10)
    modeldata=FinEtools.AlgoDeforLinearModule.exportmode(modeldata)
    @async run(`"paraview.exe" $(modeldata["postprocessing"]["file"]*"1.vtk")`)
    
    true
    
end # alum_cyl_modes_algo


function alum_cyl_modes_export()
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
    d = d - OmegaShift;
    fs = real(sqrt.(complex(d)))/(2*pi)
    println("Eigenvalues: $fs [Hz]")
    
    mode = 7
    scattersysvec!(u, v[:,mode])
    File =  "alum_cyl_modes.vtk"
    vtkexportmesh(File, fens, fes; vectors=[("mode$mode", u.values)])
    @async run(`"paraview.exe" $File`)
    
    AE = AbaqusExporter("alum_cyl_modes_h20");
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
    
end # alum_cyl_modes_export


function alum_cyl_modes_msh8_algo()
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
    
    modeldata["postprocessing"] = FDataDict("file"=>"alum_cyl_mode",
    "mode"=>10)
    modeldata=FinEtools.AlgoDeforLinearModule.exportmode(modeldata)
    @async run(`"paraview.exe" $(modeldata["postprocessing"]["file"]*"1.vtk")`)
    
    true
    
end # alum_cyl_modes_msh8_algo

function allrun()
    println("#####################################################") 
    println("# alum_cyl_mode_nice_t4 ")
    alum_cyl_mode_nice_t4()
    println("#####################################################") 
    println("# alum_cyl_mode_esnice_t4 ")
    alum_cyl_mode_esnice_t4()
    return true
end # function allrun

end # module alum_cyl_mode_examples
