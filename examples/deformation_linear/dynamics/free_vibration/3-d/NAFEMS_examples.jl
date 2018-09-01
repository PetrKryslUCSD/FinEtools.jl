module NAFEMS_examples
using FinEtools
using FinEtools
using FinEtools.AlgoDeforLinearModule: ssit
using LinearAlgebra
using Arpack

function NAFEMS_FV32_algo()
    println("""
    FV32: Cantilevered tapered membrane
    This is a test recommended by the National Agency for Finite Element Methods and
    Standards (U.K.): Test FV32 from NAFEMS publication TNSB, Rev. 3, “The Standard
    NAFEMS Benchmarks,” October 1990.
    
    Reference solution: 44.623	130.03	162.70	246.05	379.90	391.44 for the first
    six modes.
    """)
    
    t0 = time()
    
    
    E = 200*phun("GPA");
    nu = 0.3;
    rho= 8000*phun("KG/M^3");
    L = 10*phun("M");
    W0 = 5*phun("M");
    WL = 1*phun("M");
    H = 0.05*phun("M");
    nL, nW, nH = 8, 4, 1;# How many element edges per side?
    neigvs=20                   # how many eigenvalues
    Reffs = [44.623	130.03	162.70	246.05	379.90	391.44]
    
    fens,fes =H20block(1.0, 2.0, 1.0, nL, nW, nH)
    for i = 1:count(fens)
        xi, eta, theta = fens.xyz[i,:];
        eta = eta - 1.0
        fens.xyz[i,:] = [xi*L eta*(1.0 - 0.8*xi)*W0/2 theta*H/2];
    end
    # File =  "mesh.vtk"
    # vtkexportmesh(File, fens, fes)
    # @async run(`"paraview.exe" $File`)
    
    # Make the region
    MR = DeforModelRed3D
    material = MatDeforElastIso(MR, rho, E, nu, 0.0)
    region1 = FDataDict("femm"=>FEMMDeforLinear(MR, IntegData(fes, GaussRule(3,2)),
    material), "femm_mass"=>FEMMDeforLinear(MR, IntegData(fes, GaussRule(3,3)),
    material))
    
    nl1 = selectnode(fens; plane=[1.0 0.0 0.0 0.0], thickness=H/1.0e4)
    ebc1 = FDataDict("node_list"=>nl1, "component"=>1, "displacement"=>0.0)
    ebc2 = FDataDict("node_list"=>nl1, "component"=>2, "displacement"=>0.0)
    ebc3 = FDataDict("node_list"=>nl1, "component"=>3, "displacement"=>0.0)
    
    nl4 = selectnode(fens; plane=[0.0 0.0 1.0 0.0], thickness=H/1.0e4)
    ebc4 = FDataDict("node_list"=>nl4, "component"=>3, "displacement"=>0.0)
    
    # Make model data
    modeldata =  FDataDict(
    "fens"=> fens, "regions"=>  [region1], "essential_bcs"=>[ebc1 ebc2 ebc3 ebc4],
    "neigvs"=>neigvs)
    
    # Solve
    modeldata = FinEtools.AlgoDeforLinearModule.modal(modeldata)
    
    fs = modeldata["omega"]/(2*pi)
    println("Eigenvalues: $fs [Hz]")
    println("Percentage frequency errors: $((vec(fs[1:6]) - vec(Reffs))./vec(Reffs)*100)")
    
    modeldata["postprocessing"] = FDataDict("file"=>"FV32-modes", "mode"=>1:10)
    modeldata=FinEtools.AlgoDeforLinearModule.exportmode(modeldata)
    @async run(`"paraview.exe" $(modeldata["postprocessing"]["file"]*"1.vtk")`)
    
    true
    
end # NAFEMS_FV32_algo
    
    
function NAFEMS_TEST13H_vib()
    # Harmonic forced vibration problem is solved for a homogeneous square plate,
    # simply-supported on the circumference.
    # This is the TEST 13H from the Abaqus v 6.12 Benchmarks manual.
    # The test is recommended by the National Agency for Finite Element Methods and Standards (U.K.):
    # Test 13 from NAFEMS “Selected Benchmarks for Forced Vibration,” R0016, March 1993.
    #
    #
    # The plate is discretized with hexahedral solid elements. The simple support
    # condition is approximated by distributed rollers on the boundary.
    # Because only the out of plane displacements are prevented, the structure
    # has three rigid body modes in the plane of the plate.
    #
    #
    # The nonzero benchmark frequencies are (in hertz): 2.377, 5.961, 5.961,
    # 9.483, 12.133, 12.133, 15.468, 15.468 [Hz].
    
    println("""
    Homogeneous square plate, simply-supported on the circumference from
    the test 13 from NAFEMS “Selected Benchmarks for Forced Vibration,” R0016, March 1993.
    The nonzero benchmark frequencies are (in hertz): 2.377, 5.961, 5.961,
    9.483, 12.133, 12.133, 15.468, 15.468 [Hz].
    """)
    
    # t0 = time()
    
    E = 200*phun("GPa");# Young's modulus
    nu = 0.3;# Poisson ratio
    rho = 8000*phun("KG*M^-3");# mass density
    L = 10.0*phun("M"); # side of the square plate
    t = 0.05*phun("M"); # thickness of the square plate
    nL = 8; nt = 4;
    tolerance = t/nt/100;
    neigvs = 11;
    OmegaShift = (2*pi*0.5) ^ 2; # to resolve rigid body modes
    
    MR = DeforModelRed3D
    fens,fes  = H8block(L, L, t, nL, nL, nt)
    
    geom = NodalField(fens.xyz)
    u = NodalField(zeros(size(fens.xyz,1),3)) # displacement field
    nl = selectnode(fens, box=[0.0 0.0 -Inf Inf -Inf Inf], inflate=tolerance)
    setebc!(u, nl, true, 3)
    nl = selectnode(fens, box=[L L -Inf Inf -Inf Inf], inflate=tolerance)
    setebc!(u, nl, true, 3)
    nl = selectnode(fens, box=[-Inf Inf 0.0 0.0 -Inf Inf], inflate=tolerance)
    setebc!(u, nl, true, 3)
    nl = selectnode(fens, box=[-Inf Inf L L -Inf Inf], inflate=tolerance)
    setebc!(u, nl, true, 3)
    applyebc!(u)
    numberdofs!(u)
    println("nfreedofs = $(u.nfreedofs)")
    
    material=MatDeforElastIso(MR, rho, E, nu, 0.0)
    
    femm = FEMMDeforLinearMSH8(MR, IntegData(fes, GaussRule(3,2)), material)
    femm = associategeometry!(femm, geom)
    K = stiffness(femm, geom, u)
    femm = FEMMDeforLinear(MR, IntegData(fes, GaussRule(3,3)), material)
    M = mass(femm, geom, u)
    
    if true
        t0 = time()
        d,v,nev,nconv = eigs(K+OmegaShift*M, M; nev=neigvs, which=:SM)
        d = d .- OmegaShift;
        fs = real(sqrt.(complex(d)))/(2*pi)
        println("Reference Eigenvalues: $fs [Hz]")
        println("eigs solution ($(time() - t0) sec)")
    end
    true
    
end # NAFEMS_TEST13H_vib

function allrun()
    println("#####################################################") 
    println("# NAFEMS_FV32_algo ")
    NAFEMS_FV32_algo()
    println("#####################################################") 
    println("# NAFEMS_TEST13H_vib ")
    NAFEMS_TEST13H_vib()
    return true
end # function allrun

end # module NAFEMS_examples
