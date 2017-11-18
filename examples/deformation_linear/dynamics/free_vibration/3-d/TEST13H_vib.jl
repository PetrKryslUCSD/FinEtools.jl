module mTEST13H_vib
using FinEtools
using FinEtools.AlgoDeforLinearModule: ssit
function test()
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
        d = d - OmegaShift;
        fs = real(sqrt.(complex(d)))/(2*pi)
        println("Reference Eigenvalues: $fs [Hz]")
        println("eigs solution ($(time() - t0) sec)")
    end
    true
end
end
using mTEST13H_vib
mTEST13H_vib.test()
