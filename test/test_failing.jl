module mocylpull14nnn
using FinEtools
using Test
using InteractiveUtils
function test()
    # Cylinder  compressed by enforced displacement, axially symmetric model


    # Parameters:
    E1=1.0;
    E2=1.0;
    E3=3.0;
    nu12=0.29;
    nu13=0.29;
    nu23=0.19;
    G12=0.3;
    G13=0.3;
    G23=0.3;
    p= 0.15;
    rin=1.;
    rex =1.2;
    Length = 1*rex
    ua = -0.05*Length
    tolerance=rin/1000.

    ##
    # Note that the FinEtools objects needs to be created with the proper
    # model-dimension reduction at hand.  In this case that is the axial symmetry
    # assumption.
    MR = DeforModelRed2DAxisymm

    fens,fes = Q4block(rex-rin,Length,5,20);
    fens.xyz[:, 1] = fens.xyz[:, 1] .+ rin
    bdryfes = meshboundary(fes);

    # now we create the geometry and displacement fields
    geom = NodalField(fens.xyz)
    u = NodalField(zeros(size(fens.xyz,1),2)) # displacement field

    # # the symmetry plane
    # l1 =selectnode(fens; box=[0 rex 0 0], inflate = tolerance)
    # setebc!(u,l1,true, 2, 0.0)
    # # The other end
    # l1 =selectnode(fens; box=[0 rex Length Length], inflate = tolerance)
    # setebc!(u,l1,true, 2, ua)

    # applyebc!(u)
    # numberdofs!(u)
    # # println("Number of degrees of freedom = $(u.nfreedofs)")
    # @test u.nfreedofs == 240

    # Property and material
    material = MatDeforElastIso(MR, 00.0, E1, nu23, 0.0)
    
    femm = FEMMDeforLinear(MR, IntegData(fes, GaussRule(2, 2), true), material, true)
    println("========== With printing ==========")
    @show @code_lowered FEMMDeforLinear(MR, IntegData(fes, GaussRule(2, 2), true), material, true)
    println("========== Original ==========")
    @show @code_lowered FEMMDeforLinear(MR, IntegData(fes, GaussRule(2, 2), true), material)
    femm = FEMMDeforLinear(MR, IntegData(fes, GaussRule(2, 2), true), material)
    
    true
end
end
using .mocylpull14nnn
mocylpull14nnn.test()
