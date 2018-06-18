module mocylpull14
using FinEtools
using Compat.Test
function test()
    E1=1.0;
    nu23=0.19;
    rin=1.;
    rex =1.2;
    Length = 1*rex
    tolerance=rin/1000.

    MR = DeforModelRed2DAxisymm

    fens,fes = Q4block(rex-rin,Length,5,20);
    fens.xyz[:, 1] = fens.xyz[:, 1] .+ rin
    bdryfes = meshboundary(fes);

    # now we create the geometry and displacement fields
    geom = NodalField(fens.xyz)
    u = NodalField(zeros(size(fens.xyz,1),2)) # displacement field

    # the symmetry plane
    l1 =selectnode(fens; box=[0 rex 0 0], inflate = tolerance)
    setebc!(u,l1,true, 2, 0.0)
    # The other end
    l1 =selectnode(fens; box=[0 rex Length Length], inflate = tolerance)
    setebc!(u,l1,true, 2, 0.0)

    applyebc!(u)
    numberdofs!(u)
    # println("Number of degrees of freedom = $(u.nfreedofs)")
    @test u.nfreedofs == 240

    material=MatDeforElastIso(MR, 00.0, E1, nu23, 0.0)
    @show mr 
    @show material.mr

    femm = FEMMDeforLinear(MR, IntegData(fes, GaussRule(2, 2), true), material)

    true
end
end
using .mocylpull14
mocylpull14.test()

