module mophysun13
using FinEtools
using Test
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

    geom = NodalField(fens.xyz)
    u = NodalField(zeros(size(fens.xyz,1),2)) 
    l1 =selectnode(fens; box=[0 rex 0 0], inflate = tolerance)
    setebc!(u,l1,true, 2, 0.0)
    l1 =selectnode(fens; box=[0 rex Length Length], inflate = tolerance)
    setebc!(u,l1,true, 2, 0.0)
    applyebc!(u)
    numberdofs!(u)
    @test u.nfreedofs == 240

    material=MatDeforElastIso(MR, 00.0, E1, nu23, 0.0)
    println("success? ")
    # @code_llvm FEMMDeforLinear(MR, IntegData(fes, GaussRule(2, 2), true), material, true)
    femm = FEMMDeforLinear(MR, IntegData(fes, GaussRule(2, 2), true), material, true)
    println("failure? ")
    # @code_llvm FEMMDeforLinear(MR, IntegData(fes, GaussRule(2, 2), true), material)
    femm = FEMMDeforLinear(MR, IntegData(fes, GaussRule(2, 2), true), material)

    true
end
end
using .mophysun13
mophysun13.test()

module mmMeasurement_1
using FinEtools
using Test
function test()
    W = 1.1;
    L = 12.;
    t =  3.32;
    nl, nt, nw = 5, 3, 4;

    println("New segmentation fault?")
    for orientation in [:a :b :ca :cb]
        fens,fes  = T4block(L,W,t, nl,nw,nt, orientation)
        geom  =  NodalField(fens.xyz)

        femm  =  FEMMBase(IntegData(fes, TetRule(5)))
        V = integratefunction(femm, geom, (x) ->  1.0)
        @test abs(V - W*L*t)/V < 1.0e-5
    end

end
end
using .mmMeasurement_1
mmMeasurement_1.test()

