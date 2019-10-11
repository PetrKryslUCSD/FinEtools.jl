
module mtfieldtm1
using FinEtools
using Test
using BenchmarkTools
function test()
    W = 1.1;
    L = 12.;
    t =  0.32;
    nl, nt, nw = 2, 3, 4;

    fens,fes  = H8block(L,W,t, nl,nw,nt)
    geom  =  NodalField(fens.xyz)

ecoords = fill(zero(FFlt), nodesperelem(fes), ndofs(geom)); # array of Element coordinates
@btime gathervalues_asmat!(geom, ecoords, fes.conn[1]);


    femm  =  FEMMBase(IntegDomain(fes, GaussRule(3, 2)))
    V = integratefunction(femm, geom, (x) ->  1.0)
    @test abs(V - W*L*t)/V < 1.0e-5
end
end
using .mtfieldtm1
mtfieldtm1.test()
