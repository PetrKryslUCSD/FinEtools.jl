
module mmmeasurementm4
using FinEtools
using Base.Test
function test()
    W = 1.1;
    L = 12.;
    t =  0.32;
    nl, nt, nw = 2, 3, 4;

    fens,fes  = L2block(L,nl)
    geom  =  NodalField(fens.xyz)
    bfes = meshboundary(fes)
    femm  =  FEMMBase(GeoD(bfes, PointRule()))
    S = integratefunction(femm, geom, (x) ->  1.0)
    @test abs(S - 2)/S < 1.0e-5
end
end
using mmmeasurementm4
mmmeasurementm4.test()
