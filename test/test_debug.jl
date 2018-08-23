module mmMeasurement_3a
using FinEtools
using Test
function test()
    W = 1.1;
    L = 12.;
    t =  3.32;
    nl, nt, nw = 5, 3, 4;
    Ea =  210000*phun("MEGA*Pa")
    nua =  0.3;

    # println("New segmentation fault?")
    for orientation in [:a :b :ca :cb]
        fens,fes  = T4block(L,W,t, nl,nw,nt, orientation)
        geom  =  NodalField(fens.xyz)
        MR  =  DeforModelRed3D
        material = MatDeforElastIso(MR, 0.0, Ea, nua, 0.0)

        femm  =  FEMMDeforLinearNICET4(MR, IntegData(fes, NodalSimplexRule(3)), material)
        V = integratefunction(femm, geom, (x) ->  1.0)
        @test abs(V - W*L*t)/V < 1.0e-5
    end

end
end
using .mmMeasurement_3a
mmMeasurement_3a.test()
