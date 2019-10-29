
module mxmeasure1
using FinEtools
using Test
function test()
    W = 4.1;
    L = 12.;
    t =  3.32;
    nl, nt, nw = 2, 3, 4;

    fens,fes  = H27block(L,W,t, nl,nw,nt)
    geom  =  NodalField(fens.xyz)

    femm  =  FEMMBase(IntegDomain(fes, GaussRule(3, 3)))

    # Test the calculation of the volume
    V = integratefunction(femm, geom, (x) ->  1.0)
    @test abs(V - W*L*t)/V < 1.0e-5

    # Test the calculation of the center of gravity
    Sx = integratefunction(femm, geom, (x) ->  x[1])
    @test Sx / V ≈ L / 2

    # Test the calculation of the moments of inertia
    # The block is translated to be centered at the origin
    fens.xyz[:, 1] .-= L / 2
    fens.xyz[:, 2] .-= W / 2
    fens.xyz[:, 3] .-= t / 2
    geom  =  NodalField(fens.xyz)
    Sx = integratefunction(femm, geom, (x) ->  x[1])
    @test abs(Sx / V) / L <= eps(1.0)
    Sy = integratefunction(femm, geom, (x) ->  x[2])
    @test abs(Sy / V) / W <= eps(1.0)
    Sz = integratefunction(femm, geom, (x) ->  x[3])
    @test abs(Sz / V) / t <= eps(1.0)
    Ixx = integratefunction(femm, geom, (x) ->  x[2]^2 + x[3]^2)
    @test abs(Ixx - V * (W^2 + t^2) / 12) / V <= eps(V)
    Iyy = integratefunction(femm, geom, (x) ->  x[1]^2 + x[3]^2)
    @test abs(Iyy - V * (L^2 + t^2) / 12) / V <= eps(V)
    Ixy = integratefunction(femm, geom, (x) ->  (-x[1] * x[2]))
    @test abs(Ixy) / V <= eps(V)
    true
end
end
using .mxmeasure1
mxmeasure1.test()
