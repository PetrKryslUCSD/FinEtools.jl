module mocylpull14
using FinEtools
using Compat.Test
function test()
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

    MR = DeforModelRed2DAxisymm

    fens,fes = Q4block(rex-rin,Length,5,20);
    fens.xyz[:, 1] = fens.xyz[:, 1] .+ rin
    bdryfes = meshboundary(fes);

    # now we create the geometry and displacement fields
    geom = NodalField(fens.xyz)
    u = NodalField(zeros(size(fens.xyz,1),2)) # displacement field

    # Property and material
    material=MatDeforElastIso(MR, 00.0, E1, nu23, 0.0)
    # @show mr 
    # @show material.mr

    femm = FEMMDeforLinear(MR, IntegData(fes, GaussRule(2, 2), true), material)

    # K =stiffness(femm, geom, u)
    # F = nzebcloadsstiffness(femm, geom, u)
    # U=  K\(F)
    # scattersysvec!(u,U[:])

    # fld= fieldfromintegpoints(femm, geom, u, :princCauchy, 1)
    # # println("Minimum/maximum = $(minimum(fld.values))/$(maximum(fld.values))")
    # @test abs(minimum(fld.values) - 0.0) < 1.0e-5
    # @test abs(maximum(fld.values) - 0.0) < 1.0e-5
    # fld= fieldfromintegpoints(femm, geom, u, :princCauchy, 2)
    # # println("Minimum/maximum = $(minimum(fld.values))/$(maximum(fld.values))")
    # @test abs(minimum(fld.values) - 0.0) < 1.0e-5
    # @test abs(maximum(fld.values) - 0.0) < 1.0e-5
    # fld= fieldfromintegpoints(femm, geom, u, :princCauchy, 3)
    # # println("Minimum/maximum = $(minimum(fld.values))/$(maximum(fld.values))")
    # @test abs(minimum(fld.values) - -0.050) < 1.0e-5
    # @test abs(maximum(fld.values) - -0.04999999999999919) < 1.0e-5
end
end
using .mocylpull14
mocylpull14.test()

