module m111ocylpull14n1 # From miscellaneous
using FinEtools
using Test
using InteractiveUtils
function test()
    E1=1.0;
    nu23=0.19;
    rin=1.;
    rex =1.2;
    Length = 1*rex

    MR = DeforModelRed2DAxisymm
    fens,fes = Q4block(rex-rin,Length,5,20);
    material = MatDeforElastIso(MR, 00.0, E1, nu23, 0.0)
    
    femm = FEMMDeforLinear(MR, IntegData(fes, GaussRule(2, 2), true), material, true)
    println("isaisaisa= With printing isaisaisa=")
    @show @code_lowered FEMMDeforLinear(MR, IntegData(fes, GaussRule(2, 2), true), material, true)
    println("isaisaisa= Original isaisaisa=")
    @show @code_lowered FEMMDeforLinear(MR, IntegData(fes, GaussRule(2, 2), true), material)
    femm = FEMMDeforLinear(MR, IntegData(fes, GaussRule(2, 2), true), material)
    
    true
end
end
using .m111ocylpull14n1
m111ocylpull14n1.test()

module m111ocylpull14nnn # From miscellaneous
using FinEtools
using Test
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

    # the symmetry plane
    l1 =selectnode(fens; box=[0 rex 0 0], inflate = tolerance)
    setebc!(u,l1,true, 2, 0.0)
    # The other end
    l1 =selectnode(fens; box=[0 rex Length Length], inflate = tolerance)
    setebc!(u,l1,true, 2, ua)

    applyebc!(u)
    numberdofs!(u)
    # println("Number of degrees of freedom = $(u.nfreedofs)")
    @test u.nfreedofs == 240

    # Property and material
    material = MatDeforElastIso(MR, 00.0, E1, nu23, 0.0)
    # display(material)
    # println("$(material.D)")
    # @show MR
    
    femm = FEMMDeforLinear(MR, IntegData(fes, GaussRule(2, 2), true), material)

    K =stiffness(femm, geom, u)
    F = nzebcloadsstiffness(femm, geom, u)
    U=  K\(F)
    scattersysvec!(u,U[:])

    fld= fieldfromintegpoints(femm, geom, u, :princCauchy, 1)
    # println("Minimum/maximum = $(minimum(fld.values))/$(maximum(fld.values))")
    @test abs(minimum(fld.values) - 0.0) < 1.0e-5
    @test abs(maximum(fld.values) - 0.0) < 1.0e-5
    fld= fieldfromintegpoints(femm, geom, u, :princCauchy, 2)
    # println("Minimum/maximum = $(minimum(fld.values))/$(maximum(fld.values))")
    @test abs(minimum(fld.values) - 0.0) < 1.0e-5
    @test abs(maximum(fld.values) - 0.0) < 1.0e-5
    fld= fieldfromintegpoints(femm, geom, u, :princCauchy, 3)
    # println("Minimum/maximum = $(minimum(fld.values))/$(maximum(fld.values))")
    @test abs(minimum(fld.values) - -0.050) < 1.0e-5
    @test abs(maximum(fld.values) - -0.04999999999999919) < 1.0e-5

    # File =  "mocylpull14.vtk"
    # vtkexportmesh(File, fens, fes; scalars=[("sigmaz", fld.values)],
    #               vectors=[("u", u.values)])
    # @async run(`"paraview.exe" $File`)
end
end
using .m111ocylpull14nnn
m111ocylpull14nnn.test()

module mophysun13 # From miscellaneous
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
    # println("success? ")
    # @code_llvm FEMMDeforLinear(MR, IntegData(fes, GaussRule(2, 2), true), material, true)
    # femm = FEMMDeforLinear(MR, IntegData(fes, GaussRule(2, 2), true), material, true)
    # println("failure? ")
    # @code_llvm FEMMDeforLinear(MR, IntegData(fes, GaussRule(2, 2), true), material)
    femm = FEMMDeforLinear(MR, IntegData(fes, GaussRule(2, 2), true), material)

    true
end
end
using .mophysun13
mophysun13.test()

module mocylpull14 # From linear deformation
using FinEtools
using Test
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

    # the symmetry plane
    l1 =selectnode(fens; box=[0 rex 0 0], inflate = tolerance)
    setebc!(u,l1,true, 2, 0.0)
    # The other end
    l1 =selectnode(fens; box=[0 rex Length Length], inflate = tolerance)
    setebc!(u,l1,true, 2, ua)

    applyebc!(u)
    numberdofs!(u)
    # println("Number of degrees of freedom = $(u.nfreedofs)")
    @test u.nfreedofs == 240

    # Property and material
    material=MatDeforElastIso(MR, 00.0, E1, nu23, 0.0)
    # display(material)
    # println("$(material.D)")

    femm = FEMMDeforLinear(MR, IntegData(fes, GaussRule(2, 2), true), material)

    K =stiffness(femm, geom, u)
    F = nzebcloadsstiffness(femm, geom, u)
    U=  K\(F)
    scattersysvec!(u,U[:])

    fld= fieldfromintegpoints(femm, geom, u, :princCauchy, 1)
    # println("Minimum/maximum = $(minimum(fld.values))/$(maximum(fld.values))")
    @test abs(minimum(fld.values) - 0.0) < 1.0e-5
    @test abs(maximum(fld.values) - 0.0) < 1.0e-5
    fld= fieldfromintegpoints(femm, geom, u, :princCauchy, 2)
    # println("Minimum/maximum = $(minimum(fld.values))/$(maximum(fld.values))")
    @test abs(minimum(fld.values) - 0.0) < 1.0e-5
    @test abs(maximum(fld.values) - 0.0) < 1.0e-5
    fld= fieldfromintegpoints(femm, geom, u, :princCauchy, 3)
    # println("Minimum/maximum = $(minimum(fld.values))/$(maximum(fld.values))")
    @test abs(minimum(fld.values) - -0.050) < 1.0e-5
    @test abs(maximum(fld.values) - -0.04999999999999919) < 1.0e-5

    # File =  "mocylpull14.vtk"
    # vtkexportmesh(File, fens, fes; scalars=[("sigmaz", fld.values)],
    #               vectors=[("u", u.values)])
    # @async run(`"paraview.exe" $File`)
end
end
using .mocylpull14
mocylpull14.test()


module mocylpull1 # From deformation
using FinEtools
using Test
function test()
    # Cylinder  pulled by enforced displacement, axially symmetric model


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


    # Create the mesh and initialize the geometry.  First we are going
    # to construct the block of elements with the first coordinate
    # corresponding to the thickness in the radial direction, and the second
    # coordinate is the thickness in the axial direction.

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
    setebc!(u,l1,true, 2, ua)

    applyebc!(u)
    numberdofs!(u)
    # println("Number of degrees of freedom = $(u.nfreedofs)")
    @test u.nfreedofs == 240

    # Property and material
    material=MatDeforElastOrtho(MR, E1,E2,E3,nu12,nu13,nu23,G12,G13,G23)

    femm = FEMMDeforLinear(MR, IntegData(fes, GaussRule(2, 2), true), material)

    K =stiffness(femm, geom, u)
    F = nzebcloadsstiffness(femm, geom, u)
    U=  K\(F)
    scattersysvec!(u,U[:])

    # Produce a plot of the radial stress component in the cylindrical
    # coordinate system. Note that this is the usual representation of
    # stress using nodal stress field.

    fld= fieldfromintegpoints(femm, geom, u, :Cauchy, 2)
    # println("Minimum/maximum = $(minimum(fld.values))/$(maximum(fld.values))")
    @test abs(minimum(fld.values) - -0.050318853446676436) < 1.0e-5
    @test abs(maximum(fld.values) - -0.0497395167360893) < 1.0e-5
    # File =  "orthoballoon_sigmaz.vtk"
    # vtkexportmesh(File, fens, fes; scalars=[("sigmaz", fld.values)],
    #               vectors=[("u", u.values)])
    # @async run(`"paraview.exe" $File`)
end
end
using .mocylpull1
mocylpull1.test()

module mocylpull13
using FinEtools
using Test
function test()
    # Cylinder  compressed by enforced displacement, axially symmetric model


    # Parameters: make sure this represents isotropic material
    E1=1.0;
    E2=E1;
    E3=E1;
    nu12=0.29;
    nu13=nu12;
    nu23=nu12;
    G12=E1/2.0/(1.0 + nu12);
    G13=G12;
    G23=G12;
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

    # the symmetry plane
    l1 =selectnode(fens; box=[0 rex 0 0], inflate = tolerance)
    setebc!(u,l1,true, 2, 0.0)
    # The other end
    l1 =selectnode(fens; box=[0 rex Length Length], inflate = tolerance)
    setebc!(u,l1,true, 2, ua)

    applyebc!(u)
    numberdofs!(u)
    # println("Number of degrees of freedom = $(u.nfreedofs)")
    @test u.nfreedofs == 240

    # Property and material
    material=MatDeforElastOrtho(MR, E1,E2,E3,nu12,nu13,nu23,G12,G13,G23)
    # display(material)
    # println("$(material.D)")

    femm = FEMMDeforLinear(MR, IntegData(fes, GaussRule(2, 2), true), material)

    K =stiffness(femm, geom, u)
    F = nzebcloadsstiffness(femm, geom, u)
    U=  K\(F)
    scattersysvec!(u,U[:])

    fld= fieldfromintegpoints(femm, geom, u, :princCauchy, 1)
    # println("Minimum/maximum = $(minimum(fld.values))/$(maximum(fld.values))")
    @test abs(minimum(fld.values) - 0.0) < 1.0e-5
    @test abs(maximum(fld.values) - 0.0) < 1.0e-5
    fld= fieldfromintegpoints(femm, geom, u, :princCauchy, 2)
    # println("Minimum/maximum = $(minimum(fld.values))/$(maximum(fld.values))")
    @test abs(minimum(fld.values) - 0.0) < 1.0e-5
    @test abs(maximum(fld.values) - 0.0) < 1.0e-7
    fld= fieldfromintegpoints(femm, geom, u, :princCauchy, 3)
    # println("Minimum/maximum = $(minimum(fld.values))/$(maximum(fld.values))")
    @test abs(minimum(fld.values) - -0.050) < 1.0e-5
    @test abs(maximum(fld.values) - -0.04999999999999919) < 1.0e-5

    # File =  "mocylpull13.vtk"
    # vtkexportmesh(File, fens, fes; scalars=[("sigmaz", fld.values)],
    #               vectors=[("u", u.values)])
    # @async run(`"paraview.exe" $File`)
end
end
using .mocylpull13
mocylpull13.test()