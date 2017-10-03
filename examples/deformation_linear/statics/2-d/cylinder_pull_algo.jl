module mocylpullFun
using FinEtools
using FinEtools.AlgoDeforLinearModule
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
    fens.xyz[:, 1] += rin
    bdryfes = meshboundary(fes);

    # the symmetry plane
    la1 =selectnode(fens; box=[0 rex 0 0], inflate = tolerance)
    # The other end
    la2 =selectnode(fens; box=[0 rex Length Length], inflate = tolerance)

    e1 = FDataDict("node_list"=>la1, "component"=>2, "displacement"=>x -> 0.0)
    e2 = FDataDict("node_list"=>la2, "component"=>2, "displacement"=>x -> ua)

    # Property and material
    material=MatDeforElastOrtho(MR, E1,E2,E3,nu12,nu13,nu23,G12,G13,G23)

    femm = FEMMDeforLinear(MR, GeoD(fes, GaussRule(2, 2), true), material)

    # Make region
    region = FDataDict("femm"=>femm)

    # Make model data
    modeldata =  FDataDict(
    "fens"=> fens, "regions"=>  [region],
    "essential_bcs"=>[e1, e2])

    # Call the solver
    modeldata = AlgoDeforLinearModule.linearstatics(modeldata)
    geom = modeldata["geom"]
    u = modeldata["u"]

    # Produce a plot of the radial stress component in the cylindrical
    # coordinate system. Note that this is the usual representation of
    # stress using nodal stress field.

    fld= fieldfromintegpoints(femm, geom, u, :Cauchy, 2)
    println("Minimum/maximum = $(minimum(fld.values))/$(maximum(fld.values))")

    File =  "orthoballoon_sigmaz.vtk"
    vtkexportmesh(File, fens, fes; scalars=[("sigmaz", fld.values)],
                  vectors=[("u", u.values)])
    @async run(`"paraview.exe" $File`)
end
end
using mocylpullFun
mocylpullFun.test()
