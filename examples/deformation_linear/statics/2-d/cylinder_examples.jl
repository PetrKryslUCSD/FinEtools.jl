module cylinder_examples
using FinEtools
using FinEtools.AlgoDeforLinearModule


function cylinder_bend()
    # Cylinder  bent by edge by enforced displacement, axially symmetric model
    println("Cylinder  bent by edge by enforced displacement, axially symmetric model")
    
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
    fens.xyz[:, 1] .+= rin
    bdryfes = meshboundary(fes);
    
    # now we create the geometry and displacement fields
    geom = NodalField(fens.xyz)
    u = NodalField(zeros(size(fens.xyz,1),2)) # displacement field
    
    # the symmetry plane
    l1 =selectnode(fens; box=[0 rex 0 0], inflate = tolerance)
    setebc!(u,l1,true, 2, 0.0)
    # The other end
    l1 =selectnode(fens; box=[0 rex Length Length], inflate = tolerance)
    for ixxxx = 1:length(l1)
        r = fens.xyz[l1[ixxxx], 1]
        setebc!(u,[l1[ixxxx]],true, 2, (r-(rex+rin)/2)/((rex+rin)/2)*ua)
    end
    
    
    applyebc!(u)
    numberdofs!(u)
    println("Number of degrees of freedom = $(u.nfreedofs)")
    
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
    println("Minimum/maximum = $(minimum(fld.values))/$(maximum(fld.values))")
    
    File =  "cylinder_bend_sigmaz.vtk"
    vtkexportmesh(File, fens, fes; scalars=[("sigmaz", fld.values)],
    vectors=[("u", u.values)])
    @async run(`"paraview.exe" $File`)
    
end # cylinder_bend


function cylinder_pull()
    # Cylinder  pulled by enforced displacement, axially symmetric model
    println("Cylinder  pulled by enforced displacement, axially symmetric model")
    
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
    fens.xyz[:, 1] .+= rin
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
    println("Number of degrees of freedom = $(u.nfreedofs)")
    
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
    println("Minimum/maximum = $(minimum(fld.values))/$(maximum(fld.values))")
    
    File =  "orthoballoon_sigmaz.vtk"
    vtkexportmesh(File, fens, fes; scalars=[("sigmaz", fld.values)],
    vectors=[("u", u.values)])
    @async run(`"paraview.exe" $File`)
    
end # cylinder_pull


function cylinder_pull_algo()
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
    fens.xyz[:, 1] .+= rin
    bdryfes = meshboundary(fes);
    
    # the symmetry plane
    la1 =selectnode(fens; box=[0 rex 0 0], inflate = tolerance)
    # The other end
    la2 =selectnode(fens; box=[0 rex Length Length], inflate = tolerance)
    
    e1 = FDataDict("node_list"=>la1, "component"=>2, "displacement"=>x -> 0.0)
    e2 = FDataDict("node_list"=>la2, "component"=>2, "displacement"=>x -> ua)
    
    # Property and material
    material=MatDeforElastOrtho(MR, E1,E2,E3,nu12,nu13,nu23,G12,G13,G23)
    
    femm = FEMMDeforLinear(MR, IntegData(fes, GaussRule(2, 2), true), material)
    
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
end # cylinder_pull_algo

function allrun()
    println("#####################################################") 
    println("# cylinder_bend ")
    cylinder_bend()
    println("#####################################################") 
    println("# cylinder_pull ")
    cylinder_pull()
    println("#####################################################") 
    println("# cylinder_pull_algo ")
    cylinder_pull_algo()
    return true
end # function allrun

end # module cylinder_examples
