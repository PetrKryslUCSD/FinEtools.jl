module orthoballoon_examples
using FinEtools
using LinearAlgebra

function orthoballoon()
    # Orthotropic balloon inflation, axially symmetric model
    
    
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
    
    fens,fes = Q4block(rex-rin,pi/2,5,20);
    bdryfes = meshboundary(fes);
    
    icl = selectelem(fens, bdryfes, box=[0.,0.,0.,pi/2],inflate=tolerance);
    for i=1:count(fens)
        r=rin+fens.xyz[i,1]; a=fens.xyz[i,2];
        fens.xyz[i,:]=[r*cos(a) r*sin(a)];
    end
    
    
    # now we create the geometry and displacement fields
    geom = NodalField(fens.xyz)
    u = NodalField(zeros(size(fens.xyz,1),2)) # displacement field
    
    # the symmetry plane
    l1 =selectnode(fens; box=[0 rex 0 0], inflate = tolerance)
    setebc!(u,l1,true, 2, 0.0)
    # the axis of symmetry
    l1 =selectnode(fens; box=[0 0 0 rex], inflate = tolerance)
    setebc!(u,l1,true, 1, 0.0)
    
    applyebc!(u)
    numberdofs!(u)
    println("Number of degrees of freedom = $(u.nfreedofs)")
    
    # The traction boundary condition is applied in the radial
    # direction.
    
    el1femm =  FEMMBase(IntegData(subset(bdryfes,icl), GaussRule(1, 3), true))
    function pressureloading!(forceout::FFltVec, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)
        copyto!(forceout, XYZ/norm(XYZ)*p)
        return forceout
    end
    fi = ForceIntensity(FFlt, 2, pressureloading!); # pressure normal to the internal cylindrical surface
    F2= distribloads(el1femm, geom, u, fi, 2);
    
    # Property and material
    material=MatDeforElastOrtho(MR, E1,E2,E3,nu12,nu13,nu23,G12,G13,G23)
    
    femm = FEMMDeforLinear(MR, IntegData(fes, GaussRule(2, 2), true), material)
    
    K =stiffness(femm, geom, u)
    U=  K\(F2)
    scattersysvec!(u,U[:])
    
    # Produce a plot of the radial stress component in the cylindrical
    # coordinate system. Note that this is the usual representation of
    # stress using nodal stress field.
    
    fld= fieldfromintegpoints(femm, geom, u, :Cauchy, 3)
    println("Minimum/maximum = $(minimum(fld.values))/$(maximum(fld.values))")
    
    File =  "orthoballoon_sigmaz.vtk"
    vtkexportmesh(File, fens, fes; scalars=[("sigmaz", fld.values)],
    vectors=[("u", u.values)])
    @async run(`"paraview.exe" $File`)
#pub_thick_pipe_axi()

end # orthoballoon


function orthoballoon_penalty()
    # Orthotropic balloon inflation, axially symmetric model
    
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
    tolerance=rin/1000.
    
    MR = DeforModelRed2DAxisymm
    
    fens,fes = Q4block(rex-rin,pi/2,5,20);
    bdryfes = meshboundary(fes);
    
    icl = selectelem(fens, bdryfes, box=[0.,0.,0.,pi/2], inflate=tolerance);
    for i=1:count(fens)
        r=rin+fens.xyz[i,1]; a=fens.xyz[i,2];
        fens.xyz[i,:]=[r*cos(a) r*sin(a)];
    end
    
    # now we create the geometry and displacement fields
    geom = NodalField(fens.xyz)
    u = NodalField(zeros(size(fens.xyz,1),2)) # displacement field
    
    # the symmetry plane
    ly =selectelem(fens, bdryfes; box=[0 rex 0 0], inflate = tolerance)
    # the axis of symmetry
    lx =selectelem(fens, bdryfes; box=[0 0 0 rex], inflate = tolerance)
    
    # No EBC
    applyebc!(u)
    numberdofs!(u)
    println("Number of degrees of freedom = $(u.nfreedofs)")
    
    # The traction boundary condition is applied in the radial
    # direction.
    
    el1femm =  FEMMBase(IntegData(subset(bdryfes,icl), GaussRule(1, 3), true))
    function pressureloading!(forceout::FFltVec, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)
        copyto!(forceout, XYZ/norm(XYZ)*p)
        return forceout
    end
    fi = ForceIntensity(FFlt, 2, pressureloading!); # pressure normal to the internal cylindrical surface
    F2= distribloads(el1femm, geom, u, fi, 2);
    
    # Property and material
    material=MatDeforElastOrtho(MR, E1,E2,E3,nu12,nu13,nu23,G12,G13,G23)
    
    femm = FEMMDeforLinear(MR, IntegData(fes, GaussRule(2, 2), true), material)
    
    ##
    # The restraints of the nodes on the bounding cross-sections in the direction
    # of the normal to the plane of the cross-section  in the
    # circumferential direction are introduced using a penalty formulation.
    # For that purpose we introduce  a finite element model machine for the
    # surface  finite elements on the cross-sections.
    springcoefficient =1.0e9/ (abs(p)/E1)
    xsfemm = FEMMDeforWinkler(IntegData(subset(bdryfes, lx), GaussRule(1, 3), true))
    ysfemm = FEMMDeforWinkler(IntegData(subset(bdryfes, ly), GaussRule(1, 3), true))
    H = surfacenormalspringstiffness(xsfemm,  geom, u, springcoefficient, SurfaceNormal(3)) +
    surfacenormalspringstiffness(ysfemm,  geom, u, springcoefficient, SurfaceNormal(3))
    K =stiffness(femm, geom, u)
    U=  (K + H)\(F2)
    scattersysvec!(u,U[:])
    
    # Produce a plot of the radial stress component in the cylindrical
    # coordinate system. Note that this is the usual representation of
    # stress using nodal stress field.
    
    fld= fieldfromintegpoints(femm, geom, u, :Cauchy, 3)
    println("Minimum/maximum = $(minimum(fld.values))/$(maximum(fld.values))")
    
    File =  "orthoballoon_penalty_sigmaz.vtk"
    vtkexportmesh(File, fens, fes; scalars=[("sigmaz", fld.values)],
    vectors=[("u", u.values)])
    @async run(`"paraview.exe" $File`)
#pub_thick_pipe_axi()

end # orthoballoon_penalty

function allrun()
    println("#####################################################") 
    println("# orthoballoon ")
    orthoballoon()
    println("#####################################################") 
    println("# orthoballoon_penalty ")
    orthoballoon_penalty()
    return true
end # function allrun

end # module orthoballoon_examples
