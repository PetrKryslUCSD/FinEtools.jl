module LE11NAFEMS_examples
using FinEtools
using Statistics

function LE11NAFEMS_H20()
    ## Solid cylinder/taper/sphere—-temperature loading; quadratic brick mesh
    #
    
    ## Description
    #
    # The solid cylinder/taper/sphere axially-symmetric part represented in
    # Figure 1 is exposed to linearly varying temperature in the plane of the
    # cross-section. The temperature in the coordinates $r$  (the coordinate)
    # and $z$ (the axial ccoordinate)  is given as $T=r+z$. The goal is to find
    # the mechanical stress at the point A induced by the thermal expansion.
    #
    
    ##
    # The part is constrained against axial expansion along the faces of HIH'I'
    # and ABA'B'. The Young's modulus is 210 GPa, the Poisson's ratio is .3,
    # and the coefficient of thermal expansion is 2.3e-4/degree Celsius.
    
    ##
    # This is a test recommended by the National Agency for Finite Element
    # Methods and Standards (U.K.): Test LE11 from NAFEMS Publication TNSB,
    # Rev. 3, “The Standard NAFEMS Benchmarks,” October 1990.
    #
    # Target solution: Compressive  axial stress $\sigma_z$  = –105 MPa along
    # the circle passing through point A.
    
    ##
    # The toolkit has a helpful physical-units facility.  The function phun()
    # allows use of basic  units and basic
    # multipliers (for instance, mega).
    
    ##
    # Set the material properties.
    Ea = 210000*phun("MEGA*PA");# Young's modulus
    nua = 0.3;# Poisson ratio
    alphaa = 2.3e-4;# coefficient of thermal expansion
    
    ##
    # This is the target stress value.
    sigmaA = -105*phun("MEGA*PA");
    
    ##
    # The mesh  will be created in a very coarse representation from the
    # key points in the drawing. The first coordinate is radial, the second coordinate is axial.
    rz=[1.     0.;#A
    1.4    0.;#B
    0.995184726672197   0.098017140329561;
    1.393258617341076 0.137223996461385;
    0.980785  0.195090;#
    1.37309939 0.27312645;
    0.956940335732209   0.290284677254462
    1.339716470025092 0.406398548156247
    0.9238795  0.38268;#C
    1.2124  0.7;#D
    0.7071  0.7071;#E
    1.1062  1.045;#F
    0.7071  (0.7071+1.79)/2;#(E+H)/2
    1.      1.39;#G
    0.7071  1.79;#H
    1.      1.79;#I
    ]*phun("M")
    tolerance =1.e-6*phun("M")
    
    ##
    # Note that the material object needs to be created with the proper
    # model-dimension reduction in mind.  In this case that is the fully three-dimensional solid.
    MR = DeforModelRed3D
    
    # This is the quadrilateral mesh of the cross-section.   It will be modified and
    # refined as  we go.
    fens = FENodeSet(rz);
    fes = FESetQ4([1 2 4 3; 3 4 6 5; 5 6 8 7; 7 8 10 9; 9 10 12 11; 11 12 14 13; 13 14 16 15]);
    
    ##
    # If needed, the initial mesh  can be refined by bisection.  Just set
    # `nref` greater than zero.  Note that  the nodes located along the
    # edges are moved onto the  spherical surface when they _should be_ on
    # the spherical surface.  This is important in order to ensure
    # convergence to the proper value of the stress.  Just refining  the
    # initial mesh without repositioning of the nodes onto the spherical surface would mean that the
    # refinement would preserve a concave corner where in reality there is
    # none.  The stress would be artificially raised and convergence would
    # not be guaranteed.
    
    nref = 0;
    for ref = 1:nref
        fens,fes = Q4refine(fens,fes);
        list = selectnode(fens, distance=1.0+0.1/2^nref, from=[0. 0.], inflate=tolerance);
        fens.xyz[list,:] = FinEtools.MeshUtilModule.ontosphere(fens.xyz[list,:],1.0);
    end
    
    ##
    # The mesh is extruded by sweeping around the axis of symmetry.
    # Only a single layer of elements is generated of internal angle
    # |angslice|.
    nLayers = 7;
    angslice  = 5*pi/16;
    
    ##
    # First the mesh is extruded to a block whose third dimension
    # represents the angular coordinate.
    fens,fes = H8extrudeQ4(fens, fes, nLayers,
    (rz,k)->[rz[1],rz[2],0.0]-(k)/nLayers*[0.,0.,angslice]);
    
    ##
    # The mesh is now converted to the serendipity 20-node elements.
    # We will reposition the nodes later.
    fens,fes = H8toH20(fens,fes);
    
    ##
    # The boundary of the block is extracted and the faces of the mesh on
    # the bounding cross-sections are identified. Recall that this is just
    # about the topology (connectivity), the geometry does not matter at
    # this point.
    bfes = meshboundary(fes);
    f1l = selectelem(fens, bfes, box=[-Inf,Inf,-Inf,Inf,0.0,0.0], inflate=tolerance);
    f2l = selectelem(fens, bfes, box=[-Inf,Inf,-Inf,Inf,-angslice,-angslice],
    inflate=tolerance);
    
    ##
    # The block is now converted  to the axially symmetric geometry by using the
    # third (angular) coordinate  to sweep out  an axially symmetric domain. The
    # ccoordinates of the nodes at this point are |rza|,  radial distance,
    # Z-coordinate, angle.
    sweep(rza) = [-rza[1]*sin(rza[3]+angslice/2.0), rza[1]*cos(rza[3]+angslice/2.0), rza[2]]
    for j=1:size(fens.xyz,1)
        fens.xyz[j,:] = sweep(fens.xyz[j,:])
    end
    
    
    ##
    # The nodes within the radial distance of 1.0 of the origin (i. e.
    # those on the spherical surface)  are repositioned one more time to be
    # located on the spherical surface for sure. (Recall  that we have
    # inserted additional nodes at the midpoints of the edges when the mesh
    # was converted to quadratic elements.)
    list = selectnode(fens,distance=1.0+0.1/2^nref,
    from=[0. 0. 0.], inflate=tolerance);
    fens.xyz[list,:]= FinEtools.MeshUtilModule.ontosphere(fens.xyz[list,:], 1.0);
    
    ##
    # We are ready to create the  finite element model machine and to use
    # it to construct  the global system for the displacements.
    ##
    # The material is created from the property object.  Note that the
    # |alpha| attribute is the thermal expansion coefficient.
    
    # Create isotropic elastic material
    material = MatDeforElastIso(MR, 1.0, Ea, nua, alphaa)
    
    ##
    # The finite element  model machine puts together the material, the
    # finite elements,  and the integration rule. The Gauss quadrature with
    # 3x3x3 points  gives good accuracy in this case. Compare it with 2x2x2
    # quadrature to appreciate the difference.
    
    femm = FEMMDeforLinear(MR, IntegData(fes, GaussRule(3, 3)), material)
    
    ##
    # The geometry nodal field is created from the node set.   The
    # displacement field is created by cloning the geometry and then
    # zeroing out the nodal parameters.
    geom = NodalField(fens.xyz)
    u = NodalField(zeros(size(fens.xyz,1),3)) # displacement field
    nnodes(geom)
    ##
    # The EBCs are applied  next.  Only the axial (Z) degrees of freedom at
    # the bottom and top are fixed to zero.
    l1 = selectnode(fens, box=[-Inf Inf -Inf Inf 0.0 0.0], inflate=tolerance)
    setebc!(u, l1, true, 3, zeros(size(l1)))
    l1 = selectnode(fens, box=[-Inf Inf -Inf Inf 1.79  1.79], inflate=tolerance)
    setebc!(u, l1, true, 3, zeros(size(l1)))
    applyebc!(u)
    numberdofs!(u)
    
    
    ##
    # The restraints of the nodes on the bounding cross-sections in the direction
    # of the normal to the plane of the cross-section  in the
    # circumferential direction are introduced using a penalty formulation.
    # For that purpose we introduce  a finite element model machine for the
    # surface  finite elements on the cross-sections.
    springcoefficient =1.0 / ((abs(sigmaA)/1.0e12)/Ea)
    fl = vcat(f1l, f2l)
    xsfemm = FEMMDeforWinkler(IntegData(subset(bfes,fl), GaussRule(2, 3)))
    
    ##
    # We create the temperature field using the formula $T=r+z$.
    dT = NodalField(reshape(sqrt.(fens.xyz[:,1].^2+fens.xyz[:,2].^2)+fens.xyz[:,3],size(fens.xyz,1),1));
    
    ##
    # And we are ready to assemble the system matrix. Both the elastic stiffness of
    # the hexahedral elements ...
    K = stiffness(femm, geom, u)
    # ...  and the elastic stiffness    of the springs on the contact surfaces of the cross-sections.
    H = surfacenormalspringstiffness(xsfemm,  geom, u, springcoefficient, SurfaceNormal(3))
    
    ##
    # The mechanical loads are computed from the thermal strains.
    F = thermalstrainloads(femm, geom, u, dT)
    
    ##
    # And  the solution for the free degrees of freedom is obtained.
    U=  (K+H)\F
    scattersysvec!(u, U[:])
    
    
    ##
    # The stress  is recovered from the stress calculated at the
    # integration points.
    
    fld= fieldfromintegpoints(femm, geom, u, dT, :Cauchy, 3)
    
    
    ##
    # Now that we have the nodal field  for the axial stress, we can plot
    # the axial stress painted on the deformed geometry.
    
    
    File =  "LE11NAFEMS_H20_sigmaz.vtk"
    vtkexportmesh(File, fens, fes;
    scalars=[("sigmaz", fld.values)], vectors=[("u", u.values)])
    @async run(`"paraview.exe" $File`)
    # File =  "LE11NAFEMS_H20_dT.vtk"
    # vtkexportmesh(File, fens, fes; scalars=dT.values,scalars_name ="dT", vectors=u.values,vectors_name="u")
    
    ##
    # The  computed stress at the node that is located at the point A  is
    # going to be now extracted from the nodal field for the stress.
    # Nodes at level Z=0.0
    l1 =selectnode(fens,box=FFlt[-Inf  Inf -Inf  Inf 0.0 0.0],inflate=tolerance);
    l2 =selectnode(fens,distance=1.0+0.1/2^nref,from=FFlt[0.0 0.0 0.0],inflate=tolerance);
    nA=intersect(l1,l2);
    sA = mean(fld.values[nA])/phun("MEGA*Pa")
    sAn = mean(fld.values[nA])/sigmaA
    println("Stress at point A: $(sA) i. e.  $( sAn*100  )% of reference value")
    
    
    ## Discussion
    #
    ##
    # The 3-D solution corresponds well to the 2-D axially symmetric model.
    # We also see good correspondence to other published solutions for
    # comparable finite element models.  For instance, Abaqus 6.11
    # Benchmark manual lists very similar numbers.
    
    
end # LE11NAFEMS_H20

function allrun()
    println("#####################################################") 
    println("# LE11NAFEMS_H20 ")
    LE11NAFEMS_H20()
    return true
end # function allrun

end # module LE11NAFEMS_examples
