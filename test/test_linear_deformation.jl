
module mmLE11NAFEMSQ8algo2
using FinEtools
using FinEtools.AlgoDeforLinearModule: linearstatics, exportdeformation,
exportstress, exportstresselementwise
using Test
import LinearAlgebra: norm, cholesky, cross
function test()
    # NAFEMS LE11 benchmark with Q8 elements.
    # # This is a test recommended by the National Agency for Finite Element
    # # Methods and Standards (U.K.): Test LE11 from NAFEMS Publication TNSB,
    # # Rev. 3, “The Standard NAFEMS Benchmarks,” October 1990.
    # #
    # # Target solution: Direct stress,   =  –105 MPa at point A.
    #function  LE11NAFEMS()
    # Parameters:
    Ea =  210000*phun("MEGA*Pa")
    nua =  0.3;
    alphaa = 2.3e-4;              # thermal expansion coefficient
    sigmaA = -105*phun("MEGA*Pa")
    nref =  1;                        # how many times should we refine the mesh?
    X = [1.     0.;#A
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
    tolerance  = 1.e-6*phun("M")
    ##
    # Note that the material object needs to be created with the proper
    # model-dimension reduction in mind.  In this case that is the axial symmetry
    # assumption.
    MR  =  DeforModelRed2DAxisymm

    fens = FENodeSet(X);
    fes = FESetQ4([1 2 4 3; 3 4 6 5; 5 6 8 7; 7 8 10 9; 9 10 12 11; 11 12 14 13; 13 14 16 15]);
    for ref = 1:nref
        fens,fes = Q4refine(fens,fes);
        list = selectnode(fens,distance = 1.0+0.1/2^nref, from = [0. 0.], inflate = tolerance);
        fens.xyz[list,:] =  FinEtools.MeshUtilModule.ontosphere(fens.xyz[list,:],1.0);
    end
    fens,fes = Q4toQ8(fens,fes)
    list = selectnode(fens,distance = 1.0+0.1/2^nref, from = [0. 0.], inflate = tolerance);
    fens.xyz[list,:] =  FinEtools.MeshUtilModule.ontosphere(fens.xyz[list,:],1.0);


    # EBC's
    l1 = selectnode(fens,box = [-Inf Inf 0 0],inflate = tolerance)
    e1 = FDataDict("node_list"=>l1, "component"=>2, "displacement"=>x -> 0.0)
    l1 = selectnode(fens,box = [-Inf Inf 1.79  1.79],inflate = tolerance)
    e2 = FDataDict("node_list"=>l1, "component"=>2, "displacement"=>x -> 0.0)

    # Temperature field
    dtemp = FDataDict("temperature"=>x -> x[1] + x[2])

    # Property and material
    material = MatDeforElastIso(MR, 0.0, Ea, nua, alphaa)

    femm  =  FEMMDeforLinear(MR, IntegDomain(fes, GaussRule(2, 3), true), material)

    # Make region 1
    region = FDataDict("femm"=>femm)
    # Make model data
    modeldata =  FDataDict(
    "fens"=> fens, "regions"=>  [region],
    "essential_bcs"=>[e1, e2], "temperature_change"=>dtemp)

    # Call the solver
    modeldata = linearstatics(modeldata)
    geom = modeldata["geom"]
    u = modeldata["u"]
    dT = modeldata["temp"]

    modeldata["postprocessing"] = FDataDict("boundary_only"=> true,
    "file"=>"LE11NAFEMS_Q8_deformation.vtk")
    modeldata = exportdeformation(modeldata)
    # @async run(`"paraview.exe" $(modeldata["postprocessing"]["exported"][1]["file"])`)
    try rm(modeldata["postprocessing"]["exported"][1]["file"]) catch end

    nA  = selectnode(fens,box = FFlt[1.0  1.0 0.0 0.0], inflate = tolerance);

    modeldata["postprocessing"] = FDataDict("boundary_only"=> true,
    "file"=>"LE11NAFEMS_Q8_sigmay.vtk", "quantity"=>:Cauchy,
    "component"=>2)
    modeldata = exportstress(modeldata)
    modeldata["postprocessing"] = FDataDict("boundary_only"=> false,
    "file"=>"LE11NAFEMS_Q8_sigmay.vtk", "quantity"=>:Cauchy,
    "component"=>2)
    modeldata = exportstress(modeldata)
    # @async run(`"paraview.exe" $(modeldata["postprocessing"]["exported"][1]["file"])`)
    try rm(modeldata["postprocessing"]["exported"][1]["file"]) catch end
    fld =  modeldata["postprocessing"]["exported"][1]["field"]

    sA  =  fld.values[nA]/phun("MEGA*Pa")
    sAn  =  fld.values[nA]/sigmaA
    # println("Stress at point A: $(sA) i. e.  $( sAn*100  )% of reference value")
    @test norm(sA - [-93.8569]) < 1.0e-2

    # Loop over only those elements that share the node nA
    fen2fe  = FENodeToFEMap(connasarray(fes), nnodes(geom))
    function inspector(idat, elnum, conn, xe,  out,  xq)
        # println("loc = $(  xq  ) : $(  transpose(out)/phun("MEGA*Pa")  )")
        return idat
    end

    inspectintegpoints(femm, geom, u, dT,  fen2fe.map[nA[1]],
    inspector, []; quantity = :Cauchy)

    modeldata["postprocessing"] = FDataDict("boundary_only"=> false,
    "file"=>"LE11NAFEMS_Q8_sigmay_ew.vtk", "quantity"=>:Cauchy,
    "component"=>2)
    modeldata = exportstresselementwise(modeldata)
    # @async run(`"paraview.exe" $(modeldata["postprocessing"]["exported"][1]["file"])`)
    try rm(modeldata["postprocessing"]["exported"][1]["file"]) catch end
end
end
using .mmLE11NAFEMSQ8algo2
mmLE11NAFEMSQ8algo2.test()


module sscratch_06112017
using FinEtools
using Test
import Statistics: mean
function test()
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

  femm = FEMMDeforLinear(MR, IntegDomain(fes, GaussRule(3, 3)), material)

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
  xsfemm = FEMMDeforWinkler(IntegDomain(subset(bfes,fl), GaussRule(2, 3)))

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


  # File =  "LE11NAFEMS_H20_sigmaz.vtk"
  # vtkexportmesh(File, fens, fes;
  # scalars=[("sigmaz", fld.values)], vectors=[("u", u.values)])
  # @async run(`"paraview.exe" $File`)
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
  # println("Stress at point A: $(sA) i. e.  $( sAn*100  )% of reference value")

  @test abs(sA-(-83.7322285847101))/(-83.7322285847101) < 1.e-3
  ## Discussion
  #
  ##
  # The 3-D solution corresponds well to the 2-D axially symmetric model.
  # We also see good correspondence to other published solutions for
  # comparable finite element models.  For instance, Abaqus 6.11
  # Benchmark manual lists very similar numbers.
end
end
using .sscratch_06112017
sscratch_06112017.test()

module cookstress_1
using Test
using FinEtools
using FinEtools.MeshExportModule
import LinearAlgebra: norm, cholesky, cross
function test()
  # println("Cook membrane problem,  plane stress."        )
  t0 = time()

  E = 1.0;
  nu = 1.0/3;
  width = 48.0; height = 44.0; thickness  = 1.0;
  free_height  = 16.0;
  Mid_edge  = [48.0, 52.0];# Location of tracked  deflection
  magn = 1.0/free_height;# Magnitude of applied load
  convutip = 23.97;
  n = 32;#*int(round(sqrt(170.)/2.)); # number of elements per side
  tolerance = minimum([width, height])/n/1000.;#Geometrical tolerance

  fens, fes = T3block(width, height,  n,  n)

  # Reshape into a trapezoidal panel
  for i=1:count(fens)
    fens.xyz[i, 2]=fens.xyz[i, 2]+(fens.xyz[i, 1]/width)*(height -fens.xyz[i, 2]/height*(height-free_height));
  end

  geom = NodalField(fens.xyz)
  u = NodalField(zeros(size(fens.xyz, 1), 2)) # displacement field

  l1 = selectnode(fens; box=[0, 0, -Inf,  Inf],  inflate = tolerance)
  setebc!(u, l1, 1, val=0.0)
  setebc!(u, l1, 2, val=0.0)
  applyebc!(u)
  numberdofs!(u)

  boundaryfes =  meshboundary(fes);
  Toplist = selectelem(fens, boundaryfes,  box= [width,  width,  -Inf,  Inf ],  inflate=  tolerance);
  el1femm =  FEMMBase(IntegDomain(subset(boundaryfes, Toplist),  GaussRule(1, 2)))
  fi = ForceIntensity([0.0, +magn]);
  F2 = distribloads(el1femm,  geom,  u,  fi,  2);


  MR = DeforModelRed2DStress
  material = MatDeforElastIso(MR,  0.0, E, nu, 0.0)

  femm = FEMMDeforLinear(MR, IntegDomain(fes,  TriRule(1)),  material)

  K = stiffness(femm,  geom,  u)
  K = cholesky(K)
  U=  K\(F2)
  scattersysvec!(u, U[:])

  nl = selectnode(fens,  box=[Mid_edge[1], Mid_edge[1], Mid_edge[2], Mid_edge[2]], inflate=tolerance);
  theutip = zeros(FFlt, 1, 2)
  gathervalues_asmat!(u, theutip, nl);
  # println("$(time()-t0) [s];  displacement =$(theutip[2]) as compared to converged $convutip")

  # File =  "a.vtk"
  # vtkexportmesh(File,  fes.conn,  geom.values+u.values,
  # FinEtools.MeshExportModule.T3; vectors=[("u", u.values)])

  @test abs(theutip[2]-23.8155)/23.8155 < 1.e-3 # FinEALE solution
end
end
using .cookstress_1
cookstress_1.test()


module scratch1_06092017

using FinEtools
using Test

mutable struct MyIData
  c::FInt
  r::FFltVec
  s::FFltVec
end

function test()

  # println("Thick pipe with internal pressure: axially symmetric model")
  #=
  This is a simple modification of the full three-dimensional simulation of
  the tutorial pub_thick_pipe that implements the axially-symmetric model
  reduction procedure.

  An infinitely long thick walled cylindrical pipe
  with inner boundary radius of 3 mm and outer boundary radius of 9 mm is
  subjected to an internal pressure of 1.0 MPa. A wedge   with thickness of
  2 mm and a 90-degree angle sector is considered for the finite element
  analysis. The material properties are taken as  isotropic linear elastic
  with $E=1000$ MPa and $\nu=0.4999$ to represent nearly incompressible
  behavior. This problem has been proposed to by MacNeal and Harder as a
  test of an element's ability to represent the  response of a nearly
  incompressible material. The plane-strain condition is assumed in the
  axial direction of the pipe which together with the radial symmetry
  confines the material in all but the radial direction and therefore
  amplifies the numerical difficulties associated with the confinement of
  the nearly incompressible material.

  There is an analytical solution to this problem. Timoshenko and Goodier
  presented the original solution of Lame in their textbook. We are going
  to compare with  both the stress distribution (radial and hoop stresses)
  and the displacement of the inner  cylindrical surface.

  References:
  - Macneal RH, Harder RL (1985) A proposed standard set of problems to test
  finite element accuracy. Finite Elements in Analysis and Design 1: 3-20.
  - Timoshenko S. and Goodier J. N., Theory of Elasticity, McGraw-Hill, 2nd ed., 1951.

  =#

  # Internal radius of the pipe.
  a=3*phun("MM");
  ##
  # External radius of the pipe.
  b = 9*phun("MM");
  ##
  # Thickness of the slice.
  t = 2*phun("MM");

  ##
  # Geometrical tolerance.
  tolerance   = a/10000.;
  ##
  # Young's modulus and Poisson's ratio.
  E = 1000*phun("MEGA*PA");
  nu = 0.499;
  ##
  # Applied pressure on the internal surface.
  press =   1.0*phun("MEGA*PA");

  ##
  # Analytical solutions.   Radial stress:
  radial_stress(r)  = press*a.^2/(b^2-a^2).*(1-(b^2)./r.^2);
  ##
  # Circumferential (hoop) stress:
  hoop_stress(r) = press*a.^2/(b^2-a^2).*(1+(b^2)./r.^2);

  ##
  # Radial displacement:
  radial_displacement(r) = press*a^2*(1+nu)*(b^2+r.^2*(1-2*nu))/(E*(b^2-a^2).*r);;

  ##
  # Therefore the radial displacement of the loaded surface will be:
  urex  =  radial_displacement(a);


  ##
  # The mesh parameters: The numbers of element edges axially,
  # and through the thickness of the pipe wall (radially).

  na = 1; nt = 10;

  ##
  # Note that the material object needs to be created with the proper
  # model-dimension reduction in effect.  In this case that is the axial symmetry
  # assumption.
  MR = DeforModelRed2DAxisymm
  axisymmetric = true

  # Create the mesh and initialize the geometry.  First we are going
  # to construct the block of elements with the first coordinate
  # corresponding to the thickness in the radial direction, and the second
  # coordinate is the thickness in the axial direction.
  fens,fes  =   Q8block(b-a, t, nt, na);

  # Extract the boundary  and mark the finite elements on the
  # interior surface.
  bdryfes = meshboundary(fes);

  bcl = selectelem(fens, bdryfes, box=[0.,0.,-Inf,Inf], inflate=tolerance);
  internal_fenids= connectednodes(subset(bdryfes,bcl));
  # Now  shape the block  into  the actual wedge piece of the pipe.
  for i=1:count(fens)
    fens.xyz[i,:] = fens.xyz[i,:] + [a; 0.0];
  end

  # now we create the geometry and displacement fields
  geom = NodalField(fens.xyz)
  u = NodalField(zeros(size(fens.xyz,1),2)) # displacement field

  # The plane-strain condition in the axial direction  is specified by selecting nodes
  # on the plane y=0 and y=t.
  l1 = selectnode(fens; box=[-Inf Inf 0.0 0.0], inflate = tolerance)
  setebc!(u, l1, true, 2, 0.0)
  l1 = selectnode(fens; box=[-Inf Inf t t], inflate = tolerance)
  setebc!(u, l1, true, 2, 0.0)

  applyebc!(u)
  numberdofs!(u)

  # The traction boundary condition is applied in the radial
  # direction.

  el1femm =  FEMMBase(IntegDomain(subset(bdryfes,bcl), GaussRule(1, 3), axisymmetric))
  fi = ForceIntensity([press; 0.0]);
  F2= distribloads(el1femm, geom, u, fi, 2);

  # Property and material
  material = MatDeforElastIso(MR,  E, nu)

  femm = FEMMDeforLinear(MR, IntegDomain(fes, GaussRule(2, 2), axisymmetric), material)

  K = stiffness(femm, geom, u)
  #K=cholesky(K)
  U =  K\(F2)
  scattersysvec!(u,U[:])

  # Transfer the solution of the displacement to the nodes on the
  # internal cylindrical surface and convert to
  # cylindrical-coordinate displacements there.
  uv=u.values[internal_fenids,:]
  # Report the  relative displacement on the internal surface:
  # println("(Approximate/true displacement) at the internal surface: $( mean(uv[:,1])/urex*100  ) %")

  # Produce a plot of the radial stress component in the cylindrical
  # coordinate system. Note that this is the usual representation of
  # stress using nodal stress field.

  fld = fieldfromintegpoints(femm, geom, u, :Cauchy, 1)


  # File =  "thick_pipe_sigmax.vtk"
  # vtkexportmesh(File, fens, fes; scalars=[("sigmax", fld.values)])

  # Produce a plot of the solution components in the cylindrical
  # coordinate system.

  function inspector(idat::MyIData, elnum, conn, xe,  out,  xq)
    push!(idat.r, xq[1])
    push!(idat.s, out[idat.c])
    return idat
  end

  idat = MyIData(1, FInt[], FInt[])
  idat = inspectintegpoints(femm, geom, u, collect(1:count(fes)),
  inspector, idat, :Cauchy)

  # using Plots
  # plotly()
  #
  # # Plot the analytical solution.
  # r = linearspace(a,b,100);
  # plot(r, radial_stress(r))
  # # Plot the computed  integration-point data
  # plot!(idat.r, idat.s, m=:circle, color=:red)
  # gui()

  @test  abs(idat.r[1]-0.003126794919243112)<1.0e-9
  @test  abs(idat.s[1]- -910911.9777008593)<1.0e-2

  ## Discussion
  #
  ##
  # The axially symmetric model is clearly very effective
  # computationally, as the size is much reduced compared to the 3-D
  # model.  In conjunction with  uniform or selective reduced integration
  # it can be very accurate as well.
end
end
using .scratch1_06092017
scratch1_06092017.test()

module scratch2_06102017

using FinEtools
using Test
import Arpack: eigs
function test()
  # println("""
  #         % Vibration modes of unit cube  of almost incompressible material.
  #         %
  #         % Reference: Puso MA, Solberg J (2006) A stabilized nodally integrated
  #         % tetrahedral. International Journal for Numerical Methods in
  #         % Engineering 67: 841-867.""")
  #         t0 = time()


  E = 1*phun("PA");
  nu = 0.499;
  rho = 1*phun("KG/M^3");
  a = 1*phun("M"); b = a; h =  a;
  n1 = 10;# How many element edges per side?
  na =  n1; nb =  n1; nh  = n1;
  neigvs = 20                   # how many eigenvalues
  OmegaShift = (0.01*2*pi)^2;

  MR = DeforModelRed3D
  fens,fes  = H20block(a,b,h, na,nb,nh)

  geom = NodalField(fens.xyz)
  u = NodalField(zeros(size(fens.xyz,1),3)) # displacement field

  numberdofs!(u)

  material=MatDeforElastIso(MR, rho, E, nu, 0.0)

  femm = FEMMDeforLinear(MR, IntegDomain(fes, GaussRule(3,2)), material)

  K =stiffness(femm, geom, u)
  femm = FEMMDeforLinear(MR, IntegDomain(fes, GaussRule(3,3)), material)
  M =mass(femm, geom, u)
  d,v,nev,nconv = eigs(K+OmegaShift*M, M; nev=neigvs, which=:SM)
  d = d .- OmegaShift;
  fs = real(sqrt.(complex(d)))/(2*pi)
  # println("Eigenvalues: $fs [Hz]")

  # mode = 17
  # scattersysvec!(u, v[:,mode])
  # File =  "unit_cube_modes.vtk"
  # vtkexportmesh(File, fens, fes; vectors=[("mode$mode", u.values)])

  @test abs(fs[7]-0.26259869196259) < 1.0e-5
end
end

using .scratch2_06102017
scratch2_06102017.test()

module mxxxx1_06102017

using FinEtools
using FinEtools.AlgoDeforLinearModule
import Statistics: mean
using Test

# println("""
# The initially twisted cantilever beam is one of the standard test
# problems for verifying the finite-element accuracy [1]. The beam is
# clamped at one end and loaded either with unit in-plane or
# unit out-of-plane force at the other. The centroidal axis of the beam is
# straight at the undeformed  configuration, while its cross-sections are
# twisted about the centroidal axis from 0 at the clamped end to pi/2 at
# the free end.
#
# Reference:
# Zupan D, Saje M (2004) On "A proposed standard set of problems to test
# finite element accuracy": the twisted beam. Finite Elements in Analysis
# and Design 40: 1445-1451.
# """)
function  Twisted_beam(dir)
  E = 0.29e8;
  nu = 0.22;
  W = 1.1;
  L = 12.;
  t =  0.32;
  nl = 2; nt = 1; nw = 1; ref = 7;
  p =   1/W/t;
  #  Loading in the Z direction
  if dir == 2
    loadv = [0;p;0]; dir = 2; uex = 0.001753248285256; # Harder: 1.754e-3;
  else
    loadv = [0;0;p]; dir = 3; uex = 0.005424534868469; # Harder: 5.424e-3;
  end

  # loadv = [0;0;p]; dir = 3; uex = 0.005424534868469; # Harder: 5.424e-3;
  #   Loading in the Y direction
  #loadv = [0;p;0]; dir = 2; uex = 0.001753248285256; # Harder: 1.754e-3;
  tolerance  = t/1000;

  fens,fes  = H20block(L,W,t, nl*ref,nw*ref,nt*ref)

  # Reshape into a twisted beam shape
  for i = 1:count(fens)
    a = fens.xyz[i,1]/L*(pi/2); y = fens.xyz[i,2]-(W/2); z = fens.xyz[i,3]-(t/2);
    fens.xyz[i,:] = [fens.xyz[i,1],y*cos(a)-z*sin(a),y*sin(a)+z*cos(a)];
  end

  # Clamped end of the beam
  l1  = selectnode(fens; box = [0 0 -100*W 100*W -100*W 100*W], inflate  =  tolerance)
  e1 = FDataDict("node_list"=>l1, "component"=>1, "displacement"=>0.0)
  e2 = FDataDict("node_list"=>l1, "component"=>2, "displacement"=>0.0)
  e3 = FDataDict("node_list"=>l1, "component"=>3, "displacement"=>0.0)

  # Traction on the opposite edge
  boundaryfes  =   meshboundary(fes);
  Toplist   = selectelem(fens,boundaryfes, box =  [L L -100*W 100*W -100*W 100*W], inflate =   tolerance);
  el1femm  = FEMMBase(IntegDomain(subset(boundaryfes,Toplist), GaussRule(2, 2)))
  flux1 = FDataDict("femm"=>el1femm, "traction_vector"=>loadv)


  # Make the region
  MR = DeforModelRed3D
  material = MatDeforElastIso(MR, 00.0, E, nu, 0.0)
  region1 = FDataDict("femm"=>FEMMDeforLinear(MR, IntegDomain(fes, GaussRule(3,2)),
  material))

  # Make model data
  modeldata =  FDataDict(
  "fens"=> fens, "regions"=>  [region1],
  "essential_bcs"=>[e1, e2, e3], "traction_bcs"=>  [flux1])

  # Call the solver
  modeldata = AlgoDeforLinearModule.linearstatics(modeldata)
  geom = modeldata["geom"]
  u = modeldata["u"]

  # Extract the solution
  nl = selectnode(fens, box = [L L -100*W 100*W -100*W 100*W],inflate = tolerance);
  nl1 = selectnode(fens, nearestto = fens.xyz[nl[1],:]+ tolerance/30*rand(size(fens.xyz, 2)));
  @test nl[1] == nl1[1]
  theutip = mean(u.values[nl,:], dims = 1)
  # println("displacement  = $(theutip[dir]) as compared to converged $uex")
  @test    abs(theutip[dir]-uex)<1.0e-5

  # # Write out mesh with displacements
  # modeldata["postprocessing"] = FDataDict("file"=>"twisted_beam")
  # modeldata = AlgoDeforLinearModule.exportdeformation(modeldata)
  #
  # # Write out mesh with stresses
  # modeldata["postprocessing"] = FDataDict("file"=>"twisted_beam",
  # "quantity"=> :Cauchy, "component"=> :xy)
  # modeldata = AlgoDeforLinearModule.exportstress(modeldata)
  #
  # # Write out mesh with stresses
  # modeldata["postprocessing"] = FDataDict("file"=>"twisted_beam",
  # "quantity"=> :Cauchy, "component"=> :xz)
  # modeldata = AlgoDeforLinearModule.exportstress(modeldata)
  #
  # # Write out mesh with von Mises stresses
  # modeldata["postprocessing"] = FDataDict("file"=>"twisted_beam",
  # "quantity"=> :vm)
  # modeldata = AlgoDeforLinearModule.exportstress(modeldata)
  #
  # # Write out mesh with von Mises stresses, elementwise
  # modeldata["postprocessing"] = FDataDict("file"=>"twisted_beam-ew",
  # "quantity"=> :vm)
  # modeldata = AlgoDeforLinearModule.exportstresselementwise(modeldata)
  #
  # # Write out mesh with von Mises stresses, elementwise
  # modeldata["postprocessing"] = FDataDict("file"=>"twisted_beam-ew",
  # "quantity"=> :Cauchy, "component"=> :xz)
  # modeldata = AlgoDeforLinearModule.exportstresselementwise(modeldata)
  #
  # println("Done")
  true

end

function test()
  Twisted_beam(2)
  Twisted_beam(3)
end

end
using .mxxxx1_06102017
mxxxx1_06102017.test()




module mx_06112017
using FinEtools
using Test
import Statistics: mean
function test()
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

  femm = FEMMDeforLinear(MR, IntegDomain(fes, GaussRule(3, 3)), material)

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
  xsfemm = FEMMDeforWinkler(IntegDomain(subset(bfes,fl), GaussRule(2, 3)))

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


  # File =  "LE11NAFEMS_H20_sigmaz.vtk"
  # vtkexportmesh(File, fens, fes;
  #      scalars=[("sigmaz", fld.values)], vectors=[("u", u.values)])
  # @async run(`"paraview.exe" $File`)
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
  # println("Stress at point A: $(sA) i. e.  $( sAn*100  )% of reference value")

  @test abs(sA-(-83.7322285847101)) < 1e-4

  ## Discussion
  #
  ##
  # The 3-D solution corresponds well to the 2-D axially symmetric model.
  # We also see good correspondence to other published solutions for
  # comparable finite element models.  For instance, Abaqus 6.11
  # Benchmark manual lists very similar numbers.
end

end
using .mx_06112017
mx_06112017.test()


module my_06112017

using FinEtools
using Test
function test()
  # NAFEMS LE11 benchmark with Q8 elements.
  # # This is a test recommended by the National Agency for Finite Element
  # # Methods and Standards (U.K.): Test LE11 from NAFEMS Publication TNSB,
  # # Rev. 3, “The Standard NAFEMS Benchmarks,” October 1990.
  # #
  # # Target solution: Direct stress,   =  –105 MPa at point A.
  #function  LE11NAFEMS()
  # Parameters:
  Ea =  210000*phun("MEGA*Pa")
  nua =  0.3;
  alphaa = 2.3e-4;              # thermal expansion coefficient
  sigmaA = -105*phun("MEGA*Pa")
  nref =  1;                        # how many times should we refine the mesh?
  X = [1.     0.;#A
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
  tolerance  = 1.e-6*phun("M")
  ##
  # Note that the material object needs to be created with the proper
  # model-dimension reduction in mind.  In this case that is the axial symmetry
  # assumption.
  MR  =  DeforModelRed2DAxisymm



  fens = FENodeSet(X);
  fes = FESetQ4([1 2 4 3; 3 4 6 5; 5 6 8 7; 7 8 10 9; 9 10 12 11; 11 12 14 13; 13 14 16 15]);
  for ref = 1:nref
    fens,fes = Q4refine(fens,fes);
    list = selectnode(fens,distance = 1.0+0.1/2^nref, from = [0. 0.], inflate = tolerance);
    fens.xyz[list,:] =  FinEtools.MeshUtilModule.ontosphere(fens.xyz[list,:],1.0);
  end
  fens,fes = Q4toQ8(fens,fes)
  list = selectnode(fens,distance = 1.0+0.1/2^nref, from = [0. 0.], inflate = tolerance);
  fens.xyz[list,:] =  FinEtools.MeshUtilModule.ontosphere(fens.xyz[list,:],1.0);

  #     File  =   "mesh.vtk"
  # vtkexportmesh(File, fens, fes)

  # now we create the geometry and displacement fields
  geom  =  NodalField(fens.xyz)
  u  =  NodalField(zeros(size(fens.xyz,1),2)) # displacement field

  # Apply EBC's
  l1 = selectnode(fens,box = [-Inf Inf 0 0],inflate = tolerance)
  setebc!(u, l1, true, 2, 00.0)
  l1 = selectnode(fens,box = [-Inf Inf 1.79  1.79],inflate = tolerance)
  setebc!(u, l1, true, 2, 00.0)
  applyebc!(u)
  numberdofs!(u)

  # Temperature field
  dT  = NodalField(reshape(fens.xyz[:,1]+fens.xyz[:,2],size(fens.xyz,1),1));


  # Property and material
  material = MatDeforElastIso(MR, 0.0, Ea, nua, alphaa)

  femm  =  FEMMDeforLinear(MR, IntegDomain(fes, GaussRule(2, 3), true), material)

  K  = stiffness(femm, geom, u)
  F  =  thermalstrainloads(femm, geom, u, dT)
  #K = cholesky(K)
  U =   K\F
  scattersysvec!(u, U[:])

  nA  = selectnode(fens,box = FFlt[1.0  1.0 0.0 0.0], inflate = tolerance);

  fld =  fieldfromintegpoints(femm, geom, u, dT, :Cauchy, 2)


  File  =   "LE11NAFEMS_Q8_sigmay.vtk"
  vtkexportmesh(File, fens, fes; scalars = [("sigmay", fld.values)],
  vectors = [("u", u.values)])
  try rm(File); catch end

  sA  =  fld.values[nA]/phun("MEGA*Pa")
  sAn  =  fld.values[nA]/sigmaA
  # println("Stress at point A: $(sA) i. e.  $( sAn*100  )% of reference value")
  @test abs(sA[1]-(-93.8569)) < 1e-3

  fen2fe  = FENodeToFEMap(connasarray(fes), nnodes(geom))
  function inspector(idat, elnum, conn, xe,  out,  xq)
    # println("loc = $(  xq  ) : $(  transpose(out)/phun("MEGA*Pa")  )")
    return idat
  end

  inspectintegpoints(femm, geom, u, dT,  fen2fe.map[nA[1]],
  inspector, []; quantity = :Cauchy)

  #finealemesh(fens,fes,"meshmfile")

  # end
  # LE11NAFEMS()
end

end
using .my_06112017
my_06112017.test()

module mmmZenkourm
using FinEtools
using FinEtools.AlgoDeforLinearModule
import Statistics: mean
using Test
function test()


  # println("""
  # % Three-dimensional Elasticity Solution for Uniformly Loaded Cross-ply
  # % Laminates and Sandwich Plates
  # % Ashraf M. Zenkour, Journal of Sandwich Structures and Materials 2007 9: 213-238
  # % DOI: 10.1177/1099636207065675
  # """)

  t0 = time()
  # Lamina material parameters
  E1s = 25.0e6*phun("psi")
  E2s = 1.0e6*phun("psi")
  E3s = E2s
  nu12s = nu13s = nu23s = 0.25
  G12s = 0.5e6*phun("psi")
  G13s = G12s
  G23s = 0.2e6*phun("psi")

  a = 200.0*phun("mm") # side of the square plate
  b = 600.0*phun("mm") # side of the square plate
  q0 = 1.0*phun("psi");
  # The below values come from Table 2
  h = a/4; wc_analytical = 3.65511/(100*E3s*h^3/a^4/q0);
  # h = a/10; wc_analytical = 1.16899/(100*E3s*h^3/a^4/q0);
  # h = a/50; wc_analytical = 0.66675/(100*E3s*h^3/a^4/q0);
  # h = a/100; wc_analytical = 0.65071/(100*E3s*h^3/a^4/q0);
  angles =[0,90,0];
  nLayers = length(angles);
  tolerance = 0.0001*h

  # Generate mesh
  na = 8 # number of elements along the side of the plate
  nb = 24 # number of elements along the side of the plate
  xs = collect(linearspace(0.0, a, na+1))
  ys = collect(linearspace(0.0, b, nb+1))
  ts = h/nLayers*ones(nLayers);# layer thicknesses
  nts= 3*ones(Int, nLayers);# number of elements per layer
  fens,fes = H8layeredplatex(xs, ys, ts, nts)
  fens,fes = H8toH20(fens,fes)

  MR = DeforModelRed3D
  laminamaterial = MatDeforElastOrtho(MR,
    0.0, E1s, E2s, E3s,
    nu12s, nu13s, nu23s,
    G12s, G13s, G23s,
    0.0, 0.0, 0.0)

  function updatecs!(csmatout::FFltMat, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)
    rotmat3!(csmatout, angles[fe_label]/180.0*pi* [0.;0.;1.]);
  end

  gr = GaussRule(3, 3)

  region = FDataDict("femm"=>FEMMDeforLinear(MR,
      IntegDomain(fes, gr), CSys(3, 3, updatecs!), laminamaterial))

  lx0 = selectnode(fens, box=[0.0 0.0 -Inf Inf -Inf Inf], inflate=tolerance)
  lxa = selectnode(fens, box=[a a -Inf Inf -Inf Inf], inflate=tolerance)
  ly0 = selectnode(fens, box=[-Inf Inf 0.0 0.0 -Inf Inf], inflate=tolerance)
  lyb = selectnode(fens, box=[-Inf Inf b b -Inf Inf], inflate=tolerance)

  ex02 = FDataDict( "displacement"=>  0.0, "component"=> 2, "node_list"=>lx0 )
  ex03 = FDataDict( "displacement"=>  0.0, "component"=> 3, "node_list"=>lx0 )
  exa2 = FDataDict( "displacement"=>  0.0, "component"=> 2, "node_list"=>lxa )
  exa3 = FDataDict( "displacement"=>  0.0, "component"=> 3, "node_list"=>lxa )
  ey01 = FDataDict( "displacement"=>  0.0, "component"=> 1, "node_list"=>ly0 )
  ey03 = FDataDict( "displacement"=>  0.0, "component"=> 3, "node_list"=>ly0 )
  eyb1 = FDataDict( "displacement"=>  0.0, "component"=> 1, "node_list"=>lyb )
  eyb3 = FDataDict( "displacement"=>  0.0, "component"=> 3, "node_list"=>lyb )

  bfes = meshboundary(fes)
  ttopl = selectelem(fens, bfes; facing=true, direction = [0.0 0.0 1.0])
  Trac = FDataDict("traction_vector"=>[0.0; 0.0; -q0],
      "femm"=>FEMMBase(IntegDomain(subset(bfes, ttopl), GaussRule(2, 3))))

  modeldata = FDataDict("fens"=>fens,
   "regions"=>[region],
   "essential_bcs"=>[ex02, ex03, exa2, exa3, ey01, ey03, eyb1, eyb3],
   "traction_bcs"=> [Trac]
   )
  modeldata = AlgoDeforLinearModule.linearstatics(modeldata)

  u = modeldata["u"]
  geom = modeldata["geom"]
  lcenter = selectnode(fens, box=[a/2 a/2  b/2 b/2 -Inf Inf], inflate=tolerance)
  cdis = abs(mean(u.values[lcenter, 3]))
  # println("")
  # println("Normalized Center deflection: $(cdis/wc_analytical)")
  @test abs(cdis/wc_analytical-1.0) < 6.0e-3
  true

end
end
using .mmmZenkourm
mmmZenkourm.test()

module mmmultimaterial_beam_xz
using FinEtools
using FinEtools.AlgoDeforLinearModule
using Test
function test()
  # println("""
  # Multi-material beam. Rubber-like and metal-like halfs,
  # clamped, with shear traction at free end.
  # """)
  E1 = 0.29e3;
  nu1 = 0.49;
  E2 = 0.4e4;
  nu2 = 0.3;
  W = 4.1;
  L = 12.;
  t = 6.5;
  nl = 2; nt = 1; nw = 1; ref = 9;
  p =   200.0/W/t;
  #  Loading in the Z direction
  loadv = [0;0;p]; dir = 3;
  tolerance  = t/1000;

  fens,fes  = H20block(L,W,t, nl*ref,nw*ref,nt*ref)

  # Clamped end of the beam
  l1  = selectnode(fens; box = [0 0 -Inf Inf -Inf Inf], inflate  =  tolerance)
  e1 = FDataDict("node_list"=>l1, "component"=>1, "displacement"=>0.0)
  e2 = FDataDict("node_list"=>l1, "component"=>2, "displacement"=>0.0)
  e3 = FDataDict("node_list"=>l1, "component"=>3, "displacement"=>0.0)

  # Traction on the opposite edge
  boundaryfes  =   meshboundary(fes);
  Toplist   = selectelem(fens,boundaryfes, box =  [L L -Inf Inf -Inf Inf], inflate =   tolerance);
  el1femm  =   FEMMBase(IntegDomain(subset(boundaryfes,Toplist), GaussRule(2, 2)))
  flux1 = FDataDict("femm"=>el1femm, "traction_vector"=>loadv)

  r1list   = selectelem(fens,fes, box =  [0 L/2. -Inf Inf -Inf Inf], inflate =   tolerance);
  r2list   = selectelem(fens,fes, box =  [L/2. L -Inf Inf -Inf Inf], inflate =   tolerance);

  # Model reduction type
  MR = DeforModelRed3D

  # Make region 1
  region1 = FDataDict("femm"=>FEMMDeforLinear(MR,
  IntegDomain(subset(fes,r1list), GaussRule(3,2)),
  MatDeforElastIso(MR, 0.0, E1, nu1, 0.0)))

  # Make region 2
  region2 = FDataDict("femm"=>FEMMDeforLinear(MR,
  IntegDomain(subset(fes,r2list), GaussRule(3,2)),
  MatDeforElastIso(MR, 0.0, E2, nu2, 0.0)))

  # Make model data
  modeldata =  FDataDict(
  "fens"=> fens, "regions"=>  [region1, region2],
  "essential_bcs"=>[e1, e2, e3], "traction_bcs"=>  [flux1])

  # Call the solver
  modeldata = AlgoDeforLinearModule.linearstatics(modeldata)
  geom = modeldata["geom"]
  u = modeldata["u"]

  # Write out mesh with displacements
  modeldata["postprocessing"] = FDataDict("file"=>"multimaterial_beam")
  modeldata = AlgoDeforLinearModule.exportdeformation(modeldata)
  for exported in modeldata["postprocessing"]["exported"]
    rm(exported["file"])
  end

  # Write out mesh with stresses
  modeldata["postprocessing"] = FDataDict("file"=>"multimaterial_beam_xy",
  "quantity"=> :Cauchy, "component"=> :xy)
  modeldata = AlgoDeforLinearModule.exportstress(modeldata)
  for exported in modeldata["postprocessing"]["exported"]
    rm(exported["file"])
  end

  # Write out mesh with stresses
  modeldata["postprocessing"] = FDataDict("file"=>"multimaterial_beam_xz",
  "quantity"=> :Cauchy, "component"=> :xz)
  modeldata = AlgoDeforLinearModule.exportstress(modeldata)
  for exported in modeldata["postprocessing"]["exported"]
    rm(exported["file"])
  end

  # Write out mesh with von Mises stresses
  modeldata["postprocessing"] = FDataDict("file"=>"multimaterial_beam_vm",
    "quantity"=> :vm)
  modeldata = AlgoDeforLinearModule.exportstress(modeldata)
  for exported in modeldata["postprocessing"]["exported"]
    rm(exported["file"])
  end

  # Write out mesh with von Mises stresses, elementwise
  modeldata["postprocessing"] = FDataDict("file"=>"multimaterial_beam_vm-ew",
   "quantity"=> :vm)
  modeldata = AlgoDeforLinearModule.exportstresselementwise(modeldata)
  for exported in modeldata["postprocessing"]["exported"]
    rm(exported["file"])
  end

  true

end
end
using .mmmultimaterial_beam_xz
mmmultimaterial_beam_xz.test()

module mmmmunitmccubemm
using FinEtools
using Test
import LinearAlgebra: norm, cholesky, cross
function test()

  # println("""
  # % Vibration modes of unit cube  of almost incompressible material.
  # %
  # % Reference: Puso MA, Solberg J (2006) A stabilized nodally integrated
  # % tetrahedral. International Journal for Numerical Methods in
  # % Engineering 67: 841-867.""")
  t0 = time()


  E = 1*phun("PA");
  nu = 0.499;
  rho= 1*phun("KG/M^3");
  a=1*phun("M"); b=a; h= a;
  n1=2;# How many element edges per side?
  na= n1; nb= n1; nh =n1;
  neigvs=20                   # how many eigenvalues
  omega_shift=(0.1*2*pi)^2;

  fens,fes =H20block(a,b,h, na,nb,nh)

  # Make the region
  MR = DeforModelRed3D
  material = MatDeforElastIso(MR, rho, E, nu, 0.0)
  region1 = FDataDict("femm"=>FEMMDeforLinear(MR, IntegDomain(fes, GaussRule(3,2)),
    material), "femm_mass"=>FEMMDeforLinear(MR, IntegDomain(fes, GaussRule(3,3)),
    material))

  # Make model data
  modeldata =  FDataDict(
    "fens"=> fens, "regions"=>  [region1],
    "omega_shift"=>omega_shift, "neigvs"=>neigvs)

  # Solve
  modeldata = FinEtools.AlgoDeforLinearModule.modal(modeldata)

  fs = modeldata["omega"]/(2*pi)
  # println("Eigenvalues: $fs [Hz]")

  @test norm(fs-[2.01585e-8, 7.50057e-8, 7.71588e-8, 8.83876e-8, 8.84035e-8, 9.99803e-8, 0.267319, 0.267319, 0.365102, 0.365102, 0.365102, 0.369647, 0.369647, 0.369647, 0.399873, 0.408512, 0.408512, 0.464131, 0.464131, 0.464131]) < 1.0e-5

  modeldata["postprocessing"] = FDataDict("file"=>"unit_cube_mode",
    "mode"=>10)
  modeldata = FinEtools.AlgoDeforLinearModule.exportmode(modeldata)
  # @async run(`"paraview.exe" $(modeldata["postprocessing"]["file"]*"1.vtk")`)
  for exported in modeldata["postprocessing"]["exported"]
    rm(exported["file"])
  end

  true
end
end
using .mmmmunitmccubemm
mmmmunitmccubemm.test()

module mmpipemmPSmmm
using FinEtools
using Test
import LinearAlgebra: norm, cholesky, cross, dot
import Statistics: mean

mutable struct MyIData
    c::FInt
    r::FFltVec
    s::FFltVec
end

function test()

## Thick pipe with internal pressure: plane strain
#

##
# Link to the  <matlab:edit('pub_thick_pipe_ps') m-file>.

## Description
##
# This is a simple modification of the full three-dimensional simulation of
# the tutorial pub_thick_pipe takes advantage of the plane-strain model
# reduction procedure.
##
# An infinitely long thick walled cylindrical pipe
# with inner boundary radius of 3 mm and outer boundary radius of 9 mm is
# subjected to an internal pressure of 1.0 MPa. A wedge   with thickness of
# 2 mm and a 90-degree angle sector is considered for the finite element
# analysis. The material properties are taken as  isotropic linear elastic
# with $E=1000$ MPa and $\nu=0.4999$ to represent nearly incompressible
# behavior. This problem has been proposed to by MacNeal and Harder as a
# test of an element's ability to represent the  response of a nearly
# incompressible material. The plane-strain condition is assumed in the
# axial direction of the pipe which together with the radial symmetry
# confines the material in all but the radial direction and therefore
# amplifies the numerical difficulties associated with the confinement of
# the nearly incompressible material.
##
# There is an analytical solution to this problem. Timoshenko and Goodier
# presented the original solution of Lame in their textbook. We are going
# to compare with  both the stress distribution (radial and hoop stresses)
# and the displacement of the inner  cylindrical surface.

##
#
# <html>
# <table border=0><tr><td>
# <img src="../docs/pub_thick_pipe_ps.png" width = "30#">
# </td></tr>
# <tr><td>Figure 1. Definition of the geometry of the internally pressurized thick pipe</td></tr>
# </table>
# </html>

##
# References:
#
# # Macneal RH, Harder RL (1985) A proposed standard set of problems to test
# finite element accuracy. Finite Elements in Analysis and Design 1: 3-20.
#
# # Timoshenko S. and Goodier J. N., Theory of Elasticity, McGraw-Hill, 2nd ed., 1951.

## Solution
#

##
# Internal radius of the pipe.
a = 3*phun("MM");
##
# External radius of the pipe.
b = 9*phun("MM");
##
# Thickness of the slice.
t = 2*phun("MM");

##
# Geometrical tolerance.
tolerance  =a/10000.;
##
# Young's modulus and Poisson's ratio.
E = 1000*phun("MEGA*PA");
nu = 0.499;
##
# Applied pressure on the internal surface.
press = 1.0*phun("MEGA*PA");

##
# Analytical solutions.   Radial stress:
radial_stress(r) =press*a.^2/(b^2-a^2).*(1-(b^2)./r.^2);
##
# Circumferential (hoop) stress:
hoop_stress(r)=press*a.^2/(b^2-a^2).*(1+(b^2)./r.^2);

##
# Radial displacement:
radial_displacement(r)=press*a^2*(1+nu)*(b^2+r.^2*(1-2*nu))/(E*(b^2-a^2).*r);;

##
# Therefore the radial displacement of the loaded surface will be:
urex = radial_displacement(a);


##
# The mesh parameters: The numbers of element edges axially,
# and through the thickness of the pipe wall (radially).

nc=3; nt=3;

##
# Note that the material object needs to be created with the proper
# model-dimension reduction in mind.  In this case that is the axial symmetry
# assumption.
MR = DeforModelRed2DStrain

# Create the mesh and initialize the geometry.  First we are going
# to construct the block of elements with the first coordinate
# corresponding to the angle, and the second
# coordinate is the thickness in the radial direction.
anglrrange = 90.0/180*pi;
fens,fes =  Q8block(anglrrange, b-a, nc, nt);

# Extract the boundary  and mark the finite elements on the
# interior surface.
bdryfes = meshboundary(fes);
bcl = selectelem(fens, bdryfes, box=[-Inf,Inf,0.,0.], inflate=tolerance);
internal_fenids = connectednodes(subset(bdryfes,bcl));
# Now  shape the block  into  the actual wedge piece of the pipe.
ayr=fens.xyz;
for i=1:count(fens)
    angl=ayr[i,1]; r=a+ayr[i,2];
    fens.xyz[i,:] = [r*sin(angl),(r*cos(angl))];
end

# now we create the geometry and displacement fields
geom = NodalField(fens.xyz)
u = NodalField(zeros(size(fens.xyz,1),2)) # displacement field

# The symmetry boundary condition  is specified by selecting nodes
# on the plane x=0.
l1 = selectnode(fens; box=[0.0 0.0 -Inf Inf], inflate = tolerance)
setebc!(u, l1, true, 1, 0.0)
# The second symmetry boundary condition is specified by selecting
# nodes on the plane y=0.
l1 = selectnode(fens; box=[-Inf Inf 0.0 0.0], inflate = tolerance)
setebc!(u, l1, true, 2, 0.0)

applyebc!(u)
numberdofs!(u)

# The traction boundary condition is applied in the radial
# direction.

el1femm =  FEMMBase(IntegDomain(subset(bdryfes,bcl), GaussRule(1, 3)))
function pressureloading!(forceout::FFltVec, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)
  copyto!(forceout, XYZ/norm(XYZ)*press)
  return forceout
end
fi = ForceIntensity(FFlt, 2, pressureloading!); # pressure normal to the internal cylindrical surface
F2 = distribloads(el1femm, geom, u, fi, 2);

# Property and material
material = MatDeforElastIso(MR, E, nu)

femm = FEMMDeforLinear(MR, IntegDomain(fes, GaussRule(2, 2)), material)

K =stiffness(femm, geom, u)
#K=cholesky(K)
U=  K\(F2)
scattersysvec!(u, U[:])

# Transfer the solution of the displacement to the nodes on the
# internal cylindrical surface and convert to
# cylindrical-coordinate displacements there.
uv = u.values[internal_fenids,:];
ur = zeros(FFlt,length(internal_fenids));

for j = 1:length(internal_fenids)
    n = fens.xyz[internal_fenids[j],:];
    n = n'/norm(n);# normal to the cylindrical internal surface
    ur[j] = dot(vec(uv[j,:]),vec(n));
end

# Report the  relative displacement on the internal surface:
# println("(Approximate/true displacement) at the internal surface: $( mean(ur)/urex*100  ) %")
@test abs(mean(ur)/urex*100 - 100) < 0.1

# Produce a plot of the radial stress component in the cylindrical
# coordinate system. Note that this is the usual representation of
# stress using nodal stress field.

fld = fieldfromintegpoints(femm, geom, u, :Cauchy, 1)


File =  "thick_pipe_sigmax.vtk"
vtkexportmesh(File, fens, fes; scalars=[("sigmax", fld.values)])
try rm(File); catch end

# Produce a plot of the solution components in the cylindrical
# coordinate system.
# Plot the analytical solution.

function inspector(idat::MyIData, elnum, conn, xe,  out,  xq)
  function outputRm(c)
    theNormal=c;
    r=norm(theNormal);# distance from the axis of symmetry
    theNormal =theNormal/r;# compute the unit normal vector
    e1p=[theNormal';0.];# local cylind. coordinate  basis vectors
    e3p=[0.,0.,1.]';# this one points along the axis of the cylinder
    e2p=cross(vec(e3p),vec(e1p));# this one is along the hoop direction
    R= [vec(e1p) vec(e2p) vec(e3p)];# transformation matrix for the stress
    return R
  end
  Rm=outputRm(xq)
  tm=zeros(FFlt,3,3)
  stress4vto3x3t!(tm, out);# stress in global XYZ
  tpm = Rm'*tm*Rm;#  stress matrix in cylindrical coordinates
  sp=zeros(FFlt,6)
  stress3x3tto6v!(sp, tpm);# stress vector in cylindr. coord.
  push!(idat.r,norm(xq))
  push!(idat.s,sp[idat.c])
  return idat
end

idat = MyIData(1, FFltVec[], FFltVec[])
idat = inspectintegpoints(femm, geom, u, collect(1:count(fes)),
 inspector, idat, :Cauchy)
# show(idat)

@test norm(idat.s - [-7.44858e5, -3.55143e5, -7.44858e5, -3.55143e5, -2.19961e5,
-1.08612e5, -2.19961e5, -1.08612e5, -58910.9, -12517.6, -58910.9, -12517.6,
-7.44858e5, -3.55143e5, -7.44858e5, -3.55143e5, -2.19961e5, -1.08612e5,
-2.19961e5, -1.08612e5, -58910.9, -12517.6, -58910.9, -12517.6, -7.44858e5,
-3.55143e5, -7.44858e5, -3.55143e5, -2.19961e5, -1.08612e5, -2.19961e5,
 -1.08612e5, -58910.9, -12517.6, -58910.9, -12517.6])/1.0e5 < 1.e-3
# using Plots
# plotly()
#
# # Plot the analytical solution.
# r = linearspace(a,b,100);
# plot(r, radial_stress(r))
# # Plot the computed  integration-point data
# scatter!(idat.r, idat.s, m=:circle, color=:red)
# gui()
end
end
using .mmpipemmPSmmm
mmpipemmPSmmm.test()

module mmOrthotropicmm
using FinEtools
using Test
import LinearAlgebra: norm, cholesky, cross
function test()

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

  # The traction boundary condition is applied in the radial
  # direction.

  el1femm =  FEMMBase(IntegDomain(subset(bdryfes,icl), GaussRule(1, 3), true))
  function pressureloading!(forceout::FFltVec, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)
    copyto!(forceout, XYZ/norm(XYZ)*p)
    return forceout
  end
  fi = ForceIntensity(FFlt, 2, pressureloading!); # pressure normal to the internal cylindrical surface
  F2= distribloads(el1femm, geom, u, fi, 2);

  # Property and material
  material=MatDeforElastOrtho(MR, E1,E2,E3,nu12,nu13,nu23,G12,G13,G23)

  femm = FEMMDeforLinear(MR, IntegDomain(fes, GaussRule(2, 2), true), material)

  K =stiffness(femm, geom, u)
  U=  K\(F2)
  scattersysvec!(u,U[:])

  # Produce a plot of the radial stress component in the cylindrical
  # coordinate system. Note that this is the usual representation of
  # stress using nodal stress field.

  fld= fieldfromintegpoints(femm, geom, u, :Cauchy, 3)


  File =  "orthoballoon_sigmaz.vtk"
  vtkexportmesh(File, fens, fes; scalars=[("sigmaz", fld.values)],
                vectors=[("u", u.values)])
  # @async run(`"paraview.exe" $File`)
  try rm(File); catch end

  @test abs(minimum(fld.values)-(-0.05139098099088048)) < 1.0e-5
  @test abs(maximum(fld.values)-0.5704453140236726) < 1.0e-5

end
end
using .mmOrthotropicmm
mmOrthotropicmm.test()

module mmCookmm
using FinEtools
using FinEtools.AlgoDeforLinearModule
using Test
import LinearAlgebra: norm, cholesky, cross
function test()
  # println("Cook plane stress, with quadratic triangles.")
  E = 1.0;
  nu = 1.0/3;
  width = 48.0; height = 44.0; thickness  = 1.0;
  free_height  = 16.0;
  Mid_edge  = [48.0, 52.0];# Location of tracked  deflection
  magn = 1.0/free_height;# Magnitude of applied load
  convutip = 23.97;
  n = 10;#*int(round(sqrt(170.)/2.)); # number of elements per side
  tolerance = minimum([width, height])/n/1000.;#Geometrical tolerance

  fens,fes = T6block(width, height, n, n)

  # Reshape into a trapezoidal panel
  for i = 1:count(fens)
      fens.xyz[i,2] = fens.xyz[i,2]+(fens.xyz[i,1]/width)*(height -fens.xyz[i,2]/height*(height-free_height));
  end

  # Clamped edge of the membrane
  l1 = selectnode(fens; box=[0.,0.,-Inf, Inf], inflate = tolerance)
  ess1 = FDataDict("displacement"=>  0.0, "component"=> 1, "node_list"=>l1)
  ess2 = FDataDict("displacement"=>  0.0, "component"=> 2, "node_list"=>l1)

  # Traction on the opposite edge
  boundaryfes =  meshboundary(fes);
  Toplist  = selectelem(fens, boundaryfes, box= [width, width, -Inf, Inf ], inflate=  tolerance);
  el1femm = FEMMBase(IntegDomain(subset(boundaryfes, Toplist), GaussRule(1, 3)))
  flux1 = FDataDict("traction_vector"=>[0.0,+magn],
      "femm"=>el1femm
      )

  # Make the region
  MR = DeforModelRed2DStress
  material = MatDeforElastIso(MR,  0.0, E, nu, 0.0)
  region1 = FDataDict("femm"=>FEMMDeforLinear(MR,
      IntegDomain(fes, TriRule(3)), material))

  modeldata = FDataDict("fens"=>fens,
   "regions"=>[region1],
   "essential_bcs"=>[ess1, ess2],
   "traction_bcs"=>[flux1]
   )

  # Call the solver
  modeldata = AlgoDeforLinearModule.linearstatics(modeldata)

  u = modeldata["u"]
  geom = modeldata["geom"]

  # Extract the solution
  nl = selectnode(fens, box=[Mid_edge[1],Mid_edge[1],Mid_edge[2],Mid_edge[2]],
            inflate=tolerance);
  theutip = u.values[nl,:]
  # println("displacement =$(theutip[2]) as compared to converged $convutip")
  @test abs(theutip[2]-23.93976266124016) < 1.e-5

  modeldata["postprocessing"] = FDataDict("file"=>"cookstress",
     "quantity"=>:Cauchy, "component"=>:xy)
  modeldata = AlgoDeforLinearModule.exportstress(modeldata)
  File = modeldata["postprocessing"]["exported"][1]["file"]
  # @async run(`"paraview.exe" $File`)
  for exported in modeldata["postprocessing"]["exported"]
    rm(exported["file"])
  end
  fld = modeldata["postprocessing"]["exported"][1]["field"]
  # println("$(minimum(fld.values)) $(maximum(fld.values))")
  @test norm([minimum(fld.values) maximum(fld.values)] -
    [-0.06292574794975273 0.12022571940768388]) < 1.0e-5
end
end
using .mmCookmm
mmCookmm.test()

module mCookmmfakeorthom
using FinEtools
using FinEtools.AlgoDeforLinearModule
using FinEtools.MeshExportModule
using Test
import LinearAlgebra: norm, cholesky, cross
function test()
  # println("Cook plane stress, with quadratic triangles. With orthotropic  material model.")
  E = 1.0;
  nu = 1.0/3;
  E1 = E2 = E3 = E;
  nu12 = nu13 = nu23 = nu;
  G12 = G13 = G23 = E/2.0/(1+nu);
  width = 48.0; height = 44.0; thickness  = 1.0;
  free_height  = 16.0;
  Mid_edge  = [48.0, 52.0];# Location of tracked  deflection
  magn = 1.0/free_height;# Magnitude of applied load
  convutip = 23.97;
  n = 10;#*int(round(sqrt(170.)/2.)); # number of elements per side
  tolerance = minimum([width, height])/n/1000.;#Geometrical tolerance

  fens,fes = T6block(width, height, n, n)

  # Reshape into a trapezoidal panel
  for i = 1:count(fens)
      fens.xyz[i,2] = fens.xyz[i,2]+(fens.xyz[i,1]/width)*(height -fens.xyz[i,2]/height*(height-free_height));
  end

  # Clamped edge of the membrane
  l1 = selectnode(fens; box=[0.,0.,-Inf, Inf], inflate = tolerance)
  ess1 = FDataDict("displacement"=>  0.0, "component"=> 1, "node_list"=>l1)
  ess2 = FDataDict("displacement"=>  0.0, "component"=> 2, "node_list"=>l1)

  # Traction on the opposite edge
  boundaryfes =  meshboundary(fes);
  Toplist  = selectelem(fens, boundaryfes, box= [width, width, -Inf, Inf ], inflate=  tolerance);
  el1femm = FEMMBase(IntegDomain(subset(boundaryfes, Toplist), GaussRule(1, 3)))
  flux1 = FDataDict("traction_vector"=>[0.0,+magn],
      "femm"=>el1femm
      )

  # Make the region
  MR = DeforModelRed2DStress
  # This material model is orthotropic,  but the input parameters correspond to an
  # isotropiic material  model..
  material=MatDeforElastOrtho(MR, E1,E2,E3,nu12,nu13,nu23,G12,G13,G23)
  region1 = FDataDict("femm"=>FEMMDeforLinear(MR,
      IntegDomain(fes, TriRule(3)), material))

  modeldata = FDataDict("fens"=>fens,
   "regions"=>[region1],
   "essential_bcs"=>[ess1, ess2],
   "traction_bcs"=>[flux1]
   )

  # Call the solver
  modeldata = AlgoDeforLinearModule.linearstatics(modeldata)

  u = modeldata["u"]
  geom = modeldata["geom"]

  # Extract the solution
  nl = selectnode(fens, box=[Mid_edge[1],Mid_edge[1],Mid_edge[2],Mid_edge[2]],
            inflate=tolerance);
  theutip = u.values[nl,:]
  @test abs(theutip[2]-23.93976266124016) < 1.e-5

  modeldata["postprocessing"] = FDataDict("file"=>"cookstress",
     "quantity"=>:Cauchy, "component"=>:xy)
  modeldata = AlgoDeforLinearModule.exportstress(modeldata)
  # @async run(`"paraview.exe" $File`)
  for exported in modeldata["postprocessing"]["exported"]
    rm(exported["file"])
  end
  fld = modeldata["postprocessing"]["exported"][1]["field"]
  # println("$(minimum(fld.values)) $(maximum(fld.values))")
  @test norm([minimum(fld.values) maximum(fld.values)] -
    [-0.06292574794975273 0.12022571940768388]) < 1.0e-5

end
end
using .mCookmmfakeorthom
mCookmmfakeorthom.test()

module mmCanttronglymm
using FinEtools
using FinEtools.AlgoDeforLinearModule
using Test
import LinearAlgebra: norm, cholesky, cross
import Statistics: mean
function test()
  # println("""
  # Cantilever example.  Strongly orthotropic material. Orientation "y".
  # @article{
  # author = {Krysl, P.},
  # title = {Mean-strain 8-node hexahedron with optimized energy-sampling stabilization},
  # journal = {Finite Elements in Analysis and Design},
  # volume = {108}, pages = {41-53}, DOI = {10.1016/j.finel.2015.09.008}, year = {2016}
  # }
  # """)

  # t0 = time()
  # # Orthotropic material
  E1s = 100000.0*phun("GPa")
  E2s = 1.0*phun("GPa")
  E3s = E2s
  nu23s = nu12s = nu13s = 0.25
  G12s = 0.2*phun("GPa")
  G23s = G13s = G12s
  CTE1 = 0.0
  CTE2 = 0.0
  CTE3 = 0.0
  # # Isotropic material
  # E = 1.0e9*phun("Pa")
  # nu = 0.25
  # CTE = 0.0

  # Reference value for  the vertical deflection of the tip
  uz_ref = -1.027498445054843e-05;

  a = 90.0*phun("mm") # length of the cantilever
  b = 10.0*phun("mm") # width of the cross-section
  t = 20.0*phun("mm") # height of the cross-section
  q0 = -1000.0*phun("Pa") # shear traction
  dT = 0*phun("K") # temperature rise

  tolerance = 0.00001*t

  # Generate mesh
  n = 4
  na = 4*n # number of elements lengthwise
  nb = n # number of elements through the wwith
  nt = n # number of elements through the thickness
  xs = collect(linearspace(0.0, a, na+1))
  ys = collect(linearspace(0.0, b, nb+1))
  ts = collect(linearspace(0.0, t, nt+1))
  fens,fes = H8blockx(xs, ys, ts)
  fens,fes = H8toH20(fens,fes)
  bfes = meshboundary(fes)
  # end cross-section surface  for the shear loading
  sshearl = selectelem(fens, bfes; facing=true, direction = [+1.0 0.0 0.0])

  MR = DeforModelRed3D
  material = MatDeforElastOrtho(MR,
    0.0, E1s, E2s, E3s,
    nu12s, nu13s, nu23s,
    G12s, G13s, G23s,
    CTE1, CTE2, CTE3)
  # material = MatDeforElastIso(MR,
  #   0.0, E, nu, CTE)

  # Material orientation matrix
  csmat = zeros(3, 3)
  rotmat3!(csmat, -45.0/180.0*pi*[0,1,0])

  function updatecs!(csmatout::FFltMat, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)
    copyto!(csmatout, csmat)
  end

  gr = GaussRule(3, 2)

  region = FDataDict("femm"=>FEMMDeforLinear(MR, IntegDomain(fes, gr), CSys(3, 3, updatecs!), material))

  lx0 = selectnode(fens, box=[0.0 0.0 -Inf Inf -Inf Inf], inflate=tolerance)

  ex01 = FDataDict( "displacement"=>  0.0, "component"=> 1, "node_list"=>lx0 )
  ex02 = FDataDict( "displacement"=>  0.0, "component"=> 2, "node_list"=>lx0 )
  ex03 = FDataDict( "displacement"=>  0.0, "component"=> 3, "node_list"=>lx0 )

  function getshr!(forceout::FFltVec, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)
    copyto!(forceout, q0*[0.0; 0.0; 1.0])
  end

  Trac = FDataDict("traction_vector"=>getshr!,
      "femm"=>FEMMBase(IntegDomain(subset(bfes, sshearl), GaussRule(2, 3))))

  modeldata = FDataDict("fens"=>fens,
   "regions"=>[region],
   "essential_bcs"=>[ex01, ex02, ex03],
   "traction_bcs"=>[Trac],
   "temperature_change"=>FDataDict("temperature"=>dT)
   )
  modeldata = AlgoDeforLinearModule.linearstatics(modeldata)

  u = modeldata["u"]
  geom = modeldata["geom"]

  Tipl = selectnode(fens, box=[a a b b 0. 0.], inflate=tolerance)
  utip = mean(u.values[Tipl, 3])
  # println("Deflection $utip, normalized: $(utip/uz_ref)")
  @test abs(abs(utip/uz_ref) - 0.8140937283389692) < 1.0e-5
  # println("Solution: $(  time()-t0 )")

  # File =  "NAFEMS-R0031-2-plate.vtk"
  # vtkexportmesh(File, fes.conn, geom.values, FinEtools.MeshExportModule.H20;
  #     scalars = [("Layer", fes.label)], vectors = [("displacement", u.values)])
  # @async run(`"paraview.exe" $File`)

  modeldata["postprocessing"] = FDataDict("file"=>"fiber_reinf_cant_yn_strong",
    "outputcsys"=>CSys(3, 3, updatecs!), "quantity"=>:Cauchy, "component"=>5)
  modeldata = AlgoDeforLinearModule.exportstress(modeldata)
  # @async run(`"paraview.exe" $File`)
  for exported in modeldata["postprocessing"]["exported"]
    rm(exported["file"])
  end
  fld = modeldata["postprocessing"]["exported"][1]["field"]
  # println("$(minimum(fld.values)) $(maximum(fld.values))")
  @test norm([minimum(fld.values) maximum(fld.values)] - [-9990.853445671826 8228.919034163288]) < 1.e-3

  # modeldata = AlgoDeforLinearModule.exportdeformation(modeldata)
  # File = modeldata["postprocessing"]["exported"][1]["file"]
  # @async run(`"paraview.exe" $File`)

  # println("Done: $(  time()-t0 )")
  true

end
end
using .mmCanttronglymm
mmCanttronglymm.test()

module mmmNAFEMS_R0031_3m
using FinEtools
using FinEtools.AlgoDeforLinearModule
using Test
import Statistics: mean
function test()
  # println("""
  # NAFEMS publication R0031/3 Composite plate test.
  # Simply supported on all four edges.  Uniform transverse  loading.
  # The modeled part is one quarter of the full plate here.
  # """)

  # This is a test recommended by the National Agency for Finite Element Methods
  # and Standards (U.K.): Test R0031/3 from NAFEMS publication R0031, “Composites
  # Benchmarks,” February 1995.
  t0 = time()
  # Skin (face) material parameters
  E1s = 1.0e7*phun("psi")
  E2s = 0.4e7*phun("psi")
  E3s = 0.4e7*phun("psi")
  nu12s = 0.3
  nu13s = 0.3
  nu23s = 0.3
  G12s = 0.1875e7*phun("psi")
  G13s = 0.1875e7*phun("psi")
  G23s = 0.1875e7*phun("psi")
  # Core material parameters
  E1c = 10.0*phun("psi")
  E2c = 10.0*phun("psi")
  E3c = 10e4.*phun("psi")
  nu12c = 0.
  nu13c = 0.
  nu23c = 0.
  G12c = 10.0*phun("psi")
  G13c = 3.0e4*phun("psi")
  G23c = 1.2e4*phun("psi")
  L = 10.0*phun("in") # side of the square plate
  nL = 8 # number of elements along the side of the plate
  tolerance = 0.0001*phun("in")
  xs = collect(linearspace(0.0, L/2, nL+1))
  ys = collect(linearspace(0.0, L/2, nL+1))
  ts = [0.028; 0.75; 0.028]*phun("in")
  nts = [2; 3;  2; ] # number of elements through the thickness
  tmag = 100*phun("psi")

  # Generate mesh
  fens,fes = H8layeredplatex(xs, ys, ts, nts)
  fens,fes = H8toH20(fens,fes)

  MR = DeforModelRed3D
  skinmaterial = MatDeforElastOrtho(MR,
    0.0, E1s, E2s, E3s,
    nu12s, nu13s, nu23s,
    G12s, G13s, G23s,
    0.0, 0.0, 0.0)
  corematerial = MatDeforElastOrtho(MR,
    0.0, E1c, E2c, E3c,
    nu12c, nu13c, nu23c,
    G12c, G13c, G23c,
    0.0, 0.0, 0.0)

  gr = GaussRule(3, 3)

  rl1 = selectelem(fens, fes, label=1)
  skinbot = FDataDict("femm"=>FEMMDeforLinear(MR,
      IntegDomain(subset(fes, rl1), gr), skinmaterial))

  rl3 = selectelem(fens, fes, label=3)
  skintop = FDataDict("femm"=>FEMMDeforLinear(MR,
      IntegDomain(subset(fes, rl3), gr), skinmaterial))

  rl2 = selectelem(fens, fes, label=2)
  core = FDataDict("femm"=>FEMMDeforLinear(MR,
      IntegDomain(subset(fes, rl2), gr), corematerial))

  lx0 = selectnode(fens, box=[0.0 0.0 -Inf Inf -Inf Inf], inflate=tolerance)
  lxL2 = selectnode(fens, box=[L/2 L/2 -Inf Inf -Inf Inf], inflate=tolerance)
  ly0 = selectnode(fens, box=[-Inf Inf 0.0 0.0 -Inf Inf], inflate=tolerance)
  lyL2 = selectnode(fens, box=[-Inf Inf L/2 L/2 -Inf Inf], inflate=tolerance)

  ex0 = FDataDict( "displacement"=>  0.0, "component"=> 3, "node_list"=>lx0 )
  exL2 = FDataDict( "displacement"=>  0.0, "component"=> 1, "node_list"=>lxL2 )
  ey0 = FDataDict( "displacement"=>  0.0, "component"=> 3, "node_list"=>ly0 )
  eyL2 = FDataDict( "displacement"=>  0.0, "component"=> 2, "node_list"=>lyL2 )

  bfes = meshboundary(fes)
  ttopl = selectelem(fens, bfes; facing=true, direction = [0.0 0.0 1.0])
  Trac = FDataDict("traction_vector"=>[0.0; 0.0; -tmag],
      "femm"=>FEMMBase(IntegDomain(subset(bfes, ttopl), GaussRule(2, 3))))

  modeldata = FDataDict("fens"=>fens,
   "regions"=>[skinbot, core, skintop], "essential_bcs"=>[ex0, exL2, ey0, eyL2],
   "traction_bcs"=> [Trac]
   )
  modeldata = AlgoDeforLinearModule.linearstatics(modeldata)

  u = modeldata["u"]
  geom = modeldata["geom"]
  lcenter = selectnode(fens, box=[L/2 L/2  L/2 L/2 -Inf Inf], inflate=tolerance)
  cdis = mean(u.values[lcenter, 3])/phun("in")
  # println("Center node displacements $(cdis) [in]; NAFEMS-R0031-3 lists –0.123	[in]")
  # println("")
  @test abs(cdis - (-0.13634800328800462)) < 1.0e-5

  File =  "NAFEMS-R0031-3-plate.vtk"
  vtkexportmesh(File, connasarray(fes), geom.values, FinEtools.MeshExportModule.H20;
      scalars = [("Layer", fes.label)], vectors = [("displacement", u.values)])
  # @async run(`"paraview.exe" $File`)
  try  rm(File); catch end

end
end
using .mmmNAFEMS_R0031_3m
mmmNAFEMS_R0031_3m.test()

module mmtwistedmsh8mmm
using FinEtools
using FinEtools.AlgoDeforLinearModule
import LinearAlgebra: norm, cholesky, cross
using Test
import Statistics: mean
function test()

  # println("""
  # The initially twisted cantilever beam is one of the standard test
  # problems for verifying the finite-element accuracy [1]. The beam is
  #   clamped at one end and loaded either with unit in-plane or
  #   unit out-of-plane force at the other. The centroidal axis of the beam is
  #   straight at the undeformed  configuration, while its cross-sections are
  #   twisted about the centroidal axis from 0 at the clamped end to pi/2 at
  #   the free end.
  #
  #   Reference:
  #   Zupan D, Saje M (2004) On "A proposed standard set of problems to test
  #   finite element accuracy": the twisted beam. Finite Elements in Analysis
  #   and Design 40: 1445-1451.
  #   """)
  E = 0.29e8;
  nu = 0.22;
  W = 1.1;
  L = 12.;
  t =  0.32;
  nl = 2; nt = 1; nw = 1; ref = 7;
  p =   1/W/t;
  #  Loading in the Z direction
  loadv = [0;0;p]; dir = 3; uex = 0.005424534868469; # Harder: 5.424e-3;
  #   Loading in the Y direction
  #loadv = [0;p;0]; dir = 2; uex = 0.001753248285256; # Harder: 1.754e-3;
  tolerance  = t/1000;

  fens,fes  = H8block(L,W,t, nl*ref,nw*ref,nt*ref)

  # Reshape into a twisted beam shape
  for i = 1:count(fens)
    a = fens.xyz[i,1]/L*(pi/2); y = fens.xyz[i,2]-(W/2); z = fens.xyz[i,3]-(t/2);
    fens.xyz[i,:] = [fens.xyz[i,1],y*cos(a)-z*sin(a),y*sin(a)+z*cos(a)];
  end

  # Clamped end of the beam
  l1  = selectnode(fens; box = [0 0 -100*W 100*W -100*W 100*W], inflate  =  tolerance)
  e1 = FDataDict("node_list"=>l1, "component"=>1, "displacement"=>0.0)
  e2 = FDataDict("node_list"=>l1, "component"=>2, "displacement"=>0.0)
  e3 = FDataDict("node_list"=>l1, "component"=>3, "displacement"=>0.0)

  # Traction on the opposite edge
  boundaryfes  =   meshboundary(fes);
  Toplist   = selectelem(fens,boundaryfes, box =  [L L -100*W 100*W -100*W 100*W], inflate =   tolerance);
  el1femm  = FEMMBase(IntegDomain(subset(boundaryfes,Toplist), GaussRule(2, 2)))
  flux1 = FDataDict("femm"=>el1femm, "traction_vector"=>loadv)


  # Make the region
  MR = DeforModelRed3D
  material = MatDeforElastIso(MR, 00.0, E, nu, 0.0)
  region1 = FDataDict("femm"=>FEMMDeforLinearMSH8(MR, IntegDomain(fes, GaussRule(3,2)),
  material))

  # Make model data
  modeldata =  FDataDict(
  "fens"=> fens, "regions"=>  [region1],
  "essential_bcs"=>[e1, e2, e3], "traction_bcs"=>  [flux1])

  # Call the solver
  modeldata = AlgoDeforLinearModule.linearstatics(modeldata)
  geom = modeldata["geom"]
  u = modeldata["u"]

  # Extract the solution
  nl = selectnode(fens, box = [L L -100*W 100*W -100*W 100*W],inflate = tolerance);
  theutip = mean(u.values[nl,:], dims = 1)
  # println("displacement  = $(theutip[dir]) as compared to converged $uex")
  @test abs(theutip[dir]-uex)/uex < 0.0012

  # Write out mesh with displacements
  modeldata["postprocessing"] = FDataDict("file"=>"twisted_beam_msh8")
  modeldata = AlgoDeforLinearModule.exportdeformation(modeldata)
  try rm(modeldata["postprocessing"]["exported"][1]["file"]); catch end

  # Write out mesh with stresses
  modeldata["postprocessing"] = FDataDict("file"=>"twisted_beam_msh8",
  "quantity"=> :Cauchy, "component"=> :xy)
  modeldata = AlgoDeforLinearModule.exportstress(modeldata)
  try rm(modeldata["postprocessing"]["exported"][1]["file"]); catch end

  # Write out mesh with stresses
  modeldata["postprocessing"] = FDataDict("file"=>"twisted_beam_msh8",
  "quantity"=> :Cauchy, "component"=> :xz)
  modeldata = AlgoDeforLinearModule.exportstress(modeldata)
  try rm(modeldata["postprocessing"]["exported"][1]["file"]); catch end

  # Write out mesh with von Mises stresses
  modeldata["postprocessing"] = FDataDict("file"=>"twisted_beam_msh8",
  "quantity"=> :vm)
  modeldata = AlgoDeforLinearModule.exportstress(modeldata)
  vm  = modeldata["postprocessing"]["exported"][1]["field"]
  # println("extremes of vm, nodal: $([minimum(vm.values),   maximum(vm.values)])")
  try rm(modeldata["postprocessing"]["exported"][1]["file"]); catch end
  @test norm([minimum(vm.values),   maximum(vm.values)]-[4.78774, 522.126]) < 0.01

  # Write out mesh with von Mises stresses, elementwise
  modeldata["postprocessing"] = FDataDict("file"=>"twisted_beam_msh8-ew",
  "quantity"=> :vm)
  modeldata = AlgoDeforLinearModule.exportstresselementwise(modeldata)
  vm  = modeldata["postprocessing"]["exported"][1]["field"]
  # println("extremes of vm, elemental: $([minimum(vm.values),   maximum(vm.values)])")
  try rm(modeldata["postprocessing"]["exported"][1]["file"]); catch end
  @test norm([minimum(vm.values),   maximum(vm.values)]-[1.85882, 522.126]) < 0.01

  # Write out mesh with von Mises stresses, elementwise
  modeldata["postprocessing"] = FDataDict("file"=>"twisted_beam_msh8-ew",
  "quantity"=> :Cauchy, "component"=> :xz)
  modeldata = AlgoDeforLinearModule.exportstresselementwise(modeldata)
  try rm(modeldata["postprocessing"]["exported"][1]["file"]); catch end

end
end
using .mmtwistedmsh8mmm
mmtwistedmsh8mmm.test()

module mmunitmmccubemmvibrationmmms
using FinEtools
using Test
import LinearAlgebra: norm, cholesky, cross
function test()

  # println("""
  # % Vibration modes of unit cube  of almost incompressible material.
  # % Mean-strain hexahedron.
  # % Reference: Puso MA, Solberg J (2006) A stabilized nodally integrated
  # % tetrahedral. International Journal for Numerical Methods in
  # % Engineering 67: 841-867.""")
  t0 = time()


  E = 1*phun("PA");
  nu = 0.499;
  rho= 1*phun("KG/M^3");
  a=1*phun("M"); b=a; h= a;
  n1=8 # How many element edges per side?
  na= n1; nb= n1; nh =n1;
  neigvs=20                   # how many eigenvalues
  omega_shift=(0.1*2*pi)^2;

  fens,fes = H8block(a,b,h, na,nb,nh)

  # Make the region
  MR = DeforModelRed3D
  material = MatDeforElastIso(MR, rho, E, nu, 0.0)
  region1 = FDataDict("femm"=>FEMMDeforLinearMSH8(MR, IntegDomain(fes, GaussRule(3,2)),
    material), "femm_mass"=>FEMMDeforLinearMSH8(MR, IntegDomain(fes, GaussRule(3,3)),
    material))

  # Make model data
  modeldata =  FDataDict(
    "fens"=> fens, "regions"=>  [region1],
    "omega_shift"=>omega_shift, "neigvs"=>neigvs)

  # Solve
  modeldata = FinEtools.AlgoDeforLinearModule.modal(modeldata)

  fs = modeldata["omega"]/(2*pi)
  # println("Eigenvalues: $fs [Hz]")
  @test norm(fs-[1.92866e-7, 2.07497e-7, 2.16105e-7, 2.31656e-7, 2.35711e-7, 2.53067e-7, 0.266016, 0.266016, 0.364001, 0.364001, 0.364001, 0.366888, 0.366888, 0.366888, 0.415044, 0.415044, 0.41703, 0.467364, 0.467364, 0.467364]) < 0.0001

  modeldata["postprocessing"] = FDataDict("file"=>"unit_cube_mode",
    "mode"=>10)
  modeldata=FinEtools.AlgoDeforLinearModule.exportmode(modeldata)
  # @async run(`"paraview.exe" $(modeldata["postprocessing"]["file"]*"1.vtk")`)
  try rm(modeldata["postprocessing"]["exported"][1]["file"]); catch end

  true

end
end
using .mmunitmmccubemmvibrationmmms
mmunitmmccubemmvibrationmmms.test()

module mmtwistedbeamisomm
using FinEtools
using Test
using FinEtools.AlgoDeforLinearModule
import LinearAlgebra: norm, cholesky, cross
import Statistics: mean
function test()
  E = 0.29e8;
  nu = 0.22;
  W = 1.1;
  L = 12.;
  t =  0.32;
  nl = 2; nt = 1; nw = 1; ref = 4;
  p =   1/W/t;
  #  Loading in the Z direction
  loadv = [0;0;p]; dir = 3; uex = 0.005424534868469; # Harder: 5.424e-3;
  #   Loading in the Y direction
  #loadv = [0;p;0]; dir = 2; uex = 0.001753248285256; # Harder: 1.754e-3;
  tolerance  = t/1000;

  fens,fes  = H20block(L,W,t, nl*ref,nw*ref,nt*ref)

  # Reshape into a twisted beam shape
  for i = 1:count(fens)
    a = fens.xyz[i,1]/L*(pi/2); y = fens.xyz[i,2]-(W/2); z = fens.xyz[i,3]-(t/2);
    fens.xyz[i,:] = [fens.xyz[i,1],y*cos(a)-z*sin(a),y*sin(a)+z*cos(a)];
  end

  # Clamped end of the beam
  l1  = selectnode(fens; box = [0 0 -100*W 100*W -100*W 100*W], inflate  =  tolerance)
  e1 = FDataDict("node_list"=>l1, "component"=>1, "displacement"=>0.0)
  e2 = FDataDict("node_list"=>l1, "component"=>2, "displacement"=>0.0)
  e3 = FDataDict("node_list"=>l1, "component"=>3, "displacement"=>0.0)

  # Traction on the opposite edge
  boundaryfes  =   meshboundary(fes);
  Toplist   = selectelem(fens,boundaryfes, box =  [L L -100*W 100*W -100*W 100*W], inflate =   tolerance);
  el1femm  = FEMMBase(IntegDomain(subset(boundaryfes,Toplist), GaussRule(2, 2)))
  flux1 = FDataDict("femm"=>el1femm, "traction_vector"=>loadv)


  # Make the region
  MR = DeforModelRed3D
  material = MatDeforElastIso(MR, 00.0, E, nu, 0.0)
  region1 = FDataDict("femm"=>FEMMDeforLinear(MR, IntegDomain(fes, GaussRule(3,2)),
            material))

  # Make model data
  modeldata =  FDataDict(
  "fens"=> fens, "regions"=>  [region1],
  "essential_bcs"=>[e1, e2, e3], "traction_bcs"=>  [flux1])

  # Call the solver
  modeldata = AlgoDeforLinearModule.linearstatics(modeldata)
  geom = modeldata["geom"]
  u = modeldata["u"]

  # Extract the solution
  nl = selectnode(fens, box = [L L -100*W 100*W -100*W 100*W],inflate = tolerance);
  theutip = mean(u.values[nl,:], dims = 1)
  # println("displacement  = $(theutip[dir]) as compared to converged $uex")
    # println("normalized displacement  = $(theutip[dir]/uex*100) %")
    @test abs(theutip[dir]/uex*100-99.85504856450584) < 1.0e-6

  # Write out mesh with displacements
  modeldata["postprocessing"] = FDataDict("file"=>"twisted_beam")
  modeldata = AlgoDeforLinearModule.exportdeformation(modeldata)
  try rm(modeldata["postprocessing"]["exported"][1]["file"]); catch end

  # Write out mesh with stresses
  modeldata["postprocessing"] = FDataDict("file"=>"twisted_beam",
  "quantity"=> :Cauchy, "component"=> :xy)
  modeldata = AlgoDeforLinearModule.exportstress(modeldata)
  try rm(modeldata["postprocessing"]["exported"][1]["file"]); catch end

  # Write out mesh with stresses
  modeldata["postprocessing"] = FDataDict("file"=>"twisted_beam",
  "quantity"=> :Cauchy, "component"=> :xz)
  modeldata = AlgoDeforLinearModule.exportstress(modeldata)
  try rm(modeldata["postprocessing"]["exported"][1]["file"]); catch end

  # Write out mesh with von Mises stresses
  modeldata["postprocessing"] = FDataDict("file"=>"twisted_beam",
  "quantity"=> :vm)
  modeldata = AlgoDeforLinearModule.exportstress(modeldata)
  try rm(modeldata["postprocessing"]["exported"][1]["file"]); catch end

  # Write out mesh with von Mises stresses, elementwise
  modeldata["postprocessing"] = FDataDict("file"=>"twisted_beam-ew",
  "quantity"=> :vm)
  modeldata = AlgoDeforLinearModule.exportstresselementwise(modeldata)
  vm  = modeldata["postprocessing"]["exported"][1]["field"]
  # println("extremes of von Mises: $([minimum(vm.values),   maximum(vm.values)])")
  @test norm([minimum(vm.values),   maximum(vm.values)] - [6.94796, 451.904]) < 1.e-3
  try rm(modeldata["postprocessing"]["exported"][1]["file"]); catch end

  # Write out mesh with von Mises stresses, elementwise
  modeldata["postprocessing"] = FDataDict("file"=>"twisted_beam-ew",
  "quantity"=> :Cauchy, "component"=> :xz)
  modeldata = AlgoDeforLinearModule.exportstresselementwise(modeldata)
  try rm(modeldata["postprocessing"]["exported"][1]["file"]); catch end

  # Write out mesh with principal stresses, elementwise
  modeldata["postprocessing"] = FDataDict("file"=>"twisted_beam-principal-1-ew",
  "quantity"=> :princCauchy, "component"=> 1)
  modeldata = AlgoDeforLinearModule.exportstresselementwise(modeldata)
  ps  = modeldata["postprocessing"]["exported"][1]["field"]
  # println("extremes of first principal stress: $([minimum(ps.values),   maximum(ps.values)])")
  @test norm([minimum(ps.values),   maximum(ps.values)] - [0.493918, 459.106]) < 1.e-3
  try rm(modeldata["postprocessing"]["exported"][1]["file"]); catch end

  # Write out mesh with principal stresses, elementwise
  modeldata["postprocessing"] = FDataDict("file"=>"twisted_beam-principal-3-ew",
  "quantity"=> :princCauchy, "component"=> 3)
  modeldata = AlgoDeforLinearModule.exportstresselementwise(modeldata)
  ps  = modeldata["postprocessing"]["exported"][1]["field"]
  # println("extremes of third principal stress: $([minimum(ps.values),   maximum(ps.values)])")
  @test norm([minimum(ps.values),   maximum(ps.values)] - [-459.106, -0.493918]) < 1.e-3
  try rm(modeldata["postprocessing"]["exported"][1]["file"]); catch end


  # Write out mesh with pressure, elementwise
  modeldata["postprocessing"] = FDataDict("file"=>"twisted_beam-press-ew",
  "quantity"=> :pressure, "component"=> 1)
  modeldata = AlgoDeforLinearModule.exportstresselementwise(modeldata)
  ps  = modeldata["postprocessing"]["exported"][1]["field"]
  # println("extremes of pressure: $([minimum(ps.values),   maximum(ps.values)])")
  @test norm([minimum(ps.values),   maximum(ps.values)] - [-160.396, 160.396]) < 1.e-3
  try rm(modeldata["postprocessing"]["exported"][1]["file"]); catch end

end
end
using .mmtwistedbeamisomm
mmtwistedbeamisomm.test()

module mmtwistedbeamoorthomm
using FinEtools
using Test
import LinearAlgebra: norm, cholesky, cross
using FinEtools.AlgoDeforLinearModule
import Statistics: mean
function test()
  E = 0.29e8;
  nu = 0.22;
  W = 1.1;
  L = 12.;
  t =  0.32;
  nl = 2; nt = 1; nw = 1; ref = 4;
  p =   1/W/t;
  #  Loading in the Z direction
  loadv = [0;0;p]; dir = 3; uex = 0.005424534868469; # Harder: 5.424e-3;
  #   Loading in the Y direction
  #loadv = [0;p;0]; dir = 2; uex = 0.001753248285256; # Harder: 1.754e-3;
  tolerance  = t/1000;

  fens,fes  = H20block(L,W,t, nl*ref,nw*ref,nt*ref)

  # Reshape into a twisted beam shape
  for i = 1:count(fens)
    a = fens.xyz[i,1]/L*(pi/2); y = fens.xyz[i,2]-(W/2); z = fens.xyz[i,3]-(t/2);
    fens.xyz[i,:] = [fens.xyz[i,1],y*cos(a)-z*sin(a),y*sin(a)+z*cos(a)];
  end

  # Clamped end of the beam
  l1  = selectnode(fens; box = [0 0 -100*W 100*W -100*W 100*W], inflate  =  tolerance)
  e1 = FDataDict("node_list"=>l1, "component"=>1, "displacement"=>0.0)
  e2 = FDataDict("node_list"=>l1, "component"=>2, "displacement"=>0.0)
  e3 = FDataDict("node_list"=>l1, "component"=>3, "displacement"=>0.0)

  # Traction on the opposite edge
  boundaryfes  =   meshboundary(fes);
  Toplist   = selectelem(fens,boundaryfes, box =  [L L -100*W 100*W -100*W 100*W], inflate =   tolerance);
  el1femm  = FEMMBase(IntegDomain(subset(boundaryfes,Toplist), GaussRule(2, 2)))
  flux1 = FDataDict("femm"=>el1femm, "traction_vector"=>loadv)


  # Make the region
  MR = DeforModelRed3D
  material = MatDeforElastOrtho(MR, E, nu)
  region1 = FDataDict("femm"=>FEMMDeforLinear(MR, IntegDomain(fes, GaussRule(3,2)),
            material))

  # Make model data
  modeldata =  FDataDict(
  "fens"=> fens, "regions"=>  [region1],
  "essential_bcs"=>[e1, e2, e3], "traction_bcs"=>  [flux1])

  # Call the solver
  modeldata = AlgoDeforLinearModule.linearstatics(modeldata)
  geom = modeldata["geom"]
  u = modeldata["u"]

  # Extract the solution
  nl = selectnode(fens, box = [L L -100*W 100*W -100*W 100*W],inflate = tolerance);
  theutip = mean(u.values[nl,:], dims = 1)
  # println("displacement  = $(theutip[dir]) as compared to converged $uex")
    # println("normalized displacement  = $(theutip[dir]/uex*100) %")
    @test abs(theutip[dir]/uex*100-99.85504856450584) < 1.0e-6

  # Write out mesh with displacements
  modeldata["postprocessing"] = FDataDict("file"=>"twisted_beam")
  modeldata = AlgoDeforLinearModule.exportdeformation(modeldata)
  try rm(modeldata["postprocessing"]["exported"][1]["file"]); catch end

  # Write out mesh with stresses
  modeldata["postprocessing"] = FDataDict("file"=>"twisted_beam",
  "quantity"=> :Cauchy, "component"=> :xy)
  modeldata = AlgoDeforLinearModule.exportstress(modeldata)
  try rm(modeldata["postprocessing"]["exported"][1]["file"]); catch end

  # Write out mesh with stresses
  modeldata["postprocessing"] = FDataDict("file"=>"twisted_beam",
  "quantity"=> :Cauchy, "component"=> :xz)
  modeldata = AlgoDeforLinearModule.exportstress(modeldata)
  try rm(modeldata["postprocessing"]["exported"][1]["file"]); catch end

  # Write out mesh with von Mises stresses
  modeldata["postprocessing"] = FDataDict("file"=>"twisted_beam",
  "quantity"=> :vm)
  modeldata = AlgoDeforLinearModule.exportstress(modeldata)
  try rm(modeldata["postprocessing"]["exported"][1]["file"]); catch end

  # Write out mesh with von Mises stresses, elementwise
  modeldata["postprocessing"] = FDataDict("file"=>"twisted_beam-ew",
  "quantity"=> :vm)
  modeldata = AlgoDeforLinearModule.exportstresselementwise(modeldata)
  vm  = modeldata["postprocessing"]["exported"][1]["field"]
  # println("extremes of von Mises: $([minimum(vm.values),   maximum(vm.values)])")
  @test norm([minimum(vm.values),   maximum(vm.values)] - [6.94796, 451.904]) < 1.e-3
  try rm(modeldata["postprocessing"]["exported"][1]["file"]); catch end

  # Write out mesh with von Mises stresses, elementwise
  modeldata["postprocessing"] = FDataDict("file"=>"twisted_beam-ew",
  "quantity"=> :Cauchy, "component"=> :xz)
  modeldata = AlgoDeforLinearModule.exportstresselementwise(modeldata)
  try rm(modeldata["postprocessing"]["exported"][1]["file"]); catch end

  # Write out mesh with principal stresses, elementwise
  modeldata["postprocessing"] = FDataDict("file"=>"twisted_beam-principal-1-ew",
  "quantity"=> :princCauchy, "component"=> 1)
  modeldata = AlgoDeforLinearModule.exportstresselementwise(modeldata)
  ps  = modeldata["postprocessing"]["exported"][1]["field"]
  # println("extremes of first principal stress: $([minimum(ps.values),   maximum(ps.values)])")
  @test norm([minimum(ps.values),   maximum(ps.values)] - [0.493918, 459.106]) < 1.e-3
  try rm(modeldata["postprocessing"]["exported"][1]["file"]); catch end

  # Write out mesh with principal stresses, elementwise
  modeldata["postprocessing"] = FDataDict("file"=>"twisted_beam-principal-3-ew",
  "quantity"=> :princCauchy, "component"=> 3)
  modeldata = AlgoDeforLinearModule.exportstresselementwise(modeldata)
  ps  = modeldata["postprocessing"]["exported"][1]["field"]
  # println("extremes of third principal stress: $([minimum(ps.values),   maximum(ps.values)])")
  @test norm([minimum(ps.values),   maximum(ps.values)] - [-459.106, -0.493918]) < 1.e-3
  try rm(modeldata["postprocessing"]["exported"][1]["file"]); catch end


  # Write out mesh with pressure, elementwise
  modeldata["postprocessing"] = FDataDict("file"=>"twisted_beam-press-ew",
  "quantity"=> :pressure, "component"=> 1)
  modeldata = AlgoDeforLinearModule.exportstresselementwise(modeldata)
  ps  = modeldata["postprocessing"]["exported"][1]["field"]
  # println("extremes of pressure: $([minimum(ps.values),   maximum(ps.values)])")
  @test norm([minimum(ps.values),   maximum(ps.values)] - [-160.396, 160.396]) < 1.e-3
  try rm(modeldata["postprocessing"]["exported"][1]["file"]); catch end

end
end
using .mmtwistedbeamoorthomm
mmtwistedbeamoorthomm.test()

module muunit_cube_modes_exportmmm
using FinEtools
using FinEtools.MeshExportModule
using Test
import Arpack: eigs
import LinearAlgebra: norm, cholesky, cross
function test()


  # println("""
  # Vibration modes of unit cube  of almost incompressible material.
  #
  # This example EXPORTS the model to Abaqus.
  #
  # Reference: Puso MA, Solberg J (2006) A stabilized nodally integrated
  # tetrahedral. International Journal for Numerical Methods in
  # Engineering 67: 841-867.
  # """)
  t0 = time()


  E = 1*phun("PA");
  nu = 0.499;
  rho = 1*phun("KG/M^3");
  a = 1*phun("M"); b = a; h =  a;
  n1 = 5;# How many element edges per side?
  na =  n1; nb =  n1; nh  = n1;
  neigvs = 20                   # how many eigenvalues
  OmegaShift = (0.01*2*pi)^2;

  MR = DeforModelRed3D
  fens,fes  = H20block(a,b,h, na,nb,nh)

  geom = NodalField(fens.xyz)
  u = NodalField(zeros(size(fens.xyz,1),3)) # displacement field

  numberdofs!(u)

  material=MatDeforElastIso(MR, rho, E, nu, 0.0)

  femm = FEMMDeforLinear(MR, IntegDomain(fes, GaussRule(3,2)), material)

  K =stiffness(femm, geom, u)
  femm = FEMMDeforLinear(MR, IntegDomain(fes, GaussRule(3,3)), material)
  M =mass(femm, geom, u)
  d,v,nev,nconv = eigs(K+OmegaShift*M, M; nev=neigvs, which=:SM)
  d = broadcast(-, d, OmegaShift);
  fs = real(sqrt.(complex(d)))/(2*pi)
  # println("Eigenvalues: $fs [Hz]")
  @test norm(fs - [2.73674e-7, 3.00469e-7, 3.14245e-7, 3.19946e-7, 3.42634e-7, 3.56347e-7, 0.262723, 0.262723, 0.357791, 0.357791, 0.357791, 0.36088, 0.36088, 0.36088, 0.408199, 0.408397, 0.408397, 0.461756, 0.461756, 0.461756]) < 1.0e-3

  mode = 7
  scattersysvec!(u, v[:,mode])
  File =  "unit_cube_modes.vtk"
  vtkexportmesh(File, fens, fes; vectors=[("mode$mode", u.values)])
  # @async run(`"paraview.exe" $File`)
  try rm(File); catch end


  AE = AbaqusExporter("unit_cube_modes_h20");
  # AE.ios = STDOUT;
  HEADING(AE, "Vibration modes of unit cube  of almost incompressible material.");
  COMMENT(AE, "The  first six frequencies are rigid body modes.");
  COMMENT(AE, "The  first nonzero frequency (7) should be around 0.26 Hz");
  PART(AE, "part1");
  END_PART(AE);
  ASSEMBLY(AE, "ASSEM1");
  INSTANCE(AE, "INSTNC1", "PART1");
  NODE(AE, fens.xyz);
  COMMENT(AE, "The hybrid form of the serendipity hexahedron is chosen because");
  COMMENT(AE, "the material is  nearly incompressible.");
  ELEMENT(AE, "c3d20rh", "AllElements", 1, connasarray(fes))
  ORIENTATION(AE, "GlobalOrientation", vec([1. 0 0]), vec([0 1. 0]));
  SOLID_SECTION(AE, "elasticity", "GlobalOrientation", "AllElements");
  END_INSTANCE(AE);
  END_ASSEMBLY(AE);
  MATERIAL(AE, "elasticity")
  ELASTIC(AE, E, nu)
  DENSITY(AE, rho)
  STEP_FREQUENCY(AE, neigvs)
  END_STEP(AE)
  close(AE)
  nlines = 0
  open("unit_cube_modes_h20.inp") do f
    s = readlines(f)
    nlines = length(s)
  end
  @test nlines == 1409
  rm("unit_cube_modes_h20.inp")
end
end
using .muunit_cube_modes_exportmmm
muunit_cube_modes_exportmmm.test()

module mmpipemmPSmorthom
using FinEtools
using Test
import LinearAlgebra: norm, cholesky, cross, dot
import Statistics: mean

mutable struct MyIData
    c::FInt
    r::FFltVec
    s::FFltVec
end

function test()

## Thick pipe with internal pressure: plane strain
#

##
# Link to the  <matlab:edit('pub_thick_pipe_ps') m-file>.

## Description
##
# This is a simple modification of the full three-dimensional simulation of
# the tutorial pub_thick_pipe takes advantage of the plane-strain model
# reduction procedure.
##
# An infinitely long thick walled cylindrical pipe
# with inner boundary radius of 3 mm and outer boundary radius of 9 mm is
# subjected to an internal pressure of 1.0 MPa. A wedge   with thickness of
# 2 mm and a 90-degree angle sector is considered for the finite element
# analysis. The material properties are taken as  isotropic linear elastic
# with $E=1000$ MPa and $\nu=0.4999$ to represent nearly incompressible
# behavior. This problem has been proposed to by MacNeal and Harder as a
# test of an element's ability to represent the  response of a nearly
# incompressible material. The plane-strain condition is assumed in the
# axial direction of the pipe which together with the radial symmetry
# confines the material in all but the radial direction and therefore
# amplifies the numerical difficulties associated with the confinement of
# the nearly incompressible material.
##
# There is an analytical solution to this problem. Timoshenko and Goodier
# presented the original solution of Lame in their textbook. We are going
# to compare with  both the stress distribution (radial and hoop stresses)
# and the displacement of the inner  cylindrical surface.

##
#
# <html>
# <table border=0><tr><td>
# <img src="../docs/pub_thick_pipe_ps.png" width = "30#">
# </td></tr>
# <tr><td>Figure 1. Definition of the geometry of the internally pressurized thick pipe</td></tr>
# </table>
# </html>

##
# References:
#
# # Macneal RH, Harder RL (1985) A proposed standard set of problems to test
# finite element accuracy. Finite Elements in Analysis and Design 1: 3-20.
#
# # Timoshenko S. and Goodier J. N., Theory of Elasticity, McGraw-Hill, 2nd ed., 1951.

## Solution
#

##
# Internal radius of the pipe.
a = 3*phun("MM");
##
# External radius of the pipe.
b = 9*phun("MM");
##
# Thickness of the slice.
t = 2*phun("MM");

##
# Geometrical tolerance.
tolerance  =a/10000.;
##
# Young's modulus and Poisson's ratio.
E = 1000*phun("MEGA*PA");
nu = 0.499;
##
# Applied pressure on the internal surface.
press = 1.0*phun("MEGA*PA");

##
# Analytical solutions.   Radial stress:
radial_stress(r) =press*a.^2/(b^2-a^2).*(1-(b^2)./r.^2);
##
# Circumferential (hoop) stress:
hoop_stress(r)=press*a.^2/(b^2-a^2).*(1+(b^2)./r.^2);

##
# Radial displacement:
radial_displacement(r)=press*a^2*(1+nu)*(b^2+r.^2*(1-2*nu))/(E*(b^2-a^2).*r);;

##
# Therefore the radial displacement of the loaded surface will be:
urex = radial_displacement(a);


##
# The mesh parameters: The numbers of element edges axially,
# and through the thickness of the pipe wall (radially).

nc=3; nt=3;

##
# Note that the material object needs to be created with the proper
# model-dimension reduction in mind.  In this case that is the axial symmetry
# assumption.
MR = DeforModelRed2DStrain

# Create the mesh and initialize the geometry.  First we are going
# to construct the block of elements with the first coordinate
# corresponding to the angle, and the second
# coordinate is the thickness in the radial direction.
anglrrange = 90.0/180*pi;
fens,fes =  Q8block(anglrrange, b-a, nc, nt);

# Extract the boundary  and mark the finite elements on the
# interior surface.
bdryfes = meshboundary(fes);
bcl = selectelem(fens, bdryfes, box=[-Inf,Inf,0.,0.], inflate=tolerance);
internal_fenids = connectednodes(subset(bdryfes,bcl));
# Now  shape the block  into  the actual wedge piece of the pipe.
ayr=fens.xyz;
for i=1:count(fens)
    angl=ayr[i,1]; r=a+ayr[i,2];
    fens.xyz[i,:] = [r*sin(angl),(r*cos(angl))];
end

# now we create the geometry and displacement fields
geom = NodalField(fens.xyz)
u = NodalField(zeros(size(fens.xyz,1),2)) # displacement field

# The symmetry boundary condition  is specified by selecting nodes
# on the plane x=0.
l1 = selectnode(fens; box=[0.0 0.0 -Inf Inf], inflate = tolerance)
setebc!(u, l1, true, 1, 0.0)
# The second symmetry boundary condition is specified by selecting
# nodes on the plane y=0.
l1 = selectnode(fens; box=[-Inf Inf 0.0 0.0], inflate = tolerance)
setebc!(u, l1, true, 2, 0.0)

applyebc!(u)
numberdofs!(u)

# The traction boundary condition is applied in the radial
# direction.

el1femm =  FEMMBase(IntegDomain(subset(bdryfes,bcl), GaussRule(1, 3)))
function pressureloading!(forceout::FFltVec, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)
  copyto!(forceout, XYZ/norm(XYZ)*press)
  return forceout
end
fi = ForceIntensity(FFlt, 2, pressureloading!); # pressure normal to the internal cylindrical surface
F2 = distribloads(el1femm, geom, u, fi, 2);

# Property and material
material = MatDeforElastOrtho(MR, E, nu)

femm = FEMMDeforLinear(MR, IntegDomain(fes, GaussRule(2, 2)), material)

K =stiffness(femm, geom, u)
#K=cholesky(K)
U=  K\(F2)
scattersysvec!(u, U[:])

# Transfer the solution of the displacement to the nodes on the
# internal cylindrical surface and convert to
# cylindrical-coordinate displacements there.
uv = u.values[internal_fenids,:];
ur = zeros(FFlt,length(internal_fenids));

for j = 1:length(internal_fenids)
    n = fens.xyz[internal_fenids[j],:];
    n = n'/norm(n);# normal to the cylindrical internal surface
    ur[j] = dot(vec(uv[j,:]),vec(n));
end

# Report the  relative displacement on the internal surface:
# println("(Approximate/true displacement) at the internal surface: $( mean(ur)/urex*100  ) %")
@test abs(mean(ur)/urex*100 - 100) < 0.1

# Produce a plot of the radial stress component in the cylindrical
# coordinate system. Note that this is the usual representation of
# stress using nodal stress field.

fld = fieldfromintegpoints(femm, geom, u, :Cauchy, 1)


File =  "thick_pipe_sigmax.vtk"
vtkexportmesh(File, fens, fes; scalars=[("sigmax", fld.values)])
try rm(File); catch end

# Produce a plot of the solution components in the cylindrical
# coordinate system.
# Plot the analytical solution.

function inspector(idat::MyIData, elnum, conn, xe,  out,  xq)
  function outputRm(c)
    theNormal=c;
    r=norm(theNormal);# distance from the axis of symmetry
    theNormal =theNormal/r;# compute the unit normal vector
    e1p=[theNormal';0.];# local cylind. coordinate  basis vectors
    e3p=[0.,0.,1.]';# this one points along the axis of the cylinder
    e2p=cross(vec(e3p),vec(e1p));# this one is along the hoop direction
    R= [vec(e1p) vec(e2p) vec(e3p)];# transformation matrix for the stress
    return R
  end
  Rm=outputRm(xq)
  tm=zeros(FFlt,3,3)
  stress4vto3x3t!(tm, out);# stress in global XYZ
  tpm = Rm'*tm*Rm;#  stress matrix in cylindrical coordinates
  sp=zeros(FFlt,6)
  stress3x3tto6v!(sp, tpm);# stress vector in cylindr. coord.
  push!(idat.r,norm(xq))
  push!(idat.s,sp[idat.c])
  return idat
end

idat = MyIData(1, FFltVec[], FFltVec[])
idat = inspectintegpoints(femm, geom, u, collect(1:count(fes)),
 inspector, idat, :Cauchy)
# show(idat)

@test norm(idat.s - [-7.44858e5, -3.55143e5, -7.44858e5, -3.55143e5, -2.19961e5,
-1.08612e5, -2.19961e5, -1.08612e5, -58910.9, -12517.6, -58910.9, -12517.6,
-7.44858e5, -3.55143e5, -7.44858e5, -3.55143e5, -2.19961e5, -1.08612e5,
-2.19961e5, -1.08612e5, -58910.9, -12517.6, -58910.9, -12517.6, -7.44858e5,
-3.55143e5, -7.44858e5, -3.55143e5, -2.19961e5, -1.08612e5, -2.19961e5,
 -1.08612e5, -58910.9, -12517.6, -58910.9, -12517.6])/1.0e5 < 1.e-3
# using Plots
# plotly()
#
# # Plot the analytical solution.
# r = linearspace(a,b,100);
# plot(r, radial_stress(r))
# # Plot the computed  integration-point data
# scatter!(idat.r, idat.s, m=:circle, color=:red)
# gui()
end
end
using .mmpipemmPSmorthom
mmpipemmPSmorthom.test()

module scratch1_06092017_ortho
using FinEtools
using Test

mutable struct MyIData
  c::FInt
  r::FFltVec
  s::FFltVec
end

function test()

  # println("Thick pipe with internal pressure: axially symmetric model")
  #=
  This is a simple modification of the full three-dimensional simulation of
  the tutorial pub_thick_pipe that implements the axially-symmetric model
  reduction procedure.

  An infinitely long thick walled cylindrical pipe
  with inner boundary radius of 3 mm and outer boundary radius of 9 mm is
  subjected to an internal pressure of 1.0 MPa. A wedge   with thickness of
  2 mm and a 90-degree angle sector is considered for the finite element
  analysis. The material properties are taken as  isotropic linear elastic
  with $E=1000$ MPa and $\nu=0.4999$ to represent nearly incompressible
  behavior. This problem has been proposed to by MacNeal and Harder as a
  test of an element's ability to represent the  response of a nearly
  incompressible material. The plane-strain condition is assumed in the
  axial direction of the pipe which together with the radial symmetry
  confines the material in all but the radial direction and therefore
  amplifies the numerical difficulties associated with the confinement of
  the nearly incompressible material.

  There is an analytical solution to this problem. Timoshenko and Goodier
  presented the original solution of Lame in their textbook. We are going
  to compare with  both the stress distribution (radial and hoop stresses)
  and the displacement of the inner  cylindrical surface.

  References:
  - Macneal RH, Harder RL (1985) A proposed standard set of problems to test
  finite element accuracy. Finite Elements in Analysis and Design 1: 3-20.
  - Timoshenko S. and Goodier J. N., Theory of Elasticity, McGraw-Hill, 2nd ed., 1951.

  =#

  # Internal radius of the pipe.
  a=3*phun("MM");
  ##
  # External radius of the pipe.
  b = 9*phun("MM");
  ##
  # Thickness of the slice.
  t = 2*phun("MM");

  ##
  # Geometrical tolerance.
  tolerance   = a/10000.;
  ##
  # Young's modulus and Poisson's ratio.
  E = 1000*phun("MEGA*PA");
  nu = 0.499;
  ##
  # Applied pressure on the internal surface.
  press =   1.0*phun("MEGA*PA");

  ##
  # Analytical solutions.   Radial stress:
  radial_stress(r)  = press*a.^2/(b^2-a^2).*(1-(b^2)./r.^2);
  ##
  # Circumferential (hoop) stress:
  hoop_stress(r) = press*a.^2/(b^2-a^2).*(1+(b^2)./r.^2);

  ##
  # Radial displacement:
  radial_displacement(r) = press*a^2*(1+nu)*(b^2+r.^2*(1-2*nu))/(E*(b^2-a^2).*r);;

  ##
  # Therefore the radial displacement of the loaded surface will be:
  urex  =  radial_displacement(a);


  ##
  # The mesh parameters: The numbers of element edges axially,
  # and through the thickness of the pipe wall (radially).

  na = 1; nt = 10;

  ##
  # Note that the material object needs to be created with the proper
  # model-dimension reduction in effect.  In this case that is the axial symmetry
  # assumption.
  MR = DeforModelRed2DAxisymm
  axisymmetric = true

  # Create the mesh and initialize the geometry.  First we are going
  # to construct the block of elements with the first coordinate
  # corresponding to the thickness in the radial direction, and the second
  # coordinate is the thickness in the axial direction.
  fens,fes  =   Q8block(b-a, t, nt, na);

  # Extract the boundary  and mark the finite elements on the
  # interior surface.
  bdryfes = meshboundary(fes);

  bcl = selectelem(fens, bdryfes, box=[0.,0.,-Inf,Inf], inflate=tolerance);
  internal_fenids= connectednodes(subset(bdryfes,bcl));
  # Now  shape the block  into  the actual wedge piece of the pipe.
  for i=1:count(fens)
    fens.xyz[i,:] = fens.xyz[i,:] + [a; 0.0];
  end

  # now we create the geometry and displacement fields
  geom = NodalField(fens.xyz)
  u = NodalField(zeros(size(fens.xyz,1),2)) # displacement field

  # The plane-strain condition in the axial direction  is specified by selecting nodes
  # on the plane y=0 and y=t.
  l1 = selectnode(fens; box=[-Inf Inf 0.0 0.0], inflate = tolerance)
  setebc!(u, l1, true, 2, 0.0)
  l1 = selectnode(fens; box=[-Inf Inf t t], inflate = tolerance)
  setebc!(u, l1, true, 2, 0.0)

  applyebc!(u)
  numberdofs!(u)

  # The traction boundary condition is applied in the radial
  # direction.

  el1femm =  FEMMBase(IntegDomain(subset(bdryfes,bcl), GaussRule(1, 3), axisymmetric))
  fi = ForceIntensity([press; 0.0]);
  F2= distribloads(el1femm, geom, u, fi, 2);

  # Property and material
  material = MatDeforElastOrtho(MR,  E, nu)

  femm = FEMMDeforLinear(MR, IntegDomain(fes, GaussRule(2, 2), axisymmetric), material)

  K = stiffness(femm, geom, u)
  #K=cholesky(K)
  U =  K\(F2)
  scattersysvec!(u,U[:])

  # Transfer the solution of the displacement to the nodes on the
  # internal cylindrical surface and convert to
  # cylindrical-coordinate displacements there.
  uv=u.values[internal_fenids,:]
  # Report the  relative displacement on the internal surface:
  # println("(Approximate/true displacement) at the internal surface: $( mean(uv[:,1])/urex*100  ) %")

  # Produce a plot of the radial stress component in the cylindrical
  # coordinate system. Note that this is the usual representation of
  # stress using nodal stress field.

  fld = fieldfromintegpoints(femm, geom, u, :Cauchy, 1)


  # File =  "thick_pipe_sigmax.vtk"
  # vtkexportmesh(File, fens, fes; scalars=[("sigmax", fld.values)])

  # Produce a plot of the solution components in the cylindrical
  # coordinate system.

  function inspector(idat::MyIData, elnum, conn, xe,  out,  xq)
    push!(idat.r, xq[1])
    push!(idat.s, out[idat.c])
    return idat
  end

  idat = MyIData(1, FInt[], FInt[])
  idat = inspectintegpoints(femm, geom, u, collect(1:count(fes)),
  inspector, idat, :Cauchy)

  # using Plots
  # plotly()
  #
  # # Plot the analytical solution.
  # r = linearspace(a,b,100);
  # plot(r, radial_stress(r))
  # # Plot the computed  integration-point data
  # plot!(idat.r, idat.s, m=:circle, color=:red)
  # gui()

  @test  abs(idat.r[1]-0.003126794919243112)<1.0e-9
  @test  abs(idat.s[1]- -910911.9777008593)<1.0e-2

  ## Discussion
  #
  ##
  # The axially symmetric model is clearly very effective
  # computationally, as the size is much reduced compared to the 3-D
  # model.  In conjunction with  uniform or selective reduced integration
  # it can be very accurate as well.
end
end
using .scratch1_06092017_ortho
scratch1_06092017_ortho.test()

module mmLE11Q8mm
using FinEtools
using Test
import LinearAlgebra: norm, cholesky, cross
function test()

  # NAFEMS LE11 benchmark with Q8 elements.
  # # This is a test recommended by the National Agency for Finite Element
  # # Methods and Standards (U.K.): Test LE11 from NAFEMS Publication TNSB,
  # # Rev. 3, “The Standard NAFEMS Benchmarks,” October 1990.
  # #
  # # Target solution: Direct stress,   =  –105 MPa at point A.
  #function  LE11NAFEMS()
  # Parameters:
  Ea =  210000*phun("MEGA*Pa")
  nua =  0.3;
  alphaa = 2.3e-4;              # thermal expansion coefficient
  sigmaA = -105*phun("MEGA*Pa")
  nref =  1;                        # how many times should we refine the mesh?
  X = [1.     0.;#A
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
  tolerance  = 1.e-6*phun("M")
  ##
  # Note that the material object needs to be created with the proper
  # model-dimension reduction in mind.  In this case that is the axial symmetry
  # assumption.
  MR  =  DeforModelRed2DAxisymm



  fens = FENodeSet(X);
  fes = FESetQ4([1 2 4 3; 3 4 6 5; 5 6 8 7; 7 8 10 9; 9 10 12 11; 11 12 14 13; 13 14 16 15]);
  for ref = 1:nref
    fens,fes = Q4refine(fens,fes);
    list = selectnode(fens,distance = 1.0+0.1/2^nref, from = [0. 0.], inflate = tolerance);
    fens.xyz[list,:] =  FinEtools.MeshUtilModule.ontosphere(fens.xyz[list,:],1.0);
  end
  fens,fes = Q4toQ8(fens,fes)
  list = selectnode(fens,distance = 1.0+0.1/2^nref, from = [0. 0.], inflate = tolerance);
  fens.xyz[list,:] =  FinEtools.MeshUtilModule.ontosphere(fens.xyz[list,:],1.0);

  #     File  =   "mesh.vtk"
  # vtkexportmesh(File, fens, fes)

  # now we create the geometry and displacement fields
  geom  =  NodalField(fens.xyz)
  u  =  NodalField(zeros(size(fens.xyz,1),2)) # displacement field

  # Apply EBC's
  l1 = selectnode(fens,box = [-Inf Inf 0 0],inflate = tolerance)
  setebc!(u, l1, true, 2, 00.0)
  l1 = selectnode(fens,box = [-Inf Inf 1.79  1.79],inflate = tolerance)
  setebc!(u, l1, true, 2, 00.0)
  applyebc!(u)
  numberdofs!(u)

  # Temperature field
  dT  = NodalField(reshape(fens.xyz[:,1]+fens.xyz[:,2],size(fens.xyz,1),1));


  # Property and material
  material = MatDeforElastIso(MR, 0.0, Ea, nua, alphaa)

  femm  =  FEMMDeforLinear(MR, IntegDomain(fes, GaussRule(2, 3), true), material)

  K  = stiffness(femm, geom, u)
  F  =  thermalstrainloads(femm, geom, u, dT)
  #K = cholesky(K)
  U =   K\F
  scattersysvec!(u, U[:])

  nA  = selectnode(fens,box = FFlt[1.0  1.0 0.0 0.0], inflate = tolerance);

  fld =  fieldfromintegpoints(femm, geom, u, dT, :Cauchy, 2)


  File  =   "LE11NAFEMS_Q8_sigmay.vtk"
  vtkexportmesh(File, fens, fes; scalars = [("sigmay", fld.values)],
  vectors = [("u", u.values)])
  # println("range of  sigmay = $((minimum(fld.values), maximum(fld.values)))")
  @test norm([minimum(fld.values), maximum(fld.values)] - [-1.443052182185006e8, -1.4106181545272605e7]) < 1.0e-2
  # @async run(`"paraview.exe" $File`)
  try rm(File); catch end

  sA  =  fld.values[nA]/phun("MEGA*Pa")
  sAn  =  fld.values[nA]/sigmaA
  # println("Stress at point A: $(sA) i. e.  $( sAn*100  )% of reference value")
  @test abs(sA[1] - (-93.8569)) < 1.0e-3

  fen2fe  = FENodeToFEMap(connasarray(fes), nnodes(geom))
  function inspector(idat, elnum, conn, xe,  out,  xq)
    println("loc = $(  xq  ) : $(  transpose(out)/phun("MEGA*Pa")  )")
    return idat
  end

  # inspectintegpoints(femm, geom, u, dT,  fen2fe.map[nA[1]],
  # inspector, []; quantity = :Cauchy)

end
end
using .mmLE11Q8mm
mmLE11Q8mm.test()

module mmLE11Q8mmortho
using FinEtools
using Test
import LinearAlgebra: norm, cholesky, cross
function test()

  # NAFEMS LE11 benchmark with Q8 elements.
  # # This is a test recommended by the National Agency for Finite Element
  # # Methods and Standards (U.K.): Test LE11 from NAFEMS Publication TNSB,
  # # Rev. 3, “The Standard NAFEMS Benchmarks,” October 1990.
  # #
  # # Target solution: Direct stress,   =  –105 MPa at point A.
  #function  LE11NAFEMS()
  # Parameters:
  Ea =  210000*phun("MEGA*Pa")
  nua =  0.3;
  alphaa = 2.3e-4;              # thermal expansion coefficient
  sigmaA = -105*phun("MEGA*Pa")
  nref =  1;                        # how many times should we refine the mesh?
  X = [1.     0.;#A
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
  tolerance  = 1.e-6*phun("M")
  ##
  # Note that the material object needs to be created with the proper
  # model-dimension reduction in mind.  In this case that is the axial symmetry
  # assumption.
  MR  =  DeforModelRed2DAxisymm



  fens = FENodeSet(X);
  fes = FESetQ4([1 2 4 3; 3 4 6 5; 5 6 8 7; 7 8 10 9; 9 10 12 11; 11 12 14 13; 13 14 16 15]);
  for ref = 1:nref
    fens,fes = Q4refine(fens,fes);
    list = selectnode(fens,distance = 1.0+0.1/2^nref, from = [0. 0.], inflate = tolerance);
    fens.xyz[list,:] =  FinEtools.MeshUtilModule.ontosphere(fens.xyz[list,:],1.0);
  end
  fens,fes = Q4toQ8(fens,fes)
  list = selectnode(fens,distance = 1.0+0.1/2^nref, from = [0. 0.], inflate = tolerance);
  fens.xyz[list,:] =  FinEtools.MeshUtilModule.ontosphere(fens.xyz[list,:],1.0);

  #     File  =   "mesh.vtk"
  # vtkexportmesh(File, fens, fes)

  # now we create the geometry and displacement fields
  geom  =  NodalField(fens.xyz)
  u  =  NodalField(zeros(size(fens.xyz,1),2)) # displacement field

  # Apply EBC's
  l1 = selectnode(fens,box = [-Inf Inf 0 0],inflate = tolerance)
  setebc!(u, l1, true, 2, 00.0)
  l1 = selectnode(fens,box = [-Inf Inf 1.79  1.79],inflate = tolerance)
  setebc!(u, l1, true, 2, 00.0)
  applyebc!(u)
  numberdofs!(u)

  # Temperature field
  dT  = NodalField(reshape(fens.xyz[:,1]+fens.xyz[:,2],size(fens.xyz,1),1));


  # Property and material
  material = MatDeforElastOrtho(MR, 0.0, Ea, nua, alphaa)
  # display(material )

  femm  =  FEMMDeforLinear(MR, IntegDomain(fes, GaussRule(2, 3), true), material)

  K  = stiffness(femm, geom, u)
  F  =  thermalstrainloads(femm, geom, u, dT)
  #K = cholesky(K)
  U =   K\F
  scattersysvec!(u, U[:])

  nA  = selectnode(fens,box = FFlt[1.0  1.0 0.0 0.0], inflate = tolerance);

  fld =  fieldfromintegpoints(femm, geom, u, dT, :Cauchy, 2)


  File  =   "LE11NAFEMS_Q8_sigmay.vtk"
  vtkexportmesh(File, fens, fes; scalars = [("sigmay", fld.values)],
  vectors = [("u", u.values)])
  # println("range of  sigmay = $((minimum(fld.values), maximum(fld.values)))")
  @test norm([minimum(fld.values), maximum(fld.values)] - [-1.443052182185006e8, -1.4106181545272605e7]) < 1.0e-2
  # @async run(`"paraview.exe" $File`)
  try rm(File); catch end

  sA  =  fld.values[nA]/phun("MEGA*Pa")
  sAn  =  fld.values[nA]/sigmaA
  # println("Stress at point A: $(sA) i. e.  $( sAn*100  )% of reference value")
  @test abs(sA[1] - (-93.8569)) < 1.0e-3

  fen2fe  = FENodeToFEMap(connasarray(fes), nnodes(geom))
  function inspector(idat, elnum, conn, xe,  out,  xq)
    println("loc = $(  xq  ) : $(  transpose(out)/phun("MEGA*Pa")  )")
    return idat
  end

  # inspectintegpoints(femm, geom, u, dT,  fen2fe.map[nA[1]],
  # inspector, []; quantity = :Cauchy)

end
end
using .mmLE11Q8mmortho
mmLE11Q8mmortho.test()

module mLE11Q8aximmm
using FinEtools
using Test
import LinearAlgebra: norm, cholesky, cross
function test()
  # NAFEMS LE11 benchmark with Q8 elements.
  # # This is a test recommended by the National Agency for Finite Element
  # # Methods and Standards (U.K.): Test LE11 from NAFEMS Publication TNSB,
  # # Rev. 3, “The Standard NAFEMS Benchmarks,” October 1990.
  # #
  # # Target solution: Direct stress,   =  –105 MPa at point A.
  #function  LE11NAFEMS()
  # Parameters:
  Ea =  210000*phun("MEGA*Pa")
  nua =  0.3;
  alphaa = 2.3e-4;              # thermal expansion coefficient
  sigmaA = -105*phun("MEGA*Pa")
  nref = 2;                        # how many times should we refine the mesh?
  X = [1.     0.;#A
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
  tolerance  = 1.e-6*phun("M")
  ##
  # Note that the material object needs to be created with the proper
  # model-dimension reduction in mind.  In this case that is the axial symmetry
  # assumption.
  MR  =  DeforModelRed2DAxisymm



  fens = FENodeSet(X);
  fes = FESetQ4([1 2 4 3; 3 4 6 5; 5 6 8 7; 7 8 10 9; 9 10 12 11; 11 12 14 13; 13 14 16 15]);
  for ref = 1:nref
    fens,fes = Q4refine(fens,fes);
    list = selectnode(fens,distance = 1.0+0.1/2^nref, from = [0. 0.], inflate = tolerance);
    fens.xyz[list,:] =  FinEtools.MeshUtilModule.ontosphere(fens.xyz[list,:],1.0);
  end
  fens,fes = Q4toQ8(fens,fes)
  list = selectnode(fens,distance = 1.0+0.1/2^nref, from = [0. 0.], inflate = tolerance);
  fens.xyz[list,:] =  FinEtools.MeshUtilModule.ontosphere(fens.xyz[list,:],1.0);

  #     File  =   "mesh.vtk"
  # vtkexportmesh(File, fens, fes)

  # now we create the geometry and displacement fields
  geom  =  NodalField(fens.xyz)
  u  =  NodalField(zeros(size(fens.xyz,1),2)) # displacement field

  # Apply EBC's
  lbottom = selectnode(fens,box = [-Inf Inf 0 0],inflate = tolerance)
  setebc!(u, lbottom, true, 2, 00.0)
  ltop = selectnode(fens,box = [-Inf Inf 1.79  1.79],inflate = tolerance)
  setebc!(u, ltop, true, 2, 00.0)
  applyebc!(u)
  numberdofs!(u)

  # Temperature field
  dT  = NodalField(reshape(fens.xyz[:,1]+fens.xyz[:,2],size(fens.xyz,1),1));


  # Property and material
  material = MatDeforElastIso(MR, 0.0, Ea, nua, alphaa)

  femm  =  FEMMDeforLinear(MR, IntegDomain(fes, GaussRule(2, 3), true), material)

  K  = stiffness(femm, geom, u)
  F  =  thermalstrainloads(femm, geom, u, dT)
  #K = cholesky(K)
  U =   K\F
  scattersysvec!(u, U[:])

  nA  = selectnode(fens,box = FFlt[1.0  1.0 0.0 0.0], inflate = tolerance);

  fld =  fieldfromintegpoints(femm, geom, u, dT, :Cauchy, 2)


  File  =   "LE11NAFEMS_Q8_sigmay.vtk"
  vtkexportmesh(File, fens, fes; scalars = [("sigmay", fld.values)],
  vectors = [("u", u.values)])
  # println("range of  sigmay = $((minimum(fld.values), maximum(fld.values)))")
  @test norm([minimum(fld.values), maximum(fld.values)] - [-1.6338426447540134e8, -4.961956343464769e6]) < 1.e-1
  # @async run(`"paraview.exe" $File`)
  try rm(File); catch end


  sA  =  fld.values[nA]/phun("MEGA*Pa")
  sAn  =  fld.values[nA]/sigmaA
  # println("Stress at point A: $(sA) i. e.  $( sAn*100  )% of reference value")

  fen2fe  = FENodeToFEMap(connasarray(fes), nnodes(geom))
  function inspector(idat, elnum, conn, xe,  out,  xq)
    # println("loc = $(  xq  ) : $(  transpose(out)/phun("MEGA*Pa")  )")
    return idat
  end

  inspectintegpoints(femm, geom, u, dT,  fen2fe.map[nA[1]],
  inspector, []; quantity = :Cauchy)

  fld =  fieldfromintegpoints(femm, geom, u, dT, :Pressure, 1)
  File  =   "LE11NAFEMS_Q8_pressure.vtk"
  vtkexportmesh(File, fens, fes; scalars = [("pressure", fld.values)],
  vectors = [("u", u.values)])
  # println("range of  pressure = $((minimum(fld.values), maximum(fld.values)))")
  @test norm([minimum(fld.values), maximum(fld.values)] - [-1.1881819144904878e7, 7.555030948761216e7]) < 1.e-1
  # @async run(`"paraview.exe" $File`)
  try rm(File); catch end

  fld =  fieldfromintegpoints(femm, geom, u, dT, :vm, 1)
  File  =   "LE11NAFEMS_Q8_vm.vtk"
  vtkexportmesh(File, fens, fes; scalars = [("pressure", fld.values)],
  vectors = [("u", u.values)])
  # println("range of von Mises = $((minimum(fld.values), maximum(fld.values)))")
  @test norm([minimum(fld.values), maximum(fld.values)] - [3.221370699152578e7, 1.4437590830351183e8]) < 1.e-1
  # @async run(`"paraview.exe" $File`)
  try rm(File); catch end

  AE = AbaqusExporter("LE11NAFEMS_Q8_export_stress");
  HEADING(AE, "NAFEMS LE11 benchmark with Q8 elements.");
  COMMENT(AE, "sigmaA = -105 MPa ");
  PART(AE, "part1");
  END_PART(AE);
  ASSEMBLY(AE, "ASSEM1");
  INSTANCE(AE, "INSTNC1", "PART1");
  NODE(AE, fens.xyz);
  ELEMENT(AE, "cax8", "AllElements", 1, connasarray(fes))
  NSET_NSET(AE, "ltop", ltop)
  NSET_NSET(AE, "lbottom", lbottom)
  ORIENTATION(AE, "GlobalOrientation", vec([1. 0 0]), vec([0 1. 0]));
  SOLID_SECTION(AE, "elasticity", "GlobalOrientation", "AllElements");
  SURFACE_SECTION(AE, "TractionElements")
  END_INSTANCE(AE);
  END_ASSEMBLY(AE);
  MATERIAL(AE, "elasticity")
  ELASTIC(AE, Ea, nua)
  EXPANSION(AE, alphaa)
  STEP_PERTURBATION_STATIC(AE)
  BOUNDARY(AE, "ASSEM1.INSTNC1.ltop", 2)
  BOUNDARY(AE, "ASSEM1.INSTNC1.lbottom", 2)
  TEMPERATURE(AE, "ASSEM1.INSTNC1.", collect(1:count(fens)), vec(dT.values))
  END_STEP(AE)
  close(AE)
  nlines = 0
  open(AE.filename) do f
    s = readlines(f)
    nlines = length(s)
  end
  @test nlines == 963
  try rm(AE.filename); catch end

end
end
using .mLE11Q8aximmm
mLE11Q8aximmm.test()


module mLE11Q8aximorthom
using FinEtools
using Test
import LinearAlgebra: norm, cholesky, cross
function test()
  # NAFEMS LE11 benchmark with Q8 elements.
  # # This is a test recommended by the National Agency for Finite Element
  # # Methods and Standards (U.K.): Test LE11 from NAFEMS Publication TNSB,
  # # Rev. 3, “The Standard NAFEMS Benchmarks,” October 1990.
  # #
  # # Target solution: Direct stress,   =  –105 MPa at point A.
  #function  LE11NAFEMS()
  # Parameters:
  Ea =  210000*phun("MEGA*Pa")
  nua =  0.3;
  alphaa = 2.3e-4;              # thermal expansion coefficient
  sigmaA = -105*phun("MEGA*Pa")
  nref = 2;                        # how many times should we refine the mesh?
  X = [1.     0.;#A
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
  tolerance  = 1.e-6*phun("M")
  ##
  # Note that the material object needs to be created with the proper
  # model-dimension reduction in mind.  In this case that is the axial symmetry
  # assumption.
  MR  =  DeforModelRed2DAxisymm



  fens = FENodeSet(X);
  fes = FESetQ4([1 2 4 3; 3 4 6 5; 5 6 8 7; 7 8 10 9; 9 10 12 11; 11 12 14 13; 13 14 16 15]);
  for ref = 1:nref
    fens,fes = Q4refine(fens,fes);
    list = selectnode(fens,distance = 1.0+0.1/2^nref, from = [0. 0.], inflate = tolerance);
    fens.xyz[list,:] =  FinEtools.MeshUtilModule.ontosphere(fens.xyz[list,:],1.0);
  end
  fens,fes = Q4toQ8(fens,fes)
  list = selectnode(fens,distance = 1.0+0.1/2^nref, from = [0. 0.], inflate = tolerance);
  fens.xyz[list,:] =  FinEtools.MeshUtilModule.ontosphere(fens.xyz[list,:],1.0);

  #     File  =   "mesh.vtk"
  # vtkexportmesh(File, fens, fes)

  # now we create the geometry and displacement fields
  geom  =  NodalField(fens.xyz)
  u  =  NodalField(zeros(size(fens.xyz,1),2)) # displacement field

  # Apply EBC's
  lbottom = selectnode(fens,box = [-Inf Inf 0 0],inflate = tolerance)
  setebc!(u, lbottom, true, 2, 00.0)
  ltop = selectnode(fens,box = [-Inf Inf 1.79  1.79],inflate = tolerance)
  setebc!(u, ltop, true, 2, 00.0)
  applyebc!(u)
  numberdofs!(u)

  # Temperature field
  dT  = NodalField(reshape(fens.xyz[:,1]+fens.xyz[:,2],size(fens.xyz,1),1));


  # Property and material
  material = MatDeforElastOrtho(MR, 0.0, Ea, nua, alphaa)

  femm  =  FEMMDeforLinear(MR, IntegDomain(fes, GaussRule(2, 3), true), material)

  K  = stiffness(femm, geom, u)
  F  =  thermalstrainloads(femm, geom, u, dT)
  #K = cholesky(K)
  U =   K\F
  scattersysvec!(u, U[:])

  nA  = selectnode(fens,box = FFlt[1.0  1.0 0.0 0.0], inflate = tolerance);

  fld =  fieldfromintegpoints(femm, geom, u, dT, :Cauchy, 2)


  File  =   "LE11NAFEMS_Q8_sigmay.vtk"
  vtkexportmesh(File, fens, fes; scalars = [("sigmay", fld.values)],
  vectors = [("u", u.values)])
  # println("range of  sigmay = $((minimum(fld.values), maximum(fld.values)))")
  @test norm([minimum(fld.values), maximum(fld.values)] - [-1.6338426447540134e8, -4.961956343464769e6]) < 1.e-1
  # @async run(`"paraview.exe" $File`)
  try rm(File); catch end


  sA  =  fld.values[nA]/phun("MEGA*Pa")
  sAn  =  fld.values[nA]/sigmaA
  # println("Stress at point A: $(sA) i. e.  $( sAn*100  )% of reference value")

  fen2fe  = FENodeToFEMap(connasarray(fes), nnodes(geom))
  function inspector(idat, elnum, conn, xe,  out,  xq)
    # println("loc = $(  xq  ) : $(  transpose(out)/phun("MEGA*Pa")  )")
    return idat
  end

  inspectintegpoints(femm, geom, u, dT,  fen2fe.map[nA[1]],
  inspector, []; quantity = :Cauchy)

  fld =  fieldfromintegpoints(femm, geom, u, dT, :Pressure, 1)
  File  =   "LE11NAFEMS_Q8_pressure.vtk"
  vtkexportmesh(File, fens, fes; scalars = [("pressure", fld.values)],
  vectors = [("u", u.values)])
  # println("range of  pressure = $((minimum(fld.values), maximum(fld.values)))")
  @test norm([minimum(fld.values), maximum(fld.values)] - [-1.1881819144904878e7, 7.555030948761216e7]) < 1.e-1
  # @async run(`"paraview.exe" $File`)
  try rm(File); catch end

  fld =  fieldfromintegpoints(femm, geom, u, dT, :vm, 1)
  File  =   "LE11NAFEMS_Q8_vm.vtk"
  vtkexportmesh(File, fens, fes; scalars = [("pressure", fld.values)],
  vectors = [("u", u.values)])
  # println("range of von Mises = $((minimum(fld.values), maximum(fld.values)))")
  @test norm([minimum(fld.values), maximum(fld.values)] - [3.221370699152578e7, 1.4437590830351183e8]) < 1.e-1
  # @async run(`"paraview.exe" $File`)
  try rm(File); catch end

  AE = AbaqusExporter("LE11NAFEMS_Q8_export_stress");
  HEADING(AE, "NAFEMS LE11 benchmark with Q8 elements.");
  COMMENT(AE, "sigmaA = -105 MPa ");
  PART(AE, "part1");
  END_PART(AE);
  ASSEMBLY(AE, "ASSEM1");
  INSTANCE(AE, "INSTNC1", "PART1");
  NODE(AE, fens.xyz);
  ELEMENT(AE, "cax8", "AllElements", 1, connasarray(fes))
  NSET_NSET(AE, "ltop", ltop)
  NSET_NSET(AE, "lbottom", lbottom)
  ORIENTATION(AE, "GlobalOrientation", vec([1. 0 0]), vec([0 1. 0]));
  SOLID_SECTION(AE, "elasticity", "GlobalOrientation", "AllElements");
  SURFACE_SECTION(AE, "TractionElements")
  END_INSTANCE(AE);
  END_ASSEMBLY(AE);
  MATERIAL(AE, "elasticity")
  ELASTIC(AE, Ea, nua)
  EXPANSION(AE, alphaa)
  STEP_PERTURBATION_STATIC(AE)
  BOUNDARY(AE, "ASSEM1.INSTNC1.ltop", 2)
  BOUNDARY(AE, "ASSEM1.INSTNC1.lbottom", 2)
  TEMPERATURE(AE, "ASSEM1.INSTNC1.", collect(1:count(fens)), vec(dT.values))
  END_STEP(AE)
  close(AE)
  nlines = 0
  open(AE.filename) do f
    s = readlines(f)
    nlines = length(s)
  end
  @test nlines == 963
  try rm(AE.filename); catch end

end
end
using .mLE11Q8aximorthom
mLE11Q8aximorthom.test()

module mmmCookmmstrainmmisommm
using FinEtools
using FinEtools.AlgoDeforLinearModule
using FinEtools.MeshExportModule
using Test
import LinearAlgebra: norm, cholesky, cross
function test()


  E = 1.0;
  nu = 1.0/3;
  width = 48.0; height = 44.0; thickness  = 1.0;
  free_height  = 16.0;
  Mid_edge  = [48.0, 52.0];# Location of tracked  deflection
  magn = 1.0/(free_height*thickness);# Density of applied load
  convutip = 23.97;
  n = 30;# number of elements per side
  tolerance = minimum([width, height])/n/1000.;#Geometrical tolerance

  fens,fes = T3block(width, height, n, n)

  # Reshape into a trapezoidal panel
  for i = 1:count(fens)
      fens.xyz[i,2] = fens.xyz[i,2]+(fens.xyz[i,1]/width)*(height -fens.xyz[i,2]/height*(height-free_height));
  end

  # Clamped edge of the membrane
  l1 = selectnode(fens; box=[0.,0.,-Inf, Inf], inflate = tolerance)
  ess1 = FDataDict("displacement"=>  0.0, "component"=> 1, "node_list"=>l1)
  ess2 = FDataDict("displacement"=>  0.0, "component"=> 2, "node_list"=>l1)

  # Traction on the opposite edge
  boundaryfes =  meshboundary(fes);
  Toplist  = selectelem(fens, boundaryfes, box= [width, width, -Inf, Inf ], inflate=  tolerance);
  el1femm = FEMMBase(IntegDomain(subset(boundaryfes, Toplist), GaussRule(1, 2), thickness))
  flux1 = FDataDict("traction_vector"=>[0.0,+magn],
      "femm"=>el1femm
      )

  # Make the region
  MR = DeforModelRed2DStrain
  material = MatDeforElastIso(MR,  0.0, E, nu, 0.0)
  region1 = FDataDict("femm"=>FEMMDeforLinear(MR,
      IntegDomain(fes, TriRule(1), thickness), material))

  modeldata = FDataDict("fens"=>fens,
   "regions"=>[region1],
   "essential_bcs"=>[ess1, ess2],
   "traction_bcs"=>[flux1]
   )

  # Call the solver
  modeldata = AlgoDeforLinearModule.linearstatics(modeldata)

  u = modeldata["u"]
  geom = modeldata["geom"]

  # Extract the solution
  nl = selectnode(fens, box=[Mid_edge[1],Mid_edge[1],Mid_edge[2],Mid_edge[2]],
            inflate=tolerance);
  theutip = u.values[nl,:]
  # println("displacement =$(theutip[2]) as compared to converged $convutip")
  @test abs(theutip[2]-21.35955642390279) < 1.0e-3

  modeldata["postprocessing"] = FDataDict("file"=>"cookstress-ew",
     "quantity"=>:Cauchy, "component"=>:xy)
  modeldata = AlgoDeforLinearModule.exportstresselementwise(modeldata)
  fld = modeldata["postprocessing"]["exported"][1]["field"]
  # println("range of Cauchy_xy = $((minimum(fld.values), maximum(fld.values)))")
  @test norm([minimum(fld.values), maximum(fld.values)] - [-0.07028875155067116, 0.1301698279821655]) < 1.0e-5
  # File = modeldata["postprocessing"]["exported"][1]["file"]
  # @async run(`"paraview.exe" $File`)
  try rm(modeldata["postprocessing"]["exported"][1]["file"]); catch end

  modeldata["postprocessing"] = FDataDict("file"=>"cookstress-ew-vm",
     "quantity"=>:vm, "component"=>1)
  modeldata = AlgoDeforLinearModule.exportstresselementwise(modeldata)
  fld = modeldata["postprocessing"]["exported"][1]["field"]
  # println("range of vm = $((minimum(fld.values), maximum(fld.values)))")
  @test norm([minimum(fld.values), maximum(fld.values)] - [0.007503804468283987, 0.33798754356331173]) < 1.0e-5
  # File = modeldata["postprocessing"]["exported"][1]["file"]
  # @async run(`"paraview.exe" $File`)
  try rm(modeldata["postprocessing"]["exported"][1]["file"]); catch end

  modeldata["postprocessing"] = FDataDict("file"=>"cookstress-ew-pressure",
     "quantity"=>:pressure, "component"=>1)
  modeldata = AlgoDeforLinearModule.exportstresselementwise(modeldata)
  fld = modeldata["postprocessing"]["exported"][1]["field"]
  # println("range of pressure = $((minimum(fld.values), maximum(fld.values)))")
  @test norm([minimum(fld.values), maximum(fld.values)] - [-0.11777749217431245, 0.23457099031101358]) < 1.0e-5
  # File = modeldata["postprocessing"]["exported"][1]["file"]
  # @async run(`"paraview.exe" $File`)
  try rm(modeldata["postprocessing"]["exported"][1]["file"]); catch end

  modeldata["postprocessing"] = FDataDict("file"=>"cookstress-ew-princ1",
     "quantity"=>:princCauchy, "component"=>1)
  modeldata = AlgoDeforLinearModule.exportstresselementwise(modeldata)
  fld = modeldata["postprocessing"]["exported"][1]["field"]
  # println("range of princCauchy Max = $((minimum(fld.values), maximum(fld.values)))")
  @test norm([minimum(fld.values), maximum(fld.values)] - [-0.16098150217425994, 0.24838761904231466]) < 1.0e-5
  # File = modeldata["postprocessing"]["exported"][1]["file"]
  # @async run(`"paraview.exe" $File`)
  try rm(modeldata["postprocessing"]["exported"][1]["file"]); catch end

  modeldata["postprocessing"] = FDataDict("file"=>"cookstress-ew-princ3",
     "quantity"=>:princCauchy, "component"=>3)
  modeldata = AlgoDeforLinearModule.exportstresselementwise(modeldata)
  fld = modeldata["postprocessing"]["exported"][1]["field"]
  # println("range of princCauchy Min = $((minimum(fld.values), maximum(fld.values)))")
  @test norm([minimum(fld.values), maximum(fld.values)] - [-0.4523049370669106, 0.02980811337406548]) < 1.0e-5
  # File = modeldata["postprocessing"]["exported"][1]["file"]
  # @async run(`"paraview.exe" $File`)
  try rm(modeldata["postprocessing"]["exported"][1]["file"]); catch end

  AE = AbaqusExporter("Cookstress_algo_stress");
  HEADING(AE, "Cook trapezoidal panel, plane stress");
  COMMENT(AE, "Converged free mid-edge displacement = 23.97");
  PART(AE, "part1");
  END_PART(AE);
  ASSEMBLY(AE, "ASSEM1");
  INSTANCE(AE, "INSTNC1", "PART1");
  NODE(AE, fens.xyz);
  COMMENT(AE, "We are assuming three node triangles in plane-stress");
  COMMENT(AE, "CPE3 are pretty poor-accuracy elements, but here we don't care about it.");
  @test nodesperelem(modeldata["regions"][1]["femm"].integdomain.fes) == 3
  ELEMENT(AE, "CPE3", "AllElements", connasarray(modeldata["regions"][1]["femm"].integdomain.fes))
  NSET_NSET(AE, "clamped", modeldata["essential_bcs"][1]["node_list"])
  ORIENTATION(AE, "GlobalOrientation", vec([1. 0 0]), vec([0 1. 0]));
  SOLID_SECTION(AE, "elasticity", "GlobalOrientation", "AllElements", thickness);
  END_INSTANCE(AE);
  END_ASSEMBLY(AE);
  MATERIAL(AE, "elasticity")
  ELASTIC(AE, E, nu)
  STEP_PERTURBATION_STATIC(AE)
  BOUNDARY(AE, "ASSEM1.INSTNC1.clamped", 1)
  BOUNDARY(AE, "ASSEM1.INSTNC1.clamped", 2)
  bfes = modeldata["traction_bcs"][1]["femm"].integdomain.fes
  COMMENT(AE, "Concentrated loads: we are assuming that the elements on the boundary");
  COMMENT(AE, "have two nodes each and also that they are the same length.");
  COMMENT(AE, "Then the concentrated loads below will be correctly lumped.");
  nl = connectednodes(bfes)
  F = zeros(count(modeldata["fens"]))
  for ix = 1:count(bfes)
    for jx = 1:2
      F[bfes.conn[ix][jx]] += 1.0/n/2/thickness
    end
  end
  for ixxxx = 1:length(F)
    if F[ixxxx] != 0.0
      CLOAD(AE, "ASSEM1.INSTNC1.$(ixxxx)", 2, F[ixxxx])
    end
  end
  END_STEP(AE)
  close(AE)

  nlines = 0
  open(AE.filename) do f
    s = readlines(f)
    nlines = length(s)
  end
  @test nlines == 2886
  rm(AE.filename)
end
end
using .mmmCookmmstrainmmisommm
mmmCookmmstrainmmisommm.test()

module mmmCookmmstrainmorthommm
using FinEtools
using FinEtools.AlgoDeforLinearModule
using FinEtools.MeshExportModule
using Test
import LinearAlgebra: norm, cholesky, cross
function test()


  E = 1.0;
  nu = 1.0/3;
  width = 48.0; height = 44.0; thickness  = 1.0;
  free_height  = 16.0;
  Mid_edge  = [48.0, 52.0];# Location of tracked  deflection
  magn = 1.0/(free_height*thickness);# Density of applied load
  convutip = 23.97;
  n = 30;# number of elements per side
  tolerance = minimum([width, height])/n/1000.;#Geometrical tolerance

  fens,fes = T3block(width, height, n, n)

  # Reshape into a trapezoidal panel
  for i = 1:count(fens)
      fens.xyz[i,2] = fens.xyz[i,2]+(fens.xyz[i,1]/width)*(height -fens.xyz[i,2]/height*(height-free_height));
  end

  # Clamped edge of the membrane
  l1 = selectnode(fens; box=[0.,0.,-Inf, Inf], inflate = tolerance)
  ess1 = FDataDict("displacement"=>  0.0, "component"=> 1, "node_list"=>l1)
  ess2 = FDataDict("displacement"=>  0.0, "component"=> 2, "node_list"=>l1)

  # Traction on the opposite edge
  boundaryfes =  meshboundary(fes);
  Toplist  = selectelem(fens, boundaryfes, box= [width, width, -Inf, Inf ], inflate=  tolerance);
  el1femm = FEMMBase(IntegDomain(subset(boundaryfes, Toplist), GaussRule(1, 2), thickness))
  flux1 = FDataDict("traction_vector"=>[0.0,+magn],
      "femm"=>el1femm
      )

  # Make the region
  MR = DeforModelRed2DStrain
  material = MatDeforElastOrtho(MR,  0.0, E, nu, 0.0)
  region1 = FDataDict("femm"=>FEMMDeforLinear(MR,
      IntegDomain(fes, TriRule(1), thickness), material))

  modeldata = FDataDict("fens"=>fens,
   "regions"=>[region1],
   "essential_bcs"=>[ess1, ess2],
   "traction_bcs"=>[flux1]
   )

  # Call the solver
  modeldata = AlgoDeforLinearModule.linearstatics(modeldata)

  u = modeldata["u"]
  geom = modeldata["geom"]

  # Extract the solution
  nl = selectnode(fens, box=[Mid_edge[1],Mid_edge[1],Mid_edge[2],Mid_edge[2]],
            inflate=tolerance);
  theutip = u.values[nl,:]
  # println("displacement =$(theutip[2]) as compared to converged $convutip")
  @test abs(theutip[2]-21.35955642390279) < 1.0e-3

  modeldata["postprocessing"] = FDataDict("file"=>"cookstress-ew",
     "quantity"=>:Cauchy, "component"=>:xy)
  modeldata = AlgoDeforLinearModule.exportstresselementwise(modeldata)
  fld = modeldata["postprocessing"]["exported"][1]["field"]
  # println("range of Cauchy_xy = $((minimum(fld.values), maximum(fld.values)))")
  @test norm([minimum(fld.values), maximum(fld.values)] - [-0.07028875155067116, 0.1301698279821655]) < 1.0e-5
  # File = modeldata["postprocessing"]["exported"][1]["file"]
  # @async run(`"paraview.exe" $File`)
  try rm(modeldata["postprocessing"]["exported"][1]["file"]); catch end

  modeldata["postprocessing"] = FDataDict("file"=>"cookstress-ew-vm",
     "quantity"=>:vm, "component"=>1)
  modeldata = AlgoDeforLinearModule.exportstresselementwise(modeldata)
  fld = modeldata["postprocessing"]["exported"][1]["field"]
  # println("range of vm = $((minimum(fld.values), maximum(fld.values)))")
  @test norm([minimum(fld.values), maximum(fld.values)] - [0.007503804468283987, 0.33798754356331173]) < 1.0e-5
  # File = modeldata["postprocessing"]["exported"][1]["file"]
  # @async run(`"paraview.exe" $File`)
  try rm(modeldata["postprocessing"]["exported"][1]["file"]); catch end

  modeldata["postprocessing"] = FDataDict("file"=>"cookstress-ew-pressure",
     "quantity"=>:pressure, "component"=>1)
  modeldata = AlgoDeforLinearModule.exportstresselementwise(modeldata)
  fld = modeldata["postprocessing"]["exported"][1]["field"]
  # println("range of pressure = $((minimum(fld.values), maximum(fld.values)))")
  @test norm([minimum(fld.values), maximum(fld.values)] - [-0.11777749217431245, 0.23457099031101358]) < 1.0e-5
  # File = modeldata["postprocessing"]["exported"][1]["file"]
  # @async run(`"paraview.exe" $File`)
  try rm(modeldata["postprocessing"]["exported"][1]["file"]); catch end

  modeldata["postprocessing"] = FDataDict("file"=>"cookstress-ew-princ1",
     "quantity"=>:princCauchy, "component"=>1)
  modeldata = AlgoDeforLinearModule.exportstresselementwise(modeldata)
  fld = modeldata["postprocessing"]["exported"][1]["field"]
  # println("range of princCauchy Max = $((minimum(fld.values), maximum(fld.values)))")
  @test norm([minimum(fld.values), maximum(fld.values)] - [-0.16098150217425994, 0.24838761904231466]) < 1.0e-5
  # File = modeldata["postprocessing"]["exported"][1]["file"]
  # @async run(`"paraview.exe" $File`)
  try rm(modeldata["postprocessing"]["exported"][1]["file"]); catch end

  modeldata["postprocessing"] = FDataDict("file"=>"cookstress-ew-princ3",
     "quantity"=>:princCauchy, "component"=>3)
  modeldata = AlgoDeforLinearModule.exportstresselementwise(modeldata)
  fld = modeldata["postprocessing"]["exported"][1]["field"]
  # println("range of princCauchy Min = $((minimum(fld.values), maximum(fld.values)))")
  @test norm([minimum(fld.values), maximum(fld.values)] - [-0.4523049370669106, 0.02980811337406548]) < 1.0e-5
  # File = modeldata["postprocessing"]["exported"][1]["file"]
  # @async run(`"paraview.exe" $File`)
  try rm(modeldata["postprocessing"]["exported"][1]["file"]); catch end

  AE = AbaqusExporter("Cookstress_algo_stress");
  HEADING(AE, "Cook trapezoidal panel, plane stress");
  COMMENT(AE, "Converged free mid-edge displacement = 23.97");
  PART(AE, "part1");
  END_PART(AE);
  ASSEMBLY(AE, "ASSEM1");
  INSTANCE(AE, "INSTNC1", "PART1");
  NODE(AE, fens.xyz);
  COMMENT(AE, "We are assuming three node triangles in plane-stress");
  COMMENT(AE, "CPE3 are pretty poor-accuracy elements, but here we don't care about it.");
@test  nodesperelem(modeldata["regions"][1]["femm"].integdomain.fes) == 3
  ELEMENT(AE, "CPE3", "AllElements", connasarray(modeldata["regions"][1]["femm"].integdomain.fes))
  NSET_NSET(AE, "clamped", modeldata["essential_bcs"][1]["node_list"])
  ORIENTATION(AE, "GlobalOrientation", vec([1. 0 0]), vec([0 1. 0]));
  SOLID_SECTION(AE, "elasticity", "GlobalOrientation", "AllElements", thickness);
  END_INSTANCE(AE);
  END_ASSEMBLY(AE);
  MATERIAL(AE, "elasticity")
  ELASTIC(AE, E, nu)
  STEP_PERTURBATION_STATIC(AE)
  BOUNDARY(AE, "ASSEM1.INSTNC1.clamped", 1)
  BOUNDARY(AE, "ASSEM1.INSTNC1.clamped", 2)
  bfes = modeldata["traction_bcs"][1]["femm"].integdomain.fes
  COMMENT(AE, "Concentrated loads: we are assuming that the elements on the boundary");
  COMMENT(AE, "have two nodes each and also that they are the same length.");
  COMMENT(AE, "Then the concentrated loads below will be correctly lumped.");
  nl = connectednodes(bfes)
  F = zeros(count(modeldata["fens"]))
  for ix = 1:count(bfes)
    for jx = 1:2
      F[bfes.conn[ix][jx]] += 1.0/n/2/thickness
    end
  end
  for ixxxx = 1:length(F)
    if F[ixxxx] != 0.0
      CLOAD(AE, "ASSEM1.INSTNC1.$(ixxxx)", 2, F[ixxxx])
    end
  end
  END_STEP(AE)
  close(AE)

  nlines = 0
  open(AE.filename) do f
    s = readlines(f)
    nlines = length(s)
  end
  @test nlines == 2886
  rm(AE.filename)
end
end
using .mmmCookmmstrainmorthommm
mmmCookmmstrainmorthommm.test()


module mmmCookmstressmorthommm
using FinEtools
using FinEtools.AlgoDeforLinearModule
using FinEtools.MeshExportModule
using Test
import LinearAlgebra: norm, cholesky, cross
function test()


  E = 1.0;
  nu = 1.0/3;
  width = 48.0; height = 44.0; thickness  = 1.0;
  free_height  = 16.0;
  Mid_edge  = [48.0, 52.0];# Location of tracked  deflection
  magn = 1.0/(free_height*thickness);# Density of applied load
  convutip = 23.97;
  n = 30;# number of elements per side
  tolerance = minimum([width, height])/n/1000.;#Geometrical tolerance

  fens,fes = T3block(width, height, n, n)

  # Reshape into a trapezoidal panel
  for i = 1:count(fens)
      fens.xyz[i,2] = fens.xyz[i,2]+(fens.xyz[i,1]/width)*(height -fens.xyz[i,2]/height*(height-free_height));
  end

  # Clamped edge of the membrane
  l1 = selectnode(fens; box=[0.,0.,-Inf, Inf], inflate = tolerance)
  ess1 = FDataDict("displacement"=>  0.0, "component"=> 1, "node_list"=>l1)
  ess2 = FDataDict("displacement"=>  0.0, "component"=> 2, "node_list"=>l1)

  # Traction on the opposite edge
  boundaryfes =  meshboundary(fes);
  Toplist  = selectelem(fens, boundaryfes, box= [width, width, -Inf, Inf ], inflate=  tolerance);
  el1femm = FEMMBase(IntegDomain(subset(boundaryfes, Toplist), GaussRule(1, 2), thickness))
  flux1 = FDataDict("traction_vector"=>[0.0,+magn],
      "femm"=>el1femm
      )

  # Make the region
  MR = DeforModelRed2DStress
  material = MatDeforElastOrtho(MR,  0.0, E, nu, 0.0)
  region1 = FDataDict("femm"=>FEMMDeforLinear(MR,
      IntegDomain(fes, TriRule(1), thickness), material))

  modeldata = FDataDict("fens"=>fens,
   "regions"=>[region1],
   "essential_bcs"=>[ess1, ess2],
   "traction_bcs"=>[flux1]
   )

  # Call the solver
  modeldata = AlgoDeforLinearModule.linearstatics(modeldata)

  u = modeldata["u"]
  geom = modeldata["geom"]

  # Extract the solution
  nl = selectnode(fens, box=[Mid_edge[1],Mid_edge[1],Mid_edge[2],Mid_edge[2]],
            inflate=tolerance);
  theutip = u.values[nl,:]
  # println("displacement =$(theutip[2]) as compared to converged $convutip")
  @test abs(theutip[2]-23.79623002934245) < 1.0e-3

  modeldata["postprocessing"] = FDataDict("file"=>"cookstress-ew",
     "quantity"=>:Cauchy, "component"=>:xy)
  modeldata = AlgoDeforLinearModule.exportstresselementwise(modeldata)
  fld = modeldata["postprocessing"]["exported"][1]["field"]
  # println("range of Cauchy_xy = $((minimum(fld.values), maximum(fld.values)))")
  @test norm([minimum(fld.values), maximum(fld.values)] - [-0.07275125002229098, 0.1309644027374448]) < 1.0e-5
  # File = modeldata["postprocessing"]["exported"][1]["file"]
  # @async run(`"paraview.exe" $File`)
  try rm(modeldata["postprocessing"]["exported"][1]["file"]); catch end

  modeldata["postprocessing"] = FDataDict("file"=>"cookstress-ew-vm",
     "quantity"=>:vm, "component"=>1)
  modeldata = AlgoDeforLinearModule.exportstresselementwise(modeldata)
  fld = modeldata["postprocessing"]["exported"][1]["field"]
  # println("range of vm = $((minimum(fld.values), maximum(fld.values)))")
  @test norm([minimum(fld.values), maximum(fld.values)] - [0.014136291979824203, 0.4181381117075297]) < 1.0e-5
  # File = modeldata["postprocessing"]["exported"][1]["file"]
  # @async run(`"paraview.exe" $File`)
  try rm(modeldata["postprocessing"]["exported"][1]["file"]); catch end

  modeldata["postprocessing"] = FDataDict("file"=>"cookstress-ew-pressure",
     "quantity"=>:pressure, "component"=>1)
  modeldata = AlgoDeforLinearModule.exportstresselementwise(modeldata)
  fld = modeldata["postprocessing"]["exported"][1]["field"]
  # println("range of pressure = $((minimum(fld.values), maximum(fld.values)))")
  # display(norm([minimum(fld.values), maximum(fld.values)] - [-0.12996180159464202, 0.2436183072544255]))
  @test norm([minimum(fld.values), maximum(fld.values)] - [-0.08664120106309468, 0.16241220483628366]) < 1.0e-5
  # File = modeldata["postprocessing"]["exported"][1]["file"]
  # @async run(`"paraview.exe" $File`)
  try rm(modeldata["postprocessing"]["exported"][1]["file"]); catch end

  modeldata["postprocessing"] = FDataDict("file"=>"cookstress-ew-princ1",
     "quantity"=>:princCauchy, "component"=>1)
  modeldata = AlgoDeforLinearModule.exportstresselementwise(modeldata)
  fld = modeldata["postprocessing"]["exported"][1]["field"]
  # println("range of princCauchy Max = $((minimum(fld.values), maximum(fld.values)))")
  @test norm([minimum(fld.values), maximum(fld.values)] - [-0.10280794467415574, 0.24831137383158813]) < 1.0e-5
  # File = modeldata["postprocessing"]["exported"][1]["file"]
  # @async run(`"paraview.exe" $File`)
  try rm(modeldata["postprocessing"]["exported"][1]["file"]); catch end

  modeldata["postprocessing"] = FDataDict("file"=>"cookstress-ew-princ3",
     "quantity"=>:princCauchy, "component"=>2)
  modeldata = AlgoDeforLinearModule.exportstresselementwise(modeldata)
  fld = modeldata["postprocessing"]["exported"][1]["field"]
  # println("range of princCauchy Min = $((minimum(fld.values), maximum(fld.values)))")
  @test norm([minimum(fld.values), maximum(fld.values)] - [-0.4398236425818378, 0.022575693063719465]) < 1.0e-5
  # File = modeldata["postprocessing"]["exported"][1]["file"]
  # @async run(`"paraview.exe" $File`)
  try rm(modeldata["postprocessing"]["exported"][1]["file"]); catch end

  AE = AbaqusExporter("Cookstress_algo_stress");
  HEADING(AE, "Cook trapezoidal panel, plane stress");
  COMMENT(AE, "Converged free mid-edge displacement = 23.97");
  PART(AE, "part1");
  END_PART(AE);
  ASSEMBLY(AE, "ASSEM1");
  INSTANCE(AE, "INSTNC1", "PART1");
  NODE(AE, fens.xyz);
  COMMENT(AE, "We are assuming three node triangles in plane-stress");
  COMMENT(AE, "CPS3 are pretty poor-accuracy elements, but here we don't care about it.");
@test nodesperelem(modeldata["regions"][1]["femm"].integdomain.fes) == 3
  ELEMENT(AE, "CPS3", "AllElements", connasarray(modeldata["regions"][1]["femm"].integdomain.fes))
  NSET_NSET(AE, "clamped", modeldata["essential_bcs"][1]["node_list"])
  ORIENTATION(AE, "GlobalOrientation", vec([1. 0 0]), vec([0 1. 0]));
  SOLID_SECTION(AE, "elasticity", "GlobalOrientation", "AllElements", thickness);
  END_INSTANCE(AE);
  END_ASSEMBLY(AE);
  MATERIAL(AE, "elasticity")
  ELASTIC(AE, E, nu)
  STEP_PERTURBATION_STATIC(AE)
  BOUNDARY(AE, "ASSEM1.INSTNC1.clamped", 1)
  BOUNDARY(AE, "ASSEM1.INSTNC1.clamped", 2)
  bfes = modeldata["traction_bcs"][1]["femm"].integdomain.fes
  COMMENT(AE, "Concentrated loads: we are assuming that the elements on the boundary");
  COMMENT(AE, "have two nodes each and also that they are the same length.");
  COMMENT(AE, "Then the concentrated loads below will be correctly lumped.");
  nl = connectednodes(bfes)
  F = zeros(count(modeldata["fens"]))
  for ix = 1:count(bfes)
    for jx = 1:2
      F[bfes.conn[ix][jx]] += 1.0/n/2/thickness
    end
  end
  for ixxxx = 1:length(F)
    if F[ixxxx] != 0.0
      CLOAD(AE, "ASSEM1.INSTNC1.$(ixxxx)", 2, F[ixxxx])
    end
  end
  END_STEP(AE)
  close(AE)

  nlines = 0
  open(AE.filename) do f
    s = readlines(f)
    nlines = length(s)
  end
  @test nlines == 2886
  rm(AE.filename)
end
end
using .mmmCookmstressmorthommm
mmmCookmstressmorthommm.test()

module mmmCookmstressisommm
using FinEtools
using FinEtools.AlgoDeforLinearModule
using FinEtools.MeshExportModule
using Test
import LinearAlgebra: norm, cholesky, cross
function test()


  E = 1.0;
  nu = 1.0/3;
  width = 48.0; height = 44.0; thickness  = 1.0;
  free_height  = 16.0;
  Mid_edge  = [48.0, 52.0];# Location of tracked  deflection
  magn = 1.0/(free_height*thickness);# Density of applied load
  convutip = 23.97;
  n = 30;# number of elements per side
  tolerance = minimum([width, height])/n/1000.;#Geometrical tolerance

  fens,fes = T3block(width, height, n, n)

  # Reshape into a trapezoidal panel
  for i = 1:count(fens)
      fens.xyz[i,2] = fens.xyz[i,2]+(fens.xyz[i,1]/width)*(height -fens.xyz[i,2]/height*(height-free_height));
  end

  # Clamped edge of the membrane
  l1 = selectnode(fens; box=[0.,0.,-Inf, Inf], inflate = tolerance)
  ess1 = FDataDict("displacement"=>  0.0, "component"=> 1, "node_list"=>l1)
  ess2 = FDataDict("displacement"=>  0.0, "component"=> 2, "node_list"=>l1)

  # Traction on the opposite edge
  boundaryfes =  meshboundary(fes);
  Toplist  = selectelem(fens, boundaryfes, box= [width, width, -Inf, Inf ], inflate=  tolerance);
  el1femm = FEMMBase(IntegDomain(subset(boundaryfes, Toplist), GaussRule(1, 2), thickness))
  flux1 = FDataDict("traction_vector"=>[0.0,+magn],
      "femm"=>el1femm
      )

  # Make the region
  MR = DeforModelRed2DStress
  material = MatDeforElastIso(MR,  0.0, E, nu, 0.0)
  region1 = FDataDict("femm"=>FEMMDeforLinear(MR,
      IntegDomain(fes, TriRule(1), thickness), material))

  modeldata = FDataDict("fens"=>fens,
   "regions"=>[region1],
   "essential_bcs"=>[ess1, ess2],
   "traction_bcs"=>[flux1]
   )

  # Call the solver
  modeldata = AlgoDeforLinearModule.linearstatics(modeldata)

  u = modeldata["u"]
  geom = modeldata["geom"]

  # Extract the solution
  nl = selectnode(fens, box=[Mid_edge[1],Mid_edge[1],Mid_edge[2],Mid_edge[2]],
            inflate=tolerance);
  theutip = u.values[nl,:]
  # println("displacement =$(theutip[2]) as compared to converged $convutip")
  @test abs(theutip[2]-23.79623002934245) < 1.0e-3

  modeldata["postprocessing"] = FDataDict("file"=>"cookstress-ew",
     "quantity"=>:Cauchy, "component"=>:xy)
  modeldata = AlgoDeforLinearModule.exportstresselementwise(modeldata)
  fld = modeldata["postprocessing"]["exported"][1]["field"]
  # println("range of Cauchy_xy = $((minimum(fld.values), maximum(fld.values)))")
  @test norm([minimum(fld.values), maximum(fld.values)] - [-0.07275125002229098, 0.1309644027374448]) < 1.0e-5
  # File = modeldata["postprocessing"]["exported"][1]["file"]
  # @async run(`"paraview.exe" $File`)
  try rm(modeldata["postprocessing"]["exported"][1]["file"]); catch end

  modeldata["postprocessing"] = FDataDict("file"=>"cookstress-ew-vm",
     "quantity"=>:vm, "component"=>1)
  modeldata = AlgoDeforLinearModule.exportstresselementwise(modeldata)
  fld = modeldata["postprocessing"]["exported"][1]["field"]
  # println("range of vm = $((minimum(fld.values), maximum(fld.values)))")
  @test norm([minimum(fld.values), maximum(fld.values)] - [0.014136291979824203, 0.4181381117075297]) < 1.0e-5
  # File = modeldata["postprocessing"]["exported"][1]["file"]
  # @async run(`"paraview.exe" $File`)
  try rm(modeldata["postprocessing"]["exported"][1]["file"]); catch end

  modeldata["postprocessing"] = FDataDict("file"=>"cookstress-ew-pressure",
     "quantity"=>:pressure, "component"=>1)
  modeldata = AlgoDeforLinearModule.exportstresselementwise(modeldata)
  fld = modeldata["postprocessing"]["exported"][1]["field"]
  # println("range of pressure = $((minimum(fld.values), maximum(fld.values)))")
  @test norm([minimum(fld.values), maximum(fld.values)] - [-0.08664120106309468, 0.16241220483628366]) < 1.0e-5
  # File = modeldata["postprocessing"]["exported"][1]["file"]
  # @async run(`"paraview.exe" $File`)
  try rm(modeldata["postprocessing"]["exported"][1]["file"]); catch end

  modeldata["postprocessing"] = FDataDict("file"=>"cookstress-ew-princ1",
     "quantity"=>:princCauchy, "component"=>1)
  modeldata = AlgoDeforLinearModule.exportstresselementwise(modeldata)
  fld = modeldata["postprocessing"]["exported"][1]["field"]
  # println("range of princCauchy Max = $((minimum(fld.values), maximum(fld.values)))")
  @test norm([minimum(fld.values), maximum(fld.values)] - [-0.10280794467415574, 0.24831137383158813]) < 1.0e-5
  # File = modeldata["postprocessing"]["exported"][1]["file"]
  # @async run(`"paraview.exe" $File`)
  try rm(modeldata["postprocessing"]["exported"][1]["file"]); catch end

  modeldata["postprocessing"] = FDataDict("file"=>"cookstress-ew-princ3",
     "quantity"=>:princCauchy, "component"=>2)
  modeldata = AlgoDeforLinearModule.exportstresselementwise(modeldata)
  fld = modeldata["postprocessing"]["exported"][1]["field"]
  # println("range of princCauchy Min = $((minimum(fld.values), maximum(fld.values)))")
  @test norm([minimum(fld.values), maximum(fld.values)] - [-0.4398236425818378, 0.022575693063719465]) < 1.0e-5
  # File = modeldata["postprocessing"]["exported"][1]["file"]
  # @async run(`"paraview.exe" $File`)
  try rm(modeldata["postprocessing"]["exported"][1]["file"]); catch end

  AE = AbaqusExporter("Cookstress_algo_stress");
  HEADING(AE, "Cook trapezoidal panel, plane stress");
  COMMENT(AE, "Converged free mid-edge displacement = 23.97");
  PART(AE, "part1");
  END_PART(AE);
  ASSEMBLY(AE, "ASSEM1");
  INSTANCE(AE, "INSTNC1", "PART1");
  NODE(AE, fens.xyz);
  COMMENT(AE, "We are assuming three node triangles in plane-stress");
  COMMENT(AE, "CPS3 are pretty poor-accuracy elements, but here we don't care about it.");
@test  nodesperelem(modeldata["regions"][1]["femm"].integdomain.fes) == 3
  ELEMENT(AE, "CPS3", "AllElements", connasarray(modeldata["regions"][1]["femm"].integdomain.fes))
  NSET_NSET(AE, "clamped", modeldata["essential_bcs"][1]["node_list"])
  ORIENTATION(AE, "GlobalOrientation", vec([1. 0 0]), vec([0 1. 0]));
  SOLID_SECTION(AE, "elasticity", "GlobalOrientation", "AllElements", thickness);
  END_INSTANCE(AE);
  END_ASSEMBLY(AE);
  MATERIAL(AE, "elasticity")
  ELASTIC(AE, E, nu)
  STEP_PERTURBATION_STATIC(AE)
  BOUNDARY(AE, "ASSEM1.INSTNC1.clamped", 1)
  BOUNDARY(AE, "ASSEM1.INSTNC1.clamped", 2)
  bfes = modeldata["traction_bcs"][1]["femm"].integdomain.fes
  COMMENT(AE, "Concentrated loads: we are assuming that the elements on the boundary");
  COMMENT(AE, "have two nodes each and also that they are the same length.");
  COMMENT(AE, "Then the concentrated loads below will be correctly lumped.");
  nl = connectednodes(bfes)
  F = zeros(count(modeldata["fens"]))
  for ix = 1:count(bfes)
    for jx = 1:2
      F[bfes.conn[ix][jx]] += 1.0/n/2/thickness
    end
  end
  for ixxxx = 1:length(F)
    if F[ixxxx] != 0.0
      CLOAD(AE, "ASSEM1.INSTNC1.$(ixxxx)", 2, F[ixxxx])
    end
  end
  END_STEP(AE)
  close(AE)

  nlines = 0
  open(AE.filename) do f
    s = readlines(f)
    nlines = length(s)
  end
  @test nlines == 2886
  rm(AE.filename)
end
end
using .mmmCookmstressisommm
mmmCookmstressisommm.test()

module mmLE10expimpmm
using FinEtools
using FinEtools.MeshExportModule
using FinEtools.MeshImportModule
using Test
import LinearAlgebra: norm, cholesky, cross
function test()


    # Thick elliptical plate with an elliptical hole is clamped on its exterior
    # boundary and is loaded with transverse  pressure.
    # This is a NAFEMS Benchmark, Test No. LE10.
    # The plate is discretized with solid elements.
    # Because of the symmetries of the geometry and load, only quarter of the plate is modeled.
    # The $\sigma_y=\sigma_2$ at the point $P$ is to be determined. Since the
    # target point is on the boundary of the domain it will not be an
    # integration node as we use Gauss quadrature. The reference value is -5.38 MPa.

    # println("LE10NAFEMS: Transverse deflection of elliptical plate with elliptical hole."        )
    t0 = time()

    E = 210e3*phun("MEGA*PA");# 210e3 MPa
    nu = 0.3;
    qmagn = 1.0*phun("MEGA*PA");# transverse pressure
    sigma_yP = -5.38*phun("MEGA*PA");# tensile stress at [2.0, 0.0] meters
    Ae =3.25*phun("m"); # Major radius of the exterior ellipse
    Be =2.75*phun("m"); # Minor radius of the exterior ellipse
    Ai =2.0*phun("m"); # Major radius of the interior ellipse
    Bi =1.0*phun("m"); # Minor radius of the interior ellipse
    Thickness = 0.6*phun("m")
    tolerance = Thickness/1000.; # Geometrical tolerance

    INP_file = """
    *HEADING
    NAFEMS TEST NLE10, COARSE MESH, C3D10M ELEMENTS
    *NODE
           1,     2.83277,       1.348
           2,     2.48116,     1.04837
           3,     2.12955,    0.748738
           4,     3.14146,    0.704779
           5,     2.33893,    0.399071
           6,     2.68977,    0.374369
           7,        3.25,          0.
           8,      2.8335,          0.
           9,       2.417,          0.
          10,     2.83277,       1.348,        0.15
          11,     2.48116,     1.04837,        0.15
          12,     2.12955,    0.748738,        0.15
          13,     2.33841,    0.400381,        0.15
          14,     3.14128,     0.70533,        0.15
          15,        3.25,          0.,        0.15
          16,      2.8335,          0.,        0.15
          17,       2.417,          0.,        0.15
          18,     2.83277,       1.348,         0.3
          19,     2.48116,     1.04837,         0.3
          20,     2.12955,    0.748738,         0.3
          21,     2.62488,       0.674,         0.3
          22,     2.33893,    0.399071,         0.3
          23,     3.14146,    0.704779,         0.3
          24,        3.25,          0.,         0.3
          25,      2.8335,          0.,         0.3
          26,       2.417,          0.,         0.3
          27,     2.83277,       1.348,        0.45
          28,     2.48116,     1.04837,        0.45
          29,     2.12955,    0.748738,        0.45
          30,     2.33841,    0.400381,        0.45
          31,     3.14128,     0.70533,        0.45
          32,        3.25,          0.,        0.45
          33,      2.8335,          0.,        0.45
          34,       2.417,          0.,        0.45
          35,     2.83277,       1.348,         0.6
          36,     2.48116,     1.04837,         0.6
          37,     2.12955,    0.748738,         0.6
          38,     3.14146,    0.704779,         0.6
          39,     2.33893,    0.399071,         0.6
          40,     2.68977,    0.374369,         0.6
          41,        3.25,          0.,         0.6
          42,      2.8335,          0.,         0.6
          43,       2.417,          0.,         0.6
          45,     1.95628,    0.600869
          46,     1.78302,       0.453
          47,     2.06477,    0.374369
          48,     1.93715,    0.248725
          51,      2.2085,          0.
          52,          2.,          0.
          54,     1.95628,    0.600869,        0.15
          55,     1.78302,       0.453,        0.15
          56,     1.93661,    0.249767,        0.15
          59,      2.2085,          0.,        0.15
          60,          2.,          0.,        0.15
          62,     1.95628,    0.600869,         0.3
          63,     1.78302,       0.453,         0.3
          64,     1.93715,    0.248725,         0.3
          66,     2.10001,      0.2265,         0.3
          68,      2.2085,          0.,         0.3
          69,          2.,          0.,         0.3
          71,     1.95628,    0.600869,        0.45
          72,     1.78302,       0.453,        0.45
          74,     1.93661,    0.249767,        0.45
          76,      2.2085,          0.,        0.45
          77,          2.,          0.,        0.45
          79,     1.95628,    0.600869,         0.6
          80,     1.78302,       0.453,         0.6
          81,     2.06477,    0.374369,         0.6
          83,     1.93715,    0.248725,         0.6
          85,      2.2085,          0.,         0.6
          86,          2.,          0.,         0.6
          87,       1.783,     2.29921
          88,     1.57618,     1.80182
          89,     1.36937,     1.30443
          90,     1.95627,     1.52397
          91,     2.36495,     1.88628
          92,     1.78146,     1.06985
          96,       1.783,     2.29921,        0.15
          97,     1.57618,     1.80182,        0.15
          98,     1.36937,     1.30443,        0.15
          99,     1.78071,     1.07038,        0.15
         100,     2.36449,     1.88669,        0.15
         104,       1.783,     2.29921,         0.3
         105,     1.57618,     1.80182,         0.3
         106,     1.36937,     1.30443,         0.3
         107,     2.36495,     1.88628,         0.3
         108,     1.78146,     1.06985,         0.3
         109,     2.10107,     1.32621,         0.3
         113,       1.783,     2.29921,        0.45
         114,     1.57618,     1.80182,        0.45
         115,     1.36937,     1.30443,        0.45
         116,     2.36449,     1.88669,        0.45
         117,     1.78071,     1.07038,        0.45
         121,       1.783,     2.29921,         0.6
         122,     1.57618,     1.80182,         0.6
         123,     1.36937,     1.30443,         0.6
         124,     1.95627,     1.52397,         0.6
         125,     2.36495,     1.88628,         0.6
         126,     1.78146,     1.06985,         0.6
         131,     1.26718,     1.05863
         132,       1.165,    0.812831
         134,     1.49321,    0.665266
         135,     1.64727,    0.780784
         140,     1.26718,     1.05863,        0.15
         141,       1.165,    0.812831,        0.15
         143,      1.4924,    0.665723,        0.15
         148,     1.26718,     1.05863,         0.3
         149,       1.165,    0.812831,         0.3
         150,     1.57619,    0.878714,         0.3
         152,     1.49321,    0.665266,         0.3
         157,     1.26718,     1.05863,        0.45
         158,       1.165,    0.812831,        0.45
         160,      1.4924,    0.665723,        0.45
         165,     1.26718,     1.05863,         0.6
         166,       1.165,    0.812831,         0.6
         168,     1.49321,    0.665266,         0.6
         169,     1.64727,    0.780784,         0.6
         173,          0.,       1.583
         174,          0.,      1.2915
         175,          0.,          1.
         176,      0.5825,     1.19792
         177,    0.699442,     1.51527
         178,    0.590531,    0.955415
         182,          0.,       1.583,        0.15
         183,          0.,      1.2915,        0.15
         184,          0.,          1.,        0.15
         185,    0.698861,     1.51538,        0.15
         186,    0.590076,    0.955486,        0.15
         190,          0.,       1.583,         0.3
         191,          0.,      1.2915,         0.3
         192,          0.,          1.,         0.3
         193,    0.699442,     1.51527,         0.3
         194,    0.590531,    0.955415,         0.3
         195,    0.684684,     1.15221,         0.3
         199,          0.,       1.583,        0.45
         200,          0.,      1.2915,        0.45
         201,          0.,          1.,        0.45
         202,    0.698861,     1.51538,        0.45
         203,    0.590076,    0.955486,        0.45
         207,          0.,       1.583,         0.6
         208,          0.,      1.2915,         0.6
         209,          0.,          1.,         0.6
         210,      0.5825,     1.19792,         0.6
         211,    0.699442,     1.51527,         0.6
         212,    0.590531,    0.955415,         0.6
         216,          0.,        2.75
         217,          0.,      2.1665
         219,    0.920251,     2.63745
         221,      0.8915,      1.9411
         225,          0.,        2.75,        0.15
         226,          0.,      2.1665,        0.15
         228,    0.919707,     2.63759,        0.15
         233,          0.,        2.75,         0.3
         234,          0.,      2.1665,         0.3
         236,    0.684684,     2.02721,         0.3
         237,    0.920251,     2.63745,         0.3
         242,          0.,        2.75,        0.45
         243,          0.,      2.1665,        0.45
         245,    0.919707,     2.63759,        0.45
         250,          0.,        2.75,         0.6
         251,          0.,      2.1665,         0.6
         253,    0.920251,     2.63745,         0.6
         255,      0.8915,      1.9411,         0.6
    **
    **
    *ELEMENT, TYPE=C3D10M, ELSET=EALL
           1,       1,       7,      18,       3,       4,      14,      10,
           2,       6,      11
           2,      24,      18,       7,      26,      23,      14,      15,
          25,      21,      16
           3,       9,       3,      26,       7,       5,      13,      17,
           8,       6,      16
           4,      20,      26,       3,      18,      22,      13,      12,
          19,      21,      11
           5,       3,       7,      18,      26,       6,      14,      11,
          13,      16,      21
           6,      24,      41,      18,      26,      32,      31,      23,
          25,      33,      21
           7,      35,      18,      41,      37,      27,      31,      38,
          36,      28,      40
           8,      20,      37,      26,      18,      29,      30,      22,
          19,      28,      21
           9,      43,      26,      37,      41,      34,      30,      39,
          42,      33,      40
          10,      18,      26,      41,      37,      21,      33,      31,
          28,      30,      40
          11,       9,      26,       3,      52,      17,      13,       5,
          51,      59,      47
          12,      20,       3,      26,      63,      12,      13,      22,
          62,      54,      66
          13,      46,      63,      52,       3,      55,      56,      48,
          45,      54,      47
          14,      69,      52,      63,      26,      60,      56,      64,
          68,      59,      66
          15,       3,      52,      26,      63,      47,      59,      13,
          54,      56,      66
          16,      20,      26,      37,      63,      22,      30,      29,
          62,      66,      71
          17,      43,      37,      26,      86,      39,      30,      34,
          85,      81,      76
          18,      69,      63,      86,      26,      64,      74,      77,
          68,      66,      76
          19,      80,      86,      63,      37,      83,      74,      72,
          79,      81,      71
          20,      63,      26,      37,      86,      66,      30,      71,
          74,      76,      81
          21,       1,      18,      87,       3,      10,     100,      91,
           2,      11,      90
          22,     104,      87,      18,     106,      96,     100,     107,
         105,      97,     109
          23,      89,     106,       3,      87,      98,      99,      92,
          88,      97,      90
          24,      20,       3,     106,      18,      12,      99,     108,
          19,      11,     109
          25,      87,       3,      18,     106,      90,      11,     100,
          97,      99,     109
          26,     104,      18,     121,     106,     107,     116,     113,
         105,     109,     114
          27,      35,     121,      18,      37,     125,     116,      27,
          36,     124,      28
          28,      20,     106,      37,      18,     108,     117,      29,
          19,     109,      28
          29,     123,      37,     106,     121,     126,     117,     115,
         122,     124,     114
          30,     106,      18,     121,      37,     109,     116,     114,
         117,      28,     124
          31,      89,       3,     106,     132,      92,      99,      98,
         131,     135,     140
          32,      20,     106,       3,      63,     108,      99,      12,
          62,     150,      54
          33,      46,     132,      63,       3,     134,     143,      55,
          45,     135,      54
          34,     149,      63,     132,     106,     152,     143,     141,
         148,     150,     140
          35,     132,       3,     106,      63,     135,      99,     140,
         143,      54,     150
          36,      20,      37,     106,      63,      29,     117,     108,
          62,      71,     150
          37,     123,     106,      37,     166,     115,     117,     126,
         165,     157,     169
          38,     149,     166,      63,     106,     158,     160,     152,
         148,     157,     150
          39,      80,      63,     166,      37,      72,     160,     168,
          79,      71,     169
          40,     106,      63,      37,     166,     150,      71,     117,
         157,     160,     169
          41,      89,     106,     173,     132,      98,     185,     177,
         131,     140,     176
          42,     190,     173,     106,     192,     182,     185,     193,
         191,     183,     195
          43,     175,     192,     132,     173,     184,     186,     178,
         174,     183,     176
          44,     149,     132,     192,     106,     141,     186,     194,
         148,     140,     195
          45,     173,     132,     106,     192,     176,     140,     185,
         183,     186,     195
          46,     190,     106,     207,     192,     193,     202,     199,
         191,     195,     200
          47,     123,     207,     106,     166,     211,     202,     115,
         165,     210,     157
          48,     149,     192,     166,     106,     194,     203,     158,
         148,     195,     157
          49,     209,     166,     192,     207,     212,     203,     201,
         208,     210,     200
          50,     192,     106,     207,     166,     195,     202,     200,
         203,     157,     210
          51,     216,      87,     233,     173,     219,     228,     225,
         217,     221,     226
          52,     104,     233,      87,     106,     237,     228,      96,
         105,     236,      97
          53,      89,     173,     106,      87,     177,     185,      98,
          88,     221,      97
          54,     190,     106,     173,     233,     193,     185,     182,
         234,     236,     226
          55,     173,      87,     233,     106,     221,     228,     226,
         185,      97,     236
          56,     104,     121,     233,     106,     113,     245,     237,
         105,     114,     236
          57,     250,     233,     121,     207,     242,     245,     253,
         251,     243,     255
          58,     190,     207,     106,     233,     199,     202,     193,
         234,     243,     236
          59,     123,     106,     207,     121,     115,     202,     211,
         122,     114,     255
          60,     233,     106,     121,     207,     236,     114,     245,
         243,     202,     255
    """
    write("NLE10.inp", INP_file)

    output = MeshImportModule.import_ABAQUS("NLE10.inp")
    fens, fes = output["fens"], output["fesets"][1]

    try rm("NLE10.inp"); catch end

    # Select the  boundary faces, on the boundary that is clamped,  and on the part
    # of the boundary that is loaded with the transverse pressure
    bdryfes = meshboundary(fes);
    exteriorbfl = selectelem(fens, bdryfes, facing=true,
        direction=[1.0, 1.0, 0.0], dotmin= 0.001);
    topbfl = selectelem(fens, bdryfes, box=[0.0, Inf, 0.0, Inf, Thickness, Thickness], inflate=tolerance);

    geom = NodalField(fens.xyz)
    u = NodalField(zeros(size(fens.xyz,1),3)) # displacement field

    L12 =connectednodes(subset(bdryfes, exteriorbfl)) # external boundary
    setebc!(u, L12, true, 1, 0.0)
    setebc!(u, L12, true, 2, 0.0)
    LL = selectnode(fens; box=[0.0, Inf, 0.0, Inf, Thickness/2.0, Thickness/2.0], inflate = tolerance)
    L3 = intersect(LL, connectednodes(subset(bdryfes, exteriorbfl)))
    setebc!(u, L3, true, 3, 0.0)
    L1 =selectnode(fens; box=[0.0, 0.0, 0.0, Inf, 0.0, Thickness], inflate = tolerance)
    setebc!(u,L1,true, 1, 0.0) # symmetry plane X = 0
    L2 =selectnode(fens; box=[0.0, Inf, 0.0, 0.0, 0.0, Thickness], inflate = tolerance)
    setebc!(u,L2,true, 2, 0.0) # symmetry plane Y = 0

    applyebc!(u)
    numberdofs!(u)

    eL1femm =  FEMMBase(IntegDomain(subset(bdryfes,topbfl), TriRule(3)))
    function pfun(forceout::FVec{T}, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt) where {T}
        forceout .=  [0.0, 0.0, -qmagn]
        return forceout
    end
    fi = ForceIntensity(FFlt, 3, pfun);
    F2 = distribloads(eL1femm, geom, u, fi, 2);

    # Note that the material object needs to be created with the proper
    # model-dimension reduction in mind.  In this case that is the fully three-dimensional solid.
    MR = DeforModelRed3D

    material = MatDeforElastIso(MR, E, nu)

    femm = FEMMDeforLinearMST10(MR, IntegDomain(fes, TetRule(4)), material)

    # The geometry field now needs to be associated with the FEMM
    femm = associategeometry!(femm, geom)

    K = stiffness(femm, geom, u)
    K = cholesky(K)
    U = K\(F2)
    scattersysvec!(u, U[:])

    nl = selectnode(fens, box=[Ai,Ai,0,0,Thickness,Thickness],inflate=tolerance);
    thecorneru = zeros(FFlt,1,3)
    gathervalues_asmat!(u, thecorneru, nl);
    thecorneru = thecorneru/phun("mm")
    # println("displacement =$(thecorneru) [MM] as compared to reference [-0.030939, 0, -0.10488] [MM]")
    @test norm(thecorneru - [-0.0268854 0.0 -0.0919955]) < 1.0e-5

    fld= fieldfromintegpoints(femm, geom, u, :Cauchy, 2; nodevalmethod = :invdistance, reportat = :meanonly)#
    # println("Sigma_y =$(fld.values[nl,1][1]/phun("MPa")) as compared to reference sigma_yP = $(sigma_yP/phun("MPa")) [MPa]")
    @test abs(fld.values[nl,1][1]/phun("MPa") - -2.2616980060965024) < 1.0e-3

    fld= fieldfromintegpoints(femm, geom, u, :Cauchy, 2; nodevalmethod = :averaging, reportat = :extrapmean)#
    # println("Sigma_y =$(fld.values[nl,1][1]/phun("MPa")) as compared to reference sigma_yP = $(sigma_yP/phun("MPa")) [MPa]")
    @test abs(fld.values[nl,1][1]/phun("MPa") - -2.382478776709117) < 1.0e-3

    fld= fieldfromintegpoints(femm, geom, u, :Cauchy, 2; nodevalmethod = :averaging, reportat = :extraptrend)#
    # println("Sigma_y =$(fld.values[nl,1][1]/phun("MPa")) as compared to reference sigma_yP = $(sigma_yP/phun("MPa")) [MPa]")
# println("$(fld.values[nl,1][1]/phun("MPa"))")
    @test abs(fld.values[nl,1][1]/phun("MPa") - -5.470291697493607) < 1.0e-3

    File =  "LE10NAFEMS_MST10_sigmay.vtk"
    vtkexportmesh(File, connasarray(fes), geom.values,
                   FinEtools.MeshExportModule.T10; vectors=[("u", u.values)],
                   scalars=[("sigmay", fld.values)])
    # @async run(`"paraview.exe" $File`)
    try rm(File) catch end

    AE = AbaqusExporter("LE10NAFEMS_MST10");
    HEADING(AE, "LE10NAFEMS: Transverse deflection of elliptical plate with elliptical hole.");
    PART(AE, "part1");
    END_PART(AE);
    ASSEMBLY(AE, "ASSEM1");
    INSTANCE(AE, "INSTNC1", "PART1");
    NODE(AE, fens.xyz);
    ELEMENT(AE, "c3d10", "AllElements", 1, connasarray(femm.integdomain.fes))
    ELEMENT(AE, "SFM3D6", "TractionElements",
    1+count(femm.integdomain.fes), connasarray(eL1femm.integdomain.fes))
    NSET_NSET(AE, "L1", L1)
    NSET_NSET(AE, "L2", L2)
    NSET_NSET(AE, "L3", L3)
    NSET_NSET(AE, "L12", L12)
    ORIENTATION(AE, "GlobalOrientation", vec([1. 0 0]), vec([0 1. 0]));
    SOLID_SECTION(AE, "elasticity", "GlobalOrientation", "AllElements", "Hourglassctl");
    SURFACE_SECTION(AE, "TractionElements")
    END_INSTANCE(AE);
    END_ASSEMBLY(AE);
    MATERIAL(AE, "elasticity")
    ELASTIC(AE, E, nu)
    SECTION_CONTROLS(AE, "Hourglassctl", "HOURGLASS=ENHANCED")
    STEP_PERTURBATION_STATIC(AE)
    BOUNDARY(AE, "ASSEM1.INSTNC1.L1", 1)
    BOUNDARY(AE, "ASSEM1.INSTNC1.L2", 2)
    BOUNDARY(AE, "ASSEM1.INSTNC1.L3", 3)
    BOUNDARY(AE, "ASSEM1.INSTNC1.L12", 1)
    BOUNDARY(AE, "ASSEM1.INSTNC1.L12", 2)
    DLOAD(AE, "ASSEM1.INSTNC1.TractionElements", vec([0.0, 0.0, -qmagn]))
    END_STEP(AE)
    close(AE)

    output = MeshImportModule.import_ABAQUS(AE.filename)
    fens, fes = output["fens"], output["fesets"][1]
try rm(AE.filename) catch end

    # Select the  boundary faces, on the boundary that is clamped,  and on the part
    # of the boundary that is loaded with the transverse pressure
    bdryfes = meshboundary(fes);
    exteriorbfl = selectelem(fens, bdryfes, facing=true, direction=[1.0, 1.0, 0.0], dotmin= 0.001);
    topbfl = selectelem(fens, bdryfes, box=[0.0, Inf, 0.0, Inf, Thickness, Thickness], inflate=tolerance);

    geom = NodalField(fens.xyz)
    u = NodalField(zeros(size(fens.xyz,1),3)) # displacement field

    L12 =connectednodes(subset(bdryfes, exteriorbfl)) # external boundary
    setebc!(u, L12, true, 1, 0.0)
    setebc!(u, L12, true, 2, 0.0)
    LL = selectnode(fens; box=[0.0, Inf, 0.0, Inf, Thickness/2.0, Thickness/2.0], inflate = tolerance)
    L3 = intersect(LL, connectednodes(subset(bdryfes, exteriorbfl)))
    setebc!(u, L3, true, 3, 0.0)
    L1 =selectnode(fens; box=[0.0, 0.0, 0.0, Inf, 0.0, Thickness], inflate = tolerance)
    setebc!(u,L1,true, 1, 0.0) # symmetry plane X = 0
    L2 =selectnode(fens; box=[0.0, Inf, 0.0, 0.0, 0.0, Thickness], inflate = tolerance)
    setebc!(u,L2,true, 2, 0.0) # symmetry plane Y = 0

    applyebc!(u)
    numberdofs!(u)

    eL1femm =  FEMMBase(IntegDomain(subset(bdryfes,topbfl), TriRule(3)))
    # function pfun(forceout::FVec{T}, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt) where {T}
    #     forceout .=  [0.0, 0.0, -qmagn]
    #     return forceout
    # end
    fi = ForceIntensity(FFlt, 3, pfun);
    F2 = distribloads(eL1femm, geom, u, fi, 2);

    # Note that the material object needs to be created with the proper
    # model-dimension reduction in mind.  In this case that is the fully three-dimensional solid.
    MR = DeforModelRed3D

    material = MatDeforElastIso(MR, E, nu)

    femm = FEMMDeforLinearMST10(MR, IntegDomain(fes, TetRule(4)), material)

    # The geometry field now needs to be associated with the FEMM
    femm = associategeometry!(femm, geom)

    K = stiffness(femm, geom, u)
    K = cholesky(K)
    U = K\(F2)
    scattersysvec!(u, U[:])

    nl = selectnode(fens, box=[Ai,Ai,0,0,Thickness,Thickness],inflate=tolerance);
    thecorneru = zeros(FFlt,1,3)
    gathervalues_asmat!(u, thecorneru, nl);
    thecorneru = thecorneru/phun("mm")
    # println("displacement =$(thecorneru) [MM] as compared to reference [-0.030939, 0, -0.10488] [MM]")
    @test norm(thecorneru - [-0.0268854 0.0 -0.0919955]) < 1.0e-5

    fld= fieldfromintegpoints(femm, geom, u, :Cauchy, 2; reportat = :meanonly)#
    # println("Sigma_y =$(fld.values[nl,1][1]/phun("MPa")) as compared to reference sigma_yP = $(sigma_yP/phun("MPa")) [MPa]")
    @test abs(fld.values[nl,1][1]/phun("MPa") - -2.2616980060965024) < 1.0e-3

    fld= fieldfromintegpoints(femm, geom, u, :Cauchy, 2; nodevalmethod = :averaging, reportat = :extrapmean)#
    # println("Sigma_y =$(fld.values[nl,1][1]/phun("MPa")) as compared to reference sigma_yP = $(sigma_yP/phun("MPa")) [MPa]")
    @test abs(fld.values[nl,1][1]/phun("MPa") - -2.382478776709117) < 1.0e-3

    fld= fieldfromintegpoints(femm, geom, u, :Cauchy, 2; nodevalmethod = :averaging, reportat = :extraptrend)#
    # println("Sigma_y =$(fld.values[nl,1][1]/phun("MPa")) as compared to reference sigma_yP = $(sigma_yP/phun("MPa")) [MPa]")
# println("$(fld.values[nl,1][1]/phun("MPa"))")
    @test abs(fld.values[nl,1][1]/phun("MPa") - -5.470291697493607) < 1.0e-3

    File =  "LE10NAFEMS_MST10_sigmay.vtk"
    vtkexportmesh(File, connasarray(fes), geom.values,
                   FinEtools.MeshExportModule.T10; vectors=[("u", u.values)],
                   scalars=[("sigmay", fld.values)])
    # @async run(`"paraview.exe" $File`)
    try rm(File) catch end
end
end
using .mmLE10expimpmm
mmLE10expimpmm.test()


module mmtruncatedmfreem1
using FinEtools
using FinEtools.AlgoDeforLinearModule
using Test
import Arpack: eigs
import LinearAlgebra: norm, cholesky, cross
function test()
    # println("""
    # Vibration modes of truncated cylindrical shell.
    # """)

    # t0 = time()

    E = 205000*phun("MPa");# Young's modulus
    nu = 0.3;# Poisson ratio
    rho = 7850*phun("KG*M^-3");# mass density
    OmegaShift = (2*pi*100) ^ 2; # to resolve rigid body modes
    h = 0.05*phun("M");
    l = 10*h;
    Rmed = h/0.2;
    psi   = 0;    # Cylinder
    nh = 5; nl  = 12; nc = 40;
    tolerance = h/nh/100;
    neigvs = 20;

    MR = DeforModelRed3D
    fens,fes  = H8block(h,l,2*pi,nh,nl,nc)
    # Shape into a cylinder
    R = zeros(3, 3)
    for i = 1:count(fens)
        x, y, z = fens.xyz[i,:];
        rotmat3!(R, [0, z, 0])
        Q = [cos(psi*pi/180) sin(psi*pi/180) 0;
            -sin(psi*pi/180) cos(psi*pi/180) 0;
            0 0 1]
        fens.xyz[i,:] = reshape([x+Rmed-h/2, y-l/2, 0], 1, 3)*Q*R;
    end
    candidates = selectnode(fens, plane = [0.0 0.0 1.0 0.0], thickness = h/1000)
    fens,fes = mergenodes(fens, fes,  tolerance, candidates);

    geom = NodalField(fens.xyz)
    u = NodalField(zeros(size(fens.xyz,1),3)) # displacement field

    numberdofs!(u)

    material=MatDeforElastIso(MR, rho, E, nu, 0.0)

    femm = FEMMDeforLinearMSH8(MR, IntegDomain(fes, GaussRule(3,2)), material)
    femm = associategeometry!(femm, geom)
    K =stiffness(femm, geom, u)
    femm = FEMMDeforLinear(MR, IntegDomain(fes, GaussRule(3,3)), material)
    M =mass(femm, geom, u)


    # eigs returns the nev requested eigenvalues in d, the corresponding Ritz vectors
    # v (only if ritzvec=true), the number of converged eigenvalues nconv, the number
    # of iterations niter and the number of matrix vector multiplications nmult, as
    # well as the final residual vector resid.

    if true
        d,v,nev,nconv = eigs(K+OmegaShift*M, M; nev=neigvs, which=:SM)
        d[:] = d .- OmegaShift;
        fs = real(sqrt.(complex(d)))/(2*pi)
        # println("Eigenvalues: $fs [Hz]")
        @test norm(sort(fs)[1:8] - [0.0, 0.0, 0.0, 0.0, 0.000166835, 0.000182134, 517.147, 517.147]) < 2.0e-2
        # mode = 7
        # scattersysvec!(u, v[:,mode])
        # File =  "unit_cube_modes.vtk"
        # vtkexportmesh(File, fens, fes; vectors=[("mode$mode", u.values)])
        # @async run(`"paraview.exe" $File`)
    end

    if true
        solver = AlgoDeforLinearModule.ssit
        v0 = [i==j ? one(FFlt) : zero(FFlt) for i=1:size(K,1), j=1:2*neigvs]
        tol = 1.0e-2
        maxiter = 20
        lamb, v, nconv, niter, nmult, lamberr =
            solver(K+OmegaShift*M, M; nev=neigvs, v0=v0, tol=tol, maxiter=maxiter)
        @test nconv != neigvs
        # if nconv < neigvs
            # println("NOT converged")
        # end
        broadcast!(+, lamb, lamb, - OmegaShift);
        fs = real(sqrt.(complex(lamb)))/(2*pi)
        # println("Eigenvalues: $fs [Hz]")
        # println("$(sort(fs))")
        @test norm(sort(fs)[1:8] - [0.0, 0.0, 0.0, 0.0, 7.9048e-5, 0.0, 517.147, 517.147]) < 2.0e-2
        # println("Eigenvalue errors: $lamberr [ND]")
        # mode = 7
        # scattersysvec!(u, v[:,mode])
        # File =  "unit_cube_modes.vtk"
        # vtkexportmesh(File, fens, fes; vectors=[("mode$mode", u.values)])
        # @async run(`"paraview.exe" $File`)
    end

end
end
using .mmtruncatedmfreem1
mmtruncatedmfreem1.test()

module mmFV32mm1
using FinEtools
using Test
import LinearAlgebra: norm, cholesky, cross
function test()
    # println("""
    # FV32: Cantilevered tapered membrane
    # This is a test recommended by the National Agency for Finite Element Methods and
    # Standards (U.K.): Test FV32 from NAFEMS publication TNSB, Rev. 3, “The Standard
    # NAFEMS Benchmarks,” October 1990.
    #
    # Reference solution: 44.623	130.03	162.70	246.05	379.90	391.44 for the first
    # six modes.
    # """)

    t0 = time()


    E = 200*phun("GPA");
    nu = 0.3;
    rho= 8000*phun("KG/M^3");
    L = 10*phun("M");
    W0 = 5*phun("M");
    WL = 1*phun("M");
    H = 0.05*phun("M");
    nL, nW, nH = 8, 4, 1;# How many element edges per side?
    neigvs=20                   # how many eigenvalues
    Reffs = [44.623	130.03	162.70	246.05	379.90	391.44]

    fens,fes =H20block(1.0, 2.0, 1.0, nL, nW, nH)
    for i = 1:count(fens)
        xi, eta, theta = fens.xyz[i,:];
        eta = eta - 1.0
        fens.xyz[i,:] = [xi*L eta*(1.0 - 0.8*xi)*W0/2 theta*H/2];
    end
    # File =  "mesh.vtk"
    # vtkexportmesh(File, fens, fes)
    # @async run(`"paraview.exe" $File`)

    # Make the region
    MR = DeforModelRed3D
    material = MatDeforElastIso(MR, rho, E, nu, 0.0)
    region1 = FDataDict("femm"=>FEMMDeforLinear(MR, IntegDomain(fes, GaussRule(3,2)),
      material), "femm_mass"=>FEMMDeforLinear(MR, IntegDomain(fes, GaussRule(3,3)),
      material))

    nl1 = selectnode(fens; plane=[1.0 0.0 0.0 0.0], thickness=H/1.0e4)
    ebc1 = FDataDict("node_list"=>nl1, "component"=>1, "displacement"=>0.0)
    ebc2 = FDataDict("node_list"=>nl1, "component"=>2, "displacement"=>0.0)
    ebc3 = FDataDict("node_list"=>nl1, "component"=>3, "displacement"=>0.0)

    nl4 = selectnode(fens; plane=[0.0 0.0 1.0 0.0], thickness=H/1.0e4)
    ebc4 = FDataDict("node_list"=>nl4, "component"=>3, "displacement"=>0.0)

    # Make model data
    modeldata =  FDataDict(
      "fens"=> fens, "regions"=>  [region1], "essential_bcs"=>[ebc1 ebc2 ebc3 ebc4],
      "neigvs"=>neigvs)

    # Solve
    modeldata = FinEtools.AlgoDeforLinearModule.modal(modeldata)

    fs = modeldata["omega"]/(2*pi)
    # println("Eigenvalues: $fs [Hz]")
    # println("Percentage frequency errors: $((vec(fs[1:6]) - vec(Reffs))./vec(Reffs)*100)")
    @test norm((vec(fs[1:6]) - vec(Reffs))./vec(Reffs)*100 - [0.0162775, 0.0623384, 0.00799148, 0.151669, 0.376663, 0.0191388]) < 1.0e-6

    modeldata["postprocessing"] = FDataDict("file"=>"FV32-modes", "mode"=>1:10)
    modeldata=FinEtools.AlgoDeforLinearModule.exportmode(modeldata)
    # @async run(`"paraview.exe" $(modeldata["postprocessing"]["file"]*"1.vtk")`)
    try rm(modeldata["postprocessing"]["file"]*"1.vtk") catch end
end
end
using .mmFV32mm1
mmFV32mm1.test()

module mmtruncatedmfreem2
using FinEtools
using FinEtools.AlgoDeforLinearModule
using Test
import Arpack: eigs
import LinearAlgebra: norm, cholesky, cross
function test()
    # println("""
    # Vibration modes of truncated cylindrical shell.
    # """)

    # t0 = time()

    E = 205000*phun("MPa");# Young's modulus
    nu = 0.3;# Poisson ratio
    rho = 7850*phun("KG*M^-3");# mass density
    OmegaShift = (2*pi*100) ^ 2; # to resolve rigid body modes
    h = 0.05*phun("M");
    l = 10*h;
    Rmed = h/0.2;
    psi   = 0;    # Cylinder
    nh = 5; nl  = 12; nc = 40;
    tolerance = h/nh/100;
    neigvs = 20;

    MR = DeforModelRed3D
    fens,fes  = H8block(h,l,2*pi,nh,nl,nc)
    # Shape into a cylinder
    R = zeros(3, 3)
    for i = 1:count(fens)
        x, y, z = fens.xyz[i,:];
        rotmat3!(R, [0, z, 0])
        Q = [cos(psi*pi/180) sin(psi*pi/180) 0;
            -sin(psi*pi/180) cos(psi*pi/180) 0;
            0 0 1]
        fens.xyz[i,:] = reshape([x+Rmed-h/2, y-l/2, 0], 1, 3)*Q*R;
    end
    candidates = selectnode(fens, plane = [0.0 0.0 1.0 0.0], thickness = h/1000)
    fens,fes = mergenodes(fens, fes,  tolerance, candidates);

    geom = NodalField(fens.xyz)
    u = NodalField(zeros(size(fens.xyz,1),3)) # displacement field

    numberdofs!(u)

    material=MatDeforElastIso(MR, rho, E, nu, 0.0)

    femm = FEMMDeforLinearMSH8(MR, IntegDomain(fes, GaussRule(3,2)), material)
    femm = associategeometry!(femm, geom)
    K =stiffness(femm, geom, u)
    femm = FEMMDeforLinear(MR, IntegDomain(fes, GaussRule(3,3)), material)
    M =mass(femm, geom, u)


    # eigs returns the nev requested eigenvalues in d, the corresponding Ritz vectors
    # v (only if ritzvec=true), the number of converged eigenvalues nconv, the number
    # of iterations niter and the number of matrix vector multiplications nmult, as
    # well as the final residual vector resid.

    if true
        d,v,nev,nconv = eigs(K+OmegaShift*M, M; nev=neigvs, which=:SM)
        broadcast!(+, d, d, - OmegaShift);
        fs = real(sqrt.(complex(d)))/(2*pi)
        # println("Eigenvalues: $fs [Hz]")
        @test norm(sort(fs)[1:8] - [0.0, 0.0, 0.0, 0.0, 0.000166835, 0.000182134, 517.147, 517.147]) < 2.0e-2
        # mode = 7
        # scattersysvec!(u, v[:,mode])
        # File =  "unit_cube_modes.vtk"
        # vtkexportmesh(File, fens, fes; vectors=[("mode$mode", u.values)])
        # @async run(`"paraview.exe" $File`)
    end

    if true
        solver = AlgoDeforLinearModule.ssit
        v0 = [i==j ? one(FFlt) : zero(FFlt) for i=1:size(K,1), j=1:2*neigvs]
        tol = 1.0e-2
        maxiter = 20
        lamb, v, nconv, niter, nmult, lamberr =
            solver(K+OmegaShift*M, M; nev=neigvs, v0=v0, tol=tol, maxiter=maxiter, withrr=true)
        @test nconv == neigvs
        # if nconv < neigvs
            # println("NOT converged")
        # end
        lamb = lamb .- OmegaShift;
        fs = real(sqrt.(complex(lamb)))/(2*pi)
        # println("Eigenvalues: $fs [Hz]")
        # println("$(sort(fs))")
        @test norm(sort(fs)[1:8] - [0.0, 0.0, 0.0, 0.0, 7.9048e-5, 0.0, 517.147, 517.147]) < 2.0e-2
        # println("Eigenvalue errors: $lamberr [ND]")
        # mode = 7
        # scattersysvec!(u, v[:,mode])
        # File =  "unit_cube_modes.vtk"
        # vtkexportmesh(File, fens, fes; vectors=[("mode$mode", u.values)])
        # @async run(`"paraview.exe" $File`)
    end

end
end
using .mmtruncatedmfreem2
mmtruncatedmfreem2.test()

module mfiber_reinf_cant_yn_strong_Abaqus
using FinEtools
using FinEtools.AlgoDeforLinearModule
using Test
import LinearAlgebra: Symmetric, cholesky
import Statistics: mean
function test()


# println("""
# Cantilever example.  Strongly orthotropic material. Orientation "y".
# @article{
# author = {Krysl, P.},
# title = {Mean-strain 8-node hexahedron with optimized energy-sampling stabilization},
# journal = {Finite Elements in Analysis and Design},
# volume = {108}, pages = {41-53}, DOI = {10.1016/j.finel.2015.09.008}, year = {2016}
# }
# """)

t0 = time()
pu = ustring -> phun(ustring; system_of_units = :SIMM)

# # Orthotropic material
E1s = 100000.0*pu("GPa")
E2s = 1.0*pu("GPa")
E3s = E2s
nu23s = nu12s = nu13s = 0.25
G12s = 0.2*pu("GPa")
G23s = G13s = G12s
CTE1 = 0.0
CTE2 = 0.0
CTE3 = 0.0
# # Isotropic material
# E = 1.0e9*pu("Pa")
# nu = 0.25
# CTE = 0.0

# Reference value for  the vertical deflection of the tip
uz_ref = -1.027498445054843e-05*pu("m");

a = 90.0*pu("mm") # length of the cantilever
b = 10.0*pu("mm") # width of the cross-section
t = 20.0*pu("mm") # height of the cross-section
q0 = -1000.0*pu("Pa") # shear traction
dT = 0*pu("K") # temperature rise

tolerance = 0.00001*t

# Generate mesh
n = 2
na = 4*n # number of elements lengthwise
nb = 2*n # number of elements through the depth
nt = n # number of elements through the thickness
xs = collect(linearspace(0.0, a, na+1))
ys = collect(linearspace(0.0, b, nb+1))
ts = collect(linearspace(0.0, t, nt+1))
# println("Mesh generation")
fens,fes = H8blockx(xs, ys, ts)
fens,fes = H8toH20(fens,fes)

bfes = meshboundary(fes)
# end cross-section surface  for the shear loading
sshearl = selectelem(fens, bfes; facing=true, direction = [+1.0 0.0 0.0])

MR = DeforModelRed3D
material = MatDeforElastOrtho(MR,
  0.0, E1s, E2s, E3s,
  nu12s, nu13s, nu23s,
  G12s, G13s, G23s,
  CTE1, CTE2, CTE3)
# material = MatDeforElastIso(MR,
#   0.0, E, nu, CTE)

# Material orientation matrix
csmat = zeros(3, 3)
rotmat3!(csmat, -45.0/180.0*pi*[0,1,0])

function updatecs!(csmatout::FFltMat, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)
  copyto!(csmatout, csmat)
end

femm = FEMMDeforLinear(MR, IntegDomain(fes, GaussRule(3, 2)), CSys(3, 3, updatecs!), material)

lx0 = selectnode(fens, box=[0.0 0.0 -Inf Inf -Inf Inf], inflate=tolerance)

geom = NodalField(fens.xyz)
u = NodalField(zeros(size(fens.xyz,1),3)) # displacement field
nnodes(geom)

setebc!(u, lx0, true, 1, zeros(size(lx0)))
setebc!(u, lx0, true, 2, zeros(size(lx0)))
setebc!(u, lx0, true, 3, zeros(size(lx0)))
applyebc!(u)

# S = connectionmatrix(femm.femmbase, nnodes(geom))
numberdofs!(u)

function getshr!(forceout::FFltVec, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)
  copyto!(forceout, q0*[0.0; 0.0; 1.0])
end

Tracfemm = FEMMBase(IntegDomain(subset(bfes, sshearl), GaussRule(2, 3)))

# println("Stiffness")
K = stiffness(femm, geom, u)
fi = ForceIntensity(FFlt, 3, getshr!);
# println("Traction loads")
F =  distribloads(Tracfemm, geom, u, fi, 2);

# println("Factorization")
K = (K + K')/2;
K = cholesky(Symmetric(K))
# println("U = K\\F")
U = K\F
# # println("U = cg(K, F; tol=1e-3, maxiter=2000)")
# U = cg(K, F; tol=1e-3, maxiter=2000)
scattersysvec!(u, U[:])

Tipl = selectnode(fens, box=[a a b b 0. 0.], inflate=tolerance)
utip = mean(u.values[Tipl, 3])
# println("Deflection $utip, normalized: $(utip/uz_ref)")
@test abs(utip - -0.00653888266072445)/abs(utip) < 1.0e-6
# println("Solution: $(  time()-t0 )")

AE = AbaqusExporter("fiber_reinf_cant_yn_strong_Abaqus");
HEADING(AE, "fiber_reinf_cant_yn_strong_Abaqus.");
PART(AE, "part1");
END_PART(AE);
ASSEMBLY(AE, "ASSEM1");
INSTANCE(AE, "INSTNC1", "PART1");
NODE(AE, fens.xyz);
ELEMENT(AE, "c3d20r", "AllElements", connasarray(femm.integdomain.fes))
ELEMENT(AE, "SFM3D8", "TractionElements", connasarray(Tracfemm.integdomain.fes))
NSET_NSET(AE, "L1", lx0)
NSET_NSET(AE, "L2", lx0)
NSET_NSET(AE, "L3", lx0)
NSET_NSET(AE, "tip", Tipl)
ORIENTATION(AE, "GlobalOrientation", vec([1. 0 0]), vec([0 1. 0]));
ORIENTATION(AE, "MaterialOrientation", vec(csmat[:,1]), vec(csmat[:,2]));
SOLID_SECTION(AE, "elasticity", "MaterialOrientation", "AllElements", "Hourglassctl");
SURFACE_SECTION(AE, "TractionElements")
END_INSTANCE(AE);
END_ASSEMBLY(AE);
MATERIAL(AE, "elasticity")
ELASTIC(AE, E1s, E2s, E3s, nu12s, nu13s, nu23s, G12s, G13s, G23s)
SECTION_CONTROLS(AE, "Hourglassctl", "HOURGLASS=ENHANCED")
STEP_PERTURBATION_STATIC(AE)
BOUNDARY(AE, "ASSEM1.INSTNC1.L1", 1)
BOUNDARY(AE, "ASSEM1.INSTNC1.L2", 2)
BOUNDARY(AE, "ASSEM1.INSTNC1.L3", 3)
DLOAD(AE, "ASSEM1.INSTNC1.TractionElements", vec([0.0, 0.0, q0]))
NODE_PRINT(AE, "ASSEM1.INSTNC1.tip")
ENERGY_PRINT(AE)
END_STEP(AE)
close(AE)
try rm(AE.filename); catch end

# println("Done: $(  time()-t0 )")
true
end
end
using .mfiber_reinf_cant_yn_strong_Abaqus
mfiber_reinf_cant_yn_strong_Abaqus.test()

module mmorthoballoonpenaltymm
using FinEtools
using Test
import LinearAlgebra: norm, cholesky, cross
function test()

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

    icl = selectelem(fens, bdryfes, box=[0.,0.,0.,pi/2],inflate=tolerance);
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
    # println("Number of degrees of freedom = $(u.nfreedofs)")

    # The traction boundary condition is applied in the radial
    # direction.

    el1femm =  FEMMBase(IntegDomain(subset(bdryfes,icl), GaussRule(1, 3), true))
    function pressureloading!(forceout::FFltVec, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)
      copyto!(forceout, XYZ/norm(XYZ)*p)
      return forceout
    end
    fi = ForceIntensity(FFlt, 2, pressureloading!); # pressure normal to the internal cylindrical surface
    F2= distribloads(el1femm, geom, u, fi, 2);

    # Property and material
    material=MatDeforElastOrtho(MR, E1,E2,E3,nu12,nu13,nu23,G12,G13,G23)

    femm = FEMMDeforLinear(MR, IntegDomain(fes, GaussRule(2, 2), true), material)

    ##
    # The restraints of the nodes on the bounding cross-sections in the direction
    # of the normal to the plane of the cross-section  in the
    # circumferential direction are introduced using a penalty formulation.
    # For that purpose we introduce  a finite element model machine for the
    # surface  finite elements on the cross-sections.
    springcoefficient =1.0e9/ (abs(p)/E1)
    xsfemm = FEMMDeforWinkler(IntegDomain(subset(bdryfes, lx), GaussRule(1, 3), true))
    ysfemm = FEMMDeforWinkler(IntegDomain(subset(bdryfes, ly), GaussRule(1, 3), true))
    H = surfacenormalspringstiffness(xsfemm,  geom, u, springcoefficient, SurfaceNormal(3)) +
        surfacenormalspringstiffness(ysfemm,  geom, u, springcoefficient, SurfaceNormal(3))
    K =stiffness(femm, geom, u)
    U=  (K + H)\(F2)
    scattersysvec!(u,U[:])

    # Produce a plot of the radial stress component in the cylindrical
    # coordinate system. Note that this is the usual representation of
    # stress using nodal stress field.

    fld= fieldfromintegpoints(femm, geom, u, :Cauchy, 3)
    # println("Minimum/maximum = $(minimum(fld.values))/$(maximum(fld.values))")
@test abs(minimum(fld.values) - (-0.04635309320688638)) < 1.0e-5
@test abs(maximum(fld.values) - (0.5708149883384825)) < 1.0e-5
    # File =  "orthoballoon_penalty_sigmaz.vtk"
    # vtkexportmesh(File, fens, fes; scalars=[("sigmaz", fld.values)],
    #               vectors=[("u", u.values)])
    # @async run(`"paraview.exe" $File`)

end
end
using .mmorthoballoonpenaltymm
mmorthoballoonpenaltymm.test()

module mbar1
using FinEtools
using FinEtools.FENodeSetModule
using FinEtools.MeshExportModule
using Test
import LinearAlgebra: norm, cholesky, cross
function test()
    Area = 2.0*phun("in^2")
    E = 30e6*phun("psi") # Young's modulus
    nu = 0.0
    alpha = 5.5e-6*phun("in")/phun("in")/phun("F")
    fens = FENodeSetModule.FENodeSet([10. -5 20;
    30 25 -15]*phun("in") )
    fes = FESetL2(reshape([1,2],1,2))

    # function otherdimensionfu(loc::FFltMat,
    #   conn::CC, N::FFltMat)::FFlt where {CC<:AbstractArray{FInt}}
    #   return otherdimension::FFlt
    # end
    integdomain = IntegDomain(fes, GaussRule(1, 2), (loc, conn, N) -> Area, false)
    # display(IntegDomain)

    MR = DeforModelRed1D
    material = MatDeforElastIso(MR,  0.0, E, nu, alpha)
    # display(material )

    geom = NodalField(fens.xyz)
    u = NodalField(zeros(size(fens.xyz, 1), 3)) # displacement field
    numberdofs!(u)

    femm = FEMMDeforLinear(MR, integdomain, CSys(3, 1), material)
    K = stiffness(femm,  geom,  u)
    # println("K = $(K/(phun("lbf")/phun("in")))")
    ref_K=   1.0e+05 *[
    1.8916    2.8373   -3.3102   -1.8916   -2.8373    3.3102
    2.8373    4.2560   -4.9653   -2.8373   -4.2560    4.9653
   -3.3102   -4.9653    5.7929    3.3102    4.9653   -5.7929
   -1.8916   -2.8373    3.3102    1.8916    2.8373   -3.3102
   -2.8373   -4.2560    4.9653    2.8373    4.2560   -4.9653
    3.3102    4.9653   -5.7929   -3.3102   -4.9653    5.7929];
    @test norm(K/(phun("lbf")/phun("in")) - ref_K)/1.0e5 < 1.0e-3

    dT = NodalField(broadcast(+, zeros(size(fens.xyz, 1), 1),
        100*phun("F"))) # temperature field
    # display(dT)

    F2 = thermalstrainloads(femm, geom, u, dT)
    # println("F2 = $(F2/(phun("lbf")))")
    ref_F=  1.0e+04 *[
    -1.313449091077187
    -1.970173636615779
    2.298535909385076
    1.313449091077187
    1.970173636615779
    -2.298535909385076]
    @test norm(F2/(phun("lbf")) - ref_F) < 1.0e-2

   # K = cholesky(K)
   # U=  K\(F2)
   # scattersysvec!(u, U[:])

    # File = "playground.vtk"
    # MeshExportModule.vtkexportmesh(File, fens, fes)
    # @async run(`"paraview.exe" $File`)
    # try rm(File) catch end

end
end
using .mbar1
mbar1.test()

module mbar2
using FinEtools
using FinEtools.FENodeSetModule
using FinEtools.MeshExportModule
using Test
import LinearAlgebra: norm, cholesky, cross
function test()
    Area = 1.5
    E = 1.0e7 # Young's modulus
    nu = 0.0
    alpha = 0.0
    fens = FENodeSetModule.FENodeSet([0.0 0;
    0 40;
    40 0;
    40 40;
    80 0;
    80 40] )
    fes = FESetL2([1     3
     1     4
     2     4
     3     4
     3     5
     5     4
     6     4
     5     6])


    MR = DeforModelRed1D
    material = MatDeforElastIso(MR,  0.0, E, nu, alpha)
    # display(material )

    geom = NodalField(fens.xyz)
    u = NodalField(zeros(size(fens.xyz, 1), 2)) # displacement field
    setebc!(u, 1)
    setebc!(u, 2)
    applyebc!(u)
    numberdofs!(u)
    # display(u)

    femm = FEMMDeforLinear(MR, IntegDomain(fes, GaussRule(1, 1),  (loc, conn, N) -> Area, false), CSys(2, 1), material)
    K = stiffness(femm,  geom,  u)

    fi = ForceIntensity(vec([0 -2000.0]))
    lfemm = FEMMBase(IntegDomain(FESetP1(reshape([3], 1,1)), PointRule()))
    F = distribloads(lfemm,  geom,  u,  fi,  3);
    fi = ForceIntensity(vec([+2000.0 0]))
    lfemm = FEMMBase(IntegDomain(FESetP1(reshape([5], 1,1)), PointRule()))
    F = F + distribloads(lfemm,  geom,  u,  fi,  3);
    fi = ForceIntensity(vec([+4000.0 +6000.0]))
    lfemm = FEMMBase(IntegDomain(FESetP1(reshape([6], 1,1)), PointRule()))
    F = F + distribloads(lfemm,  geom,  u,  fi,  3);

    K = cholesky(K)
    U=  K\F
    scattersysvec!(u, U[:])
    @test norm(u.values  - [ 0.0         0.0
                              0.0         0.0
                              0.0213333   0.0408366
                             -0.016       0.0461699
                              0.0426667   0.150091
                             -0.00533333  0.166091 ]) < 1.0e-4

    sfld =  elemfieldfromintegpoints(femm, geom, u, :Cauchy, 1)
    # display(sfld)
    # println("Cauchy = $(sfld.values)")
    @test norm(sfld.values - [5333.33; 3771.24; -4000.0; 1333.33; 5333.33; -5656.85; 2666.67; 4000.0]) < 1.0e-2
    vfld =  elemfieldfromintegpoints(femm, geom, u, :vm, 1)
    # display(vfld)

    File = "Planar_truss.vtk"
    MeshExportModule.vtkexportmesh(File, fens, fes;
    scalars=[("sx", sfld.values), ("vm", vfld.values)])
    # @async run(`"paraview.exe" $File`)
    try rm(File) catch end

end
end
using .mbar2
mbar2.test()

module mmmLE10expiAbaqus2mmmm
using FinEtools
using FinEtools.MeshExportModule
using Test
import LinearAlgebra: norm, cholesky, cross
function test()


    # Thick elliptical plate with an elliptical hole is clamped on its exterior
    # boundary and is loaded with transverse  pressure.
    # This is a NAFEMS Benchmark, Test No. LE10.
    # The plate is discretized with solid elements.
    # Because of the symmetries of the geometry and load, only quarter of the plate is modeled.
    # The $\sigma_y=\sigma_2$ at the point $P$ is to be determined. Since the
    # target point is on the boundary of the domain it will not be an
    # integration node as we use Gauss quadrature. The reference value is -5.38 MPa.

    # println("LE10NAFEMS: Transverse deflection of elliptical plate with elliptical hole."        )
    t0 = time()

    E = 210e3*phun("MEGA*PA");# 210e3 MPa
    nu = 0.3;
    qmagn = 1.0*phun("MEGA*PA");# transverse pressure
    sigma_yP = -5.38*phun("MEGA*PA");# tensile stress at [2.0, 0.0] meters
    Ae =3.25*phun("m"); # Major radius of the exterior ellipse
    Be =2.75*phun("m"); # Minor radius of the exterior ellipse
    Ai =2.0*phun("m"); # Major radius of the interior ellipse
    Bi =1.0*phun("m"); # Minor radius of the interior ellipse
    Thickness = 0.6*phun("m")
    nc = 6; # number of elements per side
    nr = 5; # number of elements per side
    nt = 2; # number of elements through the thickness
    # nc = 26; # number of elements per side
    # nr = 25; # number of elements per side
    # nt = 18; # number of elements through the thickness
    tolerance = Thickness/nt/1000.; # Geometrical tolerance

    fens,fes = Q4block(1.0, pi/2, nr, nc)
    #
    @test  nt % 2 == 0
    fens,fes  = H8extrudeQ4(fens, fes,
      nt, (xyz, layer)->[xyz[1], xyz[2], (layer)/nt*Thickness]);

    # Select the  boundary faces, on the boundary that is clamped,  and on the part
    # of the boundary that is loaded with the transverse pressure
    bdryfes = meshboundary(fes);
    exteriorbfl = selectelem(fens, bdryfes, box=[1.0, 1.0, 0.0, pi/2, 0.0, Thickness], inflate=tolerance);
    topbfl = selectelem(fens, bdryfes, box=[0.0, 1.0, 0.0, pi/2, Thickness, Thickness], inflate=tolerance);

    # Reshape the generated block into the elliptical plate
    for i=1:count(fens)
        r=fens.xyz[i,1]; a=fens.xyz[i,2]; z=fens.xyz[i,3]
        fens.xyz[i,:]=[(r*Ae+(1-r)*Ai)*cos(a) (r*Be+(1-r)*Bi)*sin(a) z];
    end


    geom = NodalField(fens.xyz)
    u = NodalField(zeros(size(fens.xyz,1),3)) # displacement field

    l12 =connectednodes(subset(bdryfes, exteriorbfl)) # external boundary
    setebc!(u, l12, true, 1, 0.0)
    setebc!(u, l12, true, 2, 0.0)
    ll = selectnode(fens; box=[0.0, Inf, 0.0, Inf, Thickness/2.0, Thickness/2.0], inflate = tolerance)
    l3 = intersect(ll, connectednodes(subset(bdryfes, exteriorbfl)))
    setebc!(u, l3, true, 3, 0.0)
    l1 =selectnode(fens; box=[0.0, 0.0, 0.0, Inf, 0.0, Thickness], inflate = tolerance)
    setebc!(u,l1,true, 1, 0.0) # symmetry plane X = 0
    l2 =selectnode(fens; box=[0.0, Inf, 0.0, 0.0, 0.0, Thickness], inflate = tolerance)
    setebc!(u,l2,true, 2, 0.0) # symmetry plane Y = 0

    applyebc!(u)
    numberdofs!(u)

    el1femm =  FEMMBase(IntegDomain(subset(bdryfes,topbfl), GaussRule(2, 2)))
    function pfun(forceout::FVec{T}, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt) where {T}
        forceout .=  [0.0, 0.0, -qmagn]
        return forceout
    end
    fi = ForceIntensity(FFlt, 3, pfun);
    F2 = distribloads(el1femm, geom, u, fi, 2);

    # Note that the material object needs to be created with the proper
    # model-dimension reduction in mind.  In this case that is the fully three-dimensional solid.
    MR = DeforModelRed3D

    material = MatDeforElastIso(MR, E, nu)

    femm = FEMMDeforLinearMSH8(MR, IntegDomain(fes, GaussRule(3, 2)), material)

    # The geometry field now needs to be associated with the FEMM
    femm = associategeometry!(femm, geom)

    K = stiffness(femm, geom, u)
    K = cholesky(K)
    U = K\(F2)
    scattersysvec!(u, U[:])

    nl = selectnode(fens, box=[Ai,Ai,0,0,Thickness,Thickness],inflate=tolerance);
    thecorneru = zeros(FFlt,1,3)
    gathervalues_asmat!(u, thecorneru, nl);
    thecorneru = thecorneru/phun("mm")
    # println("displacement =$(thecorneru) [MM] as compared to reference [-0.030939, 0, -0.10488] [MM]")


    fld= fieldfromintegpoints(femm, geom, u, :Cauchy, 2; nodevalmethod = :averaging, reportat = :extraptrend)#
    # println("Sigma_y =$(fld.values[nl,1][1]/phun("MPa")) as compared to reference sigma_yP = $(sigma_yP/phun("MPa")) [MPa]")

    # println("$((nc, nr, nt)), $(fld.values[nl,1][1]/phun("MPa"))")
# println("$(fld.values[nl,1][1]/phun("MPa"))")
    @test abs(fld.values[nl,1][1]/phun("MPa") - -4.627214556813842) < 1.0e-3

    # File =  "LE10NAFEMS_sigmay.vtk"
    # vtkexportmesh(File, fes.conn, geom.values,
    #                FinEtools.MeshExportModule.H8; vectors=[("u", u.values)],
    #                scalars=[("sigmay", fld.values)])
    # @async run(`"paraview.exe" $File`)
    # true


    AE = AbaqusExporter("LE10NAFEMS_H8");
    HEADING(AE, "LE10NAFEMS: Transverse deflection of elliptical plate with elliptical hole.");
    PART(AE, "part1");
    END_PART(AE);
    ASSEMBLY(AE, "ASSEM1");
    INSTANCE(AE, "INSTNC1", "PART1");
    NODE(AE, fens.xyz);
    ELEMENT(AE, "c3d8rh", "AllElements", 1, connasarray(femm.integdomain.fes))
    ELEMENT(AE, "SFM3D4", "TractionElements",
    1+count(femm.integdomain.fes), connasarray(el1femm.integdomain.fes))
    NSET_NSET(AE, "l1", l1)
    NSET_NSET(AE, "l2", l2)
    NSET_NSET(AE, "l3", l3)
    NSET_NSET(AE, "l12", l12)
    ORIENTATION(AE, "GlobalOrientation", vec([1. 0 0]), vec([0 1. 0]));
    SOLID_SECTION(AE, "elasticity", "GlobalOrientation", "AllElements", "Hourglassctl");
    SURFACE_SECTION(AE, "TractionElements")
    END_INSTANCE(AE);
    END_ASSEMBLY(AE);
    MATERIAL(AE, "elasticity")
    ELASTIC(AE, E, nu)
    SECTION_CONTROLS(AE, "Hourglassctl", "HOURGLASS=ENHANCED")
    STEP_PERTURBATION_STATIC(AE)
    BOUNDARY(AE, "ASSEM1.INSTNC1.", u.is_fixed, u.fixed_values)
    DLOAD(AE, "ASSEM1.INSTNC1.TractionElements", vec([0.0, 0.0, -qmagn]))
    END_STEP(AE)
    close(AE)
    lines = read(AE.filename)
    @test length(lines) - 10270 == 0
    try rm(AE.filename) catch end

end
end
using .mmmLE10expiAbaqus2mmmm
mmmLE10expiAbaqus2mmmm.test()

module mplate_w_hole_RECT_MSH8m
using FinEtools
using FinEtools.MeshExportModule
# using DataFrames
# using CSV
using Test
import LinearAlgebra: norm, cholesky, cross
import Statistics: mean
function test()
    E = 210000*phun("MEGA*PA");# 210e3 MPa
    nu = 0.3;
    Ri= 0.1*phun("M"); # hole radius
    Re = 2*Ri; # outer radius
    H = 0.01*phun("M") # thickness of the plate
    nRadial, nCircumferential=6, 3;
    sigma0=1*phun("MEGA*PA");

    function sigmaxx(x)
        local r = norm(x[1:2]);
        local th = atan(x[2],x[1]);
        return sigma0*(1-Ri^2/r^2*(3/2*cos(2*th)+cos(4*th))+3/2*Ri^4/r^4*cos(4*th));
    end
    function sigmayy(x)
        local r=norm(x[1:2]);
        local th = atan(x[2],x[1]);
        return -sigma0*(Ri^2/r^2*(1/2*cos(2*th)-cos(4*th))+3/2*Ri^4/r^4*cos(4*th));
    end
    function sigmaxy(x)
        local r=norm(x[1:2]);
        local th = atan(x[2],x[1]);
        return -sigma0*(Ri^2/r^2*(1/2*sin(2*th)+sin(4*th))-3/2*Ri^4/r^4*sin(4*th));
    end
    function sigmarr(x)
        local r = norm(x[1:2]);
        local th = atan(x[2],x[1]);
        return sigma0/2*(1-Ri^2/r^2) + sigma0/2*(1-4*Ri^2/r^2+3*Ri^4/r^4)*cos(2*th)
    end
    function sigmatt(x)
        local r = norm(x[1:2]);
        local th = atan(x[2],x[1]);
        return sigma0/2*(1+Ri^2/r^2) - sigma0/2*(1+3*Ri^4/r^4)*cos(2*th)
    end
    function sigmart(x)
        local r = norm(x[1:2]);
        local th = atan(x[2],x[1]);
        return -sigma0/2*(1+2*Ri^2/r^2-3*Ri^4/r^4)*sin(2*th)
    end

    sigyderrs = Dict{Symbol, FFltVec}()

    nelems = []
    for extrapolation in [:extrapmean]
        sigyderrs[extrapolation] = FFltVec[]
        nelems = []
        for ref in [1]
            # Thickness = H
            Thickness = H/2^ref
            tolerance = Thickness/2^ref/1000.; # Geometrical tolerance

            fens,fes = H8elliphole(Ri, Ri, Re, Re, Thickness,
                2^ref*nCircumferential, 2^ref*nCircumferential, 2^ref*nRadial, 1)
            @test count(fes) == 144
            # File =  "a.vtk"
            # vtkexportmesh(File, fes.conn, fens.xyz,
            #     FinEtools.MeshExportModule.H8)
            # @async run(`"paraview.exe" $File`)

            geom = NodalField(fens.xyz)
            u = NodalField(zeros(size(fens.xyz,1),3)) # displacement field

            l1 =selectnode(fens; box=[0.0, Inf, 0.0, 0.0, 0.0, Thickness], inflate = tolerance)
            setebc!(u,l1,true, 2, 0.0)
            l1 =selectnode(fens; box=[0.0, 0.0, 0.0, Inf, 0.0, Thickness], inflate = tolerance)
            setebc!(u,l1,true, 1, 0.0)
            l1 =selectnode(fens; box=[0.0, Inf, 0.0, Inf, 0.0, 0.0], inflate = tolerance)
            setebc!(u,l1,true, 3, 0.0)
            # l1 =selectnode(fens; box=[0.0, Inf, 0.0, Inf, Thickness, Thickness], inflate = tolerance)
            # setebc!(u,l1,true, 3, 0.0)

            applyebc!(u)
            numberdofs!(u)


            bdryfes = meshboundary(fes);
            ixl = selectelem(fens, bdryfes, plane=[1.0, 0.0, 0.0, Re], thickness=tolerance);
            elxfemm =  FEMMBase(IntegDomain(subset(bdryfes,ixl), GaussRule(2, 2)))
            function pfunx(forceout::FVec{T}, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt) where {T}
                forceout[1] = sigmaxx(XYZ)
                forceout[2] = sigmaxy(XYZ)
                forceout[3] = 0.0
                return forceout
            end
            fi = ForceIntensity(FFlt, 3, pfunx);
            Fx = distribloads(elxfemm, geom, u, fi, 2);
            iyl = selectelem(fens, bdryfes, plane=[0.0, 1.0, 0.0, Re], thickness=tolerance);
            elyfemm =  FEMMBase(IntegDomain(subset(bdryfes,iyl), GaussRule(2, 2)))
            function pfuny(forceout::FVec{T}, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt) where {T}
                forceout[1] = -sigmaxy(XYZ)
                forceout[2] = sigmayy(XYZ)
                forceout[3] = 0.0
                return forceout
            end
            fi = ForceIntensity(FFlt, 3, pfuny);
            Fy = distribloads(elyfemm, geom, u, fi, 2);

            MR = DeforModelRed3D

            material = MatDeforElastIso(MR, E, nu)

            femm = FEMMDeforLinearMSH8(MR, IntegDomain(fes, GaussRule(3, 2)), material)

            # The geometry field now needs to be associated with the FEMM
            femm = associategeometry!(femm, geom)

            K = stiffness(femm, geom, u)
            K = cholesky(K)
            U = K\(Fx + Fy)
            scattersysvec!(u, U[:])

            nlA = selectnode(fens, box=[Ri, Ri, 0.0, 0.0, 0.0, Thickness], inflate=tolerance);
            nlB = selectnode(fens, box=[0.0, 0.0, Ri, Ri, 0.0, Thickness], inflate=tolerance);
            # thecorneru = zeros(FFlt,length(nlA),3)
            # gathervalues_asmat!(u, thecorneru, nl);
            # thecorneru = mean(thecorneru, 1)[1]/phun("mm")
            # println("displacement = $(thecorneru) vs -0.10215 [MM]")

            # println("Extrapolation: $( extrapolation )")
            sigx = fieldfromintegpoints(femm, geom, u, :Cauchy, 1;
                nodevalmethod = :averaging, reportat = extrapolation)
            sigy = fieldfromintegpoints(femm, geom, u, :Cauchy, 2;
                nodevalmethod = :averaging, reportat = extrapolation)
            sigyA = mean(sigy.values[nlA,1], dims = 1)[1]
            sigyAtrue = sigmatt([Ri, 0.0, 0.0])
            # println("sig_y@A =$(sigyA/phun("MPa")) vs $(sigyAtrue/phun("MPa")) [MPa]")
            @test abs(sigyA/phun("MPa") - -0.8521990950600441) < 1.0e-3
            sigxB = mean(sigx.values[nlB,1], dims = 1)[1]
            sigxBtrue = sigmatt([0.0, Ri, 0.0])
            # println("sig_x@B =$(sigxB/phun("MPa")) vs $(sigxBtrue/phun("MPa")) [MPa]")
            @test abs(sigxB/phun("MPa") - 2.7749827820003374) < 1.0e-3
            # println("$extrapolation, $(count(fes)), $(sigyd/phun("MPa"))")
            # push!(nelems, count(fes))
            # push!(sigyderrs[extrapolation], abs(sigyd/sigma_yD - 1.0))
            # File =  "a.vtk"
            # vtkexportmesh(File, fes.conn, geom.values,
            # FinEtools.MeshExportModule.H8; vectors=[("u", u.values)],
            # scalars=[("sigmax", sigx.values/phun("MEGA*PA")),
            # ("sigmay", sigy.values/phun("MEGA*PA"))])
            # @async run(`"paraview.exe" $File`)
        end
    end

    # df = DataFrame(nelems=vec(nelems),
    #     sigyderrtrend=vec(sigyderrs[:extraptrend]),
    #     sigyderrdefault=vec(sigyderrs[:extrapmean]))
    # File = "LE1NAFEMS_MSH8_convergence.CSV"
    # CSV.write(File, df)
    # @async run(`"paraview.exe" $File`)

end
end
using .mplate_w_hole_RECT_MSH8m
mplate_w_hole_RECT_MSH8m.test()

module mplate_w_hole_RECT_H20m
using FinEtools
using FinEtools.MeshExportModule
using FinEtools.MeshImportModule: import_ABAQUS
# using DataFrames
# using CSV
using Test
import LinearAlgebra: norm, cholesky, cross
import Statistics: mean
function test()
    E = 210000*phun("MEGA*PA");# 210e3 MPa
    nu = 0.3;
    Ri= 0.15*phun("M"); # hole radius
    Re = 2*Ri; # outer radius
    H = 0.01*phun("M") # thickness of the plate
    nRadial, nCircumferential=6, 3;
    sigma0=1*phun("MEGA*PA");

    function sigmaxx(x)
        local r = norm(vec(x[1:2]));
        local th = atan(x[2],x[1]);
        return sigma0*(1-Ri^2/r^2*(3/2*cos(2*th)+cos(4*th))+3/2*Ri^4/r^4*cos(4*th));
    end
    function sigmayy(x)
        local r = norm(vec(x[1:2]));
        local th = atan(x[2],x[1]);
        return -sigma0*(Ri^2/r^2*(1/2*cos(2*th)-cos(4*th))+3/2*Ri^4/r^4*cos(4*th));
    end
    function sigmaxy(x)
        local r = norm(vec(x[1:2]));
        local th = atan(x[2],x[1]);
        return -sigma0*(Ri^2/r^2*(1/2*sin(2*th)+sin(4*th))-3/2*Ri^4/r^4*sin(4*th));
    end
    function sigmarr(x)
        local r = norm(vec(x[1:2]));
        local th = atan(x[2],x[1]);
        return sigma0/2*(1-Ri^2/r^2) + sigma0/2*(1-4*Ri^2/r^2+3*Ri^4/r^4)*cos(2*th)
    end
    function sigmatt(x)
        local r = norm(vec(x[1:2]));
        local th = atan(x[2],x[1]);
        return sigma0/2*(1+Ri^2/r^2) - sigma0/2*(1+3*Ri^4/r^4)*cos(2*th)
    end
    function sigmart(x)
        local r = norm(vec(x[1:2]));
        local th = atan(x[2],x[1]);
        return -sigma0/2*(1+2*Ri^2/r^2-3*Ri^4/r^4)*sin(2*th)
    end

    sigyderrs = Dict{Symbol, FFltVec}()

    nelems = []
    for extrapolation in [:extrapmean]
        sigyderrs[extrapolation] = FFltVec[]
        nelems = []
        for ref in [1]
            Thickness = H
            # Thickness = H/2^ref
            tolerance = Thickness/2^ref/1000.; # Geometrical tolerance

            fens,fes = H8elliphole(Ri, Ri, Re, Re, Thickness,
            2^ref*nCircumferential, 2^ref*nCircumferential, 2^ref*nRadial, 1)
            fens,fes = H8toH20(fens,fes)
            # File =  "a.vtk"
            # vtkexportmesh(File, fes.conn, fens.xyz,
            #     FinEtools.MeshExportModule.H20)
            # @async run(`"paraview.cexe" $File`)

            # println("My mesh=>$((count(fens), count(fes)))")
            @test count(fens) == 1131
            @test count(fes) == 144
            #
            # output = import_ABAQUS("plane_w_hole_m_debug.inp")
            # fens1,fes1 = output["fens"], output["fesets"][1]
            # println("Matlab mesh=>$((count(fens1), count(fes1[1])))")
            #
            #  fens3, newfes1, fes2 = mergemeshes(fens,fes, fens1,fes1[1], tolerance)
            #  fes3 = cat(2, newfes1)
            #  println("Merged mesh=>$((count(fens3), count(fes3)))")

            geom = NodalField(fens.xyz)
            u = NodalField(zeros(size(fens.xyz,1),3)) # displacement field

            l1 =selectnode(fens; box=[0.0, Inf, 0.0, 0.0, 0.0, Thickness], inflate = tolerance)
            setebc!(u,l1,true, 2, 0.0)
            l1 =selectnode(fens; box=[0.0, 0.0, 0.0, Inf, 0.0, Thickness], inflate = tolerance)
            setebc!(u,l1,true, 1, 0.0)
            l1 =selectnode(fens; box=[0.0, Inf, 0.0, Inf, 0.0, 0.0], inflate = tolerance)
            setebc!(u,l1,true, 3, 0.0)
            # l1 =selectnode(fens; box=[0.0, Inf, 0.0, Inf, Thickness, Thickness], inflate = tolerance)
            # setebc!(u,l1,true, 3, 0.0)

            applyebc!(u)
            numberdofs!(u)


            bdryfes = meshboundary(fes);
            # ixl = selectelem(fens, bdryfes, plane=[1.0, 0.0, 0.0, Re], thickness=tolerance);
            ixl = selectelem(fens, bdryfes, box=[Re, Re, -Inf, +Inf, -Inf, +Inf], inflate = tolerance);
            elxfemm =  FEMMBase(IntegDomain(subset(bdryfes,ixl), GaussRule(2, 2)))
            function pfunx(forceout::FVec{T}, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt) where {T}
                forceout[1] = sigmaxx(XYZ)
                forceout[2] = sigmaxy(XYZ)
                forceout[3] = 0.0
                return forceout
            end
            fi = ForceIntensity(FFlt, 3, pfunx);
            Fx = distribloads(elxfemm, geom, u, fi, 2);
            # iyl = selectelem(fens, bdryfes, plane=[0.0, 1.0, 0.0, Re], thickness=tolerance);
            iyl = selectelem(fens, bdryfes, box=[-Inf, +Inf, Re, Re, -Inf, +Inf], inflate = tolerance);
            elyfemm =  FEMMBase(IntegDomain(subset(bdryfes,iyl), GaussRule(2, 2)))
            function pfuny(forceout::FVec{T}, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt) where {T}
                forceout[1] = sigmaxy(XYZ)
                forceout[2] = sigmayy(XYZ)
                forceout[3] = 0.0
                return forceout
            end
            fi = ForceIntensity(FFlt, 3, pfuny);
            Fy = distribloads(elyfemm, geom, u, fi, 2);

            MR = DeforModelRed3D

            material = MatDeforElastIso(MR, E, nu)

            femm = FEMMDeforLinear(MR, IntegDomain(fes, GaussRule(3, 2)), material)

            # The geometry field now needs to be associated with the FEMM
            femm = associategeometry!(femm, geom)

            K = stiffness(femm, geom, u)
            K = cholesky(K)
            U = K\(Fx + Fy)
            scattersysvec!(u, U[:])
            # println("oof load = $(norm(Fx + Fy, 2))")
            @test abs(norm(Fx + Fy, 2) - 883.437848042617) < 1.0e-2

            nlA = selectnode(fens, box=[Ri, Ri, 0.0, 0.0, 0.0, 00.0], inflate=tolerance);
            pointu = zeros(FFlt,length(nlA),3)
            gathervalues_asmat!(u, pointu, nlA);
            # println("disp@A = $(pointu/phun("mm")) [MM]")
            @test norm(pointu/phun("mm") - [0.00213238 0.0 0.0]) < 1.0e-4
            nlB = selectnode(fens, box=[0.0, 0.0, Ri, Ri, 0.0, 0.0], inflate=tolerance);
            pointu = zeros(FFlt,length(nlB),3)
            gathervalues_asmat!(u, pointu, nlB);
            # println("disp@B = $(pointu/phun("mm")) [MM]")
            @test norm(pointu/phun("mm") - [0.0 -0.000708141 0.0]) < 1.0e-4
            nlC = selectnode(fens, box=[Re, Re, Re, Re, Thickness, Thickness], inflate=tolerance);
            pointu = zeros(FFlt,length(nlC),3)
            gathervalues_asmat!(u, pointu, nlC);
            # println("disp@C = $(pointu/phun("mm")) [MM]")
            @test norm(pointu/phun("mm") - [0.00168556 -0.000455007 -1.4286e-5]) < 1.0e-4

            nlAallz = selectnode(fens, box=[Ri, Ri, 0.0, 0.0, 0.0, Thickness], inflate=tolerance);
            nlBallz = selectnode(fens, box=[0.0, 0.0, Ri, Ri, 0.0, Thickness], inflate=tolerance);
            sigx = fieldfromintegpoints(femm, geom, u, :Cauchy, 1;
                nodevalmethod = :averaging, reportat = extrapolation)
            sigy = fieldfromintegpoints(femm, geom, u, :Cauchy, 2;
                nodevalmethod = :averaging, reportat = extrapolation)
            sigyA = mean(sigy.values[nlAallz,1], dims = 1)[1]
            sigyAtrue = sigmayy([Ri, 0.0, 0.0])
            # println("sig_y@A =$(sigyA/phun("MPa")) vs $(sigyAtrue/phun("MPa")) [MPa]")
            @test abs(sigyA/phun("MPa") - -0.8513053526935438)/(sigyAtrue/phun("MPa")) < 1.0e-4
            sigxB = mean(sigx.values[nlBallz,1], dims = 1)[1]
            sigxBtrue = sigmaxx([0.0, Ri, 0.0])
            # println("sig_x@B =$(sigxB/phun("MPa")) vs $(sigxBtrue/phun("MPa")) [MPa]")
            # @test abs(sigxB/phun("MPa") - 2.789413093796375)/3.0 < 1.0e-4
            # println("$extrapolation, $(count(fes)), $(sigyd/phun("MPa"))")
            # push!(nelems, count(fes))
            # push!(sigyderrs[extrapolation], abs(sigyd/sigma_yD - 1.0))
            File =  "a.vtk"
            vtkexportmesh(File, connasarray(fes), geom.values,
                FinEtools.MeshExportModule.H20; vectors=[("u", u.values)],
                scalars=[("sigmax", sigx.values/phun("MEGA*PA")),
                ("sigmay", sigy.values/phun("MEGA*PA"))])
            # @async run(`"paraview.exe" $File`)
            try rm(File); catch end
        end
    end

    # df = DataFrame(nelems=vec(nelems),
    #     sigyderrtrend=vec(sigyderrs[:extraptrend]),
    #     sigyderrdefault=vec(sigyderrs[:extrapmean]))
    # File = "LE1NAFEMS_MSH8_convergence.CSV"
    # CSV.write(File, df)
    # @async run(`"paraview.exe" $File`)

end
end
using .mplate_w_hole_RECT_H20m
mplate_w_hole_RECT_H20m.test()


module mplate_w_hole_MST10m
using FinEtools
using FinEtools.MeshExportModule
using Test
import LinearAlgebra: norm, cholesky, cross
import Statistics: mean
function test()
    E = 2.4*phun("MEGA*PA");# 210e3 MPa
    nu = 0.49995;
    Re = 0.3*phun("M"); # outer radius
    Ri= 0.1*phun("M"); # hole radius
    H = 0.1*phun("M") # thickness of the plate
    nRadial, nCircumferential=3, 5;
    sigma0=1*phun("MEGA*PA");

    function sigmaxx(x)
        local r = norm(x[1:2]);
        local th = atan(x[2],x[1]);
        return sigma0*(1-Ri^2/r^2*(3/2*cos(2*th)+cos(4*th))+3/2*Ri^4/r^4*cos(4*th));
    end
    function sigmayy(x)
        local r=norm(x[1:2]);
        local th = atan(x[2],x[1]);
        return -sigma0*(Ri^2/r^2*(1/2*cos(2*th)-cos(4*th))+3/2*Ri^4/r^4*cos(4*th));
    end
    function sigmaxy(x)
        local r=norm(x[1:2]);
        local th = atan(x[2],x[1]);
        return -sigma0*(Ri^2/r^2*(1/2*sin(2*th)+sin(4*th))-3/2*Ri^4/r^4*sin(4*th));
    end
    function sigmarr(x)
        local r = norm(x[1:2]);
        local th = atan(x[2],x[1]);
        return sigma0/2*(1-Ri^2/r^2) + sigma0/2*(1-4*Ri^2/r^2+3*Ri^4/r^4)*cos(2*th)
    end
    function sigmatt(x)
        local r = norm(x[1:2]);
        local th = atan(x[2],x[1]);
        return sigma0/2*(1+Ri^2/r^2) - sigma0/2*(1+3*Ri^4/r^4)*cos(2*th)
    end
    function sigmart(x)
        local r = norm(x[1:2]);
        local th = atan(x[2],x[1]);
        return -sigma0/2*(1+2*Ri^2/r^2-3*Ri^4/r^4)*sin(2*th)
    end

    sigxderrs = Dict{Symbol, FFltVec}()
    sigyderrs = Dict{Symbol, FFltVec}()
    numelements = []
    numnodes = []
    for extrapolation in [:extrapmean] # :extraptrend
        sigxderrs[extrapolation] = FFltVec[]
        sigyderrs[extrapolation] = FFltVec[]
        numelements = []
        numnodes = []
        for ref in 1:1
            # Thickness = H
            Thickness = H/2^ref
            tolerance = Thickness/2^ref/1000.; # Geometrical tolerance

            fens,fes = T10block(1.0, pi/2, Thickness, 2^ref*nRadial, 2^ref*nCircumferential, 1)

            bdryfes = meshboundary(fes);
            icl = selectelem(fens, bdryfes, box=[1.0, 1.0, 0.0, pi/2, 0.0, Thickness], inflate=tolerance);

            for i=1:count(fens)
                t=fens.xyz[i,1]; a=fens.xyz[i,2]; z=fens.xyz[i,3]
                fens.xyz[i,:] = [(t*Re+(1-t)*Ri)*cos(a), (t*Re+(1-t)*Ri)*sin(a), z];
            end

            # File =  "a.vtk"
            # vtkexportmesh(File, fes.conn, fens.xyz,
            #     FinEtools.MeshExportModule.H8)
            # @async run(`"paraview.exe" $File`)

            geom = NodalField(fens.xyz)
            u = NodalField(zeros(size(fens.xyz,1),3)) # displacement field

            l1 =selectnode(fens; box=[0.0, Inf, 0.0, 0.0, 0.0, Thickness], inflate = tolerance)
            setebc!(u,l1,true, 2, 0.0)
            l1 =selectnode(fens; box=[0.0, 0.0, 0.0, Inf, 0.0, Thickness], inflate = tolerance)
            setebc!(u,l1,true, 1, 0.0)
            # Plane-stress constraint: assume the plane z=0 is the plane of symmetry of the plate
            l1 =selectnode(fens; box=[0.0, Inf, 0.0, Inf, 0.0, 0.0], inflate = tolerance)
            setebc!(u,l1,true, 3, 0.0)
            # If this was enabled, the plane-strain  constraint would be enforced.
            # l1 =selectnode(fens; box=[0.0, Inf, 0.0, Inf, Thickness, Thickness], inflate = tolerance)
            # setebc!(u,l1,true, 3, 0.0)

            applyebc!(u)
            numberdofs!(u)
            el1femm =  FEMMBase(IntegDomain(subset(bdryfes,icl), SimplexRule(2, 3)))
            function pfun(forceout::FVec{T}, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt) where {T}
                local r = sqrt(XYZ[1]^2 + XYZ[2]^2)
                nx = XYZ[1]/r; ny = XYZ[2]/r
                # local sx, sy, txy
                # sx, sy, txy = sigmaxx(XYZ), sigmayy(XYZ), sigmaxy(XYZ)
                # sn = sx * nx^2 + sy * ny^2 + 2 * nx * ny * txy
                # tn = -(sx - sy) * nx * ny + (nx^2 - ny^2) * txy
                # forceout[1] = sn * nx - tn * ny
                # forceout[2] = sn * ny + tn * nx
                # forceout[3] = 0.0
                forceout[1] = sigmarr(XYZ) * nx - sigmart(XYZ) * ny
                forceout[2] = sigmarr(XYZ) * ny + sigmart(XYZ) * nx
                forceout[3] = 0.0
                return forceout
            end
            fi = ForceIntensity(FFlt, 3, pfun);
            F2 = distribloads(el1femm, geom, u, fi, 2);

            MR = DeforModelRed3D

            material = MatDeforElastIso(MR, E, nu)

            femm = FEMMDeforLinearMST10(MR, IntegDomain(fes, SimplexRule(3, 4)), material)

            # The geometry field now needs to be associated with the FEMM
            femm = associategeometry!(femm, geom)

            K = stiffness(femm, geom, u)
            K = cholesky(K)
            U = K\(F2)
            scattersysvec!(u, U[:])

            nlA = selectnode(fens, box=[Ri, Ri, 0.0, 0.0, 0.0, Thickness], inflate=tolerance);
            nlB = selectnode(fens, box=[0.0, 0.0, Ri, Ri, 0.0, Thickness], inflate=tolerance);
            # thecorneru = zeros(FFlt,length(nlA),3)
            # gathervalues_asmat!(u, thecorneru, nl);
            # thecorneru = mean(thecorneru, 1)[1]/phun("mm")
            # println("displacement = $(thecorneru) vs -0.10215 [MM]")

            # println("Extrapolation: $( extrapolation )")
            sigx = fieldfromintegpoints(femm, geom, u, :Cauchy, 1;
                nodevalmethod = :averaging, reportat = extrapolation)
            sigy = fieldfromintegpoints(femm, geom, u, :Cauchy, 2;
                nodevalmethod = :averaging, reportat = extrapolation)
            sigyA = mean(sigy.values[nlA,1], dims = 1)[1]
            sigyAtrue = sigmatt([Ri, 0.0, 0.0])
            # println("sig_y@A =$(sigyA/phun("MPa")) vs $(sigyAtrue/phun("MPa")) [MPa]")
            @test  abs(sigyA/phun("MPa") - -0.6705333234697736)/(sigyAtrue/phun("MPa")) < 1.0e-4
            sigxB = mean(sigx.values[nlB,1], dims = 1)[1]
            sigxBtrue = sigmatt([0.0, Ri, 0.0])
            # println("sig_x@B =$(sigxB/phun("MPa")) vs $(sigxBtrue/phun("MPa")) [MPa]")
            @test  abs(sigxB/phun("MPa") - 2.301542874107758)/3.0 < 1.0e-4
            push!(numnodes, count(fens))
            push!(numelements, count(fes))
            push!(sigxderrs[extrapolation], abs(sigxB/sigxBtrue - 1.0))
            push!(sigyderrs[extrapolation], abs(sigyA/sigyAtrue - 1.0))
            # File =  "a.vtk"
            # vtkexportmesh(File, fes.conn, geom.values,
            # FinEtools.MeshExportModule.H8; vectors=[("u", u.values)],
            # scalars=[("sigmax", sigx.values/phun("MEGA*PA"))])
            # @async run(`"paraview.exe" $File`)
        end
    end

end
end
using .mplate_w_hole_MST10m
mplate_w_hole_MST10m.test()

module mmLE1NAFEMSsstress
using FinEtools
using Test
import LinearAlgebra: norm, cholesky, cross
function test()
    E = 210e3*phun("MEGA*PA");# 210e3 MPa
    nu = 0.3;
    p = 10*phun("MEGA*PA");# 10 MPA Outward pressure on the outside ellipse
    sigma_yD = 92.7*phun("MEGA*PA");# tensile stress at [2.0, 0.0] meters
    Radius = 1.0*phun("m")
    Thickness = 0.1*phun("m")
    n = 2; # number of elements per side
    tolerance = 1.0/n/1000.; # Geometrical tolerance

    fens,fes = Q4block(1.0, pi/2, n, n*2)
    fens,fes  = H8extrudeQ4(fens, fes,
      1, (xyz, layer)->[xyz[1], xyz[2], (layer)*Thickness]);

    bdryfes = meshboundary(fes);
    icl = selectelem(fens, bdryfes, box=[1.0, 1.0, 0.0, pi/2, 0.0, Thickness], inflate=tolerance);
    for i=1:count(fens)
        t=fens.xyz[i,1]; a=fens.xyz[i,2]; z=fens.xyz[i,3]
        fens.xyz[i,:]=[(t*3.25+(1-t)*2)*cos(a), (t*2.75+(1-t)*1)*sin(a), z];
    end


    geom = NodalField(fens.xyz)
    u = NodalField(zeros(size(fens.xyz,1),3)) # displacement field

    l1 =selectnode(fens; box=[0.0, Inf, 0.0, 0.0, 0.0, Thickness], inflate = tolerance)
    setebc!(u,l1,true, 2, 0.0)
    l1 =selectnode(fens; box=[0.0, 0.0, 0.0, Inf, 0.0, Thickness], inflate = tolerance)
    setebc!(u,l1,true, 1, 0.0)
    l1 =selectnode(fens; box=[0.0, Inf, 0.0, Inf, 0.0, 0.0], inflate = tolerance)
    setebc!(u,l1,true, 3, 0.0)

    applyebc!(u)
    numberdofs!(u)


    el1femm =  FEMMBase(IntegDomain(subset(bdryfes,icl), GaussRule(2, 2)))
    function pfun(forceout::FVec{T}, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt) where {T}
        pt= [2.75/3.25*XYZ[1], 3.25/2.75*XYZ[2], 0.0]
        forceout .=    vec(p*pt/norm(pt));
        return forceout
    end
    fi = ForceIntensity(FFlt, 3, pfun);
    F2 = distribloads(el1femm, geom, u, fi, 2);

    # Note that the material object needs to be created with the proper
    # model-dimension reduction in mind.  In this case that is the fully three-dimensional solid.
    MR = DeforModelRed3D

    material = MatDeforElastIso(MR, E, nu)

    femm = FEMMDeforLinearMSH8(MR, IntegDomain(fes, GaussRule(3, 2)), material)

    # The geometry field now needs to be associated with the FEMM
    femm = associategeometry!(femm, geom)

    K = stiffness(femm, geom, u)
    K = cholesky(K)
    U = K\(F2)
    scattersysvec!(u, U[:])

    nl = selectnode(fens, box=[2.0, 2.0, 0.0, 0.0, 0.0, 0.0],inflate=tolerance);
    thecorneru = zeros(FFlt,1,3)
    gathervalues_asmat!(u, thecorneru, nl);
    thecorneru = thecorneru/phun("mm")
    # println("displacement =$(thecorneru) [MM] as compared to reference [-0.10215,0] [MM]")
    @test norm(thecorneru - [-0.107276 0.0 0.0]) < 1.0e-4


    fld= fieldfromintegpoints(femm, geom, u, :Cauchy, 2;
        nodevalmethod = :invdistance, reportat = :meanonly)
    # println("Sigma_y =$(fld.values[nl,1][1]/phun("MPa")) as compared to reference sigma_yD = $(sigma_yD/phun("MPa")) [MPa]")
    @test abs(fld.values[nl,1][1]/phun("MPa") - 42.54884174624546) < 1.0e-3
    fld= fieldfromintegpoints(femm, geom, u, :Cauchy, 2;
        nodevalmethod = :averaging, reportat = :extrapmean)
    # println("Sigma_y =$(fld.values[nl,1][1]/phun("MPa")) as compared to reference sigma_yD = $(sigma_yD/phun("MPa")) [MPa]")
    @test abs(fld.values[nl,1][1]/phun("MPa") - 42.54884174624546) < 1.0e-3
    fld= fieldfromintegpoints(femm, geom, u, :Cauchy, 2;
        nodevalmethod = :averaging, reportat = :extraptrend)
    # println("Sigma_y =$(fld.values[nl,1][1]/phun("MPa")) as compared to reference sigma_yD = $(sigma_yD/phun("MPa")) [MPa]")
    @test abs(fld.values[nl,1][1]/phun("MPa") - 45.44562958746983) < 1.0e-3

    fld= fieldfromintegpoints(femm, geom, u, :Cauchy, 2;
        nodevalmethod = :invdistance, reportat = :meanonly)
    # println("Sigma_y =$(fld.values[nl,1][1]/phun("MPa")) as compared to reference sigma_yD = $(sigma_yD/phun("MPa")) [MPa]")
    @test abs(fld.values[nl,1][1]/phun("MPa") - 42.54884174624546) < 1.0e-3

end
end
using .mmLE1NAFEMSsstress
mmLE1NAFEMSsstress.test()

module mocylpullFun
using FinEtools
using FinEtools.AlgoDeforLinearModule
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

    # the symmetry plane
    la1 =selectnode(fens; box=[0 rex 0 0], inflate = tolerance)
    # The other end
    la2 =selectnode(fens; box=[0 rex Length Length], inflate = tolerance)

    e1 = FDataDict("node_list"=>la1, "component"=>2, "displacement"=>x -> 0.0)
    e2 = FDataDict("node_list"=>la2, "component"=>2, "displacement"=>x -> ua)

    # Property and material
    material=MatDeforElastOrtho(MR, E1,E2,E3,nu12,nu13,nu23,G12,G13,G23)

    femm = FEMMDeforLinear(MR, IntegDomain(fes, GaussRule(2, 2), true), material)

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
    # println("Minimum/maximum = $(minimum(fld.values))/$(maximum(fld.values))")
    @test abs(minimum(fld.values) - -0.050318853446676574) < 1.0e-4
    @test abs(maximum(fld.values) - -0.04973951673608955) < 1.0e-4
    # File =  "orthoballoon_sigmaz.vtk"
    # vtkexportmesh(File, fens, fes; scalars=[("sigmaz", fld.values)],
    #               vectors=[("u", u.values)])
    # @async run(`"paraview.exe" $File`)
end
end
using .mocylpullFun
mocylpullFun.test()

module mmLE11malgo
using FinEtools
using FinEtools.AlgoDeforLinearModule
using Test
import LinearAlgebra: norm, cholesky, cross
function test()


    # NAFEMS LE11 benchmark with Q8 elements.
    # # This is a test recommended by the National Agency for Finite Element
    # # Methods and Standards (U.K.): Test LE11 from NAFEMS Publication TNSB,
    # # Rev. 3, “The Standard NAFEMS Benchmarks,” October 1990.
    # #
    # # Target solution: Direct stress,   =  –105 MPa at point A.
    #function  LE11NAFEMS()
    # Parameters:
    Ea =  210000*phun("MEGA*Pa")
    nua =  0.3;
    alphaa = 2.3e-4;              # thermal expansion coefficient
    sigmaA = -105*phun("MEGA*Pa")
    nref =  1;                        # how many times should we refine the mesh?
    X = [1.     0.;#A
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
    tolerance  = 1.e-6*phun("M")
    ##
    # Note that the material object needs to be created with the proper
    # model-dimension reduction in mind.  In this case that is the axial symmetry
    # assumption.
    MR  =  DeforModelRed2DAxisymm



    fens = FENodeSet(X);
    fes = FESetQ4([1 2 4 3; 3 4 6 5; 5 6 8 7; 7 8 10 9; 9 10 12 11; 11 12 14 13; 13 14 16 15]);
    for ref = 1:nref
      fens,fes = Q4refine(fens,fes);
      list = selectnode(fens,distance = 1.0+0.1/2^nref, from = [0. 0.], inflate = tolerance);
      fens.xyz[list,:] =  FinEtools.MeshUtilModule.ontosphere(fens.xyz[list,:],1.0);
    end
    fens,fes = Q4toQ8(fens,fes)
    list = selectnode(fens,distance = 1.0+0.1/2^nref, from = [0. 0.], inflate = tolerance);
    fens.xyz[list,:] =  FinEtools.MeshUtilModule.ontosphere(fens.xyz[list,:],1.0);


    # EBC's
    l1 = selectnode(fens,box = [-Inf Inf 0 0],inflate = tolerance)
    e1 = FDataDict("node_list"=>l1, "component"=>2, "displacement"=>x -> 0.0)
    l1 = selectnode(fens,box = [-Inf Inf 1.79  1.79],inflate = tolerance)
    e2 = FDataDict("node_list"=>l1, "component"=>2, "displacement"=>x -> 0.0)

    # Temperature field
    dtemp = FDataDict("temperature"=>x -> x[1] + x[2])

    # Property and material
    material = MatDeforElastIso(MR, 0.0, Ea, nua, alphaa)

    femm  =  FEMMDeforLinear(MR, IntegDomain(fes, GaussRule(2, 3), true), material)

    # Make region 1
    region = FDataDict("femm"=>femm)
    # Make model data
    modeldata =  FDataDict(
    "fens"=> fens, "regions"=>  [region],
    "essential_bcs"=>[e1, e2], "temperature_change"=>dtemp)

    # Call the solver
    modeldata = AlgoDeforLinearModule.linearstatics(modeldata)
    geom = modeldata["geom"]
    u = modeldata["u"]
    dT = modeldata["temp"]

    nA  = selectnode(fens,box = FFlt[1.0  1.0 0.0 0.0], inflate = tolerance);

    fld =  fieldfromintegpoints(femm, geom, u, dT, :Cauchy, 2)


    File  =   "LE11NAFEMS_Q8_sigmay.vtk"
    vtkexportmesh(File, fens, fes; scalars = [("sigmay", fld.values)],
        vectors = [("u", u.values)])
    # println("range of  sigmay = $((minimum(fld.values), maximum(fld.values)))")
    @test norm([minimum(fld.values), maximum(fld.values)] - [-1.443052182185007e8, -1.4106181545272522e7]) < 1.0e-1
        # @async run(`"paraview.exe" $File`)
    try rm(File)   catch end

    sA  =  fld.values[nA]/phun("MEGA*Pa")
    sAn  =  fld.values[nA]/sigmaA
    # println("Stress at point A: $(sA) i. e.  $( sAn*100  )% of reference value")
    @test norm(sA .- -93.8569) < 1.0e-4

    fen2fe  = FENodeToFEMap(connasarray(fes), nnodes(geom))
    function inspector(idat, elnum, conn, xe,  out,  xq)
    #   println("loc = $(  xq  ) : $(  transpose(out)/phun("MEGA*Pa")  )")
      return idat
    end

    inspectintegpoints(femm, geom, u, dT,  fen2fe.map[nA[1]],
      inspector, []; quantity = :Cauchy)

end
end
using .mmLE11malgo
mmLE11malgo.test()

module mmtwistedmsh8ort
using FinEtools
using FinEtools.AlgoDeforLinearModule
using Test
import LinearAlgebra: norm, cholesky, cross
import Statistics: mean
function test()

  # println("""
  # The initially twisted cantilever beam is one of the standard test
  # problems for verifying the finite-element accuracy [1]. The beam is
  #   clamped at one end and loaded either with unit in-plane or
  #   unit out-of-plane force at the other. The centroidal axis of the beam is
  #   straight at the undeformed  configuration, while its cross-sections are
  #   twisted about the centroidal axis from 0 at the clamped end to pi/2 at
  #   the free end.
  #
  #   Reference:
  #   Zupan D, Saje M (2004) On "A proposed standard set of problems to test
  #   finite element accuracy": the twisted beam. Finite Elements in Analysis
  #   and Design 40: 1445-1451.
  #   """)
  E1s =   E2s =  E3s = 0.29e8
  nu12s = nu13s = nu23s = 0.22
  G12s = G13s = G23s = E1s/2/(1+nu12s)
  # E = 0.29e8;
  # nu = 0.22;
  W = 1.1;
  L = 12.;
  t =  0.32;
  nl = 2; nt = 1; nw = 1; ref = 7;
  p =   1/W/t;
  #  Loading in the Z direction
  loadv = [0;0;p]; dir = 3; uex = 0.005424534868469; # Harder: 5.424e-3;
  #   Loading in the Y direction
  #loadv = [0;p;0]; dir = 2; uex = 0.001753248285256; # Harder: 1.754e-3;
  tolerance  = t/1000;

  fens,fes  = H8block(L,W,t, nl*ref,nw*ref,nt*ref)

  # Reshape into a twisted beam shape
  for i = 1:count(fens)
    a = fens.xyz[i,1]/L*(pi/2); y = fens.xyz[i,2]-(W/2); z = fens.xyz[i,3]-(t/2);
    fens.xyz[i,:] = [fens.xyz[i,1],y*cos(a)-z*sin(a),y*sin(a)+z*cos(a)];
  end

  # Clamped end of the beam
  l1  = selectnode(fens; box = [0 0 -100*W 100*W -100*W 100*W], inflate  =  tolerance)
  e1 = FDataDict("node_list"=>l1, "component"=>1, "displacement"=>0.0)
  e2 = FDataDict("node_list"=>l1, "component"=>2, "displacement"=>0.0)
  e3 = FDataDict("node_list"=>l1, "component"=>3, "displacement"=>0.0)

  # Traction on the opposite edge
  boundaryfes  =   meshboundary(fes);
  Toplist   = selectelem(fens,boundaryfes, box =  [L L -100*W 100*W -100*W 100*W], inflate =   tolerance);
  el1femm  = FEMMBase(IntegDomain(subset(boundaryfes,Toplist), GaussRule(2, 2)))
  flux1 = FDataDict("femm"=>el1femm, "traction_vector"=>loadv)


  # Make the region
  MR = DeforModelRed3D
  material = MatDeforElastOrtho(MR,
    0.0, E1s, E2s, E3s,
    nu12s, nu13s, nu23s,
    G12s, G13s, G23s,
    0.0, 0.0, 0.0)
  region1 = FDataDict("femm"=>FEMMDeforLinearMSH8(MR, IntegDomain(fes, GaussRule(3,2)),
  material))

  # Make model data
  modeldata =  FDataDict(
  "fens"=> fens, "regions"=>  [region1],
  "essential_bcs"=>[e1, e2, e3], "traction_bcs"=>  [flux1])

  # Call the solver
  modeldata = AlgoDeforLinearModule.linearstatics(modeldata)
  geom = modeldata["geom"]
  u = modeldata["u"]

  # Extract the solution
  nl = selectnode(fens, box = [L L -100*W 100*W -100*W 100*W],inflate = tolerance);
  theutip = mean(u.values[nl,:], dims = 1)
  # println("displacement  = $(theutip[dir]) as compared to converged $uex")
  @test abs(theutip[dir]-uex)/uex < 0.0012

  # Write out mesh with displacements
  modeldata["postprocessing"] = FDataDict("file"=>"twisted_beam_msh8")
  modeldata = AlgoDeforLinearModule.exportdeformation(modeldata)
  try rm(modeldata["postprocessing"]["exported"][1]["file"]); catch end

  # Write out mesh with stresses
  modeldata["postprocessing"] = FDataDict("file"=>"twisted_beam_msh8",
  "quantity"=> :Cauchy, "component"=> :xy)
  modeldata = AlgoDeforLinearModule.exportstress(modeldata)
  try rm(modeldata["postprocessing"]["exported"][1]["file"]); catch end

  # Write out mesh with stresses
  modeldata["postprocessing"] = FDataDict("file"=>"twisted_beam_msh8",
  "quantity"=> :Cauchy, "component"=> :xz)
  modeldata = AlgoDeforLinearModule.exportstress(modeldata)
  try rm(modeldata["postprocessing"]["exported"][1]["file"]); catch end

  # Write out mesh with von Mises stresses
  modeldata["postprocessing"] = FDataDict("file"=>"twisted_beam_msh8",
  "quantity"=> :vm)
  modeldata = AlgoDeforLinearModule.exportstress(modeldata)
  vm  = modeldata["postprocessing"]["exported"][1]["field"]
  # println("extremes of vm, nodal: $([minimum(vm.values),   maximum(vm.values)])")
  try rm(modeldata["postprocessing"]["exported"][1]["file"]); catch end
  @test norm([minimum(vm.values),   maximum(vm.values)]-[4.78774, 522.126]) < 0.01

  # Write out mesh with von Mises stresses, elementwise
  modeldata["postprocessing"] = FDataDict("file"=>"twisted_beam_msh8-ew",
  "quantity"=> :vm)
  modeldata = AlgoDeforLinearModule.exportstresselementwise(modeldata)
  vm  = modeldata["postprocessing"]["exported"][1]["field"]
  # println("extremes of vm, elemental: $([minimum(vm.values),   maximum(vm.values)])")
  try rm(modeldata["postprocessing"]["exported"][1]["file"]); catch end
  @test norm([minimum(vm.values),   maximum(vm.values)]-[1.85882, 522.126]) < 0.01

  # Write out mesh with von Mises stresses, elementwise
  modeldata["postprocessing"] = FDataDict("file"=>"twisted_beam_msh8-ew",
  "quantity"=> :Cauchy, "component"=> :xz)
  modeldata = AlgoDeforLinearModule.exportstresselementwise(modeldata)
  try rm(modeldata["postprocessing"]["exported"][1]["file"]); catch end

end
end
using .mmtwistedmsh8ort
mmtwistedmsh8ort.test()

module mmtwistedmsh9ort
using FinEtools
using FinEtools.AlgoDeforLinearModule
using Test
import LinearAlgebra: norm, cholesky, cross
import Statistics: mean
function test()

    # println("""
    # The initially twisted cantilever beam is one of the standard test
    # problems for verifying the finite-element accuracy [1]. The beam is
    #   clamped at one end and loaded either with unit in-plane or
    #   unit out-of-plane force at the other. The centroidal axis of the beam is
    #   straight at the undeformed  configuration, while its cross-sections are
    #   twisted about the centroidal axis from 0 at the clamped end to pi/2 at
    #   the free end.
    #
    #   Reference:
    #   Zupan D, Saje M (2004) On "A proposed standard set of problems to test
    #   finite element accuracy": the twisted beam. Finite Elements in Analysis
    #   and Design 40: 1445-1451.
    #   """)
    E1s =   E2s =  E3s = 0.29e8
    nu12s = nu13s = nu23s = 0.22
    G12s = G13s = G23s = E1s/2/(1+nu12s)
    # E = 0.29e8;
    # nu = 0.22;
    W = 1.1;
    L = 12.;
    t =  0.32;
    nl = 2; nt = 1; nw = 1; ref = 7;
    p =   1/W/t;
    #  Loading in the Z direction
    loadv = [0;0;p]; dir = 3; uex = 0.005424534868469; # Harder: 5.424e-3;
    #   Loading in the Y direction
    #loadv = [0;p;0]; dir = 2; uex = 0.001753248285256; # Harder: 1.754e-3;
    tolerance  = t/1000;

    fens,fes  = H8block(L,W,t, nl*ref,nw*ref,nt*ref)

    # Reshape into a twisted beam shape
    for i = 1:count(fens)
        a = fens.xyz[i,1]/L*(pi/2); y = fens.xyz[i,2]-(W/2); z = fens.xyz[i,3]-(t/2);
        fens.xyz[i,:] = [fens.xyz[i,1],y*cos(a)-z*sin(a),y*sin(a)+z*cos(a)];
    end

    # Clamped end of the beam
    l1  = selectnode(fens; box = [0 0 -100*W 100*W -100*W 100*W], inflate  =  tolerance)
    e1 = FDataDict("node_list"=>l1, "component"=>1, "displacement"=>0.0)
    e2 = FDataDict("node_list"=>l1, "component"=>2, "displacement"=>0.0)
    e3 = FDataDict("node_list"=>l1, "component"=>3, "displacement"=>0.0)

    # Traction on the opposite edge
    boundaryfes  =   meshboundary(fes);
    Toplist   = selectelem(fens,boundaryfes, box =  [L L -100*W 100*W -100*W 100*W], inflate =   tolerance);
    el1femm  = FEMMBase(IntegDomain(subset(boundaryfes,Toplist), GaussRule(2, 2)))
    flux1 = FDataDict("femm"=>el1femm, "traction_vector"=>loadv)


    # Make the region
    MR = DeforModelRed3D
    material = MatDeforElastOrtho(MR,
    0.0, E1s, E2s, E3s,
    nu12s, nu13s, nu23s,
    G12s, G13s, G23s,
    0.0, 0.0, 0.0)
    region1 = FDataDict("femm"=>FEMMDeforLinearMSH8(MR, IntegDomain(fes, GaussRule(3,2)),
    material))

    # Make model data
    modeldata =  FDataDict( "fens"=> fens, "regions"=>  [region1],
    "essential_bcs"=>[e1, e2, e3], "traction_bcs"=>  [flux1])

    # Call the solver
    modeldata = AlgoDeforLinearModule.linearstatics(modeldata)
    geom = modeldata["geom"]
    u = modeldata["u"]

    # Extract the solution
    nl = selectnode(fens, box = [L L -100*W 100*W -100*W 100*W],inflate = tolerance);
    theutip = mean(u.values[nl,:], dims = 1)
    # println("displacement  = $(theutip[dir]) as compared to converged $uex")
    @test abs(theutip[dir]-uex)/uex < 0.0012

    # Write out mesh with displacements
    modeldata["postprocessing"] = FDataDict("file"=>"twisted_beam_msh8")
    modeldata = AlgoDeforLinearModule.exportdeformation(modeldata)
    try rm(modeldata["postprocessing"]["exported"][1]["file"]); catch end

    # Write out mesh with stresses
    modeldata["postprocessing"] = FDataDict("file"=>"twisted_beam_msh8",
    "quantity"=> :Cauchy, "component"=> :xy, "outputcsys"=>CSys(3))
    modeldata = AlgoDeforLinearModule.exportstress(modeldata)
    try rm(modeldata["postprocessing"]["exported"][1]["file"]); catch end

    # Write out mesh with stresses
    modeldata["postprocessing"] = FDataDict("file"=>"twisted_beam_msh8",
    "quantity"=> :Cauchy, "component"=> :xz, "outputcsys"=>CSys(3))
    modeldata = AlgoDeforLinearModule.exportstress(modeldata)
    try rm(modeldata["postprocessing"]["exported"][1]["file"]); catch end

    # Write out mesh with von Mises stresses
    modeldata["postprocessing"] = FDataDict("file"=>"twisted_beam_msh8",
    "quantity"=> :vm, "outputcsys"=>CSys(3))
    modeldata = AlgoDeforLinearModule.exportstress(modeldata)
    vm  = modeldata["postprocessing"]["exported"][1]["field"]
    # println("extremes of vm, nodal: $([minimum(vm.values),   maximum(vm.values)])")
    try rm(modeldata["postprocessing"]["exported"][1]["file"]); catch end
    @test norm([minimum(vm.values),   maximum(vm.values)]-[4.78774, 522.126]) < 0.01

    # Write out mesh with von Mises stresses, elementwise
    modeldata["postprocessing"] = FDataDict("file"=>"twisted_beam_msh8-ew",
    "quantity"=> :vm, "outputcsys"=>CSys(3))
    modeldata = AlgoDeforLinearModule.exportstresselementwise(modeldata)
    vm  = modeldata["postprocessing"]["exported"][1]["field"]
    # println("extremes of vm, elemental: $([minimum(vm.values),   maximum(vm.values)])")
    try rm(modeldata["postprocessing"]["exported"][1]["file"]); catch end
    @test norm([minimum(vm.values),   maximum(vm.values)]-[1.85882, 522.126]) < 0.01

    # Write out mesh with von Mises stresses, elementwise
    modeldata["postprocessing"] = FDataDict("file"=>"twisted_beam_msh8-ew",
    "quantity"=> :Cauchy, "component"=> :xz, "outputcsys"=>CSys(3))
    modeldata = AlgoDeforLinearModule.exportstresselementwise(modeldata)
    try rm(modeldata["postprocessing"]["exported"][1]["file"]); catch end

end
end
using .mmtwistedmsh9ort
mmtwistedmsh9ort.test()

module mxRMSerror3a1
using FinEtools
using FinEtools.AlgoDeforLinearModule
using FinEtools.MeshUtilModule
using FinEtools.AlgoBaseModule
using Test
function test()


    elementtag = "MSH8"
    # println("""
    # Pagano_3layer_cylindrical_bending: $(elementtag)
    # """)

    # This example provides three-dimensional finite element model for  the
    # transverse shear stress calculations. The problem consists of a one-, two- or
    # three-layer plate subjected to a sinusoidal distributed load, as
    # described by Pagano (1969). The resulting transverse shear and axial
    # stresses through the thickness of the plate are compared to two existing
    # analytical solutions by Pagano (1969). The first solution is derived from
    # classical laminated plate theory (CPT), while the second is an exact
    # solution from linear elasticity theory.

    filebase = "Pagano_3layer_cylindrical_bending_$(elementtag)_convergence"

    modeldatasequence = FDataDict[]
    for Refinement = [1, 2, 4]

        # Orthotropic material for the 3 layers
        E1 = 25e6*phun("PSI"); E2 = 1e6*phun("PSI"); E3 = E2;
        G12 = 0.5e6*phun("PSI");  G13 = G12; G23 = 0.2e6*phun("PSI")
        nu12 =  0.25; nu13 =  0.25; nu23 =  0.25;
        Span_to_thickness = 4.0;
        T = 2.5*phun("in"); # total thickness of the plate
        L = Span_to_thickness*T;
        h = 1*phun("in");  # depth of the plate
        q0 = 1*phun("PSI")
        CTE1 =  CTE2 =  CTE3 = 0.0

        # Here we define the layout and the thicknesses of the layers.
        angles = vec([0.0 90.0 0.0]);
        nLayers = length(angles)
        ts = T/nLayers * ones(nLayers); # layer thicknesses

        tolerance = 0.0001*T

        # Select how find the mesh should be
        nL, nh = Refinement * 2 * 4, Refinement*1;
        nts= Refinement * 2 * ones(Int, nLayers);# number of elements per layer

        xs = collect(linearspace(0.0, L, nL+1))
        ys = collect(linearspace(0.0, h, nh+1))

        fens,fes = H8layeredplatex(xs, ys, ts, nts)
        # println("count(fens) = $(count(fens))")

        # This is the material  model
        MR = DeforModelRed3D
        skinmaterial = MatDeforElastOrtho(MR,
        0.0, E1, E2, E3,
        nu12, nu13, nu23,
        G12, G13, G23,
        CTE1, CTE2, CTE3)

        # The material coordinate system function is defined as:
        function updatecs!(csmatout::FFltMat, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)
            rotmat3!(csmatout, angles[fe_label]/180.0*pi* [0.0; 0.0; 1.0]);
        end

        # The vvolume integrals are evaluated using this rule
        gr = GaussRule(3, 2)

        # We will create 3 regions, one for each of the layers
        regions = FDataDict[]
        for layer = 1:nLayers
            rls = selectelem(fens, fes, label =  layer)
            push!(regions, FDataDict("femm"=>FEMMDeforLinearMSH8(MR,
            IntegDomain(subset(fes, rls), gr), CSys(3, 3, updatecs!), skinmaterial)))
        end

        # File =  "Meyer_Piening_sandwich-r1.vtk"
        # vtkexportmesh(File, skinregion["femm"].integdomain.fes.conn, fens.xyz, FinEtools.MeshExportModule.H8)
        # # @async run(`"paraview.exe" $File`)


        # The essential boundary conditions are applied to enforce the plane strain constraint.
        ly0 = selectnode(fens, box=[-Inf Inf 0.0 0.0 -Inf Inf], inflate=tolerance)
        lyh = selectnode(fens, box=[-Inf Inf h h -Inf Inf], inflate=tolerance)
        ey = FDataDict("displacement"=>  0.0, "component"=> 2, "node_list"=>vcat(ly0, lyh))
        # The transverse displacement is fixed at the two ends.
        lz0 = selectnode(fens, box=[0.0 0.0 -Inf Inf -Inf Inf], inflate=tolerance)
        lzL = selectnode(fens, box=[L L -Inf Inf -Inf Inf], inflate=tolerance)
        ez = FDataDict("displacement"=>  0.0, "component"=> 3, "node_list"=>vcat(lz0, lzL))
        ex = FDataDict("displacement"=>  0.0, "component"=> 1, "node_list"=>[1])

        # The traction boundary condition is applied at the top of the plate.
        bfes = meshboundary(fes)
        function pfun(forceout::FVec{T}, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt) where {T}
            forceout[1] = 0.0
            forceout[2] = 0.0
            forceout[3] = -q0*sin(pi*XYZ[1]/L)
            return forceout
        end
        # From  the entire boundary we select those quadrilaterals that lie on the plane
        # Z = thickness
        tl = selectelem(fens, bfes, box = [-Inf Inf -Inf Inf T T], inflate=tolerance)
        Trac = FDataDict("traction_vector"=>pfun,
        "femm"=>FEMMBase(IntegDomain(subset(bfes, tl), GaussRule(2, 2))))

        modeldata = FDataDict("fens"=>fens,
        "regions"=>regions,
        "essential_bcs"=>[ex, ey, ez],
        "traction_bcs"=> [Trac]
        )
        modeldata = AlgoDeforLinearModule.linearstatics(modeldata)

        modeldata["postprocessing"] = FDataDict("file"=>filebase * "-u")
        modeldata = AlgoDeforLinearModule.exportdeformation(modeldata)
        for e in modeldata["postprocessing"]["exported"]
            try rm(e["file"]) catch end
        end

        u = modeldata["u"]
        geom = modeldata["geom"]

        # The results of the displacement and stresses will be reported at
        # nodes located at the appropriate points.
        ntopcenter = selectnode(fens, box=[L/2 L/2 0.0 h T T], inflate=tolerance)
        ncenterline = selectnode(fens, box=[L/2 L/2 0.0 0.0 0.0 T], inflate=tolerance)
        nx0line = selectnode(fens, box=[0.0 0.0 0.0 0.0 0.0 T], inflate=tolerance)

        zclo = sortperm(vec(geom.values[ncenterline, 3]))
        ncenterline = ncenterline[zclo]
        centerz = geom.values[ncenterline, 3]

        # println("Top Center deflection: $(mean(u.values[ntopcenter, 3], 1)/phun("in")) [in]")

        # # extrap = :extrapmean
        extrap = :extraptrend
        nodevalmeth = :averaging
        # extrap = :default
        # nodevalmeth = :invdistance

        # Compute  all stresses
        modeldata["postprocessing"] = FDataDict("file"=>filebase * "-s",
        "quantity"=>:Cauchy, "component"=>collect(1:6), "outputcsys"=>CSys(3),
        "nodevalmethod"=>nodevalmeth, "reportat"=>extrap)
        modeldata = AlgoDeforLinearModule.exportstress(modeldata)
        for e in modeldata["postprocessing"]["exported"]
            try rm(e["file"]) catch end
        end

        modeldata["elementsize"] = 1.0/Refinement
        modeldata["geometricaltolerance"] = tolerance
        modeldata["targetfields"] = [e["field"] for e in modeldata["postprocessing"]["exported"]]
        push!(modeldatasequence, modeldata)
    end # for refinement

    elementsizes, errornorms, p = AlgoBaseModule.evalconvergencestudy(modeldatasequence)

    # println("")
    # println("Normalized Approximate Error = $(errornorms)")
    @test abs(p[1] - 1.3347513854727369)/1.3347513854727369 < 1.0e-3
    # csvFile = filebase * "_errors" * ".CSV"
    # savecsv(csvFile,
    # elementsizes=vec(elementsizes[1:end-1]),
    # elementsizes2=vec(elementsizes[1:end-1].^2),
    # elementsizes3=vec(elementsizes[1:end-1].^3),
    # errornorms=vec(errornorms)
    # )

    # @async run(`"paraview.exe" $csvFile`)

    # println("Done")

end
end
using .mxRMSerror3a1
mxRMSerror3a1.test()

module munit_cube_modes_nice_t4
using FinEtools
using Test
using Arpack
using LinearAlgebra
function test()
    # println("""
    # Vibration modes of unit cube  of almost incompressible material.
    # %
    # Reference: Puso MA, Solberg J (2006) A stabilized nodally integrated
    # tetrahedral. International Journal for Numerical Methods in
    # Engineering 67: 841-867.
    # """)
    t0 = time()

    E = 1*phun("PA");
    nu = 0.499;
    rho = 1*phun("KG/M^3");
    a = 1*phun("M"); b = a; h =  a;
    n1 = 10;# How many element edges per side?
    na =  n1; nb =  n1; nh  = n1;
    neigvs = 20                   # how many eigenvalues
    OmegaShift = (0.01*2*pi)^2;
    stabfact = 0.015
    Eigenvalues = [0.0, 5.93656e-8, 7.54751e-8, 9.80131e-8, 1.14899e-7, 1.27725e-7, 0.264544, 0.266128, 0.350568, 0.352546, 0.355279, 0.357389, 0.357701, 0.359704, 0.402389, 0.402968, 0.404977, 0.45061, 0.450974, 0.452039]

    MR = DeforModelRed3D
    fens,fes  = T4block(a,b,h, na,nb,nh)

    geom = NodalField(fens.xyz)
    u = NodalField(zeros(size(fens.xyz,1),3)) # displacement field

    numberdofs!(u)

    material = MatDeforElastIso(MR, rho, E, nu, 0.0)

    femm = FEMMDeforLinearNICET4(MR, IntegDomain(fes, NodalSimplexRule(3)), material, stabfact)
    associategeometry!(femm,  geom)
    K  = stiffness(femm, geom, u)
    M = mass(femm, geom, u)
    d,v,nev,nconv = eigs(K+OmegaShift*M, M; nev=neigvs, which=:SM)
    d = d .- OmegaShift;
    fs = real(sqrt.(complex(d)))/(2*pi)
    # # println("Eigenvalues: $fs [Hz]")
    @test norm(vec(fs) .- vec(Eigenvalues)) < 1.0e-4*maximum(vec(Eigenvalues))

    # mode = 17
    # scattersysvec!(u, v[:,mode])
    # File =  "unit_cube_modes.vtk"
    # vtkexportmesh(File, fens, fes; vectors=[("mode$mode", u.values)])

    true

end
end
using .munit_cube_modes_nice_t4
munit_cube_modes_nice_t4.test()

module malum_cyl_mode_nice_t4
using FinEtools
using Test
using Arpack
using LinearAlgebra
function test()
    # Aluminum cylinder free vibration, mesh imported from Abaqus
    # Mesh converted from quadratic tetrahedra to linear tetrahedra
    # NICE tetrahedral elements used
    E = 70000*phun("MPa");
    nu = 0.33;
    rho = 2700*phun("KG/M^3");
    radius = 0.5*phun("ft");
    neigvs = 20                   # how many eigenvalues
    OmegaShift = (10.0*2*pi)^2;
    stabfact = 0.005
    Eigenvalues = [4.54746e-5, 6.82231e-5, 8.7071e-5, 9.99708e-5, 0.000112778, 0.000116397, 2533.6, 2535.12, 2574.64, 4086.61, 4652.66, 4654.16, 5122.94, 6755.62, 6756.45, 6872.26, 6875.3, 6883.49, 6888.53, 6983.99]

    MR = DeforModelRed3D
    output = import_ABAQUS("alum_cyl.inp")
    fens, fes = output["fens"], output["fesets"][1]
    fens.xyz .*= phun("mm") # The input is provided in SI(mm) units
    fens, fes = T10toT4(fens, fes)

    geom = NodalField(fens.xyz)
    u = NodalField(zeros(size(fens.xyz,1),3)) # displacement field

    numberdofs!(u)

    material = MatDeforElastIso(MR, rho, E, nu, 0.0)

    femm = FEMMDeforLinearNICET4(MR, IntegDomain(fes, NodalSimplexRule(3)), material, stabfact)
    associategeometry!(femm,  geom)
    K  = stiffness(femm, geom, u)
    M = mass(femm, geom, u)
    d,v,nev,nconv = eigs(K+OmegaShift*M, M; nev=neigvs, which=:SM)
    d = d .- OmegaShift;
    fs = real(sqrt.(complex(d)))/(2*pi)
    # # println("Eigenvalues: $fs [Hz]")
    @test norm(vec(fs) .- vec(Eigenvalues)) < 1.0e-3*maximum(vec(Eigenvalues))

    true
end
end
using .malum_cyl_mode_nice_t4
malum_cyl_mode_nice_t4.test()

module malum_cyl_mode_esnice_t4
using FinEtools
using Test
using Arpack
using LinearAlgebra
function test()
    # Aluminum cylinder free vibration, mesh imported from Abaqus
    # Mesh converted from quadratic tetrahedra to linear tetrahedra
    # NICE tetrahedral elements used
    E = 70000*phun("MPa");
    nu = 0.33;
    rho = 2700*phun("KG/M^3");
    radius = 0.5*phun("ft");
    neigvs = 20                   # how many eigenvalues
    OmegaShift = (10.0*2*pi)^2;
    stabfact = 0.005
    Eigenvalues =   [0.0, 0.0, 0.0, 1.8846e-5, 7.35917e-5, 0.000119445, 2498.15, 2498.88, 2513.31, 4082.65, 4585.99, 4586.42, 4987.01, 6648.02, 6648.48, 6679.04, 6682.16, 6777.89, 6780.59, 6799.36]
    MR = DeforModelRed3D
    output = import_ABAQUS("alum_cyl.inp")
    fens, fes = output["fens"], output["fesets"][1]
    fens.xyz .*= phun("mm") # The input is provided in SI(mm) units
    fens, fes = T10toT4(fens, fes)

    geom = NodalField(fens.xyz)
    u = NodalField(zeros(size(fens.xyz,1),3)) # displacement field

    numberdofs!(u)

    material = MatDeforElastIso(MR, rho, E, nu, 0.0)

    femm = FEMMDeforLinearESNICET4(MR, IntegDomain(fes, NodalSimplexRule(3)), material)
    associategeometry!(femm,  geom)
    K  = stiffness(femm, geom, u)
    M = mass(femm, geom, u)
    d,v,nev,nconv = eigs(K+OmegaShift*M, M; nev=neigvs, which=:SM)
    d = d .- OmegaShift;
    fs = real(sqrt.(complex(d)))/(2*pi)
    # println("Eigenvalues: $fs [Hz]")
    @test norm(vec(fs) .- vec(Eigenvalues)) < 1.0e-3*maximum(vec(Eigenvalues))

    true
end
end
using .malum_cyl_mode_esnice_t4
malum_cyl_mode_esnice_t4.test()


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

    femm = FEMMDeforLinear(MR, IntegDomain(fes, GaussRule(2, 2), true), material)

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

    femm = FEMMDeforLinear(MR, IntegDomain(fes, GaussRule(2, 2), true), material)

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

    femm = FEMMDeforLinear(MR, IntegDomain(fes, GaussRule(2, 2), true), material)

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


module moblocknzebc13
using FinEtools
using LinearAlgebra: norm
using Test
function test()
    E = 1*phun("PA");
    nu = 0.499;
    rho = 1*phun("KG/M^3");
    a = 1*phun("M"); b = a; h =  a;
    n1 = 10;# How many element edges per side?
    na =  n1; nb =  n1; nh  = n1;
    mag = 0.001

    MR = DeforModelRed3D
    fens,fes  = H8block(a,b,h, na,nb,nh)

    bdryfes = meshboundary(fes);

    # now we create the geometry and displacement fields
    geom = NodalField(fens.xyz)
    u = NodalField(zeros(size(fens.xyz,1),3)) # displacement field

    # the constrained face
    l1 =selectnode(fens; plane = [1.0 0.0 0.0 0.0], thickness = a/1000)
    # This face rotates about the Y and Z axes
    setebc!(u,l1,true, 1, [mag*fens.xyz[i, 3]-mag*fens.xyz[i, 2] for i in l1])
    setebc!(u,l1,true, 2, 0.0)
    setebc!(u,l1,true, 3, 0.0)

    applyebc!(u)
    numberdofs!(u)

    # Property and material
    material=MatDeforElastIso(MR, E, nu)

    femm = FEMMDeforLinearMSH8(MR, IntegDomain(fes, GaussRule(3, 2)), material)
    associategeometry!(femm, geom)
    K =stiffness(femm, geom, u)
    F = nzebcloadsstiffness(femm, geom, u)
    U=  K\(F)
    scattersysvec!(u,U[:])
    @test abs(maximum(U[:])-0.001) < 1.0e-6

    fld= fieldfromintegpoints(femm, geom, u, :princCauchy, 1)
    @test norm(fld.values) < 1.0e-9 # Rigid body rotation: no stresses

    # File =  "moblocknzebc13.vtk"
    # vtkexportmesh(File, fens, fes; scalars=[("sigmaz", fld.values)],
    #               vectors=[("u", u.values)])
    # @async run(`"paraview.exe" $File`)
end
end
using .moblocknzebc13
moblocknzebc13.test()

module mmoblocknzebc13a
using FinEtools
using LinearAlgebra: norm
using Test
function test()
    E = 1*phun("PA");
    nu = 0.499;
    rho = 1*phun("KG/M^3");
    a = 1*phun("M"); b = a; h =  a;
    n1 = 10;# How many element edges per side?
    na =  n1; nb =  n1; nh  = n1;
    mag = 0.001

    MR = DeforModelRed3D
    fens,fes  = T10block(a,b,h, na,nb,nh)

    bdryfes = meshboundary(fes);

    # now we create the geometry and displacement fields
    geom = NodalField(fens.xyz)
    u = NodalField(zeros(size(fens.xyz,1),3)) # displacement field

    # the constrained face
    l1 =selectnode(fens; plane = [1.0 0.0 0.0 0.0], thickness = a/1000)
    # This face rotates about the Y and Z axes
    setebc!(u,l1,true, 1, [mag*fens.xyz[i, 3]-mag*fens.xyz[i, 2] for i in l1])
    setebc!(u,l1,true, 2, 0.0)
    setebc!(u,l1,true, 3, 0.0)

    applyebc!(u)
    numberdofs!(u)

    # Property and material
    material=MatDeforElastIso(MR, E, nu)

    femm = FEMMDeforLinearMST10(MR, IntegDomain(fes, TetRule(4)), material)
    associategeometry!(femm, geom)
    K =stiffness(femm, geom, u)
    F = nzebcloadsstiffness(femm, geom, u)
    U=  K\(F)
    scattersysvec!(u,U[:])
    @test abs(maximum(U[:])-0.001) < 1.0e-6

    fld= fieldfromintegpoints(femm, geom, u, :princCauchy, 1)
    @test norm(fld.values) < 1.0e-9 # Rigid body rotation: no stresses

    File =  "mmoblocknzebc13a.vtk"
    vtkexportmesh(File, fens, fes; scalars=[("sigmaz", fld.values)],
                  vectors=[("u", u.values)])
    # @async run(`"paraview.exe" $File`)
    try rm(File) catch end
end
end
using .mmoblocknzebc13a
mmoblocknzebc13a.test()

module trunc_cyl_shell_esnice_t4
using FinEtools
using Test
using Arpack
using LinearAlgebra
function test()
    # println("""
    # Vibration modes of truncated cylindrical shell. NASTRAN input file.
    # """)

    # t0 = time()

    E = 205000*phun("MPa");# Young's modulus
    nu = 0.3;# Poisson ratio
    rho = 7850*phun("KG*M^-3");# mass density
    OmegaShift = (2*pi*100) ^ 2; # to resolve rigid body modes
    h = 0.05*phun("M");
    l = 10*h;
    Rmed = h/0.2;
    psi   = 0;    # Cylinder
    nh = 5; nl  = 12; nc = 40;
    tolerance = h/nh/100;
    neigvs = 20;

    MR = DeforModelRed3D
    output = import_NASTRAN("trunc_cyl_shell_0.nas")
    fens, fes = output["fens"], output["fesets"][1]

    geom = NodalField(fens.xyz)
    u = NodalField(zeros(size(fens.xyz,1),3)) # displacement field

    numberdofs!(u)

    material = MatDeforElastIso(MR, rho, E, nu, 0.0)

    femm = FEMMDeforLinearESNICET4(MR, IntegDomain(fes, NodalSimplexRule(3)), material)
    femm = associategeometry!(femm, geom)
    K = stiffness(femm, geom, u)
    M = mass(femm, geom, u)


    # eigs returns the nev requested eigenvalues in d, the corresponding Ritz vectors
    # v (only if ritzvec=true), the number of converged eigenvalues nconv, the number
    # of iterations niter and the number of matrix vector multiplications nmult, as
    # well as the final residual vector resid.

    if true
        d,v,nev,nconv = eigs(K+OmegaShift*M, M; nev=neigvs, which=:SM)
        d = d .- OmegaShift;
        fs = real(sqrt.(complex(d)))/(2*pi)
        # println("Eigenvalues: $fs [Hz]")
        @test norm(fs .- [0.0, 0.0, 0.0, 1.29647e-5, 4.07636e-5, 8.90969e-5, 517.081, 557.081, 718.095, 749.177, 1423.0, 1423.01, 1690.21, 1690.23, 2451.16, 2474.49, 2475.38, 2478.54, 2479.72, 2500.88]) < 0.001 * norm(fs)
        # mode = 7
        # scattersysvec!(u, v[:,mode])
        # File =  "trunc_cyl_shell_nas.vtk"
        # vtkexportmesh(File, fens, fes; vectors=[("mode$mode", u.values)])
        # @async run(`"paraview.exe" $File`)
    end

    true
end
end
using .trunc_cyl_shell_esnice_t4
trunc_cyl_shell_esnice_t4.test()


module minfsuptest1
using FinEtools
using FinEtools.FEMMDeforLinearBaseModule: infsup_gh, infsup_sh
using Test
import LinearAlgebra: norm, cholesky, I, eigen
# using UnicodePlots
function test()
	lambdatol = sqrt(1e4*eps(1.0));
	E=1000.0;
	nu=0.24;
	parshiftmult= 0.002;
	A = [1.44 -0.741 -0.53; -0.626 1.589 -0.913; -0.55 0.43 1.756] + 1.0I;

	lambdamin = Float64[]
	h = Float64[]
	for ne = [2, 3, 4, 5, 6]
		Length::FFlt, Width::FFlt, Height::FFlt, nL::FInt, nW::FInt, nH::FInt, orientation::Symbol = (6.0, 6.0, 6.0, ne, ne, ne, :a)

		fens, fes = T4block(Length::FFlt, Width::FFlt, Height::FFlt, nL::FInt, nW::FInt, nH::FInt, orientation::Symbol)
		# fens, fes = T4toT10(fens, fes)
		# if parshiftmult != 0
		# 	ax=parshiftmult*Length;
		# 	ay=parshiftmult*Length;
		# 	az=parshiftmult*Length;
		# 	fens = transform_apply(fens,@(x, data) (x +[-ax/2*x(3)^2*x(2),ay/2*x(3)*x(1)^2,az/4*x(1)^3*x(2)+ax/4*x(3)^2*sin(0.2*x(2)^3)]), []);
		# end

		for i = 1:count(fens)
			fens.xyz[i,:] += A*fens.xyz[i,:];
		end

		# File =  "minfsuptest1.vtk"
		# vtkexportmesh(File, fens, fes)
		# try rm(File); catch end

		MR  =  DeforModelRed3D

		material = MatDeforElastIso(MR, E, nu)

		femm  =  FEMMDeforLinear(MR, IntegDomain(fes, TetRule(4)), material)

		geom = NodalField(fens.xyz)
		u = NodalField(zeros(size(fens.xyz,1), 3)) # displacement field
		bfes = meshboundary(fes)
		l1 = connectednodes(bfes)
		setebc!(u, l1, true, 1, 0.0)
		setebc!(u, l1, true, 2, 0.0)
		setebc!(u, l1, true, 3, 0.0)
		numberdofs!(u)

		Gh = infsup_gh(femm, geom, u);
		Sh = infsup_sh(femm, geom, u);

		lambda, modes = eigen(Matrix(Gh), Matrix(Sh));

		# @show lambda
		abslambda = deepcopy(lambda);
		ix = findall(y  -> y < 0.0, lambda);
		if !isempty(ix)
			abslambda[ix] .= 0;
		end

		abslambda = sqrt.(sort(abslambda));

		ix = findall(y  -> y >= lambdatol, abslambda);
		if isempty(ix)
			@error "Bad guess of the number of eigenvalues"
		end
		push!(lambdamin, abslambda[ix[1]])
		push!(h, Length / ne)
	end

# @show lambdamin
# 	a = lineplot(log.(h), log.(lambdamin), name = "infsup", xlabel = "log(Element Size)", ylabel = "log(minimum eigenvalue)", canvas = DotCanvas)
# 	display(a)
	
	@test norm(lambdamin - [0.262065, 0.1709, 0.126159, 0.100228, 0.0828139]) / norm(lambdamin) <= 1.0e-4

end
end
using .minfsuptest1
minfsuptest1.test()

module minfsuptestik1
using FinEtools
using FinEtools.FEMMDeforLinearBaseModule: infsup_gh, infsup_sh
using Test
import LinearAlgebra: norm, cholesky, I, eigen
# using UnicodePlots
function test()
	lambdatol = sqrt(1e8*eps(1.0));
	E=1000.0;
	nu=0.24;
	parshiftmult= 0.002;
	A = [1.44 -0.741 -0.53; -0.626 1.589 -0.913; -0.55 0.43 1.756] + 1.0I;

	lambdamin = Float64[]
	h = Float64[]
	for ne = [2, 3, 4]
		Length::FFlt, Width::FFlt, Height::FFlt, nL::FInt, nW::FInt, nH::FInt, orientation::Symbol = (6.0, 6.0, 6.0, ne, ne, ne, :a)

		fens, fes = T4block(Length::FFlt, Width::FFlt, Height::FFlt, nL::FInt, nW::FInt, nH::FInt, orientation::Symbol)
		fens, fes = T4toT10(fens, fes)
		# @show connasarray(fes)

		for i = 1:count(fens)
			fens.xyz[i,:] = fens.xyz[i,:] + vec(reshape(fens.xyz[i,:], 1, 3)*A);
		end
		# @show fens.xyz
	
		# File =  "minfsuptest1.vtk"
		# vtkexportmesh(File, fens, fes)
		# try rm(File); catch end

		MR  =  DeforModelRed3D

		material = MatDeforElastIso(MR, E, nu)

		
		geom = NodalField(fens.xyz)
		u = NodalField(zeros(size(fens.xyz,1), 3)) # displacement field
		bfes = meshboundary(fes)
		l1 = connectednodes(bfes)
		setebc!(u, l1, true, 1, 0.0)
		setebc!(u, l1, true, 2, 0.0)
		setebc!(u, l1, true, 3, 0.0)
		numberdofs!(u)

		femm  =  FEMMDeforLinear(MR, IntegDomain(fes, TetRule(1)), material)
		Gh = infsup_gh(femm, geom, u);
		femm  =  FEMMDeforLinear(MR, IntegDomain(fes, TetRule(4)), material)
		Sh = infsup_sh(femm, geom, u);
		
		lambda, modes = eigen(Matrix(Gh), Matrix(Sh));
		
		# @show lambda
		abslambda = real.(filter(y -> !isnan(y), lambda));
		ix = findall(y  -> y < 0.0, abslambda);
		if !isempty(ix)
			abslambda[ix] .= 0;
		end

		abslambda = sqrt.(sort(abslambda));
		ix = findall(y  -> y > 0.0, abslambda);
		# a = lineplot(1:length(abslambda[ix]), log.(abslambda[ix]), name = "infsup", xlabel = "eigenvalue", ylabel = "log(eigenvalue)", canvas = DotCanvas)
		# display(a)

		ix = findall(y  -> y >= lambdatol, abslambda);
		if isempty(ix)
			@error "Bad guess of the number of eigenvalues"
		end
		push!(lambdamin, abslambda[ix[1]])
		push!(h, 1.0/(count(fens))^(1/3))
	end

	# @show lambdamin
	# a = lineplot(log.(h), log.(lambdamin), name = "infsup", xlabel = "log(Element Size)", ylabel = "log(minimum eigenvalue)", canvas = DotCanvas)
	# display(a)
	
	@test norm(lambdamin - [0.0729658, 0.0585958, 0.0459494] ) / norm(lambdamin) <= 1.0e-4

end
end
using .minfsuptestik1
minfsuptestik1.test()


module minfsuptestik2
using FinEtools
using FinEtools.FEMMDeforLinearBaseModule: infsup_gh, infsup_sh
using Test
import LinearAlgebra: norm, cholesky, I, eigen
# using UnicodePlots
function test()
	lambdatol = sqrt(1e8*eps(1.0));
	E=1000.0;
	nu=0.24;
	parshiftmult= 0.002;
	A = [1.44 -0.741 -0.53; -0.626 1.589 -0.913; -0.55 0.43 1.756] + 1.0I;

	lambdamin = Float64[]
	h = Float64[]
	for ne = [2, 3, 4]
		Length::FFlt, Width::FFlt, Height::FFlt, nL::FInt, nW::FInt, nH::FInt, orientation::Symbol = (6.0, 6.0, 6.0, ne, ne, ne, :a)

		fens, fes = T4block(Length::FFlt, Width::FFlt, Height::FFlt, nL::FInt, nW::FInt, nH::FInt, orientation::Symbol)
		# fens, fes = T4toT10(fens, fes)
		# @show connasarray(fes)

		for i = 1:count(fens)
			fens.xyz[i,:] = fens.xyz[i,:] + vec(reshape(fens.xyz[i,:], 1, 3)*A);
		end
		# @show fens.xyz
	
		# File =  "minfsuptest1.vtk"
		# vtkexportmesh(File, fens, fes)
		# try rm(File); catch end

		MR  =  DeforModelRed3D

		material = MatDeforElastIso(MR, E, nu)

		
		geom = NodalField(fens.xyz)
		u = NodalField(zeros(size(fens.xyz,1), 3)) # displacement field
		bfes = meshboundary(fes)
		l1 = connectednodes(bfes)
		setebc!(u, l1, true, 1, 0.0)
		setebc!(u, l1, true, 2, 0.0)
		setebc!(u, l1, true, 3, 0.0)
		numberdofs!(u)

		femm  =  FEMMDeforLinear(MR, IntegDomain(fes, TetRule(1)), material)
		Gh = infsup_gh(femm, geom, u);
		femm  =  FEMMDeforLinear(MR, IntegDomain(fes, TetRule(1)), material)
		Sh = infsup_sh(femm, geom, u);

		lambda, modes = eigen(Matrix(Gh), Matrix(Sh));

		# @show lambda
		abslambda = real.(filter(y -> !isnan(y), lambda));
		ix = findall(y  -> y < 0.0, abslambda);
		if !isempty(ix)
			abslambda[ix] .= 0;
		end

		abslambda = sqrt.(sort(abslambda));
		ix = findall(y  -> y > 0.0, abslambda);
		# a = lineplot(1:length(abslambda[ix]), log.(abslambda[ix]), name = "infsup", xlabel = "eigenvalue", ylabel = "log(eigenvalue)", canvas = DotCanvas)
		# display(a)

		ix = findall(y  -> y >= lambdatol, abslambda);
		if isempty(ix)
			@error "Bad guess of the number of eigenvalues"
		end
		push!(lambdamin, abslambda[ix[1]])
		push!(h, 1.0/(count(fens))^(1/3))
	end

	# @show lambdamin
	# a = lineplot(log.(h), log.(lambdamin), name = "infsup", xlabel = "log(Element Size)", ylabel = "log(minimum eigenvalue)", canvas = DotCanvas)
	# display(a)
	
	@test norm(lambdamin - [0.270777, 0.179116, 0.132893]) / norm(lambdamin) <= 1.0e-4

end
end
using .minfsuptestik2
minfsuptestik2.test()

module minfsuptestik3
using FinEtools
using FinEtools.FEMMDeforLinearBaseModule: infsup_gh, infsup_sh
using Test
import LinearAlgebra: norm, cholesky, I, eigen
# using UnicodePlots
function test()
	lambdatol = sqrt(1e8*eps(1.0));
	E=1000.0;
	nu=0.24;
	parshiftmult= 0.002;
	A = [1.44 -0.741 -0.53; -0.626 1.589 -0.913; -0.55 0.43 1.756] + 1.0I;

	lambdamin = Float64[]
	h = Float64[]
	for ne = [2, 3, 4, 5, 6, 7, 8, 9]
		Length::FFlt, Width::FFlt, Height::FFlt, nL::FInt, nW::FInt, nH::FInt, orientation::Symbol = (6.0, 6.0, 6.0, ne, ne, ne, :a)

		fens, fes = H8block(Length::FFlt, Width::FFlt, Height::FFlt, nL::FInt, nW::FInt, nH::FInt)
		# @show connasarray(fes)

		for i = 1:count(fens)
			fens.xyz[i,:] = fens.xyz[i,:] + vec(reshape(fens.xyz[i,:], 1, 3)*A);
		end
		# @show fens.xyz
	
		# File =  "minfsuptest1.vtk"
		# vtkexportmesh(File, fens, fes)
		# try rm(File); catch end

		MR  =  DeforModelRed3D

		material = MatDeforElastIso(MR, E, nu)

		
		geom = NodalField(fens.xyz)
		u = NodalField(zeros(size(fens.xyz,1), 3)) # displacement field
		bfes = meshboundary(fes)
		l1 = connectednodes(bfes)
		setebc!(u, l1, true, 1, 0.0)
		setebc!(u, l1, true, 2, 0.0)
		setebc!(u, l1, true, 3, 0.0)
		numberdofs!(u)

		femm  =  FEMMDeforLinear(MR, IntegDomain(fes, GaussRule(3, 1)), material)
		Gh = infsup_gh(femm, geom, u);
		femm  =  FEMMDeforLinear(MR, IntegDomain(fes, GaussRule(3, 1)), material)
		Sh = infsup_sh(femm, geom, u);
		
		lambda, modes = eigen(Matrix(Gh), Matrix(Sh));
		
		# @show lambda
		abslambda = real.(filter(y -> !isnan(y), lambda));
		ix = findall(y  -> y < 0.0, abslambda);
		if !isempty(ix)
			abslambda[ix] .= 0;
		end

		abslambda = sqrt.(sort(abslambda));
		ix = findall(y  -> y > 0.0, abslambda);
		# a = lineplot(1:length(abslambda[ix]), log.(abslambda[ix]), name = "infsup", xlabel = "eigenvalue", ylabel = "log(eigenvalue)", canvas = DotCanvas)
		# display(a)

		ix = findall(y  -> y >= lambdatol, abslambda);
		if isempty(ix)
			@error "Bad guess of the number of eigenvalues"
		end
		push!(lambdamin, abslambda[ix[1]])
		push!(h, 1.0/(count(fens))^(1/3))
	end

	# @show lambdamin
	# a = lineplot(log.(h), log.(lambdamin), name = "infsup", xlabel = "log(Element Size)", ylabel = "log(minimum eigenvalue)", canvas = DotCanvas)
	# display(a)
	
	@test norm(lambdamin - [0.447936, 0.317169, 0.305056, 0.300754, 0.291477, 0.290534, 
0.285866, 0.285624]) / norm(lambdamin) <= 1.0e-4

end
end
using .minfsuptestik3
minfsuptestik3.test()


module mmLE1NAFEMSstressESNICET4
using FinEtools
using Test
import LinearAlgebra: norm, cholesky, cross
function test()
    E = 210e3*phun("MEGA*PA");# 210e3 MPa
    nu = 0.3;
    p = 10*phun("MEGA*PA");# 10 MPA Outward pressure on the outside ellipse
    sigma_yD = 92.7*phun("MEGA*PA");# tensile stress at [2.0, 0.0] meters
    Radius = 1.0*phun("m")
    Thickness = 0.1*phun("m")
    n = 8; # number of elements per side
    tolerance = 1.0/n/1000.; # Geometrical tolerance

    fens,fes = T4block(1.0, pi/2, Thickness, n, n*2, 3)

    bdryfes = meshboundary(fes);
    icl = selectelem(fens, bdryfes, box=[1.0, 1.0, 0.0, pi/2, 0.0, Thickness], inflate=tolerance);
    for i=1:count(fens)
        t=fens.xyz[i,1]; a=fens.xyz[i,2]; z=fens.xyz[i,3]
        fens.xyz[i,:]=[(t*3.25+(1-t)*2)*cos(a), (t*2.75+(1-t)*1)*sin(a), z];
    end

    
    geom = NodalField(fens.xyz)
    u = NodalField(zeros(size(fens.xyz,1),3)) # displacement field

    l1 =selectnode(fens; box=[0.0, Inf, 0.0, 0.0, 0.0, Thickness], inflate = tolerance)
    setebc!(u,l1,true, 2, 0.0)
    l1 =selectnode(fens; box=[0.0, 0.0, 0.0, Inf, 0.0, Thickness], inflate = tolerance)
    setebc!(u,l1,true, 1, 0.0)
    l1 =selectnode(fens; box=[0.0, Inf, 0.0, Inf, 0.0, 0.0], inflate = tolerance)
    setebc!(u,l1,true, 3, 0.0)

    applyebc!(u)
    numberdofs!(u)

    el1femm =  FEMMBase(IntegDomain(subset(bdryfes,icl), SimplexRule(2, 1)))
    function pfun(forceout::FVec{T}, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt) where {T}
        pt= [2.75/3.25*XYZ[1], 3.25/2.75*XYZ[2], 0.0]
        forceout .=    vec(p*pt/norm(pt));
        return forceout
    end
    fi = ForceIntensity(FFlt, 3, pfun);
    F2 = distribloads(el1femm, geom, u, fi, 2);

    # Note that the material object needs to be created with the proper
    # model-dimension reduction in mind.  In this case that is the fully three-dimensional solid.
    MR = DeforModelRed3D

    material = MatDeforElastIso(MR, E, nu)

    femm = FEMMDeforLinearESNICET4(MR, IntegDomain(fes, NodalSimplexRule(3)), material)

    # The geometry field now needs to be associated with the FEMM
    femm = associategeometry!(femm, geom)

    K = stiffness(femm, geom, u)
    K = cholesky(K)
    U = K\(F2)
    scattersysvec!(u, U[:])

    # File =  "a.vtk"
    # vtkexportmesh(File, fes.conn, fens.xyz, FinEtools.MeshExportModule.T4; vectors = [("u", u.values)])
    # @async run(`"paraview.exe" $File`)

    nl = selectnode(fens, box=[2.0, 2.0, 0.0, 0.0, 0.0, 0.0], inflate=tolerance);
    thecorneru = zeros(FFlt,1,3)
    gathervalues_asmat!(u, thecorneru, nl);
    thecorneru = thecorneru/phun("mm")
    # println("displacement =$(thecorneru) [MM] as compared to reference [-0.10215,0] [MM]")
    @test norm(thecorneru - [-0.115443 0.0 0.0]) < 1.0e-4


    fld = fieldfromintegpoints(femm, geom, u, :Cauchy, 2; nodevalmethod = :averaging)
    # println("Sigma_y =$(fld.values[nl,1][1]/phun("MPa")) as compared to reference sigma_yD = $(sigma_yD/phun("MPa")) [MPa]")
    # File =  "a.vtk"
    # vtkexportmesh(File, fes.conn, fens.xyz, FinEtools.MeshExportModule.T4; vectors = [("u", u.values)], scalars = [("sig", fld.values)])
    # @async run(`"paraview.exe" $File`)
    @test abs(fld.values[nl,1][1]/phun("MPa") - 91.07617647479451) < 1.0e-3

end
end
using .mmLE1NAFEMSstressESNICET4
mmLE1NAFEMSstressESNICET4.test()

module mmcubevibrationESNICEH8b
using FinEtools
using Test
using Arpack: eigs
import LinearAlgebra: norm, cholesky, cross
function test()
	# println("""
	#         % Vibration modes of unit cube  of almost incompressible material.
	#         %
	#         % Reference: Puso MA, Solberg J (2006) A stabilized nodally integrated
	#         % tetrahedral. International Journal for Numerical Methods in
	#         % Engineering 67: 841-867.""")
	#         t0 = time()


	E = 1*phun("PA");
	nu = 0.499;
	rho = 1*phun("KG/M^3");
	a = 1*phun("M"); b = a; h =  a;
	n1 = 5;# How many element edges per side?
	na =  n1; nb =  n1; nh  = n1;
	neigvs = 20                   # how many eigenvalues
	OmegaShift = (0.01*2*pi)^2;

	MR = DeforModelRed3D
	fens,fes  = T4block(a,b,h, na,nb,nh)
	fens,fes  = T4toH8(fens,fes)

	geom = NodalField(fens.xyz)
	u = NodalField(zeros(size(fens.xyz,1),3)) # displacement field

	numberdofs!(u)

	material=MatDeforElastIso(MR, rho, E, nu, 0.0)

	femm = FEMMDeforLinearESNICEH8(MR, IntegDomain(fes, NodalTensorProductRule(3)), material)
	femm = associategeometry!(femm, geom)

	K =stiffness(femm, geom, u)
	M =mass(femm, geom, u)
	d,v,nev,nconv = eigs(K+OmegaShift*M, M; nev=neigvs, which=:SM)
	d = d .- OmegaShift;
	fs = real(sqrt.(complex(d)))/(2*pi)
	# println("Eigenvalues: $fs [Hz]")

	# for mode = 7:neigvs
	# 	scattersysvec!(u, v[:, mode])
	# 	fld = fieldfromintegpoints(femm, geom, u, :pressure, 1; nodevalmethod = :averaging)
	# 	File =  "cube-mode$(mode).vtk"
	# 	vtkexportmesh(File, fes.conn, fens.xyz, FinEtools.MeshExportModule.H8; vectors = [("u", u.values)], scalars = [("pressure", fld.values)])
	# end
	
	@async run(`"paraview.exe" $File`)

	@test norm(fs.-[0.0, 0.0, 0.0, 1.76692e-8, 1.08952e-7, 1.40754e-7, 0.267052, 0.273879, 0.351477, 0.358597, 0.360482, 0.361265, 0.363118, 0.36379, 0.410238, 0.410963, 0.424434, 0.454963, 0.45869, 0.459053])./norm(fs) < 1.0e-5
end

end
using .mmcubevibrationESNICEH8b
mmcubevibrationESNICEH8b.test()

module minfsuptestik12
using FinEtools
using FinEtools.FEMMDeforLinearESNICEModule: infsup_gh, infsup_sh
using Test
import LinearAlgebra: norm, cholesky, I, eigen
# using UnicodePlots
function test()
	lambdatol = sqrt(1e8*eps(1.0));
	E=1000.0;
	nu=0.24;
	parshiftmult= 0.002;
	A = [1.44 -0.741 -0.53; -0.626 1.589 -0.913; -0.55 0.43 1.756] + 1.0I;

	lambdamin = Float64[]
	h = Float64[]
	for ne = [3, 4, 6]
		Length::FFlt, Width::FFlt, Height::FFlt, nL::FInt, nW::FInt, nH::FInt, orientation::Symbol = (6.0, 6.0, 6.0, ne, ne, ne, :a)

		fens, fes = T4block(Length::FFlt, Width::FFlt, Height::FFlt, nL::FInt, nW::FInt, nH::FInt, orientation::Symbol)
		
		for i = 1:count(fens)
			fens.xyz[i,:] = fens.xyz[i,:] + vec(reshape(fens.xyz[i,:], 1, 3)*A);
		end
		# @show fens.xyz
	
		# File =  "minfsuptest1.vtk"
		# vtkexportmesh(File, fens, fes)
		# try rm(File); catch end

		MR  =  DeforModelRed3D

		material = MatDeforElastIso(MR, E, nu)

		
		geom = NodalField(fens.xyz)
		u = NodalField(zeros(size(fens.xyz,1), 3)) # displacement field
		bfes = meshboundary(fes)
		l1 = connectednodes(bfes)
		setebc!(u, l1, true, 1, 0.0)
		setebc!(u, l1, true, 2, 0.0)
		setebc!(u, l1, true, 3, 0.0)
		numberdofs!(u)

		femm  =  FEMMDeforLinearESNICET4(MR, IntegDomain(fes, NodalSimplexRule(3)), material)
		femm = associategeometry!(femm, geom)
		Gh = infsup_gh(femm, geom, u);
		femm  =  FEMMDeforLinearESNICET4(MR, IntegDomain(fes, NodalSimplexRule(3)), material)
		femm = associategeometry!(femm, geom)
		Sh = infsup_sh(femm, geom, u);

		lambda, modes = eigen(Matrix(Gh), Matrix(Sh));

		# @show lambda
		abslambda = real.(filter(y -> !isnan(y), lambda));
		ix = findall(y  -> y < 0.0, abslambda);
		if !isempty(ix)
			abslambda[ix] .= 0;
		end

		abslambda = sqrt.(sort(abslambda));
		ix = findall(y  -> y > 0.0, abslambda);
		# a = lineplot(1:length(abslambda[ix]), log.(abslambda[ix]), name = "infsup", xlabel = "eigenvalue", ylabel = "log(eigenvalue)", canvas = DotCanvas)
		# display(a)

		ix = findall(y  -> y >= lambdatol, abslambda);
		if isempty(ix)
			@error "Bad guess of the number of eigenvalues"
		end
		push!(lambdamin, abslambda[ix[1]])
		push!(h, 1.0/(count(fens))^(1/3))
	end

	# @show lambdamin
	# a = lineplot(log.(h), log.(lambdamin), name = "infsup", xlabel = "log(Element Size)", ylabel = "log(minimum eigenvalue)", canvas = DotCanvas)
	# display(a)
	
	@test norm(lambdamin - [0.136112, 0.0475213, 0.0505729]) / norm(lambdamin) <= 1.0e-4

end
end
using .minfsuptestik12
minfsuptestik12.test()

module minfsuptestik13
using FinEtools
using FinEtools.FEMMDeforLinearMSModule: infsup_gh, infsup_sh
using Test
import LinearAlgebra: norm, cholesky, I, eigen
# using UnicodePlots
function test()
	lambdatol = sqrt(1e8*eps(1.0));
	E=1000.0;
	nu=0.24;
	parshiftmult= 0.002;
	A = [1.44 -0.741 -0.53; -0.626 1.589 -0.913; -0.55 0.43 1.756] + 1.0I;

	lambdamin = Float64[]
	h = Float64[]
	for ne = [3, 4, 6]
		Length::FFlt, Width::FFlt, Height::FFlt, nL::FInt, nW::FInt, nH::FInt, orientation::Symbol = (6.0, 6.0, 6.0, ne, ne, ne, :a)

		fens, fes = T4block(Length::FFlt, Width::FFlt, Height::FFlt, nL::FInt, nW::FInt, nH::FInt, orientation::Symbol)
		fens, fes = T4toT10(fens, fes)
		
		for i = 1:count(fens)
			fens.xyz[i,:] = fens.xyz[i,:] + vec(reshape(fens.xyz[i,:], 1, 3)*A);
		end
		# @show fens.xyz
	
		# File =  "minfsuptest1.vtk"
		# vtkexportmesh(File, fens, fes)
		# try rm(File); catch end

		MR  =  DeforModelRed3D

		material = MatDeforElastIso(MR, E, nu)

		
		geom = NodalField(fens.xyz)
		u = NodalField(zeros(size(fens.xyz,1), 3)) # displacement field
		bfes = meshboundary(fes)
		l1 = connectednodes(bfes)
		setebc!(u, l1, true, 1, 0.0)
		setebc!(u, l1, true, 2, 0.0)
		setebc!(u, l1, true, 3, 0.0)
		numberdofs!(u)

		femm  =  FEMMDeforLinearMST10(MR, IntegDomain(fes, TetRule(4)), material)
		femm = associategeometry!(femm, geom)
		Gh = infsup_gh(femm, geom, u);
		Sh = infsup_sh(femm, geom, u);

		lambda, modes = eigen(Matrix(Gh), Matrix(Sh));

		# @show lambda
		abslambda = real.(filter(y -> !isnan(y), lambda));
		ix = findall(y  -> y < 0.0, abslambda);
		if !isempty(ix)
			abslambda[ix] .= 0;
		end

		abslambda = sqrt.(sort(abslambda));
		ix = findall(y  -> y > 0.0, abslambda);
		# pl = lineplot(1:length(abslambda[ix]), log.(abslambda[ix]), name = "infsup", xlabel = "eigenvalue", ylabel = "log(eigenvalue)", canvas = DotCanvas)
		# display(pl)

		ix = findall(y  -> y >= lambdatol, abslambda);
		if isempty(ix)
			@error "Bad guess of the number of eigenvalues"
		end
		push!(lambdamin, abslambda[ix[1]])
		push!(h, 1.0/(count(fens))^(1/3))
	end

	# @show lambdamin
	# pl = lineplot(log.(h), log.(lambdamin), name = "infsup", xlabel = "log(Element Size)", ylabel = "log(minimum eigenvalue)", canvas = DotCanvas)
	# display(pl)
	
	# For some reason the results are quite sensitive to numerical precision
	# of the solver (even to the number of bits per floating-point number, 32
	# versus 64)!
	@test norm(lambdamin -  [0.0952635, 0.104529, 0.109738]) / norm(lambdamin) <= 5.0e-2

end
end
using .minfsuptestik13
minfsuptestik13.test()

