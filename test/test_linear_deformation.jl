module sscratch_06112017
using FinEtools
using Base.Test

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

  femm = FEMMDeforLinear(MR, GeoD(fes, GaussRule(3, 3)), material)

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
  xsfemm = FEMMDeforWinkler(GeoD(subset(bfes,fl), GaussRule(2, 3)))

  ##
  # We create the temperature field using the formula $T=r+z$.
  dT = NodalField(reshape(sqrt.(fens.xyz[:,1].^2+fens.xyz[:,2].^2)+fens.xyz[:,3],size(fens.xyz,1),1));

  ##
  # And we are ready to assemble the system matrix. Both the elastic stiffness of
  # the hexahedral elements ...
  K = stiffness(femm, geom, u)
  # ...  and the elastic stiffness    of the springs on the contact surfaces of the cross-sections.
  H = surfacenormalspringstiffness(xsfemm,  geom, u, springcoefficient)

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
using sscratch_06112017
sscratch_06112017.test()

module cookstress_1
using Base.Test
using FinEtools
using FinEtools.MeshExportModule

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
  el1femm =  FEMMBase(GeoD(subset(boundaryfes, Toplist),  GaussRule(1, 2)))
  fi = ForceIntensity([0.0, +magn]);
  F2 = distribloads(el1femm,  geom,  u,  fi,  2);


  MR = DeforModelRed2DStress
  material = MatDeforElastIso(MR,  0.0, E, nu, 0.0)

  femm = FEMMDeforLinear(MR, GeoD(fes,  TriRule(1)),  material)

  K = stiffness(femm,  geom,  u)
  K = cholfact(K)
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
using cookstress_1
cookstress_1.test()


module scratch1_06092017

using FinEtools
using Base.Test

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

  el1femm =  FEMMBase(GeoD(subset(bdryfes,bcl), GaussRule(1, 3), axisymmetric))
  fi = ForceIntensity([press; 0.0]);
  F2= distribloads(el1femm, geom, u, fi, 2);

  # Property and material
  material = MatDeforElastIso(MR,  E, nu)

  femm = FEMMDeforLinear(MR, GeoD(fes, GaussRule(2, 2), axisymmetric), material)

  K = stiffness(femm, geom, u)
  #K=cholfact(K)
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
  # r = linspace(a,b,100);
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
using scratch1_06092017
scratch1_06092017.test()

module scratch2_06102017

using FinEtools
using Base.Test

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

  femm = FEMMDeforLinear(MR, GeoD(fes, GaussRule(3,2)), material)

  K =stiffness(femm, geom, u)
  femm = FEMMDeforLinear(MR, GeoD(fes, GaussRule(3,3)), material)
  M =mass(femm, geom, u)
  d,v,nev,nconv = eigs(K+OmegaShift*M, M; nev=neigvs, which=:SM)
  d = d - OmegaShift;
  fs = real(sqrt.(complex(d)))/(2*pi)
  # println("Eigenvalues: $fs [Hz]")

  # mode = 17
  # scattersysvec!(u, v[:,mode])
  # File =  "unit_cube_modes.vtk"
  # vtkexportmesh(File, fens, fes; vectors=[("mode$mode", u.values)])

  @test abs(fs[7]-0.26259869196259) < 1.0e-5
end
end

using scratch2_06102017
scratch2_06102017.test()

module mxxxx1_06102017

using FinEtools
using FinEtools.AlgoDeforLinearModule
using Base.Test

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
  el1femm  = FEMMBase(GeoD(subset(boundaryfes,Toplist), GaussRule(2, 2)))
  flux1 = FDataDict("femm"=>el1femm, "traction_vector"=>loadv)


  # Make the region
  MR = DeforModelRed3D
  material = MatDeforElastIso(MR, 00.0, E, nu, 0.0)
  region1 = FDataDict("femm"=>FEMMDeforLinear(MR, GeoD(fes, GaussRule(3,2)),
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
  theutip = mean(u.values[nl,:],1)
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
using mxxxx1_06102017
mxxxx1_06102017.test()




module mx_06112017
using FinEtools
using Base.Test
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

  femm = FEMMDeforLinear(MR, GeoD(fes, GaussRule(3, 3)), material)

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
  xsfemm = FEMMDeforWinkler(GeoD(subset(bfes,fl), GaussRule(2, 3)))

  ##
  # We create the temperature field using the formula $T=r+z$.
  dT = NodalField(reshape(sqrt.(fens.xyz[:,1].^2+fens.xyz[:,2].^2)+fens.xyz[:,3],size(fens.xyz,1),1));

  ##
  # And we are ready to assemble the system matrix. Both the elastic stiffness of
  # the hexahedral elements ...
  K = stiffness(femm, geom, u)
  # ...  and the elastic stiffness    of the springs on the contact surfaces of the cross-sections.
  H = surfacenormalspringstiffness(xsfemm,  geom, u, springcoefficient)

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
using mx_06112017
mx_06112017.test()


module my_06112017

using FinEtools
using Base.Test
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

  femm  =  FEMMDeforLinear(MR, GeoD(fes, GaussRule(2, 3), true), material)

  K  = stiffness(femm, geom, u)
  F  =  thermalstrainloads(femm, geom, u, dT)
  #K = cholfact(K)
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

  fen2fe  = FENodeToFEMap(fes.conn, nnodes(geom))
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
using my_06112017
my_06112017.test()

module mmmmmmmmmZenkourm
using FinEtools
using FinEtools.AlgoDeforLinearModule
using Base.Test
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
  xs = collect(linspace(0.0, a, na+1))
  ys = collect(linspace(0.0, b, nb+1))
  ts = h/nLayers*ones(nLayers);# layer thicknesses
  nts= 3*ones(Int, nLayers);# number of elements per layer
  fens,fes = H8compositeplatex(xs, ys, ts, nts)
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
      GeoD(fes, gr, CSys(3, 3, updatecs!)), laminamaterial))

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
      "femm"=>FEMMBase(GeoD(subset(bfes, ttopl), GaussRule(2, 3))))

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
using mmmmmmmmmZenkourm
mmmmmmmmmZenkourm.test()

module mmmmmmultimaterial_beam_xz
using FinEtools
using FinEtools.AlgoDeforLinearModule
using Base.Test
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
  el1femm  =   FEMMBase(GeoD(subset(boundaryfes,Toplist), GaussRule(2, 2)))
  flux1 = FDataDict("femm"=>el1femm, "traction_vector"=>loadv)

  r1list   = selectelem(fens,fes, box =  [0 L/2. -Inf Inf -Inf Inf], inflate =   tolerance);
  r2list   = selectelem(fens,fes, box =  [L/2. L -Inf Inf -Inf Inf], inflate =   tolerance);

  # Model reduction type
  MR = DeforModelRed3D

  # Make region 1
  region1 = FDataDict("femm"=>FEMMDeforLinear(MR,
  GeoD(subset(fes,r1list), GaussRule(3,2)),
  MatDeforElastIso(MR, 0.0, E1, nu1, 0.0)))

  # Make region 2
  region2 = FDataDict("femm"=>FEMMDeforLinear(MR,
  GeoD(subset(fes,r2list), GaussRule(3,2)),
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
  for f in modeldata["postprocessing"]["exported_files"]
    rm(f)
  end

  # Write out mesh with stresses
  modeldata["postprocessing"] = FDataDict("file"=>"multimaterial_beam_xy",
  "quantity"=> :Cauchy, "component"=> :xy)
  modeldata = AlgoDeforLinearModule.exportstress(modeldata)
  for f in modeldata["postprocessing"]["exported_files"]
    rm(f)
  end

  # Write out mesh with stresses
  modeldata["postprocessing"] = FDataDict("file"=>"multimaterial_beam_xz",
  "quantity"=> :Cauchy, "component"=> :xz)
  modeldata = AlgoDeforLinearModule.exportstress(modeldata)
  for f in modeldata["postprocessing"]["exported_files"]
    rm(f)
  end

  # Write out mesh with von Mises stresses
  modeldata["postprocessing"] = FDataDict("file"=>"multimaterial_beam_vm",
    "quantity"=> :vm)
  modeldata = AlgoDeforLinearModule.exportstress(modeldata)
  for f in modeldata["postprocessing"]["exported_files"]
    rm(f)
  end

  # Write out mesh with von Mises stresses, elementwise
  modeldata["postprocessing"] = FDataDict("file"=>"multimaterial_beam_vm-ew",
   "quantity"=> :vm)
  modeldata = AlgoDeforLinearModule.exportstresselementwise(modeldata)
  for f in modeldata["postprocessing"]["exported_files"]
    rm(f)
  end

  true

end
end
using mmmmmmultimaterial_beam_xz
mmmmmmultimaterial_beam_xz.test()

module mmmmmmmunitmccubemm
using FinEtools
using Base.Test

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
  region1 = FDataDict("femm"=>FEMMDeforLinear(MR, GeoD(fes, GaussRule(3,2)),
    material), "femm_mass"=>FEMMDeforLinear(MR, GeoD(fes, GaussRule(3,3)),
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
  for f in modeldata["postprocessing"]["exported_files"]
    rm(f)
  end

  true
end
end
using mmmmmmmunitmccubemm
mmmmmmmunitmccubemm.test()

module mmmmmpipemmPSmmm
using FinEtools
using Base.Test

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

el1femm =  FEMMBase(GeoD(subset(bdryfes,bcl), GaussRule(1, 3)))
function pressureloading!(forceout::FFltVec, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)
  copy!(forceout, XYZ/norm(XYZ)*press)
  return forceout
end
fi = ForceIntensity(FFlt, 2, pressureloading!); # pressure normal to the internal cylindrical surface
F2 = distribloads(el1femm, geom, u, fi, 2);

# Property and material
material = MatDeforElastIso(MR, E, nu)

femm = FEMMDeforLinear(MR, GeoD(fes, GaussRule(2, 2)), material)

K =stiffness(femm, geom, u)
#K=cholfact(K)
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
# r = linspace(a,b,100);
# plot(r, radial_stress(r))
# # Plot the computed  integration-point data
# scatter!(idat.r, idat.s, m=:circle, color=:red)
# gui()
end
end
using mmmmmpipemmPSmmm
mmmmmpipemmPSmmm.test()

module mmmmmOrthotropicmmmmm
using FinEtools
using Base.Test
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

  el1femm =  FEMMBase(GeoD(subset(bdryfes,icl), GaussRule(1, 3), true))
  function pressureloading!(forceout::FFltVec, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)
    copy!(forceout, XYZ/norm(XYZ)*p)
    return forceout
  end
  fi = ForceIntensity(FFlt, 2, pressureloading!); # pressure normal to the internal cylindrical surface
  F2= distribloads(el1femm, geom, u, fi, 2);

  # Property and material
  material=MatDeforElastOrtho(MR, E1,E2,E3,nu12,nu13,nu23,G12,G13,G23)

  femm = FEMMDeforLinear(MR, GeoD(fes, GaussRule(2, 2), true), material)

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
  rm(File)

  @test abs(minimum(fld.values)-(-0.05139098099088048)) < 1.0e-5
  @test abs(maximum(fld.values)-0.5704453140236726) < 1.0e-5

end
end
using mmmmmOrthotropicmmmmm
mmmmmOrthotropicmmmmm.test()

module mmmmmmmmCookmm
using FinEtools
using FinEtools.AlgoDeforLinearModule
using Base.Test
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
  el1femm = FEMMBase(GeoD(subset(boundaryfes, Toplist), GaussRule(1, 3)))
  flux1 = FDataDict("traction_vector"=>[0.0,+magn],
      "femm"=>el1femm
      )

  # Make the region
  MR = DeforModelRed2DStress
  material = MatDeforElastIso(MR,  0.0, E, nu, 0.0)
  region1 = FDataDict("femm"=>FEMMDeforLinear(MR,
      GeoD(fes, TriRule(3)), material))

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
  File = modeldata["postprocessing"]["exported_files"][1]
  # @async run(`"paraview.exe" $File`)
  for f in modeldata["postprocessing"]["exported_files"]
    rm(f)
  end
  fld = modeldata["postprocessing"]["exported_fields"][1]
  # println("$(minimum(fld.values)) $(maximum(fld.values))")
  @test norm([minimum(fld.values) maximum(fld.values)] -
    [-0.06292574794975273 0.12022571940768388]) < 1.0e-5
end
end
using mmmmmmmmCookmm
mmmmmmmmCookmm.test()
