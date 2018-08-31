module thick_pipe_examples
using FinEtools
using PGFPlotsX
using Statistics
using LinearAlgebra

struct MyIData
    c::FInt
    r::FFltVec
    s::FFltVec
end

function thick_pipe_axi()
    println("Thick pipe with internal pressure: axially symmetric model")
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
    a = 3*phun("MM");
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
    radial_stress(r)  = press*a^2/(b^2-a^2)*(1-b^2.0/r^2);
    ##
    # Circumferential (hoop) stress:
    hoop_stress(r) = press*a^2/(b^2-a^2)*(1+b^2.0/r^2);
    
    ##
    # Radial displacement:
    radial_displacement(r) = press*a^2*(1+nu)*(b^2+r^2*(1-2*nu))/(E*(b^2-a^2)*r);;
    
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
    
    el1femm =  FEMMBase(IntegData(subset(bdryfes,bcl), GaussRule(1, 3), axisymmetric))
    fi = ForceIntensity([press; 0.0]);
    F2= distribloads(el1femm, geom, u, fi, 2);
    
    # Property and material
    material = MatDeforElastIso(MR,  E, nu)
    
    femm = FEMMDeforLinear(MR, IntegData(fes, GaussRule(2, 2), axisymmetric), material)
    
    K = stiffness(femm, geom, u)
    #K=cholesky(K)
    U =  K\(F2)
    scattersysvec!(u,U[:])
    
    # Transfer the solution of the displacement to the nodes on the
    # internal cylindrical surface and convert to
    # cylindrical-coordinate displacements there.
    uv=u.values[internal_fenids,:]
    # Report the  relative displacement on the internal surface:
    println("(Approximate/true displacement) at the internal surface: $( mean(uv[:,1])/urex*100  ) %")
    
    # Produce a plot of the radial stress component in the cylindrical
    # coordinate system. Note that this is the usual representation of
    # stress using nodal stress field.
    
    fld = fieldfromintegpoints(femm, geom, u, :Cauchy, 1)
    
    
    File =  "thick_pipe_sigmax.vtk"
    vtkexportmesh(File, fens, fes; scalars=[("sigmax", fld.values)])
    
    # Produce a plot of the solution components in the cylindrical
    # coordinate system.
    
    function inspector(idat::MyIData, elnum, conn, xe,  out,  xq)
        push!(idat.r, xq[1])
        push!(idat.s, out[idat.c])
        return idat
    end
    
    idat = MyIData(1, FFltVec[], FFltVec[])
    idat = inspectintegpoints(femm, geom, u, collect(1:count(fes)), inspector, idat, :Cauchy)
     
    # Plot the analytical solution.
    r = linearspace(a,b,100);
    @pgf a = Axis({
        xlabel = "Radial distance",
        ylabel = "Radial stress",
        grid="major",
        legend_pos  = "south east"
    },
    Plot({mark="dot"}, Table([:x => vec(r), :y => vec(radial_stress.(r))])), LegendEntry("Analytical"),
    Plot({"only marks", "red", mark="triangle"}, Table([:x => vec(idat.r), :y => vec(idat.s)])), LegendEntry("FEA"))
    display(a)
    # print_tex(a)

    # ##
    # # *Regular quadratic triangle*
    # ##
    # # We start with the workhorse of most commonly used finite element
    # # packages, the quadratic triangle. Similarly to the quadratic
    # # tetrahedron in the 3-D version of this tutorial (pub_thick_pipe), the
    # # stress is polluted with oscillations. Definitely not as bad as in the
    # # plane-strain simulations when the triangles were distorted into
    # # shapes with curved edges, but the disturbances are there.
    # description ='T6';# tetrahedron
    # mf =@T6_block;
    # femmf =@(fes)femm_deformation_linear(struct('fes',fes,...
    #     'material',mater,'integration_rule',tri_rule(struct('npts',3))));
    # surface_integration_rule=gauss_rule(struct('dim',1, 'order', 3));
    # execute_simulation (description, mf, femmf, surface_integration_rule);
    
    # ##
    # # The same remedy of selective reduced integration as in full 3-D models
    # # will also work here.  Is demonstrated by the simulation with the
    # # selective reduced integration quadratic triangle.
    # ##
    # # *Selective reduced
    # # integration quadratic triangle*
    # description ='T6-SRI';
    # mf =@T6_block;
    # femmf =@(fes)femm_deformation_linear_sri(struct('fes',fes,...
    #     'material',mater,...
    #     'integration_rule_volumetric',tri_rule(struct('npts',1)),...
    #     'integration_rule_deviatoric',tri_rule(struct('npts',3))));
    # surface_integration_rule=gauss_rule(struct('dim',1, 'order', 3));
    # execute_simulation (description, mf, femmf, surface_integration_rule);
    
    
    # ##
    # # The selective reduced integration works very well with the T6 triangle.
    # ##
    # # An element that is often used in these situations is the uniformly
    # # under integrated serendipity (8-node) quadrilateral.
    # ##
    # # *Reduced integration serendipity quadrilateral*
    # ##
    # #  The same finite
    # # element model machine as above is used, and the integration is the 2
    # # x 2 Gauss rule (one order lower than that required for full
    # # integration which would be 3 x 3).
    # description ='Q8R';
    # mf =@Q8_block;
    # femmf =@(fes)femm_deformation_linear(struct('fes',fes,...
    #     'material',mater,...
    #     'integration_rule',gauss_rule(struct('dim',2, 'order',2))));
    # surface_integration_rule=gauss_rule(struct('dim',1, 'order',3));
    # execute_simulation (description, mf, femmf, surface_integration_rule);
    
    # ##
    # # *Full integration serendipity quadrilateral*
    # ##
    # # Using the full-integration Gauss rule of 3 x 3 points clearly leads to disaster.
    # description ='Q8';
    # mf =@Q8_block;
    # femmf =@(fes)femm_deformation_linear(struct('fes',fes,...
    #     'material',mater,...
    #     'integration_rule',gauss_rule(struct('dim',2, 'order',3))));
    # surface_integration_rule=gauss_rule(struct('dim',1, 'order',3));
    # execute_simulation (description, mf, femmf, surface_integration_rule);
    
    ##
    # The stress is now totally unacceptable.
    
    
    ## Discussion
    #
    ##
    # The axially symmetric model is clearly very effective
    # computationally, as the size is much reduced compared to the 3-D
    # model.  In conjunction with  uniform or selective reduced integration
    # it can be very accurate as well.
    
    #pub_thick_pipe_axi()
    # end
    
end # thick_pipe_axi


function thick_pipe_ps()
    
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
    
    el1femm =  FEMMBase(IntegData(subset(bdryfes,bcl), GaussRule(1, 3)))
    function pressureloading!(forceout::FFltVec, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)
        copyto!(forceout, XYZ/norm(XYZ)*press)
        return forceout
    end
    fi = ForceIntensity(FFlt, 2, pressureloading!); # pressure normal to the internal cylindrical surface
    F2 = distribloads(el1femm, geom, u, fi, 2);
    
    # Property and material
    material = MatDeforElastIso(MR, E, nu)
    
    femm = FEMMDeforLinear(MR, IntegData(fes, GaussRule(2, 2)), material)
    
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
    println("(Approximate/true displacement) at the internal surface: $( mean(ur)/urex*100  ) %")
    
    # Produce a plot of the radial stress component in the cylindrical
    # coordinate system. Note that this is the usual representation of
    # stress using nodal stress field.
    
    fld = fieldfromintegpoints(femm, geom, u, :Cauchy, 1)
    
    
    File =  "thick_pipe_sigmax.vtk"
    vtkexportmesh(File, fens, fes; scalars=[("sigmax", fld.values)])
    
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
    
    # Plot the analytical solution.
    r = linearspace(a,b,100);
    @pgf a = Axis({
        xlabel = "Radial distance",
        ylabel = "Radial stress",
        grid="major",
        legend_pos  = "south east"
    },
    Plot({mark="dot"}, Table([:x => vec(r), :y => vec(radial_stress.(r))])), LegendEntry("Analytical"),
    Plot({"only marks", "red", mark="triangle"}, Table([:x => vec(idat.r), :y => vec(idat.s)])), LegendEntry("FEA"))
    display(a)
    
end # thick_pipe_ps


function thick_pipe_ps_T6()
    
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
    fens,fes =  T6block(anglrrange, b-a, nc, nt);
    
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
    
    el1femm =  FEMMBase(IntegData(subset(bdryfes,bcl), GaussRule(1, 3)))
    function pressureloading!(forceout::FFltVec, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)
        copyto!(forceout, XYZ/norm(XYZ)*press)
        return forceout
    end
    fi = ForceIntensity(FFlt, 2, pressureloading!); # pressure normal to the internal cylindrical surface
    F2 = distribloads(el1femm, geom, u, fi, 2);
    
    # Property and material
    material = MatDeforElastIso(MR, E, nu)
    
    femm = FEMMDeforLinear(MR, IntegData(fes, TriRule(3)), material)
    
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
    println("(Approximate/true displacement) at the internal surface: $( mean(ur)/urex*100  ) %")
    
    # Produce a plot of the radial stress component in the cylindrical
    # coordinate system. Note that this is the usual representation of
    # stress using nodal stress field.
    
    fld = fieldfromintegpoints(femm, geom, u, :Cauchy, 1)
    
    
    File =  "thick_pipe_sigmax.vtk"
    vtkexportmesh(File, fens, fes; scalars=[("sigmax", fld.values)])
    
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
    
    # Plot the analytical solution.
    r = linearspace(a,b,100);
    @pgf a = Axis({
        xlabel = "Radial distance",
        ylabel = "Radial stress",
        grid="major",
        legend_pos  = "south east"
    },
    Plot({mark="dot"}, Table([:x => vec(r), :y => vec(radial_stress.(r))])), LegendEntry("Analytical"),
    Plot({"only marks", "red", mark="triangle"}, Table([:x => vec(idat.r), :y => vec(idat.s)])), LegendEntry("FEA"))
    display(a)
    
    #pub_thick_pipe_ps()
    
end # thick_pipe_ps_T6

function allrun()
    println("#####################################################") 
    println("# thick_pipe_axi ")
    thick_pipe_axi()
    println("#####################################################") 
    println("# thick_pipe_ps ")
    thick_pipe_ps()
    println("#####################################################") 
    println("# thick_pipe_ps_T6 ")
    thick_pipe_ps_T6()
    return true
end # function allrun

end # module thick_pipe_examples
