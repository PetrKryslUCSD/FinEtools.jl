module sphere_examples
using FinEtools
using FinEtools.MeshExportModule


function sphere_dipole_1()
    println("The interior sphere accelerates in the positive x-direction, generating
    positive pressure ahead of it, negative pressure behind
    ")
    rho = 1.21*phun("kg/m^3");# mass density
    c  = 343.0*phun("m/s");# sound speed
    bulk =  c^2*rho;
    a_amplitude=1.*phun("mm/s^2");# amplitude of the  acceleration of the sphere
    # omega = 2000*phun("rev/s");      # frequency of the incident wave
    R = 5.0*phun("mm"); # radius of the interior sphere
    Ro = 4*R # radius of the external sphere
    P_amplitude = R*rho*a_amplitude; # pressure amplitude
    nref=2;
    nlayers=40;
    tolerance = Ro/(nlayers)/100
    frequency=2000; # Hz
    omega=2*pi*frequency;
    k = omega/c*[1.0 0.0 0.0];
    
    println("""
    Spherical domain filled with acoustic medium, no scatterer.
    Hexahedral mesh.
    """)
    
    t0  =  time()
    
    # Hexahedral mesh
    fens,fes  =  H8sphere(R,nref);
    fens1,fes1  =  mirrormesh(fens, fes, [-1.0, 0.0, 0.0], [0.0, 0.0, 0.0],
    renumb =  r(c) = c[[1, 4, 3, 2, 5, 8, 7, 6]]);
    fens,newfes1,fes2 =  mergemeshes(fens1, fes1, fens, fes, tolerance)
    fes = cat(newfes1,fes2)
    
    
    # Derive the finite element sets for the boundary
    bfes  =  meshboundary(fes)
    # Outer spherical boundary
    function dout(xyz)
        return xyz/norm(xyz)
    end
    
    louter = selectelem(fens, bfes, facing = true, direction = dout)
    outer_fes = subset(bfes, louter);
    
    fens,fes  = H8extrudeQ4(fens, outer_fes,
    nlayers, (xyz, layer)->(R+layer/nlayers*(Ro-R))*xyz/norm(xyz));
    connected = findunconnnodes(fens, fes);
    fens, new_numbering = compactnodes(fens, connected);
    fess = renumberconn!(fes, new_numbering);
    
    geom  =  NodalField(fens.xyz)
    P  =  NodalField(zeros(FCplxFlt,size(fens.xyz,1),1))
    
    bfes  =  meshboundary(fes)
    
    # File  =   "Sphere.vtk"
    # vtkexportmesh(File, bfes.conn, geom.values, FinEtools.MeshExportModule.Q4)
    # @async run(`"paraview.exe" $File`)
    
    # File  =   "Sphere.vtk"
    # vtkexportmesh (File, fes.conn, fens.xyz, MeshExportModule.H8)
    # @async run(`"C:/Program Files (x86)/ParaView 4.2.0/bin/paraview.exe" $File`)
    
    linner = selectelem(fens, bfes, distance = R, from = [0.0 0.0 0.0],
    inflate = tolerance)
    louter = selectelem(fens, bfes, facing = true, direction = dout)
    
    println("Pre-processing time elapsed  =  ",time() - t0,"s")
    
    t1  =  time()
    
    
    numberdofs!(P)
    
    material = MatAcoustFluid(bulk,rho)
    femm  =  FEMMAcoust(IntegData(fes, GaussRule(3, 2)), material)
    
    S  =  acousticstiffness(femm, geom, P);
    C  =  acousticmass(femm, geom, P);
    
    abcfemm  =  FEMMAcoustSurf(IntegData(subset(bfes, louter), GaussRule(2, 2)), material)
    D  =  acousticABC(abcfemm, geom, P);
    
    # Inner sphere pressure loading
    function dipole(dpdn, xyz, J, label)
        n = cross(J[:,1],J[:,2]);
        n = vec(n/norm(n));
        #println("$( (-1.0im)*dot(vec(k),n)*exp(-1.0im*dot(vec(k),vec(xyz))) )")
        dpdn[1] = -rho*a_amplitude*n[1]
    end
    
    fi  =  ForceIntensity(FCplxFlt, 1, dipole);
    dipfemm  =  FEMMAcoustSurf(IntegData(subset(bfes, linner), GaussRule(2, 2)), material)
    F  = distribloads(dipfemm, geom, P, fi, 2);
    
    K = lufact((1.0+0.0im)*(-omega^2*S + omega*1.0im*D + C)) # We fake a complex matrix here
    p = K\F  #
    
    scattersysvec!(P, p[:])
    
    println(" Minimum/maximum pressure= $(minimum(real(p)))/$(maximum(real(p)))")
    
    println("Computing time elapsed  =  ",time() - t1,"s")
    println("Total time elapsed  =  ",time() - t0,"s")
    
    File  =   "sphere_dipole_1.vtk"
    vtkexportmesh(File, connasarray(fes), geom.values, FinEtools.MeshExportModule.H8;
    scalars = [( "realP", real(P.values))])
    @async run(`"paraview.exe" $File`)
    
    true
    
end # sphere_dipole_1


function sphere_inc_example()
    rho = 1.21*phun("kg/m^3");# mass density
    c  = 343.0*phun("m/s");# sound speed
    bulk =  c^2*rho;
    omega = 2000*phun("rev/s");      # frequency of the incident wave
    Ro = 300.0*phun("mm"); # radius of the enclosure
    nPerradius = 18;#number of elements on the radius of the sphere
    tolerance = Ro/(2^nPerradius)/100
    pincampl = 1.0*phun("kilo*Pa");
    k = omega/c*[1.0 0.0 0.0];
    
    println("""
    Spherical domain filled with acoustic medium, no scatterer.
    Incident wave.
    Hexahedral mesh.
    """)
    
    t0  =  time()
    
    # Hexahedral mesh
    fens,fes  =  H8spheren(Ro, nPerradius);
    fens1,fes1 = mirrormesh(fens, fes, [-1.0, 0.0, 0.0], [0.0, 0.0, 0.0];
    renumb =  r(c) = c[[1, 4, 3, 2, 5, 8, 7, 6]])
    fens,newfes1,fes2 =  mergemeshes(fens1, fes1, fens, fes, tolerance)
    fes = cat(newfes1,fes2)
    
    # Derive the finite element sets for the boundary
    bfes  =  meshboundary(fes)
    # Outer spherical boundary
    function dout(xyz)
        return xyz/norm(xyz)
    end
    
    louter = selectelem(fens, bfes, facing = true, direction = dout)
    outer_fes = subset(bfes, louter);
    
    # File  =   "Sphere.vtk"
    # vtkexportmesh (File, fes.conn, fens.xyz, MeshExportModule.H8)
    # @async run(`"C:/Program Files (x86)/ParaView 4.2.0/bin/paraview.exe" $File`)
    
    println("Pre-processing time elapsed  =  ",time() - t0,"s")
    
    t1  =  time()
    
    geom  =  NodalField(fens.xyz)
    P  =  NodalField(zeros(FCplxFlt,size(fens.xyz,1),1))
    
    numberdofs!(P)
    
    pinc = deepcopy(P);
    
    material = MatAcoustFluid(bulk,rho)
    femm  =  FEMMAcoust(IntegData(fes, GaussRule(3, 2)), material)
    
    S  =  acousticstiffness(femm, geom, P);
    C  =  acousticmass(femm, geom, P);
    
    abcfemm  =  FEMMAcoustSurf(IntegData(outer_fes, GaussRule(2, 2)), material)
    D  =  acousticABC(abcfemm, geom, P);
    
    pincf(loc)  =  pincampl*exp((-1.0im*(vec(k)'*vec(loc)))[1])
    
    # Incident pressure loading
    for j = 1:size(geom.values,1)
        pinc.values[j] = pincampl[1,1]*exp((-1.0im*(vec(k)'*vec(geom.values[j,:])))[1])
    end
    vpinc = gathersysvec(pinc)
    
    F  =  0.0 * vpinc;                # Start with a zero vector
    
    #pfemm = FEMMBase(inner_fes, GaussRule(order = 2,dim = 2))
    function dpincdn(dpdn, xyz, J, label)
        n = cross(J[:,1],J[:,2]);
        n = vec(n/norm(n));
        #println("$( (-1.0im)*dot(vec(k),n)*exp(-1.0im*dot(vec(k),vec(xyz))) )")
        dpdn[1] =  pincampl*(-1.0im)*dot(vec(k),n)*exp(-1.0im*dot(vec(k),vec(xyz)))
    end
    
    fi  =  ForceIntensity(FCplxFlt, 1, dpincdn);
    #F  =  F - distribloads(pfemm, nothing, geom, P, fi, 2);
    F  =  F + distribloads(abcfemm, geom, P, fi, 2);
    
    K = lufact((1.0+0.0im)*(-omega^2*S + C)) # We fake a complex matrix here
    p = K\F  #+omega*1.0im*D
    
    scattersysvec!(P,p[:])
    
    println("Computing time elapsed  =  ",time() - t1,"s")
    println("Total time elapsed  =  ",time() - t0,"s")
    
    ptot = deepcopy(P)
    ptot.values = P.values+pinc.values
    
    # File  =   "Sphereptot.vtk"
    # vtkexportmesh (File, fes.conn, geom.values, MeshExportModule.H8; scalars = abs(ptot.values), scalars_name  = "ptot")
    # @async run(`"C:/Program Files (x86)/ParaView 4.2.0/bin/paraview.exe" $File`)
    
    File  =   "Spherepinc.vtk"
    vtkexportmesh(File, connasarray(fes), geom.values, FinEtools.MeshExportModule.H8;
    scalars = [("realpinc", real(pinc.values))])
    @async run(`"paraview.exe" $File`)
    
    File  =   "SphereP.vtk"
    vtkexportmesh(File, connasarray(fes), geom.values, FinEtools.MeshExportModule.H8;
    scalars = [( "realP", real(P.values))])
    @async run(`"paraview.exe" $File`)
    
    true
    
end # sphere_inc_example


function sphere_scatterer_example()
    rho = 1.21*phun("kg/m^3");# mass density
    c  = 343.0*phun("m/s");# sound speed
    bulk =  c^2*rho;
    R = 50.0*phun("mm");# radius of the piston
    Ro = 300.0*phun("mm"); # radius of the enclosure
    nPerradius=24;# number of elements along the radius of the scatterer
    nlayers=78;                     # number of layers of elements surrounding the piston
    tolerance = R/(nPerradius)/100
    kR =  1.0;
    k = kR/R*[1.0 0.0 0.0]; k = vec(k);
    omega = norm(k)*c;    # angular frequency of the incident wave
    pincampl = 1.0*phun("Pa");
    h = maximum([(Ro-R)/nlayers pi*Ro/(nPerradius)]);
    
    println("""
    
    Rigid sphere scattering.
    
    Frequency: $(omega/2/pi). Wave vector: $k. kR = $(norm(k)*R).
    Number of elements per wavelength: $(2*pi/(norm(k)*h)).
    
    Hexahedral mesh.
    """)
    
    t0  =  time()
    
    # Hexahedral mesh
    fens, fes  =  H8spheren(R, nPerradius);
    bfes  =  meshboundary(fes)
    l = selectelem(fens, bfes, facing = true, direction = [1.0 1.0  1.0], dotmin= 0.001)
    ex(xyz, layer) = (R+layer/nlayers*(Ro-R))*xyz/norm(xyz)
    fens,fes  =  H8extrudeQ4(fens, subset(bfes,l), nlayers, ex);
    fens1,fes1 = mirrormesh(fens, fes, [-1.0, 0.0, 0.0], [0.0, 0.0, 0.0]; renumb =  r(c) = c[[1, 4, 3, 2, 5, 8, 7, 6]])
    fens,newfes1,fes2 =  mergemeshes(fens1, fes1, fens, fes, tolerance)
    fes = cat(newfes1,fes2)
    
    # Derive the finite element sets for the boundary
    bfes  =  meshboundary(fes)
    # Outer spherical boundary
    function dout(xyz)
        return xyz/norm(xyz)
    end
    
    louter = selectelem(fens, bfes, facing = true, direction = dout)
    outer_fes = subset(bfes,louter);
    
    # Inner spherical boundary
    function din(xyz)
        return -xyz/norm(xyz)
    end
    linner = selectelem(fens, bfes, facing = true, direction = din)
    inner_fes = subset(bfes,linner);
    
    # File  =   "Sphere.vtk"
    # vtkexportmesh(File, inner_fes.conn, fens.xyz, MeshExportModule.Q4)
    # @async run(`"paraview.exe" $File`)
    println("Pre-processing time elapsed  =  ",time() - t0,"s")
    
    t1  =  time()
    
    geom  =  NodalField(fens.xyz)
    P  =  NodalField(zeros(FCplxFlt,size(fens.xyz,1),1))
    
    numberdofs!(P)
    
    pinc = deepcopy(P);
    
    material = MatAcoustFluid(bulk,rho)
    femm  =  FEMMAcoust(IntegData(fes, GaussRule(3, 3)), material)
    
    S  =  acousticstiffness(femm, geom, P);
    C  =  acousticmass(femm, geom, P);
    
    abcfemm  =  FEMMAcoustSurf(IntegData(outer_fes, GaussRule(2, 3)), material)
    D  =  acousticABC(abcfemm, geom, P);
    
    # Incident pressure loading
    for j = 1:size(geom.values,1)
        xyz = vec(geom.values[j,:]);
        pinc.values[j] = pincampl*exp((-1.0im*(k'*xyz))[1])
    end
    vpinc = gathersysvec(pinc)
    
    F  =  - (-omega^2*S + C) * vpinc;
    
    #pfemm = FEMMBase(inner_fes, GaussRule(order = 2,dim = 2))
    function dpincdn(dpdn, xyz, J, label)
        xyz = vec(xyz);
        n = cross(J[:,1],J[:,2]);
        n = vec(n/norm(n));
        dpdn[1] = pincampl*(-1.0im)*dot(vec(k),n)*exp(-1.0im*dot(vec(k),vec(xyz)))
    end
    
    fi  =  ForceIntensity(FCplxFlt, 1, dpincdn);
    F  =  F + distribloads(abcfemm, geom, P, fi, 2);
    
    K = ((-omega^2*S + omega*1.0im*D + C))
    p = K\F
    
    scattersysvec!(P,p[:])
    
    println("Computing time elapsed  =  ",time() - t1,"s")
    println("Total time elapsed  =  ",time() - t0,"s")
    
    ptot = deepcopy(P)
    ptot.values = P.values+pinc.values
    
    # File  =   "Sphereptot.vtk"
    # vtkexportmesh (File, fes.conn, geom.values, MeshExportModule.H8; scalars = abs(ptot.values), scalars_name  = "ptot")
    # @async run(`"C:/Program Files (x86)/ParaView 4.2.0/bin/paraview.exe" $File`)
    
    # File  =   "Spherepinc.vtk"
    # vtkexportmesh (File, fes.conn, geom.values, MeshExportModule.H8; scalars = real(pinc.values), scalars_name  = "pinc")
    # @async run(`"C:/Program Files (x86)/ParaView 4.2.0/bin/paraview.exe" $File`)
    
    File  =   "SphereP.vtk"
    vtkexportmesh(File, connasarray(fes), geom.values, FinEtools.MeshExportModule.H8; scalars = [("P", abs.(P.values))])
    @async run(`"paraview.exe" $File`)
    
    # File  =   "SphereP.vtk"
    # vtkexportmesh (File, cat(outer_fes,inner_fes).conn, geom.values, MeshExportModule.Q4; scalars = abs(P.values), scalars_name  = "P")
    # @async run(`"C:/Program Files (x86)/ParaView 4.2.0/bin/paraview.exe" $File`)
    
    # using Winston
    # pl  =  FramedPlot(title = "Matrix",xlabel = "x",ylabel = "Re P, Im P")
    # setattr(pl.frame, draw_grid = true)
    # add(pl, Curve([1:length(C[:])],vec(C[:]), color = "blue"))
    
    # # pl = plot(geom.values[nLx,1][ix],scalars[nLx][ix])
    # # xlabel("x")
    # # ylabel("Pressure")
    # display(pl)
    
    true
    
end # sphere_scatterer_example

function allrun()
    println("#####################################################") 
    println("# sphere_dipole_1 ")
    sphere_dipole_1()
    println("#####################################################") 
    println("# sphere_inc_example ")
    sphere_inc_example()
    println("#####################################################") 
    println("# sphere_scatterer_example ")
    sphere_scatterer_example()
end # function allrun

end # module sphere_examples
