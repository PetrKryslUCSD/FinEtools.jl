module miscellaneous_examples
using FinEtools

function circular_hole()
    
    println("""    Meshing  for an elliptical hole in a finite strip.        """)
    
    t0 = time()
    
    H = 100. # strip width
    ell = 50. # distance between holes
    a = 10. # radius of the hole
    L = 200. # length of the strip
    nL = 15
    nH = 10
    nR = 50
    fens,fes = Q4elliphole(a, a, L/2, H/2, nL, nH, nR)
    # Note that we have to number the quadrilaterals  in the mirrored mesh
    # differently in order to generate positive areas.
    fens2, fes2 = mirrormesh(fens, fes, [0.0; -1.0],
    [0.0; 0.0]; renumb = c -> c[end:-1:1])
    fens, fes1, fes2 = mergemeshes(fens, fes,  fens2, fes2, a/1000.)
    fes = cat(fes1, fes2)
    
    geom = NodalField(fens.xyz)
    
    File =  "elliptical.vtk"
    vtkexportmesh(File, fes.conn, geom.values, FinEtools.MeshExportModule.Q4)
    @async run(`"paraview.exe" $File`)
    
    println("Done")
    true
    
end # circular_hole

function L_shaped_mesh()
    println("""
    Meshing an L-shaped membrane: merging multiple meshes.
    """)
    
    t0 = time()
    
    W = 100. # width of the leg
    L = 200. # length of the leg
    nL = 15 # number of elements along the length of the leg
    nW = 10 # number of elements along the width
    tolerance = W / nW / 1.0e5 # tolerance for merging nodes
    Meshes = Array{Tuple{FENodeSet,FESet},1}()
    push!(Meshes, Q4quadrilateral([0.0 0.0; W W], nW, nW))
    push!(Meshes, Q4quadrilateral([-L 0.0; 0.0 W], nL, nW))
    push!(Meshes, Q4quadrilateral([0.0 -L; W 0.0], nW, nL))
    fens, outputfes = mergenmeshes(Meshes, tolerance);
    fes = cat(outputfes[1], cat(outputfes[2], outputfes[3]))
    
    geom = NodalField(fens.xyz)
    
    File =  "L_shape.vtk"
    vtkexportmesh(File, connasarray(fes), geom.values, FinEtools.MeshExportModule.Q4);
    @async run(`"paraview.exe" $File`)
    
    println("Done")
    true
    
end # L_shaped_mesh

function composite_plate_mesh()
    println(""" Meshing  for a layered (composite) plate. """)
    
    t0 = time()
    
    H = 100. # strip width
    L = 200. # length of the strip
    nL = 15
    nH = 10
    xs = collect(linearspace(0.0, L, nL+1))
    ys = collect(linearspace(0.0, H, nH+1))
    ts = [0.5; 2.0; 0.5; 2.0]
    nts = [1; 2; 2; 1; ]
    fens,fes = H8layeredplatex(xs, ys, ts, nts)
    
    geom = NodalField(fens.xyz)
    
    File =  "composite.vtk"
    vtkexportmesh(File, fes.conn, geom.values, FinEtools.MeshExportModule.H8;
    scalars = [( "Layer", fes.label)])
    @async run(`"paraview.exe" $File`)
    
    println("Done")
    true
    
end # composite_plate_mesh

function boundary_Q4_example()
    t0 = time()
    
    rho=1.21*1e-9;# mass density
    c =345.0*1000;# millimeters per second
    bulk= c^2*rho;
    Lx=1900.0;# length of the box, millimeters
    Ly=800.0; # length of the box, millimeters
    
    fens,fes = Q4block(Lx,Ly,3,2); # Mesh
    show(fes.conn)
    
    bfes = meshboundary(fes)
    File =  "B.vtk"
    vtkexportmesh(File, bfes.conn, fens.xyz, FinEtools.MeshExportModule.L2)
    @async run(`"paraview.exe" $File`)
    
end # boundary_Q4_example

function allrun()
    println("#####################################################") 
    println("# circular_hole ")
    circular_hole()
    println("#####################################################") 
    println("# L_shaped_mesh ")
    L_shaped_mesh()
    println("#####################################################") 
    println("# composite_plate_mesh ")
    composite_plate_mesh()
    println("#####################################################") 
    println("# boundary_Q4_example ")
    boundary_Q4_example()
    return true
end # function allrun

end # module miscellaneous_examples
