module trunc_cyl_shell_examples
using FinEtools


function trunc_cyl_shell_mesh()
    h = 0.05 * phun("M");
    l = 10 * h;
    Rmed = h / 0.2;
    psi   = 0;    # Cylinder
    nh = 12; nl  = 40; nc = 120;
    nh = 6; nl  = 20; nc = 60;
    # nh = 3; nl  = 8; nc = 30;
    tolerance = h / nh / 100;
    
    t0 = time()
    MR = DeforModelRed3D
    fens, fes  = H8block(h, l, 2.0 * pi, nh, nl, nc)
    # Shape into a cylinder
    R = zeros(3, 3)
    for i = 1:count(fens)
        x, y, z = fens.xyz[i,:];
        rotmat3!(R, [0, z, 0])
        Q = [cos(psi * pi / 180) sin(psi * pi / 180) 0;
        -sin(psi * pi / 180) cos(psi * pi / 180) 0;
        0 0 1]
        fens.xyz[i,:] = reshape([x + Rmed - h / 2, y - l / 2, 0], 1, 3) * Q * R;
    end
    println("  before merging  = $(count(fens))")
    # File =  "unit_cube_modes.vtk"
    # vtkexportmesh(File, fens, fes)
    # @async run(`"paraview.exe" $File`)
    
    candidates = selectnode(fens, box = boundingbox([Rmed - h -Inf 0.0; Rmed + h +Inf 0.0]), inflate = tolerance)
    fens, fes = mergenodes(fens, fes,  tolerance, candidates);
    
    # fens,fes = mergenodes(fens, fes,  tolerance);
    
    println("  after merging  = $(count(fens))")
    
    println("Mesh generation ($(time() - t0) sec)")
    
    File =  "trunc_cyl_shell_mesh.vtk"
    vtkexportmesh(File, fens, fes)
    @async run(`"paraview.exe" $File`)
    
end # trunc_cyl_shell_mesh


function trunc_cyl_shell_mesh_2()
    h = 0.05*phun("M");
    l = 10*h;
    Rmed = h/0.2;
    psi   = 0;    # Cylinder
    nh = 12; nl  = 40; nc = 120;
    nh = 6; nl  = 20; nc = 60;
    # nh = 3; nl  = 8; nc = 30;
    tolerance = h/nh/100;
    
    t0 = time()
    MR = DeforModelRed3D
    fens,fes  = H8block(h,l,2.0*pi,nh,nl,nc)
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
    println("  before merging  = $(count(fens))")
    # File =  "unit_cube_modes.vtk"
    # vtkexportmesh(File, fens, fes)
    # @async run(`"paraview.exe" $File`)
    
    candidates = selectnode(fens, plane = [0.0 0.0 1.0 0.0], inflate = tolerance)
    fens,fes = mergenodes(fens, fes,  tolerance, candidates);
    
    # fens,fes = mergenodes(fens, fes,  tolerance);
    
    println("  after merging  = $(count(fens))")
    
    println("Mesh generation ($(time() - t0) sec)")
    
    File =  "trunc_cyl_shell_mesh_2.vtk"
    vtkexportmesh(File, fens, fes)
    @async run(`"paraview.exe" $File`)
    
end # trunc_cyl_shell_mesh_2

function allrun()
    println("#####################################################") 
    println("# trunc_cyl_shell_mesh ")
    trunc_cyl_shell_mesh()
    println("#####################################################") 
    println("# trunc_cyl_shell_mesh_2 ")
    trunc_cyl_shell_mesh_2()
    return true
end # function allrun

end # module trunc_cyl_shell_examples
