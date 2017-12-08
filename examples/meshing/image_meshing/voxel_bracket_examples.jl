module voxel_bracket_examples
using FinEtools

function voxel_bracket_mesh()
    V = VoxelBoxVolume(Int, 6*[5,6,7], [4.0, 4.0, 5.0])
    
    b1 = solidbox((0.0, 0.0, 0.0), (1.0, 4.0, 5.0))
    b2 = solidbox((0.0, 0.0, 0.0), (4.0, 1.0, 5.0))
    h1 = solidcylinder((2.0, 2.5, 2.5), (1.0, 0.0, 0.0), 00.75)
    fillsolid!(V, differenceop(unionop(b1, b2), h1), 1)
    
    fens, fes = H8voximg(V.data, vec([voxeldims(V)...]), [1])
    fens = meshsmoothing(fens, fes; npass = 5)
    println("count(fes) = $(count(fes))")
    
    File = "voxel_bracket_mesh.vtk"
    vtkexportmesh(File, fens, fes)
    @async run(`"paraview.exe" $File`)
    
end # voxel_bracket_mesh


function voxel_bracket_mesh_tet()
    V = VoxelBoxVolume(Int, 6*[5,6,7], [4.0, 4.0, 5.0])
    
    b1 = solidbox((0.0, 0.0, 0.0), (1.0, 4.0, 5.0))
    b2 = solidbox((0.0, 0.0, 0.0), (4.0, 1.0, 5.0))
    h1 = solidcylinder((2.0, 2.5, 2.5), (1.0, 0.0, 0.0), 00.75)
    fillsolid!(V, differenceop(unionop(b1, b2), h1), 1)
    
    fens, fes = T4voximg(V.data, vec([voxeldims(V)...]), [1])
    fens = meshsmoothing(fens, fes; method = :laplace, npass = 5)
    println("count(fes) = $(count(fes))")
    
    File = "voxel_bracket_mesh_tet.vtk"
    vtkexportmesh(File, fens, fes)
    @async run(`"paraview.exe" $File`)
    
end # voxel_bracket_mesh_tet


function voxel_bracket_mesh_tet_remesh()
    V = VoxelBoxVolume(Int, 8*[5,6,7], [4.0, 4.0, 5.0])
    
    b1 = solidbox((0.0, 0.0, 0.0), (1.0, 4.0, 5.0))
    b2 = solidbox((0.0, 0.0, 0.0), (4.0, 1.0, 5.0))
    h1 = solidcylinder((2.0, 2.5, 2.5), (1.0, 0.0, 0.0), 0.75)
    fillsolid!(V, differenceop(unionop(b1, b2), h1), 1)
    
    function checkvolumes(im)
        for i = 1:size(im.t, 1)
            if FinEtools.MeshTetrahedronModule.tetv1times6(im.v, im.t[i,1], im.t[i,2], im.t[i,3], im.t[i,4]) < 0.0
                println("Negative volume = $([im.t[i,1], im.t[i,2], im.t[i,3], im.t[i,4]])")
            end 
        end
    end
    
    im = ImageMesher(V, zero(eltype(V.data)), eltype(V.data)[1])
    mesh!(im)
    checkvolumes(im)
    println("Mesh size: initial = $(size(im.t,1))")
    
    im.elementsizeweightfunctions = [ElementSizeWeightFunction(20.0, vec([0.0, 2.5, 2.5]), 1.0), ElementSizeWeightFunction(1.0, vec([0.0, 2.5, 2.5]), 3.5)]
    for i = 1:12
        println("Phase $i");
        checkvolumes(im)
        mesh!(im, 1.1)
        checkvolumes(im)
        println("Mesh size: final = $(size(im.t,1))")
        V = volumes(im)
        println("length(find(x -> x <= 0.0, V)) = $(length(find(x -> x <= 0.0, V)))")
        # open("im$(i)" * ".jls", "w") do file
        #     serialize(file, im)
        # end
    end
    
    fens = FENodeSet(im.v)
    fes = FESetT4(im.t)
    setlabel!(fes, im.tmid)
    
    # bfes = meshboundary(fes)
    # list = selectelem(fens, fes; overlappingbox = boundingbox([0.2018 2.1537 3.9064]), inflate = 0.01, allin = false)
    # File = "voxel_bracket_mesh_tet.vtk"
    # vtkexportmesh(File, fens, subset(fes, list))
    # @async run(`"paraview.exe" $File`)
    File = "voxel_bracket_mesh_tet.vtk"
    vtkexportmesh(File, fens, fes)
    @async run(`"paraview.exe" $File`)
    
    ne = NASTRANExporter("voxel_bracket_mesh_tet.nas")
    BEGIN_BULK(ne)
    for i = 1:count(fens)
        GRID(ne, i, fens.xyz[i, :])
    end
    for i = 1:count(fes)
        CTETRA(ne, i, 1, [fes.conn[i]...])
    end
    PSOLID(ne, 1, 1)
    MAT1(ne, 1, 20.0e9, 0.3, 1800.0)
    ENDDATA(ne)
    close(ne)
    
    stle = STLExporter("voxel_bracket_mesh_tet.stl")
    solid(stle) 
    bfes = meshboundary(fes)
    for i = 1:count(bfes)
        facet(stle, fens.xyz[bfes.conn[i][1], :], fens.xyz[bfes.conn[i][2], :], fens.xyz[bfes.conn[i][3], :])
    end
    endsolid(stle)
    close(stle)
end # voxel_bracket_mesh_tet_remesh

function allrun()
    println("#####################################################") 
    println("# voxel_bracket_mesh ")
    voxel_bracket_mesh()
    println("#####################################################") 
    println("# voxel_bracket_mesh_tet ")
    voxel_bracket_mesh_tet()
    println("#####################################################") 
    println("# voxel_bracket_mesh_tet_remesh ")
    voxel_bracket_mesh_tet_remesh()
    return true
end # function allrun

end # module voxel_bracket_examples
