module node_partitioning_examples
using FinEtools
using PyCall

function node_partitioning_1()
    
    H = 100. # strip width
    a = 10. # radius of the hole
    L = 200. # length of the strip
    nL = 15
    nH = 10
    nR = 50
    fens,fes = Q4elliphole(a, a, L/2, H/2, nL, nH, nR)
    
    npartitions = 4
    partitioning = nodepartitioning(fens, npartitions)
    partitionnumbers = unique(partitioning)
    
    @pyimport matplotlib.pyplot as plt
    fig = plt.figure() 
    ax = plt.axes()
    for ixxxx = 1:length(partitionnumbers)
        i1 = find(x -> x == partitionnumbers[ixxxx], partitioning)
        ax[:plot](vec(fens.xyz[i1, 1]), vec(fens.xyz[i1, 2]), linestyle="none", marker=:o)
    end
    plt.show()
    
end # node_partitioning_1


function node_partitioning_2()
    a = 10. # radius of the hole
    nC = 20
    nR = 4
    fens,fes = Q4annulus(a, 1.5*a, nR, nC, 2*pi)
    
    npartitions = 8
    partitioning = nodepartitioning(fens, npartitions)
    partitionnumbers = unique(partitioning)
    
    @pyimport matplotlib.pyplot as plt
    plt.style[:use]("seaborn-whitegrid")
    fig = plt.figure() 
    ax = plt.axes()
    for ixxxx = 1:length(partitionnumbers)
        i1 = find(x -> x == partitionnumbers[ixxxx], partitioning)
        ax[:plot](vec(fens.xyz[i1, 1]), vec(fens.xyz[i1, 2]), linestyle="none", marker=:o)
    end
    ax[:set_xlabel]("x")
    ax[:set_ylabel]("y")
    plt.axis("equal")
    plt.show()
    
    
    
    
end # node_partitioning_2


function node_partitioning_3()
    H = 100. # strip width
    a = 10. # radius of the hole
    L = 200. # length of the strip
    nL = 15
    nH = 10
    nR = 50
    fens,fes = Q4elliphole(a, a, L/2, H/2, nL, nH, nR)
    
    npartitions = 4
    partitioning = nodepartitioning(fens, npartitions)
    partitionnumbers = unique(partitioning)
    
    # Render the  points as dots ("Points" in paraview)
    Points = FESetP1(reshape(collect(1:count(fens)), count(fens), 1))
    File = "Q4elliphole_mesh.vtk"
    vtkexportmesh(File, fens, Points;
    scalars=[("Partition", vec(partitioning))])
    @async run(`"paraview.exe" $File`)
    
end # node_partitioning_3
 
function node_partitioning_4()
    H = 30. # strip width
    R = 10. # radius of the hole
    L = 20. # length of the strip
    nL = 15
    nH = 10
    nR = 5
    fens,fes = H8block(L, H, R, nL, nH, nR)

    npartitions = 16
    partitioning = nodepartitioning(fens, npartitions)
    partitionnumbers = unique(partitioning)
    @test norm(sort(partitionnumbers) - sort(1:npartitions)) < 1.e-5

    # Visualize partitioning
    for gp = partitionnumbers
      groupnodes = findall(k -> k == gp, partitioning)
      File =  "partition-nodes-$(gp).vtk"
      vtkexportmesh(File, fens, FESetP1(reshape(groupnodes, length(groupnodes), 1)))
    end 
    File =  "partition-mesh.vtk"
    vtkexportmesh(File, fens, fes)
    @async run(`"paraview.exe" $File`)
end


function allrun()
    println("#####################################################") 
    println("# node_partitioning_1 ")
    node_partitioning_1()
    println("#####################################################") 
    println("# node_partitioning_2 ")
    node_partitioning_2()
    println("#####################################################") 
    println("# node_partitioning_3 ")
    node_partitioning_3()
    println("#####################################################") 
    println("# node_partitioning_4 ")
    node_partitioning_4()
    return true
end # function allrun

end # module node_partitioning_examples
