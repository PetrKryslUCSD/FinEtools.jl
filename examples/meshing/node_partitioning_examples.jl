module node_partitioning_examples
using FinEtools
using Gaston
using Test
import LinearAlgebra: norm

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

    set(axis="normal", plotstyle="points", linewidth=2, pointsize = 1, color = "black", xlabel = "x", ylabel = "y ", grid="on", title = "")
    
    f = figure()
    for ixxxx = 1:length(partitionnumbers)
        i1 = findall(x -> x == partitionnumbers[ixxxx], partitioning)
        plot!(vec(fens.xyz[i1, 1]), vec(fens.xyz[i1, 2]))
    end
    figure(f)
    
end # node_partitioning_1


function node_partitioning_2()
    a = 10. # radius of the hole
    nC = 20
    nR = 4
    fens,fes = Q4annulus(a, 1.5*a, nR, nC, 2*pi)
    
    npartitions = 8
    partitioning = nodepartitioning(fens, npartitions)
    partitionnumbers = unique(partitioning)
    
    set(axis="normal", plotstyle="points", linewidth=2, pointsize = 2, color = "black", xlabel = "x", ylabel = "y ", grid="on", title = "")
    
    f = figure()
    for ixxxx = 1:length(partitionnumbers)
        i1 = findall(x -> x == partitionnumbers[ixxxx], partitioning)
        plot!(vec(fens.xyz[i1, 1]), vec(fens.xyz[i1, 2]))
    end
    figure(f)
     
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


function node_partitioning_5()
    H = 30. # strip width
    L = 20. # length of the strip
    nL = 15
    nH = 10
    fens,fes = Q4block(L, H, nL, nH)
    reg1 = selectelem(fens, fes; box = [0.0 L/3 0.0 H], inflate = L/nL/1000)
    fes1 = subset(fes, reg1)
    fes2 = subset(fes, setdiff(1:count(fes), reg1))

    # Partitioning of all the nodes
    partitioning = fill(0, count(fens))
    
    # Find the partitioning of the nodes in region 1
    nincluded1 = fill(false, count(fens))
    for i = connectednodes(fes1)
        nincluded1[i] = true
    end
    npartitions = 8
    partitioning1 = nodepartitioning(fens, nincluded1, npartitions)
    npartitions1 = maximum(unique(partitioning1))
    
    # Transfer the partitioning of region 1 into the overall partitioning
    partitioning[nincluded1] = partitioning1[nincluded1]

    # Find the partitioning of the nodes connected to region 2, but not region 1
    nincluded2 = fill(false, count(fens))
    for i = connectednodes(fes2)
        nincluded2[i] = (!nincluded1[i]) && true
    end
    npartitions = 4
    partitioning2 = nodepartitioning(fens, nincluded2, npartitions)
    partitioning2 = partitioning2 .+ npartitions1 # shift by the number of partitions in region 1
    partitionnumbers = unique(vcat(partitioning1, partitioning2))
    npartitions = maximum(partitionnumbers)

    # Transfer the partitioning of region 2 into the overall partitioning
    partitioning[nincluded2] = partitioning2[nincluded2]

    # Visualize partitioning
    for gp = partitionnumbers
      groupnodes = findall(k -> k == gp, partitioning)
      File =  "partition-nodes-$(gp).vtk"
      vtkexportmesh(File, fens, FESetP1(reshape(groupnodes, length(groupnodes), 1)))
    end 
    File =  "mesh-1.vtk"
    vtkexportmesh(File, fens, fes1)
    File =  "mesh-2.vtk"
    vtkexportmesh(File, fens, fes2)
    @async run(`"paraview.exe" $File`)
end

function node_partitioning_6()
    H = 30. # strip width
    L = 20. # length of the strip
    nL = 15
    nH = 10
    fens,fes = Q4block(L, H, nL, nH)
    reg1 = selectelem(fens, fes; box = [0.0 L/3 0.0 H], inflate = L/nL/1000)
    reg2 = selectelem(fens, fes; box = [L/3 L 0.0 H/2], inflate = L/nL/1000)
    reg3 = selectelem(fens, fes; box = [L/3 L H/2 H], inflate = L/nL/1000)
    fesarr = vec([subset(fes, reg1) subset(fes, reg2) subset(fes, reg3)]) 
    
    # Partitioning of all the nodes
    partitioning = nodepartitioning(fens, fesarr, vec([2 4 2]))
    partitionnumbers = unique(partitioning)

    # Visualize partitioning
    for gp = partitionnumbers
      groupnodes = findall(k -> k == gp, partitioning)
      File =  "partition-nodes-$(gp).vtk"
      vtkexportmesh(File, fens, FESetP1(reshape(groupnodes, length(groupnodes), 1)))
    end 
    for i = 1:length(fesarr)
        File =  "mesh-$(i).vtk"
        vtkexportmesh(File, fens, fesarr[i])
    end
    File =  "mesh-1.vtk"
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
    println("#####################################################") 
    println("# node_partitioning_5 ")
    node_partitioning_5()
    println("#####################################################") 
    println("# node_partitioning_6 ")
    node_partitioning_6()
    return true
end # function allrun

end # module node_partitioning_examples
