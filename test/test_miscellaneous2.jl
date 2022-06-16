
module mmptpartitioning1m
using FinEtools
using Test
import LinearAlgebra: norm
function test()
    a = 10. # radius of the hole
    nC = 20
    nR = 4
    fens,fes = Q4annulus(a, 1.5*a, nR, nC, 1.9*pi)

    npartitions = 8
    partitioning = pointpartitioning(fens.xyz, npartitions)
    partitionnumbers = unique(partitioning)
    @test norm(sort(partitionnumbers) - sort([1
    3
    2
    5
    6
    7
    8
    4])) < 1.e-5
end
test()
nothing
end



module mmptpartitioning2m
using FinEtools
using Test
import LinearAlgebra: norm
function test()
    H = 100. # strip width
    a = 10. # radius of the hole
    L = 200. # length of the strip
    nL = 15
    nH = 10
    nR = 50
    fens,fes = Q4elliphole(a, a, L/2, H/2, nL, nH, nR)
@test count(fes) == 1250
    npartitions = 4
    partitioning = pointpartitioning(fens.xyz, npartitions)
    partitionnumbers = unique(partitioning)
    @test norm(sort(partitionnumbers) - sort([1
    3
    4
    2])) < 1.e-5
end
test()
nothing
end

module mmptpartitioning3m
using FinEtools
using Test
import LinearAlgebra: norm
function test()
    H = 30. # strip width
    R = 10. # radius of the hole
    L = 20. # length of the strip
    nL = 15
    nH = 10
    nR = 5
    fens,fes = H8block(L, H, R, nL, nH, nR)

    npartitions = 16
    partitioning = pointpartitioning(fens.xyz, npartitions)
    partitionnumbers = unique(partitioning)
    @test norm(sort(partitionnumbers) - sort(1:npartitions)) < 1.e-5

    # for gp = partitionnumbers
    #   groupnodes = findall(k -> k == gp, partitioning)
    #   File =  "partition-nodes-$(gp).vtk"
    #   vtkexportmesh(File, fens, FESetP1(reshape(groupnodes, length(groupnodes), 1)))
    # end
    # File =  "partition-mesh.vtk"
    # vtkexportmesh(File, fens, fes)
    # @async run(`"paraview.exe" $File`)
end
test()
nothing
end
