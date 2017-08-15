using FinEtools

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

using Plots
plotly()
for ixxxx = 1:length(partitionnumbers)
    i1 = find(x -> x == partitionnumbers[ixxxx], partitioning)
    plot!(vec(fens.xyz[i1, 1]), vec(fens.xyz[i1, 2]), m=:o)
end
gui()
