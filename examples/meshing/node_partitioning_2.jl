using FinEtools

a = 10. # radius of the hole
nC = 20
nR = 4
fens,fes = Q4annulus(a, 1.5*a, nR, nC, 2*pi)

npartitions = 8
partitioning = nodepartitioning(fens, npartitions)
partitionnumbers = unique(partitioning)

using Plots
plotly(aspectratio = :equal)
plot(seriestype=:scatter)
for ixxxx = 1:length(partitionnumbers)
    i1 = find(x -> x == partitionnumbers[ixxxx], partitioning)
    plot!(vec(fens.xyz[i1, 1]), vec(fens.xyz[i1, 2]), seriestype=:scatter, m=:o)
end
gui()
