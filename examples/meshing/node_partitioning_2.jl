using FinEtools

a = 10. # radius of the hole
nC = 20
nR = 4
fens,fes = Q4annulus(a, 1.5*a, nR, nC, 2*pi)

npartitions = 8
partitioning = nodepartitioning(fens, npartitions)
partitionnumbers = unique(partitioning)

using PyCall
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



