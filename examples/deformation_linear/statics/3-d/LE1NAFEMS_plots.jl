using Gaston
# using PyCall
using DelimitedFiles

function loadcsv(inputcsv)
    contents = readdlm(inputcsv, ',', Float64, '\n'; header = true)
    return contents
end

function coldata(inputcsv, theset)
    contents = loadcsv(inputcsv)
    return contents[1][:, theset]
end

set(axis="loglog", plotstyle="linespoints", linewidth=2, pointsize = 2, color = "black", xlabel = "Number of nodes", ylabel = "Relative error", grid="on", title = "")
 
inputcsv = "LE1NAFEMS_MSH8_convergence.CSV"
x = coldata(inputcsv, 1)
y = coldata(inputcsv, 3)
plot(x, abs.(y), legend = "MSOE", marker = "edmd")

inputcsv = "LE1NAFEMS_MSH8_convergence.CSV"
x = coldata(inputcsv, 1)
y = coldata(inputcsv, 2)
plot!(x, abs.(y), legend="TBE", marker = "ecircle")

# inputcsv = "LE1NAFEMS_MST10_convergence.CSV"
# x = coldata(inputcsv, 1)
# y = coldata(inputcsv, 3)
# plot(x, abs.(y), legend = "MSOE", marker = "edmd")

# inputcsv = "LE1NAFEMS_MST10_convergence.CSV"
# x = coldata(inputcsv, 1)
# y = coldata(inputcsv, 2)
# plot!(x, abs.(y), legend="TBE", marker = "ecircle")



