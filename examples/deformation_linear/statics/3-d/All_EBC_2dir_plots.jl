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

set(axis="loglog", plotstyle="linespoints", linewidth=2, pointsize = 2, xrange="[0.002:0.1]", xlabel = "Element size", ylabel = "Approximate error", grid="on", title = "")

# inputcsv = "All_EBC_2dir_T10_default_Stress.CSV"
# x = coldata(inputcsv, 1)
# y = coldata(inputcsv, 4)
# plot(x, abs.(y), legend = "T10", color = "black", marker = "x")

# inputcsv = "All_EBC_2dir_MST10_extrapmean_Stress.CSV"
# x = coldata(inputcsv, 1)
# y = coldata(inputcsv, 4)
# plot!(x, abs.(y), legend="MSOE", color = "black", marker = "edmd")

# inputcsv = "All_EBC_2dir_MST10_extraptrend_Stress.CSV"
# x = coldata(inputcsv, 1)
# y = coldata(inputcsv, 4)
# plot!(x, abs.(y), legend="TBE", color = "black", marker = "ecircle")

f = figure()
inputcsv = "All_EBC_2dir_MSH8_extraptrend_trapezoidal_Stress.CSV"
x = coldata(inputcsv, 1)
y = coldata(inputcsv, 4)
plot!(x, abs.(y), legend = "Trapezoidal", color = "black", marker = "x")

inputcsv = "All_EBC_2dir_MSH8_extraptrend_Stress.CSV"
x = coldata(inputcsv, 1)
y = coldata(inputcsv, 4)
plot!(x, abs.(y), legend = "Gauss", color = "black", marker = "ecircle")

figure(f)
# inputcsv = "All_EBC_2dir_MST10_extrapmean_Stress.CSV"
# x = coldata(inputcsv, 1)
# y = coldata(inputcsv, 4)
# plot!(x, abs.(y), legend="MSOE", color = "black", marker = "edmd")

# inputcsv = "All_EBC_2dir_MST10_extraptrend_Stress.CSV"
# x = coldata(inputcsv, 1)
# y = coldata(inputcsv, 4)
# plot!(x, abs.(y), legend="TBE", color = "black", marker = "ecircle")



# julia> include("All_EBC_2dir_examples.jl"); All_EBC_2dir_examples.allrun()
# WARNING: replacing module All_EBC_2dir_examples.
# #####################################################
# # All_EBC_2dir_MST10_conv
# Fiber-reinforced block: MST10

# Mesh: na, nb, nts = 3, 3, [3]
# count(fens) = 343
# Mesh: na, nb, nts = 6, 6, [6]
# count(fens) = 2197
# Mesh: na, nb, nts = 12, 12, [12]
# count(fens) = 15625
# Mesh: na, nb, nts = 24, 24, [24]
# count(fens) = 117649

# Stress RMS error
# Normalized Approximate Error = [0.0192602, 0.00861774, 0.00354827]
# Linear log-log fit: p = [1.22022, 0.214338]
# Wrote All_EBC_2dir_MST10_extrapmean_Stress.CSV

# Displacement RMS error
# Normalized Approximate Error = [0.00508317, 0.00143038, 0.000373055]
# Linear log-log fit: p = [1.88413, 1.13915]
# Wrote All_EBC_2dir_MST10_extrapmean_Displ.CSV
# Mesh: na, nb, nts = 3, 3, [3]
# count(fens) = 343
# Mesh: na, nb, nts = 6, 6, [6]
# count(fens) = 2197
# Mesh: na, nb, nts = 12, 12, [12]
# count(fens) = 15625
# Mesh: na, nb, nts = 24, 24, [24]
# count(fens) = 117649

# Stress RMS error
# Normalized Approximate Error = [0.0167626, 0.00613563, 0.00220595]
# Linear log-log fit: p = [1.46289, 0.889953]
# Wrote All_EBC_2dir_MST10_extraptrend_Stress.CSV

# Displacement RMS error
# Normalized Approximate Error = [0.00508317, 0.00143038, 0.000373055]
# Linear log-log fit: p = [1.88413, 1.13915]
# Wrote All_EBC_2dir_MST10_extraptrend_Displ.CSV
# Mesh: na, nb, nts = 3, 3, [3]
# count(fens) = 343
# Mesh: na, nb, nts = 6, 6, [6]
# count(fens) = 2197
# Mesh: na, nb, nts = 12, 12, [12]
# count(fens) = 15625
# Mesh: na, nb, nts = 24, 24, [24]
# count(fens) = 117649

# Stress RMS error
# Normalized Approximate Error = [0.017322, 0.00778243, 0.00321]
# Linear log-log fit: p = [1.21598, 0.0942623]
# Wrote All_EBC_2dir_MST10_default_Stress.CSV

# Displacement RMS error
# Normalized Approximate Error = [0.00508317, 0.00143038, 0.000373055]
# Linear log-log fit: p = [1.88413, 1.13915]
# Wrote All_EBC_2dir_MST10_default_Displ.CSV
# Done
# #####################################################
# # All_EBC_2dir_MST10_conv
# Fiber-reinforced block: MST10

# Mesh: na, nb, nts = 3, 3, [3]
# count(fens) = 64
# Mesh: na, nb, nts = 6, 6, [6]
# count(fens) = 343
# Mesh: na, nb, nts = 12, 12, [12]
# count(fens) = 2197
# Mesh: na, nb, nts = 24, 24, [24]
# count(fens) = 15625

# Stress RMS error
# Normalized Approximate Error = [0.0288266, 0.0128894, 0.00580383]
# Linear log-log fit: p = [1.15616, 0.384708]
# Wrote All_EBC_2dir_MST10_extrapmean_Stress.CSV

# Displacement RMS error
# Normalized Approximate Error = [0.00947204, 0.00264193, 0.000695543]
# Linear log-log fit: p = [1.88373, 1.75715]
# Wrote All_EBC_2dir_MST10_extrapmean_Displ.CSV
# Mesh: na, nb, nts = 3, 3, [3]
# count(fens) = 64
# Mesh: na, nb, nts = 6, 6, [6]
# count(fens) = 343
# Mesh: na, nb, nts = 12, 12, [12]
# count(fens) = 2197
# Mesh: na, nb, nts = 24, 24, [24]
# count(fens) = 15625

# Stress RMS error
# Normalized Approximate Error = [0.0253678, 0.0103545, 0.00442978]
# Linear log-log fit: p = [1.25885, 0.599482]
# Wrote All_EBC_2dir_MST10_extraptrend_Stress.CSV

# Displacement RMS error
# Normalized Approximate Error = [0.00947204, 0.00264193, 0.000695543]
# Linear log-log fit: p = [1.88373, 1.75715]
# Wrote All_EBC_2dir_MST10_extraptrend_Displ.CSV
# Mesh: na, nb, nts = 3, 3, [3]
# count(fens) = 64
# Mesh: na, nb, nts = 6, 6, [6]
# count(fens) = 343
# Mesh: na, nb, nts = 12, 12, [12]
# count(fens) = 2197
# Mesh: na, nb, nts = 24, 24, [24]
# count(fens) = 15625

# Stress RMS error
# Normalized Approximate Error = [0.0288266, 0.0128894, 0.00580383]
# Linear log-log fit: p = [1.15616, 0.384708]
# Wrote All_EBC_2dir_MST10_default_Stress.CSV

# Displacement RMS error
# Normalized Approximate Error = [0.00947204, 0.00264193, 0.000695543]
# Linear log-log fit: p = [1.88373, 1.75715]
# Wrote All_EBC_2dir_MST10_default_Displ.CSV
# Done
# #####################################################
# # All_EBC_2dir_T10_conv
# Fiber-reinforced block: T10

# Mesh: na, nb, nts = 3, 3, [3]
# count(fens) = 343
# Mesh: na, nb, nts = 6, 6, [6]
# count(fens) = 2197
# Mesh: na, nb, nts = 12, 12, [12]
# count(fens) = 15625
# Mesh: na, nb, nts = 24, 24, [24]
# count(fens) = 117649

# Stress RMS error
# Normalized Approximate Error = [0.0129395, 0.00559644, 0.00226381]
# Linear log-log fit: p = [1.25748, -0.0593887]
# Wrote All_EBC_2dir_T10_default_Stress.CSV

# Displacement RMS error
# Normalized Approximate Error = [0.00289237, 0.000501187, 8.09591e-5]
# Linear log-log fit: p = [2.57946, 2.93927]
# Wrote All_EBC_2dir_T10_default_Displ.CSV
# Done
# #####################################################
# # All_EBC_2dir_T10_conv
# Fiber-reinforced block: T10

# Mesh: na, nb, nts = 3, 3, [3]
# count(fens) = 64
# Mesh: na, nb, nts = 6, 6, [6]
# count(fens) = 343
# Mesh: na, nb, nts = 12, 12, [12]
# count(fens) = 2197
# Mesh: na, nb, nts = 24, 24, [24]
# count(fens) = 15625

# Stress RMS error
# Normalized Approximate Error = [0.023016, 0.0101339, 0.00453253]
# Linear log-log fit: p = [1.17213, 0.212446]
# Wrote All_EBC_2dir_T10_default_Stress.CSV

# Displacement RMS error
# Normalized Approximate Error = [0.00929339, 0.00262218, 0.00069389]
# Linear log-log fit: p = [1.87171, 1.69831]
# Wrote All_EBC_2dir_T10_default_Displ.CSV
# Done
# true

# julia> include("All_EBC_2dir_plots.jl")
# 1

# julia> include("All_EBC_2dir_plots.jl")
# 1

# julia> include("All_EBC_2dir_plots.jl")
# 1