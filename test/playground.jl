module msels2
using FinEtools
using FinEtools.MeshSelectionModule: vselect
using LinearAlgebra: norm, I
# using BenchmarkTools
using Test
function test()
    box = fill(0.0, 6)
    initbox!(box, [0.0,0.0,0.0])
    xyz = [240.0 0.0 0.0; 120.00000000000001 120.0 0.0; 0.0 1.469576158976824e-14 0.0; 119.99999999999997 -120.0 0.0]                                                                                                    
    tolerance = 0.03                                                                                             
    @show clampedn = vselect(xyz; box=box, inflate=tolerance) 
    true
end
end
using .msels2
msels2.test()