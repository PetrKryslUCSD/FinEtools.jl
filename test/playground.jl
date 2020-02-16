module mrotationmatrixs1
using FinEtools
using LinearAlgebra: norm, I
# using BenchmarkTools
using Test
function test()
    for i in 1:10
        Rmout = rand(3, 3)
        a = rand(3)
        # Rmout = [0.056575544213388396 0.8235145288079604 0.4566335160952417; 0.6335215149040396 0.8852762577618685 0.803417762628938; 0.5940995203626245 0.35601280936101065 0.6825575362008556]                                                                                                 
        # a = [0.11011521831193871, 0.5097478695998647, 0.8139760429477749]  
        rotmat3!(Rmout, a)
        # Rmout = [0.5736166081143759 -0.6670428981551146 0.47541325067375245; 0.7189364978350996 0.6881251218379588 0.0980516638109532; -0.3925484670406451 
        # 0.2855478746485778 0.8742814834523946]      
        @test norm(Rmout'*Rmout - 1.0*I) <= 1.0e-6
        # @btime rotmat3!($Rmout, $a)
    end
    true
end
end
using .mrotationmatrixs1
mrotationmatrixs1.test()