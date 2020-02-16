module mrotationmatrixs2
using FinEtools
using LinearAlgebra: norm, I
# using BenchmarkTools
using Test
function rotmat3original!(Rmout, a) 
    na = norm(a);
    thetatilde = zeros(3,3);
    skewmat!(thetatilde, a);
    a = a/na;
    ca = cos(na);
    sa = sin(na);
    aa = a*a';
    copyto!(Rmout, ca * (1.0*I-aa) + sa/na*thetatilde + aa);
    return Rmout
end
function test()
    for i in 1:10
        Rmout = rand(3, 3)
        a = rand(3)
        # Rmout = [0.056575544213388396 0.8235145288079604 0.4566335160952417; 0.6335215149040396 0.8852762577618685 0.803417762628938; 0.5940995203626245 0.35601280936101065 0.6825575362008556]                                                                                                 
        # a = [0.11011521831193871, 0.5097478695998647, 0.8139760429477749]  
        rotmat3!(Rmout, a)
        Rmout2 = rand(3, 3)
        rotmat3original!(Rmout2, a) 
        # Rmout = [0.5736166081143759 -0.6670428981551146 0.47541325067375245; 0.7189364978350996 0.6881251218379588 0.0980516638109532; -0.3925484670406451 
        # 0.2855478746485778 0.8742814834523946]      
        @test norm(Rmout-Rmout2) <= 1.0e-6
        # @btime rotmat3!($Rmout, $a)
    end
    true
end
end
using .mrotationmatrixs2
mrotationmatrixs2.test()