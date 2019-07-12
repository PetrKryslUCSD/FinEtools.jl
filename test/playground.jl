module forceintensitytest1
using FinEtools
using Test
function test()
   XYZ = reshape([0.0, 0.0], 2, 1)
   tangents = reshape([0.0, 1.0], 2, 1)
   fe_label = 0
   vector = [10.0]
   fi = ForceIntensity(vector)
   v = updateforce!(fi, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)
   @test v == [10.0]
end
end
using .forceintensitytest1
forceintensitytest1.test()

module forceintensitytest2
using FinEtools
using Test
function test()
   XYZ = reshape([0.0, 0.0], 2, 1)
   tangents = reshape([0.0, 1.0], 2, 1)
   fe_label = 0
   setvector!(v, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt; time::FFlt = 0.0) = begin
      v .= [10.0]
   end
   vector = [10.0]
   fi = ForceIntensity(FFlt, length(vector), setvector!)
   v = updateforce!(fi, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)
   @test v == [10.0]
end
end
using .forceintensitytest2
forceintensitytest2.test()

module forceintensitytest3
using FinEtools
using Test
function test()
   XYZ = reshape([0.0, 0.0], 2, 1)
   tangents = reshape([0.0, 1.0], 2, 1)
   fe_label = 0
   setvector!(v, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt; time::FFlt = 0.0) = begin
       return (time < 5.0 ?  v .= [10.0] : v .= [0.0])
   end
   vector = [10.0]
   fi = ForceIntensity(FFlt, length(vector), setvector!, 0.0)
   v = updateforce!(fi, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)
   @test v == [10.0]
   settime!(fi, 6.0)
   v = updateforce!(fi, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)
   @test v == [0.0]
end
end
using .forceintensitytest3
forceintensitytest3.test()

module forceintensitytest4
using FinEtools
using Test
function test()
   XYZ = reshape([0.0, 0.0], 2, 1)
   tangents = reshape([0.0, 1.0], 2, 1)
   fe_label = 0
   setvector!(v, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt) = begin
      v .= [10.0]
      return v
   end
   vector = [10.0]
   fi = ForceIntensity(FFlt, length(vector), setvector!)
   v = updateforce!(fi, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)
   @test v == [10.0]
   settime!(fi, 6.0)
   v = updateforce!(fi, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)
   @test v == [10.0]
end
end
using .forceintensitytest4
forceintensitytest4.test()
