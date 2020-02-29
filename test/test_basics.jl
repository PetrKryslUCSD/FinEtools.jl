module mbas001
using FinEtools
using Test
function test()
   fens,fes = Q4block(1.3, 3.1,3,2); # Mesh
   @test manifdim(fes) == 2
   bfes = meshboundary(fes)
   @test manifdim(bfes) == 1
   true
end
end
using .mbas001
mbas001.test()

module mbas002
using FinEtools
using Test
function test()
   fens,fes = H8block(1.3, 3.1, 2.7, 3, 2, 4); # Mesh
   @test manifdim(fes) == 3
   bfes = meshboundary(fes)
   @test manifdim(bfes) == 2
   true
end
end
using .mbas002
mbas002.test()

module mbas003
using FinEtools
using Test
function test()
   fens,fes = T4block(1.3, 3.1, 2.7, 3, 2, 4); # Mesh
   @test manifdim(fes) == 3
   bfes = meshboundary(fes)
   @test manifdim(bfes) == 2
   true
end
end
using .mbas003
mbas003.test()


module mbas004
using FinEtools
using Test
function test()
   fens,fes = Q4block(1.3, 3.1,3,2); # Mesh
   @test nodesperelem(fes) == 4
   bfes = meshboundary(fes)
   @test nodesperelem(bfes) == 2
   true
end
end
using .mbas004
mbas004.test()

module mbas005
using FinEtools
using Test
function test()
   fens,fes = H8block(1.3, 3.1, 2.7, 3, 2, 4); # Mesh
   @test nodesperelem(fes) == 8
      bfes = meshboundary(fes)
      @test nodesperelem(bfes) == 4
   true
end
end
using .mbas005
mbas005.test()

module mbas006
using FinEtools
using Test
function test()
   fens,fes = T4block(1.3, 3.1, 2.7, 3, 2, 4); # Mesh
   @test nodesperelem(fes) == 4
      bfes = meshboundary(fes)
      @test nodesperelem(bfes) == 3
   true
end
end
using .mbas006
mbas006.test()


module mbas007
using FinEtools
using Test
function test()
   fens,fes = T10block(1.3, 3.1, 2.7, 3, 2, 4); # Mesh
   @test nodesperelem(fes) == 10
      bfes = meshboundary(fes)
      @test nodesperelem(bfes) == 6
   true
end
end
using .mbas007
mbas007.test()


module mbas008
using FinEtools
using Test
function test()
 fens,fes = L2block(1.37, 3); # Mesh
 @test nodesperelem(fes) == 2
 bfes = meshboundary(fes)
 @test nodesperelem(bfes) == 1
 true
end
end
using .mbas008
mbas008.test()


module mbas009
using FinEtools
using Test
function test()
   fens,fes = Q8block(1.3, 3.1,3,2); # Mesh
   @test nodesperelem(fes) == 8
   bfes = meshboundary(fes)
   @test nodesperelem(bfes) == 3
   true
end
end
using .mbas009
mbas009.test()
