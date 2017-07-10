module mmmmmiscellaneous1mmmmmm
using FinEtools
using Base.Test
function test()
  rho=1.21*1e-9;# mass density
  c =345.0*1000;# millimeters per second
  bulk= c^2*rho;
  Lx=1900.0;# length of the box, millimeters
  Ly=800.0; # length of the box, millimeters

  fens,fes = Q4block(Lx,Ly,3,2); # Mesh
  # show(fes.conn)

  bfes = meshboundary(fes)
  @test bfes.conn == [1 2; 5 1; 2 3; 3 4; 4 8; 9 5; 8 12; 10 9; 11 10; 12 11]
end
end
using mmmmmiscellaneous1mmmmmm
mmmmmiscellaneous1mmmmmm.test()

module mmmmstressconversionm
using FinEtools
using Base.Test
function test()
  symmtens(N) = begin t=rand(N, N); t = (t+t')/2.0; end
  t = symmtens(2)
  v = zeros(3)
  strain2x2tto3v!(v, t)
  to = zeros(2, 2)
  strain3vto2x2t!(to, v)
  @test norm(t-to) < eps(1.0)

  t = symmtens(3)
  v = zeros(6)
  strain3x3tto6v!(v, t)
  to = zeros(3, 3)
  strain6vto3x3t!(to, v)
  @test norm(t-to) < eps(1.0)

  t = symmtens(2)
  v = zeros(3)
  stress2x2to3v!(v, t)
  to = zeros(2, 2)
  stress3vto2x2t!(to, v)
  @test norm(t-to) < eps(1.0)

  v = vec([1. 2. 3.])
  t = zeros(3, 3)
  stress3vto3x3t!(t, v)
  to = [1. 3. 0; 3. 2. 0; 0 0 0]
  @test norm(t-to) < eps(1.0)

  v = vec([1. 2 3 4])
  t = zeros(3, 3)
  stress4vto3x3t!(t, v)
  to = [1. 3 0; 3 2 0; 0 0 4]
  @test norm(t-to) < eps(1.0)

  v = rand(6)
  t = zeros(3, 3)
  stress6vto3x3t!(t, v)
  vo = zeros(6)
  stress3x3tto6v!(vo, t)
  @test norm(v-vo) < eps(1.0)

  v = rand(9)
  t = zeros(3, 3)
  strain9vto3x3t!(t, v)
  t = (t + t')/2.0 # symmetrize
  strain3x3tto9v!(v, t)
  v6 = zeros(6)
  strain9vto6v!(v6, v)
  v9 = zeros(9)
  strain6vto9v!(v9, v6)
  @test norm(v-v9) < eps(1.0)

  v = vec([1. 2 3 4 4 5 5 6 6])
  v6 = zeros(6)
  stress9vto6v!(v6, v)
  v9 = zeros(9)
  stress6vto9v!(v9, v6)
  @test norm(v-v9) < eps(1.0)

end
end
using mmmmstressconversionm
mmmmstressconversionm.test()
