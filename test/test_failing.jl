
module mmmPoissonmmmt10m2m
using FinEtools
using Base.Test
function test()

  t0 = time()

  A = 1.0 # dimension of the domain (length of the side of the square)
  thermal_conductivity = eye(3, 3); # conductivity matrix
  Q = -6.0; # internal heat generation rate
  function getsource!(forceout::FFltVec, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)
    forceout[1] = Q; #heat source
  end
  tempf(x) = (1.0 + x[:,1].^2 + 2*x[:,2].^2);#the exact distribution of temperature
  N = 10;# number of subdivisions along the sides of the square domain

  # println("Mesh generation")
  fens,fes = T10block(A, A, A, N, N, N)

  # println("""
  # Heat conduction example described by Amuthan A. Ramabathiran
  # http://www.codeproject.com/Articles/579983/Finite-Element-programming-in-Julia:
  # Unit cube, with known temperature distribution along the boundary,
  # and uniform heat generation rate inside.  Mesh of regular quadratic TETRAHEDRA,
  # in a grid of $(N) x $(N) x $(N) edges ($(count(fens)) degrees of freedom).
  # Version: 07/03/2017
  # """
  # )

  geom = NodalField(fens.xyz)
  Temp = NodalField(zeros(size(fens.xyz,1),1))

  # println("Searching nodes  for BC")
  Tolerance = 1.0/N/100.0
  l1 = selectnode(fens; box=[0. 0. 0. A 0. A], inflate = Tolerance)
  l2 = selectnode(fens; box=[A A 0. A 0. A], inflate = Tolerance)
  l3 = selectnode(fens; box=[0. A 0. 0. 0. A], inflate = Tolerance)
  l4 = selectnode(fens; box=[0. A A A 0. A], inflate = Tolerance)
  l5 = selectnode(fens; box=[0. A 0. A 0. 0.], inflate = Tolerance)
  l6 = selectnode(fens; box=[0. A 0. A A A], inflate = Tolerance)
  List = vcat(l1, l2, l3, l4, l5, l6)
  setebc!(Temp, List, true, 1, tempf(geom.values[List,:])[:])
  applyebc!(Temp)
  numberdofs!(Temp)

  # println( "Number of free degrees of freedom: $(Temp.nfreedofs)")
  t1 = time()

  material = MatHeatDiff(thermal_conductivity)

  femm = FEMMHeatDiff(GeoD(fes, TetRule(5), 100.), material)


  # println("Conductivity")
  K = conductivity(femm, geom, Temp)
  # println("Nonzero EBC")
  F2 = nzebcloadsconductivity(femm, geom, Temp);
  # println("Internal heat generation")
  # fi = ForceIntensity(FFlt, getsource!);# alternative  specification
  fi = ForceIntensity(FFlt[Q]);
  F1 = distribloads(femm, geom, Temp, fi, 3);

  # println("Factorization")


  # println("Solution of the factorized system")
  U = K\(F1+F2)
  scattersysvec!(Temp,U[:])
  #
  # println("Total time elapsed = $(time() - t0) [s]")
  # println("Solution time elapsed = $(time() - t1) [s]")

  Error= 0.0
  for k=1:size(fens.xyz,1)
    Error = Error+abs.(Temp.values[k,1]-tempf(reshape(fens.xyz[k,:], (1,3))))
  end
  # println("Error =$Error")
@test abs(Error[1]) < 1.0e-3
end
end
using mmmPoissonmmmt10m2m
mmmPoissonmmmt10m2m.test()
