module scratch2_06102017

using FinEtools
using Compat.Test

function test()
  # println("""
  #         % Vibration modes of unit cube  of almost incompressible material.
  #         %
  #         % Reference: Puso MA, Solberg J (2006) A stabilized nodally integrated
  #         % tetrahedral. International Journal for Numerical Methods in
  #         % Engineering 67: 841-867.""")
  #         t0 = time()


  E = 1*phun("PA");
  nu = 0.499;
  rho = 1*phun("KG/M^3");
  a = 1*phun("M"); b = a; h =  a;
  n1 = 10;# How many element edges per side?
  na =  n1; nb =  n1; nh  = n1;
  neigvs = 20                   # how many eigenvalues
  OmegaShift = (0.01*2*pi)^2;

  MR = DeforModelRed3D
  fens,fes  = H20block(a,b,h, na,nb,nh)

  geom = NodalField(fens.xyz)
  u = NodalField(zeros(size(fens.xyz,1),3)) # displacement field

  numberdofs!(u)

  material=MatDeforElastIso(MR, rho, E, nu, 0.0)

  femm = FEMMDeforLinear(MR, IntegData(fes, GaussRule(3,2)), material)

  K =stiffness(femm, geom, u)
  femm = FEMMDeforLinear(MR, IntegData(fes, GaussRule(3,3)), material)
  M =mass(femm, geom, u)
  d,v,nev,nconv = eigs(K+OmegaShift*M, M; nev=neigvs, which=:SM)
  d = d .- OmegaShift;
  fs = real(sqrt.(complex(d)))/(2*pi)
  println("Eigenvalues: $fs [Hz]")

  # mode = 17
  # scattersysvec!(u, v[:,mode])
  # File =  "unit_cube_modes.vtk"
  # vtkexportmesh(File, fens, fes; vectors=[("mode$mode", u.values)])

  @test abs(fs[7]-0.26259869196259) < 1.0e-5
end
end
using .scratch2_06102017
scratch2_06102017.test()