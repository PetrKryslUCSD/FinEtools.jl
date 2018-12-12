module mmblock_Energy_2
using FinEtools
using FinEtools.MeshExportModule: vtkexportmesh, T3, vtkexportvectors
using Test
import LinearAlgebra: cholesky, norm, dot
function test()
  A = 1.0 # dimension of the domain (length of the side of the square)
  thermal_conductivity = [i==j ? one(FFlt) : zero(FFlt) for i=1:2, j=1:2]; # conductivity matrix
  Q = 0.0; # internal heat generation rate
  function getsource!(forceout::FFltVec, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)
    forceout[1] = Q; #heat source
  end
  gradtemp = [2.0, -0.5]
  tempf(x) = (1.0 .+ gradtemp[1] .* x[:,1] .+ gradtemp[2] .* x[:,2]);#the exact distribution of temperature
  N = 10;# number of subdivisions along the sides of the square domain
  Rm=[-0.9917568452513019 -0.12813414805267656
  -0.12813414805267656 0.9917568452513019]
  Rm=[-0.8020689950104449 -0.5972313850116512
  -0.5972313850116512 0.8020689950104447]

  fens,fes = T3block(A, A, N, N)

  geom = NodalField(fens.xyz)
  Temp = NodalField(zeros(size(fens.xyz,1),1))

  l1 = selectnode(fens; box=[0. 0. 0. A], inflate = 1.0/N/100.0)
  l2 = selectnode(fens; box=[A A 0. A], inflate = 1.0/N/100.0)
  l3 = selectnode(fens; box=[0. A 0. 0.], inflate = 1.0/N/100.0)
  l4 = selectnode(fens; box=[0. A A A], inflate = 1.0/N/100.0)
  List = vcat(l1, l2, l3, l4)
  setebc!(Temp, List, true, 1, tempf(geom.values[List,:])[:])
  applyebc!(Temp)
  numberdofs!(Temp)


  material = MatHeatDiff(thermal_conductivity)
  femm = FEMMHeatDiff(IntegDomain(fes, TriRule(1), 100.0), CSys(Rm), material)

  K = conductivity(femm, geom, Temp)
  F2 = nzebcloadsconductivity(femm, geom, Temp);
  fi = ForceIntensity(FFlt[Q]);
  F1 = distribloads(femm, geom, Temp, fi, 3);
  U = K\(F1+F2)
  scattersysvec!(Temp,U[:])

  qenergy = energy(femm, geom, Temp)
  @test abs(qenergy - A * A * 100.0 * dot(gradtemp, -thermal_conductivity * vec(gradtemp))) < 1.0e-9

  # File =  "mmblock_Energy_2-vectors.vtk"
  # vtkexportvectors(File, qplocs, [("heatflux", qpfluxes)])
  # File =  "mmblock_Energy_2.vtk"
  # vtkexportmesh(File, fes.conn, [geom.values 0.0 .* Temp.values], T3; scalars=[("Temperature", Temp.values)])

  true
end
end
using .mmblock_Energy_2
mmblock_Energy_2.test()
