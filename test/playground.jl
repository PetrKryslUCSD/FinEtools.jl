# module mmPoiss_heatflux_1
# using FinEtools
# using FinEtools.MeshExportModule: vtkexportmesh, T3, vtkexportvectors
# using Test
# import LinearAlgebra: cholesky
# function test()

#   # println("""
#   #
#   # Heat conduction example described by Amuthan A. Ramabathiran
#   # http://www.codeproject.com/Articles/579983/Finite-Element-programming-in-Julia:
#   # Unit square, with known temperature distribution along the boundary,
#   # and uniform heat generation rate inside.  Mesh of regular linear TRIANGLES,
#   # in a grid of 1000 x 1000 edges (2M triangles, 1M degrees of freedom).
#   # Version: 05/29/2017
#   # """
#   # )
#   t0 = time()

#   A = 1.0 # dimension of the domain (length of the side of the square)
#   thermal_conductivity = [i==j ? one(FFlt) : zero(FFlt) for i=1:2, j=1:2]; # conductivity matrix
#   Q = -6.0; # internal heat generation rate
#   function getsource!(forceout::FFltVec, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)
#     forceout[1] = Q; #heat source
#   end
#   tempf(x) = (1.0 .+ x[:,1].^2 .+ 2*x[:,2].^2);#the exact distribution of temperature
#   N = 10;# number of subdivisions along the sides of the square domain


#   # println("Mesh generation")
#   fens,fes = T3block(A, A, N, N)

#   geom = NodalField(fens.xyz)
#   Temp = NodalField(zeros(size(fens.xyz,1),1))

#   # println("Searching nodes  for BC")
#   l1 = selectnode(fens; box=[0. 0. 0. A], inflate = 1.0/N/100.0)
#   l2 = selectnode(fens; box=[A A 0. A], inflate = 1.0/N/100.0)
#   l3 = selectnode(fens; box=[0. A 0. 0.], inflate = 1.0/N/100.0)
#   l4 = selectnode(fens; box=[0. A A A], inflate = 1.0/N/100.0)
#   List = vcat(l1, l2, l3, l4)
#   setebc!(Temp, List, true, 1, tempf(geom.values[List,:])[:])
#   applyebc!(Temp)
#   numberdofs!(Temp)

#   t1 = time()

#   material = MatHeatDiff(thermal_conductivity)

#   femm = FEMMHeatDiff(IntegDomain(fes, TriRule(1), 100.), material)


#   # println("Conductivity")
#   K = conductivity(femm, geom, Temp)
#   # println("Nonzero EBC")
#   F2 = nzebcloadsconductivity(femm, geom, Temp);
#   # println("Internal heat generation")
#   # fi = ForceIntensity(FFlt, getsource!);# alternative  specification
#   fi = ForceIntensity(FFlt[Q]);
#   F1 = distribloads(femm, geom, Temp, fi, 3);

#   # println("Factorization")
#   K = cholesky(K)
#   # println("Solution of the factorized system")
#   U = K\(F1+F2)
#   scattersysvec!(Temp,U[:])

#   # println("Total time elapsed = $(time() - t0) [s]")
#   # println("Solution time elapsed = $(time() - t1) [s]")

#   Error= 0.0
#   for k=1:size(fens.xyz,1)
#     Error = Error.+abs.(Temp.values[k,1].-tempf(reshape(fens.xyz[k,:], (1,2))))
#   end

#   @test Error[1]<1.e-5

#   qplocs = []
#   qpfluxes = []
#   function inspector(idat, i, conn, xe, out, loc)
#     qplocs, qpfluxes = idat
#     push!(qplocs, copy(loc))
#     push!(qpfluxes, copy(out))
#     return (qplocs, qpfluxes)
#   end
#   idat = inspectintegpoints(femm, geom, Temp, collect(1:count(fes)), inspector, (qplocs, qpfluxes), :heatflux)

#   File =  "Poiss_heatflux_1-vectors.vtk"
#   vtkexportvectors(File, qplocs, [("heatflux", qpfluxes)])
#   File =  "Poiss_heatflux_1.vtk"
#   vtkexportmesh(File, fes.conn, [geom.values Temp.values], T3; scalars=[("Temperature", Temp.values)])

#   true
# end
# end
# using .mmPoiss_heatflux_1
# mmPoiss_heatflux_1.test()

