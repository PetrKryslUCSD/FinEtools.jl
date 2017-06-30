using FinEtools

println("""

Heat conduction example described by Amuthan A. Ramabathiran
http://www.codeproject.com/Articles/579983/Finite-Element-programming-in-Julia:
Unit square, with known temperature distribution along the boundary,
and uniform heat generation rate inside.  Mesh of regular linear TRIANGLES,
in a grid of 1000 x 1000 edges (2M triangles, 1M degrees of freedom).
Version: 05/29/2017
"""
)
t0 = time()

A = 1.0 # dimension of the domain (length of the side of the square)
thermal_conductivity = eye(2,2); # conductivity matrix
Q = -6.0; # internal heat generation rate
tempf(x) = (1.0 + x[:,1].^2 + 2*x[:,2].^2);#the exact distribution of temperature
N = 1000;# number of subdivisions along the sides of the square domain


println("Mesh generation")
@time fens,fes =T3block(A, A, N, N)

geom = NodalField(fens.xyz)
Temp = NodalField(zeros(size(fens.xyz,1),1))

println("Searching nodes  for BC")
@time l1 = selectnode(fens; box=[0. 0. 0. A], inflate = 1.0/N/100.0)
@time l2 = selectnode(fens; box=[A A 0. A], inflate = 1.0/N/100.0)
@time l3 = selectnode(fens; box=[0. A 0. 0.], inflate = 1.0/N/100.0)
@time l4 = selectnode(fens; box=[0. A A A], inflate = 1.0/N/100.0)
List = vcat(l1, l2, l3, l4)
@time setebc!(Temp, List, true, 1, tempf(geom.values[List,:])[:])
@time applyebc!(Temp)
@time numberdofs!(Temp)

t1 = time()

material = MatHeatDiff(thermal_conductivity)

femm = FEMMHeatDiff(GeoD(fes, TriRule(1)), material)


println("Conductivity")
@time K = conductivity(femm, geom, Temp)
println("Nonzero EBC")
@time F2 = nzebcloadsconductivity(femm, geom, Temp);
println("Internal heat generation")
# function getsource!(forceout::FFltVec, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)
#   forceout[1] = Q; #heat source
# end
# fi = ForceIntensity(FFlt, getsource!);# alternative  specification
fi = ForceIntensity(FFlt[Q]);
@time F1 = distribloads(FEMMBase(GeoD(fes, TriRule(1))), geom, Temp, fi, 3);

println("Factorization")
@time K = cholfact(K)
println("Solution of the factorized system")
@time U = K\(F1+F2)
scattersysvec!(Temp,U[:])

println("Total time elapsed = $(time() - t0) [s]")
println("Solution time elapsed = $(time() - t1) [s]")

Error= 0.0
for k=1:size(fens.xyz,1)
  Error = Error+abs.(Temp.values[k,1]-tempf(reshape(fens.xyz[k,:], (1,2))))
end
println("Error =$Error")


# File =  "a.vtk"
# MeshExportModule.vtkexportmesh (File, fes.conn, [geom.values Temp.values], MeshExportModule.T3; scalars=Temp.values, scalars_name ="Temperature")

true
