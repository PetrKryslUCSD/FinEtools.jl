using FinEtools

A= 1.0
thermal_conductivity =  [i==j ? one(FFlt) : zero(FFlt) for i=1:2, j=1:2]; # conductivity matrix
magn = -6.0; #heat source
truetempf(x)=1.0 + x[1].^2 + 2.0*x[2].^2;
N=20;

println("""
Heat conduction example described by Amuthan A. Ramabathiran
http://www.codeproject.com/Articles/579983/Finite-Element-programming-in-Julia:
Unit square, with known temperature distribution along the boundary,
and uniform heat generation rate inside.  Mesh of regular TRIANGLES,
in a grid of $N x $N edges.
This version uses the FinEtools algorithm module.
"""
)
t0 = time()

fens,fes = T3block(A, A, N, N)


# Define boundary conditions
l1  = selectnode(fens; box=[0. 0. 0. A], inflate = 1.0/N/100.0)
l2  = selectnode(fens; box=[A A 0. A], inflate = 1.0/N/100.0)
l3  = selectnode(fens; box=[0. A 0. 0.], inflate = 1.0/N/100.0)
l4  = selectnode(fens; box=[0. A A A], inflate = 1.0/N/100.0)

essential1 = FDataDict("node_list"=>vcat(l1, l2, l3, l4),
"temperature"=>truetempf);
material = MatHeatDiff(thermal_conductivity)
femm = FEMMHeatDiff(IntegData(fes, TriRule(1)), material)
region1 = FDataDict("femm"=>femm, "Q"=>magn)
# Make model data
modeldata= FDataDict("fens"=> fens,
                 "regions"=>[region1],
                 "essential_bcs"=>[essential1]);


# Call the solver
modeldata = FinEtools.AlgoHeatDiffModule.steadystate(modeldata)

println("Total time elapsed = ",time() - t0,"s")

geom=modeldata["geom"]
Temp=modeldata["temp"]
femm=modeldata["regions"][1]["femm"]
function errfh(loc,val)
    exact = truetempf(loc)
    return ((exact-val)*exact)[1]
end

E = integratefieldfunction(femm, geom, Temp, errfh, 0.0, m=3)
println("Error=$E")

# Postprocessing
# geom=modeldata["geom"]
# Temp=modeldata["temp"]
# MeshExportModule.vtkexportmesh ("a.vtk", fes.conn, [geom.values Temp.values], MeshExportModule.T3; scalars=Temp.values, scalars_name ="Temperature")
