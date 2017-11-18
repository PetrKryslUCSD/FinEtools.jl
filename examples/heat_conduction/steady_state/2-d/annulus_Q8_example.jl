using FinEtools

println("""
Annular region,  ingoing and outgoing flux. Minimum/maximum temperature ~(+/-)0.501.
Mesh of quadratic serendipity quadrilaterals.
Version: 05/29/2017
""")

t0 = time()

kappa = 0.2*[1.0 0; 0 1.0]; # conductivity matrix
magn = 0.06;# heat flux along the boundary
rin =  1.0;#internal radius

rex =  2.0; #external radius
nr = 5; nc = 40;
Angle = 2*pi;
thickness =  1.0;
tolerance = min(rin/nr,  rin/nc/2/pi)/10000;


fens, fes  =  Q8annulus(rin, rex, nr, nc, Angle)
fens, fes  =  mergenodes(fens,  fes,  tolerance);
edge_fes  =  meshboundary(fes);

geom = NodalField(fens.xyz)
Temp = NodalField(zeros(size(fens.xyz, 1), 1))


l1  = selectnode(fens; box=[0.0 0.0 -rex -rex],  inflate = tolerance)
setebc!(Temp, l1, 1; val=zero(FFlt))
applyebc!(Temp)

numberdofs!(Temp)


material = MatHeatDiff(kappa)
femm = FEMMHeatDiff(IntegData(fes,  GaussRule(2, 3)),  material)

@time K = conductivity(femm,  geom,  Temp)

l1 = selectelem(fens, edge_fes, box=[-1.1*rex -0.9*rex -0.5*rex 0.5*rex]);
el1femm = FEMMBase(IntegData(subset(edge_fes, l1),  GaussRule(1, 2)))
fi = ForceIntensity(FFlt[-magn]);#entering the domain
@time F1 = (-1.0)* distribloads(el1femm,  geom,  Temp,  fi,  2);

l1 = selectelem(fens, edge_fes, box=[0.9*rex 1.1*rex -0.5*rex 0.5*rex]);
el1femm =  FEMMBase(IntegData(subset(edge_fes, l1),  GaussRule(1, 2)))
fi = ForceIntensity(FFlt[+magn]);#leaving the domain
@time F2 = (-1.0)* distribloads(el1femm,  geom,  Temp,  fi,  2);

@time F3 = nzebcloadsconductivity(femm,  geom,  Temp);


@time K = cholfact(K)
@time U = K\(F1+F2+F3)
@time scattersysvec!(Temp, U[:])

println("Total time elapsed = ", time() - t0, "s")

File =  "annulusq8.vtk"
vtkexportmesh(File,  fes.conn,  [geom.values Temp.values],
 FinEtools.MeshExportModule.Q8; scalars=[("Temperature", Temp.values)])

println("Minimum/maximum temperature= $(minimum(Temp.values))/$(maximum(Temp.values)))")

true
