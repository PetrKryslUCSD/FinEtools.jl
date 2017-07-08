using FinEtools

# MEMS actuator.   Thermal analysis.
x0 =  0.0*phun("micro*m");
x1 = x0+5.0*phun("micro*m");
x2 = x1+10.0*phun("micro*m");
x3 = x2+10.0*phun("micro*m");
x4 = x3+10.0*phun("micro*m");
y0 = 0.0*phun("micro*m");
y4 = 250.0*phun("micro*m");
y3 = y4-10.0*phun("micro*m");
y2 = y3-10.0*phun("micro*m");
y1 = y2-10.0*phun("micro*m");
t = 2.0*phun("micro*m");
h = 0.1*phun("micro*m");
z0 = 0.0*phun("micro*m");
z3 = 2*t+h;
z2 = z3-t;
z1 = z2-h;
m1 = 2*2;
m2 = 2*2;
m3 = 2*2;
m4 = 3*2;
n1 = 20*2;
n2 = 4*2;
n3 = 2*2;
n4 = 2*2;
n5 = 7*2;
p1 = 1*2;
p2 = 1*2;
p3 = 1*2;
kappa = 157*eye(3, 3)*phun("W/m/K"); # W/m/K, conductivity matrix
DV = 5*phun("V"); # voltage drop in volt
l  = 2*(y1+y2)/2+2*(x1+x2)/2; # length of the conductor
resistivity  =  1.1e-5*phun("Ohm*m"); # Ohm m
Q = DV^2/resistivity/l^2; # rate of Joule heating, W/m^3
T_substrate = 293; # substrate temperature in degrees Kelvin

fens,fes =  H8hexahedron([x1 y0 z0; x2 y1 z1],m2,n1,p1);
fens1,fes1  =  H8hexahedron([x1 y1 z0;x2 y2 z1],m2,n2,p1);
fens,fes1,fes2  =  mergemeshes(fens1, fes1, fens, fes, eps(h));
fes =  cat(fes1,fes2);
fens1,fes1  =  H8hexahedron([x0 y1 z0;x1 y2 z1],m1,n2,p1);
fens,fes1,fes2  =  mergemeshes(fens1, fes1, fens, fes, eps(h));
fes =  cat(fes1,fes2);
fens1,fes1  =  H8hexahedron([x0 y1 z1;x1 y2 z2], m1,n2,p2);
fens,fes1,fes2  =  mergemeshes(fens1, fes1, fens, fes, eps(h));
fes =  cat(fes1,fes2);
fens1,fes1  =  H8hexahedron([x0 y1 z2;x1 y2 z3],m1,n2,p3);
fens,fes1,fes2  =  mergemeshes(fens1, fes1, fens, fes, eps(h));
fes =  cat(fes1,fes2);
fens1,fes1  =  H8hexahedron([x0 y2 z2;x1 y3 z3],m1,n3,p3);
fens,fes1,fes2  =  mergemeshes(fens1, fes1, fens, fes, eps(h));
fes =  cat(fes1,fes2);
fens1,fes1  =  H8hexahedron([x0 y3 z2;x1 y4 z3], m1,n4,p3);
fens,fes1,fes2  =  mergemeshes(fens1, fes1, fens, fes, eps(h));
fes =  cat(fes1,fes2);
fens1,fes1  =  H8hexahedron([x1 y3 z2;x3 y4 z3],m4,n4,p3);
fens,fes1,fes2  =  mergemeshes(fens1, fes1, fens, fes, eps(h));
fes =  cat(fes1,fes2);
fens1,fes1  =  H8hexahedron([x3 y3 z2;x4 y4 z3],m3,n4,p3);
fens,fes1,fes2  =  mergemeshes(fens1, fes1, fens, fes, eps(h));
fes =  cat(fes1,fes2);
fens1,fes1  =  H8hexahedron([x3 y0 z2;x4 y3 z3], m3,n5,p3);
fens,fes1,fes2  =  mergemeshes(fens1, fes1, fens, fes, eps(h));
fes =  cat(fes1,fes2);

fens,fes  =  H8toH20(fens,fes);

hotmater = MatHeatDiff(kappa)
coldmater = MatHeatDiff(kappa)
 cl =  selectelem(fens, fes, box=[x0,x2,y0,y2,z0,z1],inflate = t/100);

hotfemm  =  FEMMHeatDiff(GeoD(subset(fes,cl), GaussRule(3, 3), 0.), hotmater)
coldfemm  = FEMMHeatDiff(GeoD(subset(fes,setdiff(collect(1:count(fes)), cl)),
  GaussRule(3, 3), 0.), coldmater)
  geom = NodalField(fens.xyz)
  Temp = NodalField(zeros(size(fens.xyz,1),1))
fenids = selectnode(fens, box=[x0,x4,y0,y0,z0,z3],
    inflate=t/1000) ; # fixed temperature on substrate
setebc!(Temp, fenids, true, 1, T_substrate)
applyebc!(Temp)
numberdofs!(Temp)

K = conductivity(hotfemm, geom, Temp) +
conductivity(coldfemm, geom, Temp)
fi = ForceIntensity(FFlt[Q]);
F = distribloads(hotfemm, geom, Temp, fi, 3) +
nzebcloadsconductivity(hotfemm, geom, Temp) +
nzebcloadsconductivity(coldfemm, geom, Temp)

U = K\F
scattersysvec!(Temp,U[:])

using Plots
plotly()
nList = selectnode(fens, box=[x1,x1,y0,y1,z1,z1], inflate=t/1000)
y_i = geom.values[nList, 2]
T_i = Temp.values[nList, 1]
ix = sortperm(y_i)
plot(y_i[ix], T_i[ix], color=:red, label= "hot leg")
println("maximum(T_i) = $(maximum(T_i))")

nList = selectnode(fens, box=[x3,x3,y0,y3,z2,z2], inflate=t/1000)
y_o = geom.values[nList, 2]
T_o = Temp.values[nList, 1]
ix = sortperm(y_o)
plot!(y_o[ix], T_o[ix], color=:blue, label= "cold leg")

gui()
