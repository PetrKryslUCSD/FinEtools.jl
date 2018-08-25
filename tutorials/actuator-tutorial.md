# MEMS thermomechanical actuator

The actuator is etched from crystalline silicon, produced layer-by-layer as described in this [book](http://hogwarts.ucsd.edu/~pkrysl/femwabaquspython-book-2018). The termini are attached to contact plates which are part of the substrate, and the actuator is cantilevered from the termini.

It is actuated by thermally-generated strains.  The heat is produced by running electric current through the structure, either through the loop that consists of the inner legs, or through the loop that consists of the outer legs.
When the voltage to generate the current is applied on the termini of the inner legs, the inner legs warm up more than the rest of the structure, and since the inner legs are on a lower level then the outer legs and since they get longer, the actuator bends upwards.  If the voltage is applied to the termini of the outer legs, the outer legs warm up more than the inner legs, and since they get longer and since they are on a higher level than the inner legs the actuator bends downwards.
Given the tiny size, the thermal inertia is very small, and the actuation can be performed at the rate of hundreds of cycles per second. Mechanical inertia can also be ignored, at least in the first approximation.
Finally, we may assume that the silicon material properties do not change very much when the silicon is heated.

In this example  we solve the heat conduction problem for the actuator.

First we bring in the packages that we will use.

```julia
using FinEtools # The finite element toolkit
using LinearAlgebra  # For linear algebra
using PGFPlotsX # For plotting
```

In order to define the structured mesh we need some definitions of distances and numbers of elements in various directions.

```julia
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
```

The actuator has a plane of symmetry, which will be taken advantage of here and hence we will mesh only half of the actual geometry of the actuator. The complete mesh is generated  as a collection of meshes which are glued together. Each individual mesh is generated within a single  hexahedron volume. First we set up the geometry and  the meshing parameters. With the above parameters at hand we generate  the meshes inside the hexahedra, always merging the new mesh with the current one.

```julia
fens,fes =  H8hexahedron([x1 y0 z0; x2 y1 z1],m2,n1,p1);
fens1,fes1  =  H8hexahedron([x1 y1 z0;x2 y2 z1],m2,n2,p1);
fens,fes1,fes2  =  mergemeshes(fens1, fes1, fens, fes, 1.0e6*eps(h));
fes =  cat(fes1,fes2);
fens1,fes1  =  H8hexahedron([x0 y1 z0;x1 y2 z1],m1,n2,p1);
fens,fes1,fes2  =  mergemeshes(fens1, fes1, fens, fes, 1.0e6*eps(h));
fes =  cat(fes1,fes2);
fens1,fes1  =  H8hexahedron([x0 y1 z1;x1 y2 z2], m1,n2,p2);
fens,fes1,fes2  =  mergemeshes(fens1, fes1, fens, fes, 1.0e6*eps(h));
fes =  cat(fes1,fes2);
fens1,fes1  =  H8hexahedron([x0 y1 z2;x1 y2 z3],m1,n2,p3);
fens,fes1,fes2  =  mergemeshes(fens1, fes1, fens, fes, 1.0e6*eps(h));
fes =  cat(fes1,fes2);
fens1,fes1  =  H8hexahedron([x0 y2 z2;x1 y3 z3],m1,n3,p3);
fens,fes1,fes2  =  mergemeshes(fens1, fes1, fens, fes, 1.0e6*eps(h));
fes =  cat(fes1,fes2);
fens1,fes1  =  H8hexahedron([x0 y3 z2;x1 y4 z3], m1,n4,p3);
fens,fes1,fes2  =  mergemeshes(fens1, fes1, fens, fes, 1.0e6*eps(h));
fes =  cat(fes1,fes2);
fens1,fes1  =  H8hexahedron([x1 y3 z2;x3 y4 z3],m4,n4,p3);
fens,fes1,fes2  =  mergemeshes(fens1, fes1, fens, fes, 1.0e6*eps(h));
fes =  cat(fes1,fes2);
fens1,fes1  =  H8hexahedron([x3 y3 z2;x4 y4 z3],m3,n4,p3);
fens,fes1,fes2  =  mergemeshes(fens1, fes1, fens, fes, 1.0e6*eps(h));
fes =  cat(fes1,fes2);
fens1,fes1  =  H8hexahedron([x3 y0 z2;x4 y3 z3], m3,n5,p3);
fens,fes1,fes2  =  mergemeshes(fens1, fes1, fens, fes, 1.0e6*eps(h));
fes =  cat(fes1,fes2);
```

The eight node hexahedra are subsequently converted to the serendipity quadratic elements.

```julia
fens,fes  =  H8toH20(fens,fes);
```

We define  the other parameters of the problem, the thermal conductivity  and the thermal loading driven by the Joule (resistive) heating.

```julia
kappa = 157.0*[1 0 0; 0 1 0; 0 0 1]*phun("W/m/K"); # W/m/K, conductivity matrix
DV = 5.0*phun("V"); # voltage drop in volt
ell  = 2*(y1+y2)/2+2*(x1+x2)/2; # length of the conductor
resistivity  =  1.1e-5*phun("Ohm*m"); # Ohm m
Q = DV^2/resistivity/ell^2; # rate of Joule heating, W/m^3
T_substrate = 293; # substrate temperature in degrees Kelvin
mater = MatHeatDiff(kappa);
```

In the present tutorial we do not use an algorithm to obtain the solution.  All the steps of the solution process are spelled out.

We split the entire geometry of the actuator into  the part that is heated by the current running through the structure, and the rest.  Using a box we select the hexahedral finite elements that are  part of the hot leg of the structure.

```julia
cl =  selectelem(fens, fes, box=[x0,x2,y0,y2,z0,z1],inflate = t/100);
```

The two  FEMMs are then generated (please refer to the documentation for a detailed explanation of what the Finite Element Method Machine, FEMM, is for). One for the "hot" region (with the current running through) and the remaining "cold"  region (only affected by the heat conduction).

```julia
hotfemm  =  FEMMHeatDiff(IntegData(subset(fes,cl), GaussRule(3, 3)), mater)
coldfemm  = FEMMHeatDiff(IntegData(subset(fes,setdiff(collect(1:count(fes)), cl)),  GaussRule(3, 3)), mater);
```

We create the geometry  and the temperature  nodal fields. The nodal fields represent the locations of the nodes and the temperature at the nodes.

```julia
geom = NodalField(fens.xyz)
Temp = NodalField(zeros(size(fens.xyz,1),1));
```

We select the ends  of the actuator legs where the actuator is connected to the substrate, and we apply a fixed temperature condition at these nodes.

```julia
fenids = selectnode(fens, box=[x0,x4,y0,y0,z0,z3],
    inflate=t/1000) ; # fixed temperature on substrate
setebc!(Temp, fenids, true, 1, T_substrate)
applyebc!(Temp);
```

The degrees of freedom of the temperature field are then numbered.

```julia
numberdofs!(Temp);
```

Now we are ready to calculate the conductivity matrix.  Since the interior domain is split into two pieces, the conductivity matrix  is computed separately for each.

```julia
K = conductivity(hotfemm, geom, Temp) + conductivity(coldfemm, geom, Temp);
```

The Joule heating is applied to the heated part of the domain only.

```julia
fi = ForceIntensity(FFlt[Q]);
F = distribloads(hotfemm, geom, Temp, fi, 3);
```

The prescribed temperature condition also generates loads on the free degrees of freedom (the nonzero essential boundary condition loads).

```julia
F  = F + nzebcloadsconductivity(hotfemm, geom, Temp) + nzebcloadsconductivity(coldfemm, geom, Temp);
```

We have  constructed the linear algebra representation for the overall  discrete system. Now we solve for the free degrees of freedom and distribute the solution  into the temperature field.

```julia
U = K\F
scattersysvec!(Temp,U[:]);
```

The solution is presented as a line plot. Two series of nodes along the hot and cold legs  are selected, and their temperature is plotted against the distance along the leg.

```julia
x = range(0; stop = 2*pi, length = 100)
y = sin.(x)
x = range(0; stop = 2*pi, length = 10)
y = sin.(x)


nList = selectnode(fens, box=[x1,x1,y0,y1,z1,z1], inflate=t/1000)
y_i = geom.values[nList, 2]
T_i = Temp.values[nList, 1]
ixi = sortperm(y_i)
nList = selectnode(fens, box=[x3,x3,y0,y3,z2,z2], inflate=t/1000)
y_o = geom.values[nList, 2]
T_o = Temp.values[nList, 1]
ixo = sortperm(y_o)
@pgf p = Axis({
                xlabel = "Distance  [m]",
                ylabel = "Temperature [degree Kelvin]",
                title = "",
                legend_pos  = "south east"
            },
            Plot({
                color = "red",
                mark  = "x"
            }, Table([:x => y_i[ixi], :y => T_i[ixi]])), LegendEntry("hot leg"),
            Plot({
                color = "blue",
                mark  = "o"
            }, Table([:x => y_o[ixo], :y => T_o[ixo]])), LegendEntry("cold leg")
            );
```

PGFPlotsX.print_tex(p)#-
*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

