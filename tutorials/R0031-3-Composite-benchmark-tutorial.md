R0031/3 Composite plate test

This is a test recommended by the National Agency for Finite Element Methods and Standards (U.K.): Test R0031/3 from NAFEMS publication R0031, “Composites Benchmarks,” February 1995. It is  a composite  (sandwich) plate of square shape, simply supported along all four edges. Uniform transverse loading is applied to the top skin. The modeled part is one quarter of the full plate here. The serendipity  quadratic hexahedra  are used, with full integration.

The solution  can be compared with the benchmark results  in the Abaqus manual ["Abaqus Benchmarks Guide"](http://130.149.89.49:2080/v6.7/books/bmk/default.htm?startat=ch04s09anf83.html).

```julia
using FinEtools
using FinEtools.AlgoDeforLinearModule # Bringing the algorithm for linear statics
using PGFPlotsX # For plotting
import Statistics: mean
```

The material parameters are specified for an orthotropic material model.  The units are attached using the `phun` function which can take the specification of the units and spits out the numerical multiplier. The skin  material:

```julia
E1s = 1.0e7*phun("psi")
E2s = 0.4e7*phun("psi")
E3s = 0.4e7*phun("psi")
nu12s = 0.3
nu13s = 0.3
nu23s = 0.3
G12s = 0.1875e7*phun("psi")
G13s = 0.1875e7*phun("psi")
G23s = 0.1875e7*phun("psi");
```

The core material:

```julia
E1c = 10.0*phun("psi")
E2c = 10.0*phun("psi")
E3c = 10e4.*phun("psi")
nu12c = 0.
nu13c = 0.
nu23c = 0.
G12c = 10.0*phun("psi")
G13c = 3.0e4*phun("psi")
G23c = 1.2e4*phun("psi");
```

The magnitude  of the distributed uniform transfers loading is

```julia
tmag = 100*phun("psi");
```

Now we generate the mesh.   The sandwich plate volume is divided  into a regular Cartesian grid in the $X$ and $Y$ direction in the plane of the plate, and  in the thickness direction  it is divided  into three layers, with each layer again subdivided into multiple  elements.

```julia
L = 10.0*phun("in") # side of the square plate
nL = 8 # number of elements along the side of the plate
xs = collect(linearspace(0.0, L/2, nL+1))
ys = collect(linearspace(0.0, L/2, nL+1));;
```

The thicknesses are specified from the bottom of the plate: skin, core, and then again skin.

```julia
ts = [0.028; 0.75; 0.028]*phun("in")
nts = [2; 3;  2; ]; # number of elements through the thickness
```

The `H8layeredplatex` meshing function generates the mesh and marks the elements  with a label identifying  the layer to which they belong.  We will use the label to create separate regions, with their own separate materials.

```julia
fens,fes = H8layeredplatex(xs, ys, ts, nts)
```

The linear hexahedra are subsequently converted to serendipity (quadratic) elements.

```julia
fens,fes = H8toH20(fens,fes);
```

The model reduction  here simply says this is a fully three-dimensional model.  The two materials are created.

```julia
MR = DeforModelRed3D
skinmaterial = MatDeforElastOrtho(MR,
  0.0, E1s, E2s, E3s,
  nu12s, nu13s, nu23s,
  G12s, G13s, G23s,
  0.0, 0.0, 0.0)
corematerial = MatDeforElastOrtho(MR,
  0.0, E1c, E2c, E3c,
  nu12c, nu13c, nu23c,
  G12c, G13c, G23c,
  0.0, 0.0, 0.0);
```

Now we are ready to create three material regions:  one for the bottom skin, one for the core, and one for the top skin. The selection of elements assigned to each of the three regions is based on the label. Full Gauss quadrature  is used.

```julia
rl1 = selectelem(fens, fes, label=1)
skinbot = FDataDict("femm"=>FEMMDeforLinear(MR,
    IntegData(subset(fes, rl1), GaussRule(3, 3)), skinmaterial))

rl3 = selectelem(fens, fes, label=3)
skintop = FDataDict("femm"=>FEMMDeforLinear(MR,
    IntegData(subset(fes, rl3), GaussRule(3, 3)), skinmaterial))

rl2 = selectelem(fens, fes, label=2)
core = FDataDict("femm"=>FEMMDeforLinear(MR,
    IntegData(subset(fes, rl2), GaussRule(3, 3)), corematerial));
```

Note that since we did not specify the material coordinate system,  the default is assumed  (which is identical to the global Cartesian coordinate system).

```julia
@show skinbot["femm"].mcsys
```

Next we select the nodes to which  essential boundary conditions  will be applied.  A node is selected  if it is within the specified box  which for the purpose of the test  is inflated in all directions by `tolerance`. The  nodes on the planes of symmetry need to be selected, and also  the nodes  along the edges (faces) to be simply supported  need to be identified.

```julia
tolerance = 0.0001*phun("in")
lx0 = selectnode(fens, box=[0.0 0.0 -Inf Inf -Inf Inf], inflate=tolerance)
lxL2 = selectnode(fens, box=[L/2 L/2 -Inf Inf -Inf Inf], inflate=tolerance)
ly0 = selectnode(fens, box=[-Inf Inf 0.0 0.0 -Inf Inf], inflate=tolerance)
lyL2 = selectnode(fens, box=[-Inf Inf L/2 L/2 -Inf Inf], inflate=tolerance);
```

We have four sides  of the quarter of the plate, two on each plane of symmetry, and two  along the circumference. Hence we create  four essential boundary condition definitions.

```julia
ex0 = FDataDict( "displacement"=>  0.0, "component"=> 3, "node_list"=>lx0 )
exL2 = FDataDict( "displacement"=>  0.0, "component"=> 1, "node_list"=>lxL2 )
ey0 = FDataDict( "displacement"=>  0.0, "component"=> 3, "node_list"=>ly0 )
eyL2 = FDataDict( "displacement"=>  0.0, "component"=> 2, "node_list"=>lyL2 );
```

The traction on the top surface of the top skin is applied to the subset  of the surface mesh of the entire domain. First we find the  boundary mesh, and then from the boundary mesh we select the faces that  "face" upward (along the positive $Z$ axis).

```julia
bfes = meshboundary(fes)
ttopl = selectelem(fens, bfes; facing=true, direction = [0.0 0.0 1.0])
Trac = FDataDict("traction_vector"=>[0.0; 0.0; -tmag],
    "femm"=>FEMMBase(IntegData(subset(bfes, ttopl), GaussRule(2, 3))));
```

The model data  is composed of the  finite element nodes, an array  of the regions, an array of the essential boundary condition definitions, and  an array of  the traction (natural) boundary condition definitions.

```julia
modeldata = FDataDict("fens"=>fens,
 "regions"=>[skinbot, core, skintop], "essential_bcs"=>[ex0, exL2, ey0, eyL2],
 "traction_bcs"=> [Trac]
 );
```

With the model data assembled,  we can now call the algorithm.

```julia
modeldata = AlgoDeforLinearModule.linearstatics(modeldata);
```

The  computed solution can now be postprocessed. The displacement is reported at the center of the plate, along the line in the direction of the loading. We select all the nodes along this line.

```julia
u = modeldata["u"]
geom = modeldata["geom"]
lcenter = selectnode(fens, box=[L/2 L/2  L/2 L/2 -Inf Inf], inflate=tolerance);
```

The variation of the displacement along this line  can be plotted  as (the bottom surface of the shell is at $Z=0$):

```julia
ix = sortperm(geom.values[lcenter, 3])
@pgf p = Axis({
                xlabel = "Z coordinate [in]",
                ylabel = "Vertical displacement [in]",
                title = "",
                legend_pos  = "south east"
            },
            Plot({
                color = "red",
                mark  = "x"
            }, Table([:x => geom.values[lcenter, 3][ix], :y => u.values[lcenter, 3][ix]/phun("in")]))
            )
display(p)
```

A reasonable single number to report for the deflection at the center is the average of the displacements at the nodes at the center of the plate:

```julia
cdis = mean(u.values[lcenter, 3])/phun("in");
println("Center node displacements $(cdis) [in]; NAFEMS-R0031-3 reference: –0.123 [in]")
```

The reference displacement at the center of -0.123 [in] reported for the benchmark is due to an analytical formulation that neglects transverse  (pinching) deformation. Due to the soft core, significant pinching is observed. The solution to the benchmark  obtained in Abaqus  with incompatible hexahedral elements (with the same number of elements as in the stacked continuum shell solution) is -0.131 [in], so close to our own solution.

The deformed shape can be investigated  visually in `paraview` (uncomment the line at the bottom if you have `paraview` in your  PATH):

```julia
File =  "NAFEMS-R0031-3-plate.vtk"
vtkexportmesh(File, connasarray(fes), geom.values, FinEtools.MeshExportModule.H20;
    scalars = [("Layer", fes.label)], vectors = [("displacement", u.values)])
```

@async run(`"paraview.exe" $File`);

Note that the  VTK file will contain element labels (which can help us distinguish between the layers) as scalar field, and the displacements as a vector field.#-
*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

