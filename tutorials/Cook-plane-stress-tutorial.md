# Cook panel under plane stress

In this example we investigate the well-known benchmark of a tapered panel under plane stress conditions known under the name of Cook.  The problem has been solved many times with a variety of finite element models  and hence the solution is well-known.

The problem is solved in a script.  We begin  by `using` the top-level module `FinEtools` and also the linear-deformation algorithm module. With the algorithm modules, the problem is set up and handed off to an algorithm (in this case linear static solution).  Then for postprocessing another set of algorithms can be invoked.

```julia
using FinEtools
using FinEtools.AlgoDeforLinearModule
```

A few  input parameters are defined: the material parameters.

```julia
E = 1.0;
nu = 1.0/3;
```

The geometry of the tapered panel.

```julia
width = 48.0; height = 44.0; thickness  = 1.0;
free_height  = 16.0;
Mid_edge  = [48.0, 52.0];# Location of tracked  deflection
```

The tapered panel is loaded along the free edge with a unit force, which is here converted to loading per unit area.

```julia
magn = 1.0/free_height/thickness;# Magnitude of applied load
```

For the above input parameters the converged displacement of the tip  of the tapered panel in the direction of the applied shear load is

```julia
convutip = 23.97;
```

The mesh is generated as a rectangular block to begin with, and then the coordinates of the nodes are tweaked into the tapered panel shape. In this case we are using quadratic triangles.

```julia
n = 10; # number of elements per side
fens,fes = T6block(width, height, n, n)
```

Reshape into a trapezoidal panel

```julia
for i = 1:count(fens)
    fens.xyz[i,2] = fens.xyz[i,2]+(fens.xyz[i,1]/width)*(height -fens.xyz[i,2]/height*(height-free_height));
end
```

The  boundary conditions  are applied to selected finite element nodes.   The selection is based on the inclusion in a selection box.   The  selected nodes are then used twice,  to fix the degree of freedom  in the direction 1 and  in the direction 2.

```julia
tolerance = minimum([width, height])/n/1000.;#Geometrical tolerance
```

Clamped edge of the membrane

```julia
l1 = selectnode(fens; box=[0.,0.,-Inf, Inf], inflate = tolerance);
```

The list of nodes  is then used to construct entries  for the essential boundary conditions.  The  data is stored in  dictionaries: `ess1` and `ess2 `.  These dictionaries  are used below to compose the computational model.

```julia
ess1 = FDataDict("displacement"=>  0.0, "component"=> 1, "node_list"=>l1);
ess2 = FDataDict("displacement"=>  0.0, "component"=> 2, "node_list"=>l1);
```

The traction boundary condition is applied to the finite elements on the boundary of the panel. First we generate the three-node curve elements on the boundary.

```julia
boundaryfes =  meshboundary(fes);
```

Then from these finite elements we choose the ones that are inside the box that captures the edge of the geometry to which the traction should be applied.

```julia
Toplist  = selectelem(fens, boundaryfes, box= [width, width, -Inf, Inf ], inflate=  tolerance);
```

To apply the traction we create a finite element model machine (FEMM). For the evaluation of the traction it is sufficient to create a the base FEMM.  It consists of the geometry data `IntegData` (connectivity,  integration rule, evaluation  of the basis functions  and basis function gradients with respect to the parametric coordinates), which in turn is composed of the list of the finite elements and  an appropriate quadrature rule (Gauss rule here).

```julia
el1femm = FEMMBase(IntegData(subset(boundaryfes, Toplist), GaussRule(1, 3), thickness));
```

The traction boundary condition is specified with a constant traction vector and the FEMM that will be used to evaluate  the load vector.

```julia
flux1 = FDataDict("traction_vector"=>[0.0,+magn],
    "femm"=>el1femm
    );
```

We make the dictionary for the region (the interior of the domain).  The FEMM and the material are needed. The geometry data  now is equipped with the  triangular  three-point rule. Note the model-reduction type which is used to dispatch to appropriate specializations of the material routines and the FEMM which needs to execute different code for different reduced-dimension models.

```julia
MR = DeforModelRed2DStress
material = MatDeforElastIso(MR,  0.0, E, nu, 0.0)
region1 = FDataDict("femm"=>FEMMDeforLinear(MR, IntegData(fes, TriRule(3), thickness), material));
```

The model data is a dictionary.   In the present example it consists of the node set, the array of dictionaries for the regions, and arrays of dictionaries for each essential and natural boundary condition.

```julia
modeldata = FDataDict("fens"=>fens,
 "regions"=>[region1],
 "essential_bcs"=>[ess1, ess2],
 "traction_bcs"=>[flux1]
 );
```

When the model data is defined, we simply pass it to the algorithm.

```julia
modeldata = AlgoDeforLinearModule.linearstatics(modeldata);
```

The model data is augmented in the algorithm by the nodal field representing the geometry and the displacement field  computed by solving the system of linear algebraic equations of equilibrium.

```julia
u = modeldata["u"];
geom = modeldata["geom"];
```

The complete information returned from the algorithm  is

```julia
@show modeldata
```

Now we can extract the displacement at the mid-edge node and compare to the converged (reference) value.

```julia
nl = selectnode(fens, box=[Mid_edge[1],Mid_edge[1],Mid_edge[2],Mid_edge[2]],
          inflate=tolerance);
theutip = u.values[nl,:]
println("displacement =$(theutip[2]) as compared to converged $convutip")
```

For postprocessing  we will export a VTK file  with the displacement field (vectors)  and  one scalar field ($\sigma_{xy}$).

```julia
modeldata["postprocessing"] = FDataDict("file"=>"cookstress",
   "quantity"=>:Cauchy, "component"=>:xy);
modeldata = AlgoDeforLinearModule.exportstress(modeldata);
```

The  attribute `"postprocessing"` holds additional data computed and returned by the algorithm:

```julia
@show modeldata["postprocessing"]
```

The exported data can be digested as follows: `modeldata["postprocessing"]["exported"]` is an array of exported items.

```julia
display(modeldata["postprocessing"]["exported"])
```

Each entry of the array is a dictionary:

```julia
display(keys(modeldata["postprocessing"]["exported"][1]))
```

Provided we have  `paraview` in the PATH, we can bring it up  to display the exported data.

```julia
File = modeldata["postprocessing"]["exported"][1]["file"]
@async run(`"paraview.exe" $File`);
```

We can also extract the minimum and maximum value of the shear stress.

```julia
display(modeldata["postprocessing"]["exported"][1]["quantity"])
display(modeldata["postprocessing"]["exported"][1]["component"])
fld = modeldata["postprocessing"]["exported"][1]["field"]
println("$(minimum(fld.values)) $(maximum(fld.values))")
```

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

