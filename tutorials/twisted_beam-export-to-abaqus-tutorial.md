# Export  to Abaqus

In this example  we show how to export a model  to the finite element software Abaqus.

The model is solved in the FinEtools package, in the example `twisted_beam_algo.jl`.  Here we export the model for execution in Abaqus.

The task begins with defining the input parameters, creating the mesh, identifying the nodes  to which essential boundary conditions are to be applied,  and extracting from the boundary the surface finite elements to which the traction loading at the end of the beam is to be applied.

```julia
using FinEtools
using FinEtools.AlgoDeforLinearModule
using FinEtools.MeshExportModule
```

Define some parameters.

```julia
E = 0.29e8;
nu = 0.22;
W = 1.1;
L = 12.;
t =  0.32;
nl = 2; nt = 1; nw = 1; ref = 2;
p =   1/W/t;
```

 Loading in the Z direction

```julia
loadv = [0;0;p]; dir = 3; uex = 0.005424534868469;# Reference (Harder): 5.424e-3;
```

  Loading in the Y direction

```julia
#loadv = [0;p;0]; dir = 2; uex = 0.001753248285256; # Reference (Harder): 1.754e-3;
tolerance  = t/1000;

fens,fes  = H8block(L,W,t, nl*ref,nw*ref,nt*ref)
```

Reshape into a twisted beam shape.

```julia
for i = 1:count(fens)
    let
        a = fens.xyz[i,1]/L*(pi/2); y = fens.xyz[i,2]-(W/2); z = fens.xyz[i,3]-(t/2);
        fens.xyz[i,:] = [fens.xyz[i,1],y*cos(a)-z*sin(a),y*sin(a)+z*cos(a)];
    end
end
```

Clamped face of the beam: select all the nodes in this cross-section.

```julia
l1  = selectnode(fens; box = [0 0 -100*W 100*W -100*W 100*W], inflate  =  tolerance)
```

Traction on the opposite face

```julia
boundaryfes  =   meshboundary(fes);
Toplist   = selectelem(fens,boundaryfes, box =  [L L -100*W 100*W -100*W 100*W], inflate =   tolerance);
```

The tutorial proper begins here. We create the Abaqus exporter and start writing the .inp file.

```julia
AE = AbaqusExporter("twisted_beam");
HEADING(AE, "Twisted beam example");
```

The  part definition is trivial: all will be defined rather for the instance of the part.

```julia
PART(AE, "part1");
END_PART(AE);
```

The assembly will consist  of a single instance (of the empty part defined above).  The node set will be defined for the instance itself.

```julia
ASSEMBLY(AE, "ASSEM1");
INSTANCE(AE, "INSTNC1", "PART1");
NODE(AE, fens.xyz);
```

We export the finite elements themselves.  Note that the elements  need to have  distinct numbers.  We start numbering the hexahedra at 1. The definition of the element creates simultaneously an element set  which is used below in the section assignment (and the definition of the load).

```julia
ELEMENT(AE, "c3d8rh", "AllElements", 1, connasarray(fes))
```

The traction is applied to surface elements.  Because the elements in the Abaqus model need to have unique numbers, we need to start from an integer  which is  the number of the solid elements plus one.

```julia
ELEMENT(AE, "SFM3D4", "TractionElements", 1+count(fes), connasarray(subset(boundaryfes,Toplist)))
```

The nodes in the clamped cross-section are going to be grouped in the node set `l1`.

```julia
NSET_NSET(AE, "l1", l1)
```

We define a coordinate system  (orientation of the material  coordinate system), in this example it is the global Cartesian coordinate system. The sections are defined for the solid elements of the interior and the surface elements to which the traction is applied, and the assignment to the  elements is by element set (`AllElements` and `TractionElements`). Note that for the solid section we also define reference  to hourglass control named `Hourglassctl`.

```julia
ORIENTATION(AE, "GlobalOrientation", vec([1. 0 0]), vec([0 1. 0]));
SOLID_SECTION(AE, "elasticity", "GlobalOrientation", "AllElements", "Hourglassctl");
SURFACE_SECTION(AE, "TractionElements")
```

This concludes the definition  of the instance  and of the assembly.

```julia
END_INSTANCE(AE);
END_ASSEMBLY(AE);
```

This is the definition of the isotropic elastic material.

```julia
MATERIAL(AE, "elasticity")
ELASTIC(AE, E, nu)
```

The element properties for the interior hexahedra are controlled by the section-control.  In this case we are selecting enhanced hourglass stabilization (much preferable to the default  stiffness stabilization).

```julia
SECTION_CONTROLS(AE, "Hourglassctl", "HOURGLASS=ENHANCED")
```

The static perturbation  analysis step is defined  next.

```julia
STEP_PERTURBATION_STATIC(AE)
```

The boundary conditions are applied directly to the node set `l1`.   Since the node set is defined for the instance, we need to refer to it by the qualified name `ASSEM1.INSTNC1.l1`.

```julia
BOUNDARY(AE, "ASSEM1.INSTNC1.l1", 1)
BOUNDARY(AE, "ASSEM1.INSTNC1.l1", 2)
BOUNDARY(AE, "ASSEM1.INSTNC1.l1", 3)
```

The traction is applied to the surface  quadrilateral elements exported above.

```julia
DLOAD(AE, "ASSEM1.INSTNC1.TractionElements", vec(loadv))
```

Now we have defined  the analysis step and the definition of the model can be concluded.

```julia
END_STEP(AE)
close(AE)
```

As quick check, here is the contents of the  exported model file:

```julia
@show readlines("twisted_beam.inp")
```

What remains is to load the model into Abaqus and execute it as a job.  Alternatively Abaqus can be called on the input file to carry out the analysis at the command line as
```
abaqus job=twisted_beam.inp
```
The output database `twisted_beam.odb` can then be loaded for postprocessing, for instance from the command line as
```
abaqus viewer database=twisted_beam.odb
```

```julia
nothing#-
```

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

