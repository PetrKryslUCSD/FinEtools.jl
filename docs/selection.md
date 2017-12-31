[Table of contents](https://petrkryslucsd.github.io/FinEtools.jl)

# Selection of mesh entities

There are many instances of problem definitions where it is important to partition meshes into subsets. As an example,  consider a tube consisting of inner ABS core and  outer fiber-reinforced  laminate  layer. The mesh may consist  of hexahedra.  This mesh would then need to be partitioned into two subsets, because the materials and the  material orientation data  are different between the two regions.

As another example, consider a simple beam  of rectangular cross-section, clamped  at one end,  and  loaded with shear tractions at the  free end. The  entire boundary of the beam needs to be separated  into three subsets:  the first subset,  for the traction-free boundary, is ignored. The second subset, for the clamped cross-section, is extracted  and  its nodes  are used  to  formulate the essential boundary condition. The third subset is extracted and used to define an FEM machine to compute the load vector due to the shear traction.

There are  several  ways  in which mesh entities (nodes and finite elements) can be selected. The simplest uses element labels: some mesh-generation routines label the generated elements. For example,

```julia
fens,fes = H8layeredplatex(xs, ys, ts, nts)
```

generates a plate-like mesh where the layers are labeled. It is therefore possible to select  the bottom-most layer as

```julia
rls = selectelem(fens, fes, label = 1)
```

where `rls` is a list of integer indexes into the  set `fes`, so that we can extract a subset corresponding to this layer as

```julia
botskin = subset(fes, rls)
```

Geometrical techniques for selecting finite elements  or nodes can be based on

- the location within or overlap  with boxes;
- distance from  a given point;
- distance from a given plane;
- connectedness (selection by flooding).

Additionally, surface-like  finite elements (quadrilaterals and triangles embedded in three dimensions, or lines embedded in two dimensions) can be selected based upon the orientation of their normal (`facing`  criterion).

As an example, consider a straight duct with anechoic termination. A triangle mesh is generated as

```julia
fens,fes  =  T3block(Lx,Ly,n,2); 
```

and its boundary is extracted as

```julia
bfes  =  meshboundary(fes)
```

The finite elements from the  piece of the boundary on the left parallel to the Y axis can be extracted as

```julia
L0 = selectelem(fens,bfes,facing = true, direction = [-1.0 0.0])
```

where the numbers of the finite elements  whose normals point in the general direction of the vector [-1.0 0.0] are returned in the integer array `L0`.
