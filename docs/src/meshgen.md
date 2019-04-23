[Table of contents](https://petrkryslucsd.github.io/FinEtools.jl)

# Mesh generation

## Structured mesh generation

The simplest possible meshes can be generated in the form of one-dimensional, two-dimensional, and three-dimensional blocks. The spacing of the nodes can be either uniform (for instance `Q8block`), or the spacing can be given with an arbitrary distribution (for instance `Q4blockx`). Meshes of tetrahedra can be generated in various orientations of the "diagonals".

More complex meshes can be generated for certain element types: for instance an annulus (`Q4annulus`), quarter of a plate with a hole (`Q4elliphole`), quarter of a sphere (`H8spheren`), layered plate (`H8layeredplatex`).

Hexahedral meshes can also be created by extrusion of  quadrilateral meshes (`H8extrudeQ4`).

## Shaping

Simple meshes  such as blocks can be deformed into geometrically complex shapes, for instance  by tapering  or oother relocation of the nodes. For instance, we can generate a block  and then bend it  into one quarter  of  an annulus as

```julia
fens,fes = Q4block(rex-rin,pi/2,5,20);
for i=1:count(fens)
    r=rin+fens.xyz[i,1]; a=fens.xyz[i,2];
    fens.xyz[i,:]=[r*cos(a) r*sin(a)];
end
```

## Merging

Multiple mesh regions  can be generated and then merged together into a single mesh. Refer to the `MeshModificationModule`. Meshes can be also mirrored.

## Biomedical image mesh generation

The function `H8voximg` can generate a hexahedral mesh from a three-dimensional image (such as a CT scan). The resulting meshes can be smoothed (`meshsmoothing`).

A similar functionality  is also available for tetrahedra (`T4voximg`). A more sophisticated strategy is available in the `VoxelTetMeshingModule` module: the initial mesh can be progressively coarsened and smoothed, resulting  in a much more realistic looking geometry compared  to the initial jagged representation.

## Boundary extraction

Mesh  composed of  any element type can be passed to the function  `meshboundary`, and  the boundary of the mesh is extracted. As an example, the code

```julia
fens,fes = Q4block(rex-rin,pi/2,5,20);
bdryfes = meshboundary(fes);
```

generates a mesh of quadrilaterals in the set `fes`,  and `bdryfes = meshboundary(fes)` finds the boundary elements of the type L2 (line elements with two nodes) and stores them in the  finite element set `bdryfes`.

## Conversion  between element types

For any element shape  (line, triangle,  quadrilateral, hexahedron, tetrahedron) there is  the linear version and the quadratic version. Conversion routines are provided so tthat, for example, mesh can be generated as eight-node  hexahedra  and then converted  to twenty node hexahedra as

```julia
fens,fes = H8toH20(fens,fes)
```

Other conversion routines NEEDS  TO BE WRITTEN

## Refinement

Meshes composed of some  element types can be uniformly refined. For instance, quadrilateral meshes can be refined by bisection with `Q4refine`.
