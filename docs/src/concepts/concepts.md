[Table of contents](https://petrkryslucsd.github.io/FinEtools.jl/latest/index.html)

# Guide

## Break down into modules

The FinEtools package consists of many modules which fall into several  categories. The top-level module, `FinEtools`, includes all other modules and exports functions to constitute the public interface. The user is free to generate their own public interface, however. More details are provided in the section [Make up your own public interface](@ref).


- **Top-level**:     `FinEtools` is the  top-level module.  For interactive use
    it is enough to do `using FinEtools`, however in some  cases functions from
    modules need to be  brought into the scope individually (most importantly,
    the algorithm modules). This is the ONLY  module that EXPORTS  functions,
    none of the other modules exports a single function. The entire  public
    (i. e. exported) interface of the FinEtools package is specified  in the
    file `FinEtools.jl` (i. e. in the `FinEtools` module). The user is free to
    specify his or her own set of exported functions from the FinEtools package
    to create an [ad hoc public interface](rollyourown.html).

- **Utilities**: Refer to the modules `FTypesModule` (definition of basic
    numerical types), `PhysicalUnitModule` (for use with numbers specified
    using physical units), `AssemblyModule` (assembly of elementwise matrices
    and vectors),   `CSysModule` (coordinate system module),
    `MatrixUtilityModule` (utilities for operations on elementwise matrices),
    `BoxModule`  (support for working with bounding boxes),
    `ForceIntensityModule` (force-intensity module), `RotationUtilModule`
    (support for spatial rotations).

- **Mesh  entities**:  `FENodeSetModule`, `FESetModule` (node set and finite element set  types). 

- **Mesh Generation**:   `MeshLineModule`,  `MeshQuadrilateralModule`,
    `MeshTriangleModule`,   `MeshTetrahedronModule`, `MeshHexahedronModule`,
    `VoxelBoxModule`. 

- **Mesh manipulation**:  `MeshSelectionModule` (searching of nodes  and
    elements),  `MeshModificationModule` (mesh boundary, merging  of meshes and
    nodes, smoothing, partitioning),  `MeshUtilModule`
    (utilities), `FENodeToFEMapModule` (search structure from nodes to
    elements).

- **Mesh import/export**:  `MeshImportModule`,  `MeshExportModule`.

- **Fields**:   `FieldModule`,    `GeneralFieldModule`, `ElementalFieldModule`,
    `NodalFieldModule` (modules for representing quantities on the mesh).

- **Integration**: Support for  integration over solids, surfaces, curves, and
    points: `IntegRuleModule`,   `IntegDomainModule`.
    The package defines some common bilinear and linear forms to aid in constructing weighted residual methods.

- **General algorithms**: `AlgoBaseModule` (algorithms), `FEMMBaseModule`
    (FEM machine for general tasks).


## Arithmetic types

The FinEtools package tries to make typing arguments easier. The arithmetic
types used throughout are `FInt` for integer data, `FFlt` for floating-point
data, and `Complex{FFlt}` for applications that work with complex linear
algebra quantities.

The module `FTypesModule` defines these types, and also defines abbreviations for vectors and matrices with entries of these types.

Some algorithms expect input in the form of a *data dictionary*, `FDataDict`, and also produce output in this form.


## Physical units

The `PhysicalUnitModule` provides a simple function, `phun`, which can help with providing input numbers with the correct conversion between physical units. For instance, it is possible to specify the input data as

```julia
E = 200*phun("GPa");# Young's modulus
nu = 0.3;# Poisson ratio
rho = 8000*phun("KG*M^-3");# mass density
L = 10.0*phun("M"); # side of the square plate
t = 0.05*phun("M"); # thickness of the square plate
```

A few common sets of units are included, `:US`, `:IMPERIAL`, `:CGS`, `:SIMM` (millimeter-based SI units), and `:SI` (meter-based SI units). The resulting  values assigned to the variables are floating-point numbers, for instance

```julia
julia> E = 200*phun("GPa")
2.0e11
```

Numbers output by the simulation can also be converted  to appropriate units for printing as

```julia
julia> E/phun("MPa")
200000.0
```

It is also possible to use a macro to define physical units:

```julia
E = 200*u"GPa";# Young's modulus
nu = 0.3;# Poisson ratio
rho = 8000*u"KG*M^-3";# mass density
L = 10.0*u"M"; # side of the square plate
t = 0.05*u"M"; # thickness of the square plate
t / u"mm"
```


## Mesh entities

The mesh consists of one set of finite element nodes  and one or more sets of finite elements.

One of the  organizing principles of the  finite element collection  is that  finite elements can appear as representations of the interior  of the domain, but in a different model as parts of the boundary.  Thus  for instance  4-node  quadrilaterals  are finite elements that represent cross-sections of  axially symmetric models or surfaces  of membranes,  but they are also the boundaries of hexahedral  models.

A mesh  is generated by one of the functions specialized to a particular finite element type. Thus there are  [mesh generation](meshgen.md) functions for lines, triangles, quadrilaterals, tetrahedra, and hexahedra.

## Mesh generation

As an example,  the following code generates a hexahedral mesh of simple rectangular block.

```julia
fens, fes  = H8block(h, l, 2.0 * pi, nh, nl, nc)
```

The finite element node set and the finite element set  are returned. More complicated meshes can be constructed from such mesh parts. There are  functions for  merging  nodes  and even multiple meshes together.

The code snippet below  constructs the mesh of an  L-shaped  domain  from  the meshes of three rectangles.

```julia
W = 100. # width of the leg
L = 200. # length of the leg
nL = 15 # number of elements along the length of the leg
nW = 10 # number of elements along the width
tolerance = W / nW / 1.0e5 # tolerance for merging nodes
Meshes = Array{Tuple{FENodeSet,FESet},1}()
push!(Meshes, Q4quadrilateral([0.0 0.0; W W], nW, nW))
push!(Meshes, Q4quadrilateral([-L 0.0; 0.0 W], nL, nW))
push!(Meshes, Q4quadrilateral([0.0 -L; W 0.0], nW, nL))
fens, outputfes = mergenmeshes(Meshes, tolerance);
fes = cat(outputfes[1], cat(outputfes[2], outputfes[3]))
```

As an example of the  merging  of nodes to create  the final mesh, consider the creation of  a closed hollow tube.

```julia
fens, fes  = H8block(h, l, 2.0 * pi, nh, nl, nc) # generate a block
# Shape into a cylinder
R = zeros(3, 3)
for i = 1:count(fens)
    x, y, z = fens.xyz[i,:];
    rotmat3!(R, [0, z, 0])
    Q = [cos(psi * pi / 180) sin(psi * pi / 180) 0;
        -sin(psi * pi / 180) cos(psi * pi / 180) 0;
        0 0 1]
    fens.xyz[i,:] = reshape([x + Rmed - h / 2, y - l / 2, 0], 1, 3) * Q * R;
end
# Merge the nodes where the tube  closes up
candidates = selectnode(fens, box = boundingbox([Rmed - h -Inf 0.0; Rmed + h +Inf 0.0]), inflate = tolerance)
fens, fes = mergenodes(fens, fes,  tolerance, candidates);
```

The final mesh used for a simulation  consists of a *single  node set* and *one or more finite element sets*. The  finite elements may be  divided into separate sets  to  accommodate different material properties, different orientations of the material  coordinate systems, or different formulations  of the discrete model. The  assignment  of the finite elements to sets  may be based on geometrical proximity, topological connections, or some other characteristic. See the  "mesh selection" [discussion](selection.html) for details.

### Structured mesh generation

The simplest possible meshes can be generated in the form of one-dimensional, two-dimensional, and three-dimensional blocks. The spacing of the nodes can be either uniform (for instance `Q8block`), or the spacing can be given with an arbitrary distribution (for instance `Q4blockx`). Meshes of tetrahedra can be generated in various orientations of the "diagonals".

More complex meshes can be generated for certain element types: for instance an annulus (`Q4annulus`), quarter of a plate with a hole (`Q4elliphole`), quarter of a sphere (`H8spheren`), layered plate (`H8layeredplatex`).

Hexahedral meshes can also be created by extrusion of  quadrilateral meshes (`H8extrudeQ4`).

### Shaping

Simple meshes  such as blocks can be deformed into geometrically complex shapes, for instance  by tapering  or other relocation of the nodes. For instance, we can generate a block  and then bend it  into one quarter  of  an annulus as

```julia
fens,fes = Q4block(rex-rin,pi/2,5,20);
for i=1:count(fens)
    r=rin+fens.xyz[i,1]; a=fens.xyz[i,2];
    fens.xyz[i,:]=[r*cos(a) r*sin(a)];
end
```

### Merging

Multiple mesh regions  can be generated and then merged together into a single mesh. Refer to the `MeshModificationModule`. Meshes can be also mirrored.

### Boundary extraction

Mesh  composed of  any element type can be passed to the function  `meshboundary`, and  the boundary of the mesh is extracted. As an example, the code

```julia
fens,fes = Q4block(rex-rin,pi/2,5,20);
bdryfes = meshboundary(fes);
```

generates a mesh of quadrilaterals in the set `fes`,  and `bdryfes = meshboundary(fes)` finds the boundary elements of the type L2 (line elements with two nodes) and stores them in the  finite element set `bdryfes`.

### Conversion  between element types

For any element shape  (line, triangle,  quadrilateral, hexahedron, tetrahedron) there is  the linear version and the quadratic version. Conversion routines are provided so that, for example, mesh can be generated as eight-node  hexahedra  and then converted  to twenty-node hexahedra as

```julia
fens, fes = H8toH20(fens, fes)
```

Other conversion routines can convert triangles to quadrilaterals, tetrahedra to hexahedra, and so on.

### Refinement

Meshes composed of some  element types can be uniformly refined. For instance, quadrilateral meshes can be refined by bisection with `Q4refine`.

## Selection of mesh entities

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

Additionally, surface-like  finite elements (quadrilaterals and triangles embedded in three dimensions), or lines embedded in two dimensions, can be selected based upon the orientation of their normal (`facing`  criterion).

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

## Fields

The structure to maintain the numbering  and values of the degrees of freedom in the mesh  is the field. Consider for instance the temperature field: we write
```math
T(x) = \sum_i N_i(x) T_i
```
The understanding is that $T_i$ are the degrees of freedom, and the basis functions $N_i(x)$ are defined implicitly by the finite element mesh. (More about basis functions below.) Each element has its own set of functions, which when multiplied by the degree of freedom values describe the temperature over each individual finite element. The basis functions are implicitly associated with the nodes of the finite elements. The degrees of freedom are also (explicitly) associated with the nodes. The field may also be generalized a bit by extending the above sum simply to entities of the mesh, not only the nodes, but perhaps also the elements.

The role of the field is then to maintain the correspondence between the entities and the numbers and values of the degrees of freedom.

### Abstract  Field

The assumption is that a field has one set of degrees of freedom per node or per element. For simplicity we will refer to the nodes and elements as entities.
It assumes that concrete  subtypes of the abstract field  have the following data, one row per entity:

- `values::FMat{T}`: Array of the values of the degrees of freedom, one row  for each entity. All the arrays below have the same dimensions as this one.
- `dofnums::FIntMat`: Array  of the numbers of the free degrees of freedom. First the free degrees of freedom are numbered, then the fixed (prescribed) degrees of freedom.
- `is_fixed::Matrix{Bool}`: Array of  Boolean flags,  `true` for fixed  (prescribed) degrees of freedom, `false` otherwise.
- `fixed_values::FMat{T}`: Array  of the same size and type  as  `values`. Its entries are only relevant  for the fixed (prescribed)  degrees of freedom.
- `nalldofs::FInt`:  the total number of all degrees of freedom.
- `nfreedofs::FInt`:  the total number of free degrees of freedom.

The methods defined for the abstract field  include:

- Return the number of degrees of freedom and the number of entities.

- Gather and scatter the system vector.

- Gather elementwise  vectors or matrices of values, the degree of freedom numbers, or the fixed values of the degrees of freedom.

- Set  or clear essential boundary conditions.

- Copy a field. Clear the entries of the field.

### Nodal Field

In this case  the  abstract field is subtyped to a concrete field where the entities are nodes.

### Elemental Field

In this case  the  abstract field  is subtyped to a concrete field where the entities are the elements.


### General Field

In this case  the  abstract field  is subtyped to a concrete field where the entities are  use-case  specific.

### Numbering of the degrees of freedom

The simplest method is at the moment implemented: number all free degrees of freedom, row-by-row and column-by-column,
starting from 1 up to `nfreedofs(f)`, for the field `f`. 
Then number the prescribed degrees of freedom are numbered, up to `nalldofs(f)`.

There is also a method to supply the numbering of the nodes, perhaps  resulting from the Reverse Cuthill-McKee permutation. This may be useful when using LU or LDLT factorization as the fill-in may be minimized.

## Finite element

The  finite element set is one of the basic entities in `FinEtools`.
It is a homogeneous collection of  finite elements defined by the connectivity (collection of node numbers, listing the nodes connected by the element in  a specific order). The finite element set  provides  specialized methods  to compute the values of basis functions and the values of  the gradients of the basis functions  with respect to the parametric coordinates.

### Element types

The finite element sets are instances of concrete types. Each particular shape and order of element has its own type. There are types for  linear  and quadratic quadrilaterals, for instance, `FESetQ4` and `FESetQ8`. Each element set provides access to the number of nodes  connected by the element (`nodesperelem`),  the connectivity as the two dimensional array    `conn`,  and the  integer label vector `label`.

The concrete finite element set types are subtypes of the abstract type for elements of different manifold dimension (3, 2, 1, and 0), for instance for the quadrilaterals that would be `AbstractFESet2Manifold`. These types are in turn  subtypes of the abstract finite element set type `AbstractFESet`.

The concrete finite element set type provides specialized methods to compute the values of the basis functions, `bfun`, and methods to compute  the gradients of the basis functions with respect to the parametric coordinates, `bfundpar`. `FinEtools` at the moment supports only the so-called **nodal** basis functions: each basis function is associated with a node. And that is  true both globally (in the sense that each basis function is globally supported),  and locally over each finite element, and all such functions are  1 at its own node, and zero at all the other nodes.

### Finite element set functions

+ Methods defined for  the abstract type:

    - `nodesperelem`: Get the number of nodes  connected  by  the finite element.
    - `count`:  Get the number of individual connectivities in the FE set.
    - `setlabel!`: Set the label of the entire finite elements set.
    - `connasarray`: Retrieve  connectivity  as an integer array.
    - `fromarray!`: Set  connectivity from an integer array.
    - `subset`: Extract a subset of the finite elements from the given finite element set.
    - `cat`: Concatenate the connectivities of two FE sets.
    - `updateconn!`: Update the connectivity after the IDs of nodes changed.
    - `map2parametric`: Map a spatial location to parametric coordinates.

+ Methods dispatched based on the manifold type:

    - `manifdim`: Return the manifold dimension.
    - `Jacobian`: Evaluate the  Jacobian.
    - `gradN!`: Compute the gradient of the basis functions with the respect to the "reduced" spatial coordinates.

+ Methods dispatched on the concrete type:

    - `boundaryconn`: Get boundary connectivity.
    - `boundaryfe`: Return the constructor of the type of the boundary finite element.
    - `bfun`: Compute the values of the basis functions at a given parametric coordinate.
    - `bfundpar`: Compute the values of the basis function gradients at a given parametric coordinate.
    - `inparametric`: Are given parametric coordinates inside the element parametric domain?
    - `centroidparametric`: Return the parametric coordinates  of the centroid of the element.

## Integration

There are two kinds of integrals in the weighted-residual finite element method: integrals over the **interior**  of the domain,  and integrals over the **boundary** of the domain.

Consequently, in a typical simulation one would need  two meshes: one for the interior  of the domain,  and one for the boundary. Obviously, the mesh for the boundary will be derived from the mesh  constructed for the interior.

Often only a part  of the entire boundary   is   used:  on some parts of the boundary  the  boundary condition is implied as homogeneous (i. e. zero). For instance, a traction-free boundary. Therefore the necessary integrals are typically evaluated over a subset of the entire boundary.

### Manifold dimension

Finite elements  have  a certain manifold dimension.  Tetrahedra  and hexahedra are three-manifolds, triangles and quadrilaterals are two-manifolds, triangles and quadrilaterals are two-manifolds, lines are one-manifolds, and points are zero-manifolds.

Elements are equipped with an "other" dimension attribute which boosts the manifold dimension to produce the required dimension for  the integration. For instance,  a line element can be equipped with an "other" dimension to represent a cross-section so that a volume integral can be evaluated over a line element. Or, a line element can be given an "other" dimension as a thickness to result in a physical dimension needed to evaluate a surface integral.

The "other"  dimension  has the following meaning  for finite elements of different manifold dimensions:

| Manifold dimension        | Volume integral           | Surface integral  |
| ------------- | ------------- | ----- |
| 3     | NA | NA |
| 2    | Thickness  |  NA |
| 1 | Cross-section   |  Thickness |
| 0 | Volume   |  Cross-section |

### Integration  over the interior

The integrals are always  *volume* integrals. This means that for elements which are of  lower manifold dimension than three the "other"  dimension needs to compensate.

For  three-manifold finite elements (tetrahedra and hexahedra) the "other" dimension is always 1.0. This really means there is no "other" dimension to a volume-like element.

For  finite elements of manifold dimension  less than tthree, the  "other" dimension varies according to the model (axially symmetric versus simple  plane 2D) as shown  in the table below.

| Manifold dimension        | Axially symmetric    | Plane 2D |
| ------------- | ------------- | ----- |
| 2 | ``2\pi r``  |  Thickness |
| 1 | ``2\pi r\times``  Thickness  |  Cross-section |
| 0 | ``2\pi r\times`` Cross-section  |  Volume |  

The integral  is approximated with numerical quadrature as
```math
\int_{\Omega} f dV \approx \sum_q f(\xi_q) J(\xi_q) W_q
```    

Here ``f``  is the integrand, ``f(\xi_q)``  is the  value of the integrand  at
the quadrature point, ``J(\xi_q)``  is the  value of the Jacobian  at the
quadrature point. Importantly, the Jacobian incorporates the "other" dimension,
and therefore it is the  *volume* Jacobian. (For the interior integrals the
Jacobian  is computed by the `Jacobianvolume` method.)

### Integration  over the boundary

The integrals are always  *surface* integrals. This means that for elements
which are of  lower manifold dimension than two the "other"  dimension needs to
compensate.

For  two-manifold finite elements (triangles and quadrilaterals) the "other"
dimension is always 1.0. This really means there is no "other" dimension to a
surface-like element.

For  finite elements of manifold dimension  less than two, the  "other"
dimension varies according to the model (axially symmetric versus simple  plane
2D) as shown  in the table below.

| Manifold dimension        | Axially symmetric    | Plane 2D |
| ------------- | ------------- | ----- |
| 1 | ``2\pi r``   |  Thickness |
| 0 | ``2\pi r\times``   Thickness  |  Cross-section |

The integral  is approximated with numerical quadrature as
```math
\int_{\partial \Omega} f dS \approx \sum_q f(\xi_q) J(\xi_q) W_q  
```    

Here ``f``  is the integrand, ``f(\xi_q)`` is the  value of the integrand  at
the quadrature point, ``J(\xi_q)`` is the  value of the Jacobian  at the
quadrature point. Importantly, the Jacobian incorporates the "other" dimension,
and therefore it is the  *surface* Jacobian. (For the boundary integrals the
Jacobian  is computed by the `Jacobiansurface` method.)

#### Example: axially symmetric model, line element L2

The surface Jacobian in this case  is  equal to the curve Jacobian times `2*pi*r`.

### Integration Domain

As explained above, integrating over the interior or the boundary may mean different things based on the features of the solution domain: axially symmetric?, plane strain or plane stress?, and so forth.

The  module `IntegDomainModule` supports  the processing of  the geometry
necessary for the evaluation of the various integrals. The module data
structure  groups together  a finite element set with an appropriate
integration rule, information about the model (axially symmetric or not), and a
callback to evaluate  the "other" dimension.

### Other dimension

The  discussion of the surface and volume integrals introduces the notion  of the  "other"  dimension. In order to evaluate Jacobians of various space dimensions  the  Geometry Data module takes into account  whether or not the model is axially symmetric, and evaluates the "other" dimension based upon this information.

A finite element set is equipped with  a way of  calculating  the "other" dimension.  For instance, the line element with two nodes, L2, can be given  the "other" dimension  as a  "thickness"  so that  surface integrals  can be evaluated over the line element. However, if this line element  is used in an axially symmetric model, the same  "other" dimension  of "thickness"  will result in the integral  along the length of this line element  being a volume integral.

Thus, the way in which the "other"  dimension gets used by the integration domain methods depends on the model. As an example, consider  the  method
```julia
function Jacobianvolume(self::IntegDomain{T}, J::FFltMat, loc::FFltMat, conn::CC, N::FFltMat)::FFlt where {T<:AbstractFESet2Manifold, CC<:AbstractArray{FInt}}
    Jac = Jacobiansurface(self, J, loc, conn, N)::FFlt
    if self.axisymmetric
        return Jac*2*pi*loc[1];
    else
        return Jac*self.otherdimension(loc, conn,  N)
    end
end
```
which  evaluates the volume Jacobian  for an element  of manifold dimension  2  (surface). Note that  first  the surface Jacobian  is calculated, which is then boosted to a volume Jacobian in two different ways, depending on whether  the model is axially symmetric or not. For the axially symmetric case  the "other"  dimension is implied,

The callback function computes the "other" dimension from  two kinds of  information: (a) the physical location  of the quadrature point,  and (b) the interpolation data for the element  (connectivity and the values of the basis functions at the quadrature point).

- The approach ad (a) is suitable  when the "other" dimension is given as a function of the physical coordinates. The  simplest case is obviously  a uniform distribution of the "other" dimension. When  no  callback is explicitly provided,  the  "other"  dimension  callback is  automatically generated as the trivial
```julia
function otherdimensionunity(loc::FFltMat, conn::CC, N::FFltMat)::FFlt where {CC<:AbstractArray{FInt}}
    return 1.0
end
```
which simply returns 1.0 as the default value.

- The approach ad (b) is appropriate when the "other" dimension  is given by values given at the nodes of the  mesh. Than the connectivity  and  the array of the values of the basis functions  can be used to interpolate the "other"  dimension  to the quadrature point.



### Evaluation of integration data

Importantly, the  Integration Domain (`IntegDomain`) method `integrationdata` evaluates quantities  needed for numerical integration: locations and weights of quadrature points, and the values of basis functions and of the basis function gradients with respect to the parametric coordinates at the quadrature points.

## FEM machines

The construction of the matrices and vectors of the *discrete* form of the weighted residual equation is performed in FEM  machines. (FEM = Finite Element Method.)

As an example consider the weighted-residual form of the heat balance equation
```math
\int_{V}  \vartheta c_V\frac{\partial T}{\partial t} \; \mathrm{d} V
            +\int_{V}(\mathrm{grad}\vartheta)\; \kappa (\mathrm{grad}T
            )^T\; \mathrm{d} V
            -\int_{V}  \vartheta Q \; \mathrm{d} V  
            +\int_{S_2} \vartheta\;\overline{q}_{n}\; \mathrm{d} S+ \int_{S_3} \vartheta\;h
            (T-T_a)  \; \mathrm{d} S = 0
```    

where ``\vartheta(x) =0``  for  ``x \in{S_1}`` .

The  test function is  taken to be  one  finite element basis function at a time, ``\vartheta = N_{j}``, and the trial function is
```math
T = \sum_{i= 1} ^{N} N_{i} T_i .
```    

Here by ``N_{j}`` we mean the basis function constructed on the mesh and associated with the node where the degree of freedom ``j`` is situated. 

Now the test function and the trial function is substituted  into the  weighted residual equation.  

### Example:  internal heat generation rate term

For instance,  for the term
```math
\int_{V}  \vartheta Q \; \mathrm{d} V
```    

we obtain
```math
\int_{V} N_{j} Q \; \mathrm{d} V
```    

This integral evaluates to a number, the heat load  applied to the degree of freedom ``j``. When these numbers are evaluated for all  the free degrees of freedom,  they constitute the entries of the global heat load vector.


Evaluating integrals of this form is so common that there is a module `FEMMBaseModule` with the method `distribloads` that computes and assembles the global vector. For instance to evaluate this heat load vector  on the mesh composed of three-node triangles, for a uniform heat generation rate `Q`, we can write

```julia
fi = ForceIntensity(FFlt[Q]);
F1 = distribloads(FEMMBase(IntegDomain(fes, TriRule(1))), geom, tempr, fi, 3);
```

`IntegDomain(fes, TriRule(1))` constructs integration domain for the  finite elements `fes` using a triangular  integration rule with a single point. `FEMMBase` is the base  FEM  machine,  and all it needs at this point is the integration domain. The method  `distribloads` is defined for the  base FEM machine, the geometry field `geom`, the numbering of the degrees of freedom is taken from the field `tempr`, the internal heat generation rate is defined as the force intensity `fi`, and the integrals  are volume integrals  (3).

### Example: conductivity term

The conductivity term from the weighted residual equation
```math
\int_{V}(\mathrm{grad}\vartheta)\; \kappa (\mathrm{grad}T
            )^T\; \mathrm{d} V
```    

is rewritten with the test and trial functions as
```math
\sum_{i=1}^N \int_{V}(\mathrm{grad}N_{j})\; \kappa (\mathrm{grad}N_{i}
            )^T\; \mathrm{d} V \; T_i
```    
The sum over the degree of freedom number ``i`` should be split: some of the  coefficients ``T_i`` are for free degrees of freedom (``1 \le i \le  N_{\mathrm{f}}``, with ``N_{\mathrm{f}}`` being the total number of free degrees of freedom), while some are  fixed (prescribed) for nodes  which are located on the essential boundary condition surface ``S_1``  (``N_{\mathrm{f}} < i \le N``).

Thus the term splits into two  pieces,

```math
\sum_{i=1}^{N_{\mathrm{f}}} \int_{V}(\mathrm{grad}N_{j})\; \kappa (\mathrm{grad}N_{i}
            )^T\; \mathrm{d} V \; T_i
```    

where the  individual integrals are entries of the conductivity matrix, and

```math
\sum_{i=N_{\mathrm{f}}+1}^N \int_{V}(\mathrm{grad}N_{j})\; \kappa (\mathrm{grad}N_{i}
            )^T\; \mathrm{d} V \; T_i
```


which  will represent heat loads  due to nonzero  prescribed boundary condition.

The FEM machine  for the heat conduction problem can be created as

```julia
material = MatHeatDiff(thermal_conductivity)
femm = FEMMHeatDiff(IntegDomain(fes, TriRule(1)), material)
```

where we first create a `material` to  provide access to the thermal conductivity matrix ``\kappa``, and then  we create  the FEM  machine  from the integration domain  for a mesh  consisting of three node triangles, using one-point integration rule, and the material. This  FEM machine  can then be passed to a method, for instance the calculate the global conductivity matrix `K`

```julia
K = conductivity(femm, geom, Temp)
```

where the geometry comes from the geometry field `geom`, and the temperature field `Temp` provides the  numbering of the degrees of freedom. Note that the global conductivity matrix is square, and of size ``N_{\mathrm{f}}\timesN_{\mathrm{f}}``. In other words, it is only for the degrees of freedom that are free (actual unknowns).

The heat load term  due to the  nonzero essential boundary conditions  is evaluated with the method `nzebcloadsconductivity`

```julia
F2 = nzebcloadsconductivity(femm, geom, Temp);
```

where the geometry comes from the geometry field `geom`, and the temperature field `Temp` provides the  numbering of the degrees of freedom and the values of the prescribed (fixed) degrees of freedom. The result is a contribution to the global heat load vector. 

### Base FEM machine

The following  operations are provided  by the base FEM machine:

- Integrate  a function expressed in terms of a field. This is typically used to
  evaluate RMS discretization errors.

- Integrate a function of the position. Perhaps the evaluation of the moments of
  inertia,  or the calculation of the volume.

- Transfer field between meshes of different resolutions.

- Calculate  the distributed-load system vector.

- Construct a field  from integration-point quantities. This is typically used in the postprocessing phase, for instance to construct continuous distribution of stresses in the structure.

## Material and Material Orientation

The material response  is described in  material-point-attached coordinate system. These coordinate systems  are Cartesian, and the material coordinate system is typically chosen to make  the response particularly simple.  So for orthotropic or transversely isotropic materials the axes would be aligned with the axes of orthotropy.

The type `CSys` (module `CSysModule`) is the updater of the material coordinate system matrix. The object is equipped with a callback to store the current orientation matrix. For instance: the coordinate system for an orthotropic material wound around a cylinder could be described in the coordinate system `CSys(3, 3, updatecs!)`, where the callback `updatecs!` is defined as

```julia
function updatecs!(csmatout::FFltMat, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)
    csmatout[:, 2] = [0.0 0.0 1.0]
    csmatout[:, 3] = XYZ
    csmatout[3, 3] = 0.0
    csmatout[:, 3] = csmatout[:, 3]/norm(csmatout[:, 3])
    csmatout[:, 1] = cross(csmatout[:, 2], csmatout[:, 3])
end
```

## Algorithms

Solution procedures and other  common operations on FEM models  are expressed  in algorithms. Anything that algorithms can do,  the user of FinEtools  can do manually, but to use an algorithm is convenient.

Algorithms typically (not always) accept a single argument, `modeldata`, a dictionary of data, keyed by Strings. Algorithms  also return `modeldata`,  typically  including additional key/value pairs that represent the data computed by the algorithm.

### Base algorithms

These are not specific to the particular physics at hand. Examples of  algorithms are  Richardson extrapolation,  calculation of the norm of the field, or calculation of the norm  of the difference of two fields. These algorithms are the exceptions, they do not return `modeldata` but rather return directly computed values.

### Model data

Model data is a dictionary, with string keys, and arbitrary values.
The documentation string for each method of an algorithm lists the required input.
For instance, for the method `linearstatics` of the `AlgoDeforLinearModule`, the
`modeldata` dictionary needs to provide key-value pairs for the finite element node set, and
the regions, the boundary conditions, and so on.

The `modeldata` may be also supplemented with additional key-value pairs inside an algorithm
and returned for further processing by other algorithms.

## Queries of quadrature-point data

A number of quantities exist at integration (quadrature) points. For instance for heat conduction this data may refer to the temperature gradients and heat flux vectors. In stress analysis, such data would typically be stress invariants  or stress components.

How this data is calculated at the quadrature point obviously varies depending on the element type. Not only on the element order, but the element formulation may invoke rules other than those of simple gradient-taking: take as an example mean-strain  elements, which define strains by using averaging rules over the entire element, so not looking at a single integration point only.

For this purpose, `FinEtools` has ways of defining implementations of the function `inspectintegpoints` to take into account the particular features of the various finite element formulations. Each FEMM typically defines its own specialized method. 

## Postprocessing

One way in which quadrature-point data is postprocessed into graphical means is by constructing node-based fields. For instance, extrapolating quadrature-point data to the nodes is commonly done in finite element programs. This procedure is typically referred to as "averaging at the nodes". The name implies that not only the quadrature-point data is extrapolated to the nodes of the element, but since each element incident on a node may have predicted (extrapolated) a different value of a quantity (for example stress), these different values need to be somehow reconciled, and averaging, perhaps weighted averaging, is the usual procedure.

### Compute continuous stress fields

Individual FEMMs may have different ways of extrapolating to the nodes. These are implemented in various methods of the function `fieldfromintegpoints`. The resulting field represents quadrature-point data as a nodal field, where the degrees of freedom are extrapolated values to the nodes.

### Compute elementwise stress fields

Most finite element postprocessing softwares find it difficult to present results which are discontinuous at inter-element boundaries. Usually the only way in which data based on individual elements with no continuity across element boundaries is presented is by taking an average over the entire element and represent the values as uniform across each element. Various methods of the function `elemfieldfromintegpoints` produce elemental fields of this nature.

## Import/export

### Importing

At the moment importing is mostly limited to the mesh data (properties, boundary conditions, analysis of data, etc. are typically not imported).
The following formats of finite element input files can be handled:

- NASTRAN (`.nas` files).
- Abaqus (`.inp` files).

### Exporting

- VTK (`.vtk` so-called legacy files). Export of geometry and fields (nodal and elemental) is supported.
- Abaqus (`.inp` files). Mesh data and selected property, boundary condition, and procedure commands can be handled.
- NASTRAN (`.nas` files). Very basic mesh and select other attributes are handled.
- STL file export of surface data.
- H2Lib triangular-surface export (`.tri` files).
- CSV file export of numerical data is supported.


## Tutorials and Examples

### Tutorials

The `FinEtools` tutorials are written up in the repositories for the applications, heat diffusion, linear and nonlinear deformation and so on.

The tutorials are in the form of Julia files with markdown. These are converted
to markdown files using the [Literate]
(https://github.com/fredrikekre/Literate.jl) workflow.

### Examples

The examples of the use of the `FinEtools` package are separated in their own separate repositories, for instance 
[`FinEtoolsHeatDiff`](https://github.com/PetrKryslUCSD/FinEtoolsHeatDiff.jl.git), [`FinEtoolsAcoustics`](https://github.com/PetrKryslUCSD/FinEtoolsAcoustics.jl.git), and so on. For a complete information refer to [the list of the repositories](https://github.com/PetrKryslUCSD?tab=repositories).

The examples are in the form of  Julia files with multiple functions, where each function defines one or more related examples. Take for instance the example file `Fahy_examples.jl`. This incantation will run all the examples from the example file:

```
include("Fahy_examples.jl"); Fahy_examples.allrun()
```

This will run just a single example from this file:

```
include("Fahy_examples.jl"); Fahy_examples.fahy_H8_example()
```

The example file `Fahy_examples.jl` consists of a module (whose name matches the name of the file), and  the module defines multiple functions,  one for each example, and one to run *all* examples, `allrun`.

### Tests

Check out the numerous tests in the `test` folder. There are hundreds of tests which exercise the various functions of the library. These examples may help you understand how to extract the desired outcome.

## Make up your own public interface

Here we assume that the FinEtools package is installed. We also assume the user works in his or her own folder, which for simplicity we assume is a package folder in the same tree as the package folder for FinEtools.

The user may have his or her additions to the FinEtools library, for instance a new material implementation, or a new FEMM (finite element model machine). Additionally, the user writes some code to solve particular problems.

In order to facilitate interactive work at the command line(REPL), it is convenient to have one or two modules so that `using` them allows for the user's code to resolve function names from the FinEtools package and from the user's own code.

Here are two ways in which this can be accomplished.

1. The user exports his or her own additions from the module `add2FinEtools` (the name of this module is not obligatory, it can be anything). In addition, the public interface to the FinEtools package needs to be brought in separately.

    ```
    using FinEtools
    using add2FinEtools
    ```

2. The user may change entirely the public interface to the FinEtools package by selectively including parts of the `FinEtools.jl` file and the code to export his or her own functionality in a single module, let us say `myFinEtools` (this name is arbitrary), so that when the user invokes

    ```
    using myFinEtools
    ```

    all the functionality that the USER considers to be public is made available by exports.

Method 1 has the *advantage* that the interface definition of the FinEtools package itself does not change, which means that package code does not need to be touched. It also has a *disadvantage* that the interface to FinEtools does not change which means that if there is a conflict with one of the exported functions from FinEtools, it needs to be resolved by fiddling with other packages.

Method 2 has the advantage that when there is a conflict between one of the exported FinEtools functions and some other function, be it from another package or the user's own, the conflict can be resolved by changing the public interface to FinEtools by the USER (as opposed to  by the DEVELOPER). Also, in this method the USER has the power to define the public interface to the FinEtools package, and if the user decides that nothing should be exported for implicit resolution of functions, that is easily accomplished.

These two methods have been described by examples in the [FinEtoolsUseCase](https://github.com/PetrKryslUCSD/FinEtoolsUseCase) package. Refer to the Readme  file and to the method descriptions  in the  method 1 and 2 folders.
