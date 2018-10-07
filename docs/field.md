[Table of contents](https://petrkryslucsd.github.io/FinEtools.jl)

# Fields

The structure to maintain the numbering  and values of the degrees of freedom in the mesh  is the field.

## Abstract  Field

The assumption is that a field has one set of degrees of freedom per node or per element. For simplicity we will refer to the nodes and elements as entities.
It assumes that concrete  subtypes of the abstract field  have the following data, one row per entity:

- `values::FMat{T}`: Array of the values of the degrees of freedom, one row  for each entity. All the arrays below have the same dimensions as this one.
- `dofnums::FIntMat`: Array  of the numbers of the free degrees of freedom. If the degree of freedom is fixed (prescribed), the corresponding entry is zero.
- `is_fixed::Matrix{Bool}`: Array of  Boolean flags,  `true` for fixed  (prescribed) degrees of freedom, `false` otherwise.
- `fixed_values::FMat{T}`: Array  of the same size and type  as  `values`. Its entries are only relevant  for the fixed (prescribed)  degrees of freedom.
- `nfreedofs::FInt`:  the total number of free degrees of freedom.

The methods defined for the abstract field  include:

- Return the number of degrees of freedom and the number of entities.

- Gather and scatter the system vector.

- Gather elementwise  vectors or matrices of values, the degree of freedom numbers, or the fixed values of the degrees of freedom. 

- Set  or clear essential boundary conditions..

- Copy a field. Clear the entries of the field.

## Nodal Field

In this case  the  abstract field is subtyped to a concrete field where the entities are nodes.

## Elemental Field

In this case  the  abstract field  is subtyped to a concrete field where the entities are the elements.


## General Field

In this case  the  abstract field  is subtyped to a concrete field where the entities are  use-case  specific. 

## Numbering of the degrees of freedom

The simplest method is at the moment implemented: number all free degrees of freedom, row-by-row and column-by-column,
starting from 1 up to `f.nfreedofs`, for the field `f`.

The prescribed degrees of freedom are not numbered, and the prescribed degrees of freedom are marked with 0.
