[Table of contents](https://petrkryslucsd.github.io/FinEtools.jl)

# Finite element

The  finite element set is one of the basic entities in FinEtools.

The finite element set is a collection of  finite elements defined by the connectivity (collection of node numbers, listing the nodes connected by the element in  a specific order). The finite element set  provides  specialized methods  to compute values of basis functions and the values of  the gradients of the basis functions  with respect to the parametric coordinates.

## Types

The finite element sets are instances of concrete types. Each particular shape and order of element has its own type. There are types for  linear  and quadratic quadrilaterals, for instance, `FESetQ4` and `FESetQ8`. Each element set provides access to the number of nodes  connected by the element (`nodesperelem`),  the connectivity as the two dimensional array    `conn`,  and the  integer label vector `label`. 

The concrete finite element set types are subtypes of the abstract type for elements of different manifold dimension (3, 2, 1, and 0), for instance for the quadrilaterals that would be `FESet2Manifold`. These types are in turn  subtypes of the abstract finite element set type `FESet`.

The concrete finite element set type provides specialized methods to compute the values of the basis functions, `bfun`, and methods to compute  the gradients of the basis functions with respect to the parametric coordinates, `bfundpar`.

## Finite element set functions

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
