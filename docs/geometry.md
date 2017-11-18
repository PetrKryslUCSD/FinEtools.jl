[Table of contents](https://petrkryslucsd.github.io/FinEtools.jl)

# Geometry Data 

The  module `GeoDModule` supports  the processing of  the geometry necessary for the evaluation of the various integrals.
The module data structure  groups together  a finite element set with an appropriate integration rule, information about the model (axially symmetric or not), a callback to evaluate  the "other" dimension, and  a material orientation coordinate system.

## Other dimension

The  discussion of the surface and volume [integrals](integration.html) introduces the notion  of the  "other"  dimension. In order to evaluate Jacobians of various space dimensions  the  Geometry Data module takes into account  whether or not the model is axially symmetric, and evaluates the "other" dimension based upon this information. 

A finite element set is equipped with  a way of  calculating  the "other" dimension.  For instance, the line element with two nodes, L2, can be given  the "other" dimension  as a  "thickness"  so that  surface integrals  can be evaluated over the line element. However, if this line element  is used in an axially symmetric model, the same  "other" dimension  of "thickness"  will result in the integral  along the length of this line element  being a volume integral.

Thus, there is a difference between  the "other"  dimension as  calculated by the  Geometry Data  module method  and the  "other" dimension  evaluated  by the  callback function. The callback function computes the "other" dimension from  two kinds of  information: (a) the physical location  of the quadrature point,  and (b) the interpolation data for the element  (connectivity and the values of the basis functions at the quadrature point). 

The approach ad (a) is suitable  when the "other" dimension is given as a function of the physical coordinates. The  simplest case is obviously  a uniform distribution of the "other" dimension.

The approach ad (b) is appropriate when the "other" dimension  is given by values given at the nodes of the  mesh. Than the connectivity  and  the array of the values of the basis functions  can be used to interpolate the "other"  dimension  to the quadrature point.

## Integration data

Importantly, the  Geometry Data method `integrationdata` pre-computes  quantities  needed for numerical integration: locations and weights of quadrature points, and the values of basis functions and of the basis function gradients with respect to the parametric coordinates at the quadrature points.