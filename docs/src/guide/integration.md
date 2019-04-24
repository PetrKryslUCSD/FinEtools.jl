[Table of contents](https://petrkryslucsd.github.io/FinEtools.jl)

# Integration

There are two kinds of integrals in the weighted-residual finite element method: integrals over the interior  of the domain,  and integrals over the boundary of the domain.

Consequently, in a typical simulation one would need  two meshes: one for the interior  of the domain,  and one for the boundary. Obviously, the one for the boundary will be derived from the mesh  constructed for the interior.

Often only a part  of the entire boundary   is   used:  on some parts of the boundary  the  boundary condition is implied as homogeneous (i. e. zero). For instance, a traction-free boundary. Therefore the necessary integrals are typically evaluated over a subset of the entire boundary.

## Manifold dimension

Finite elements  have  a certain manifold dimension.  Tetrahedra  and hexahedra are three-manifolds, triangles and quadrilaterals are two-manifolds, triangles and quadrilaterals are two-manifolds, lines are one-manifolds, and points are zero-manifolds.

Elements are equipped with an "other" dimension which boosts the manifold dimension to produce the required dimension for  the integration. For instance,  a line element can be equipped with an "other" dimension to represent a cross-section so that a volume integral can be evaluated over a line element. Or, a line element can be given an "other" dimension as a thickness to result in a physical dimension needed to evaluate a surface integral.

The "other"  dimension  has the following meaning  for finite elements of different manifold dimensions:

| Manifold dimension        | Volume integral           | Surface integral  |
| ------------- | ------------- | ----- |
| 3     | NA | NA |
| 2    | Thickness  |  NA |
| 1 | Cross-section   |  Thickness |
| 0 | Volume   |  Cross-section |

## Integration  over the interior

The integrals are always  *volume* integrals. This means that for elements which are of  lower manifold dimension than three the "other"  dimension needs to compensate.

For  three-manifold finite elements (tetrahedra and hexahedra) the "other" dimension is always 1.0. This really means there is no "other" dimension to a volume-like element.

For  finite elements of manifold dimension  less than tthree, the  "other" dimension varies according to the model (axially symmetric versus simple  plane 2D) as shown  in the table below.

| Manifold dimension        | Axially symmetric    | Plane 2D |
| ------------- | ------------- | ----- |
| 2 | <img src="http://latex.codecogs.com/svg.latex? 2\pi r" border="0"/>  |  Thickness |
| 1 | <img src="http://latex.codecogs.com/svg.latex? 2\pi r\times" border="0"/> Thickness  |  Cross-section |
| 0 | <img src="http://latex.codecogs.com/svg.latex? 2\pi r\times" border="0"/> Cross-section  |  Volume |  

The integral  is approximated with numerical quadrature as

<img src="http://latex.codecogs.com/svg.latex? \int_{\Omega} f dV \approx \sum_q f(\xi_q) J(\xi_q) W_q " border="0"/>

Here <img src="http://latex.codecogs.com/svg.latex? f" border="0"/> is the integrand, <img src="http://latex.codecogs.com/svg.latex? f(\xi_q)" border="0"/> is the  value of the integrand  at the quadrature point,
<img src="http://latex.codecogs.com/svg.latex? J(\xi_q)" border="0"/> is the  value of the Jacobian  at the quadrature point.
Importantly, the Jacobian incorporates the "other" dimension,  and therefore it is the  *volume* 
Jacobian. (For the interior integrals the Jacobian  is computed by the `Jacobianvolume` method.)

## Integration  over the boundary

The integrals are always  *surface* integrals. This means that for elements which are of  lower manifold 
dimension than two the "other"  dimension needs to compensate.

For  two-manifold finite elements (triangles and quadrilaterals) the "other" dimension is always 1.0.
This really means there is no "other" dimension to a surface-like element.

For  finite elements of manifold dimension  less than two, the  "other" dimension varies according 
to the model (axially symmetric versus simple  plane 2D) as shown  in the table below.

| Manifold dimension        | Axially symmetric    | Plane 2D |
| ------------- | ------------- | ----- |
| 1 | <img src="http://latex.codecogs.com/svg.latex? 2\pi r" border="0"/>  |  Thickness |
| 0 | <img src="http://latex.codecogs.com/svg.latex? 2\pi r\times" border="0"/>  Thickness  |  Cross-section |

The integral  is approximated with numerical quadrature as

<img src="http://latex.codecogs.com/svg.latex? \int_{\partial \Omega} f dS \approx \sum_q f(\xi_q) J(\xi_q) W_q " border="0"/> 

Here <img src="http://latex.codecogs.com/svg.latex? f" border="0"/> is the integrand, <img src="http://latex.codecogs.com/svg.latex? f(\xi_q)" border="0"/> is the  value of the integrand  at the quadrature point,
<img src="http://latex.codecogs.com/svg.latex? J(\xi_q)" border="0"/> is the  value of the Jacobian  at the quadrature point. Importantly, the Jacobian incorporates the "other" dimension,  and therefore it is the  *surface* Jacobian. (For the boundary integrals the Jacobian  is computed by the `Jacobiansurface` method.)

### Example: axially symmetric model, line element L2

The surface Jacobian in this case  is  equal to the curve Jacobian times `2*pi*r`.

## Integration Data 

The  module `IntegDomainModule` supports  the processing of  the geometry necessary for the evaluation of the various integrals.
The module data structure  groups together  a finite element set with an appropriate integration rule, information about the model (axially symmetric or not), and a callback to evaluate  the "other" dimension.

## Other dimension

The  discussion of the surface and volume [integrals](integration.html) introduces the notion  of the  "other"  dimension. In order to evaluate Jacobians of various space dimensions  the  Geometry Data module takes into account  whether or not the model is axially symmetric, and evaluates the "other" dimension based upon this information. 

A finite element set is equipped with  a way of  calculating  the "other" dimension.  For instance, the line element with two nodes, L2, can be given  the "other" dimension  as a  "thickness"  so that  surface integrals  can be evaluated over the line element. However, if this line element  is used in an axially symmetric model, the same  "other" dimension  of "thickness"  will result in the integral  along the length of this line element  being a volume integral.

Thus, the way in which the "other"  dimension gets used by the Geometry Data methods depends on the model. As an example, consider  the  method
```julia
function Jacobianvolume(self::IntegDomain{T}, J::FFltMat, loc::FFltMat, conn::CC, N::FFltMat)::FFlt where {T<:FESet2Manifold, CC<:AbstractArray{FInt}}
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



## Evaluation of integration data

Importantly, the  Integration Domain (`IntegDomain`) method `integrationdata` evaluates quantities  needed for numerical integration: locations and weights of quadrature points, and the values of basis functions and of the basis function gradients with respect to the parametric coordinates at the quadrature points.
