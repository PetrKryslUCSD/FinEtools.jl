# Integration

There are two kinds of integrals: integrals over the interior  of the domain,  and integrals over the boundary of the domain.

Consequently, in a typical simulation one would need  two meshes: one for the interior  of the domain,  and one for the boundary.
Obviously, the one for the boundary will be derived from the mesh  constructed for the interior.
Often only at part  of the entire boundary  of the interior mesh  is   used:  on some parts of the boundary  the
boundary condition is boundary condition is implied as homogeneous (i. e. zero). For instance, a traction-free boundary. 
Therefore the necessary integrals are typically evaluated over a subset of the entire boundary.

## Manifold dimension

Finite elements  have  a certain manifold dimension.  Intuitively, tetrahedra  and hexahedra are three-manifolds,
triangles and quadrilaterals are two-manifolds, triangles and quadrilaterals are two-manifolds, 
lines are one-manifolds, and points are zero-manifolds.

Elements are equipped with  "other" dimension which boosts the manifold dimension to produce the required
dimension for  the integration. For instance,  a line element can be equipped dimension for  the integration. 
For instance,  a line element can be equipped with an "other" dimension to represent a cross-section 
so that a volume integral can be evaluated over a line element. Or, a line element can be given 
an "other" dimension as a thickness to result in a physical dimension needed to evaluate a surface integral.

The "other"  dimension  has the following meaning  for finite elements of different manifold dimensions:

| Manifold dimension        | Volume integral           | Surface integral  |
| ------------- | ------------- | ----- |
| 3     | NA | NA |
| 2    | Thickness  |  NA |
| 1 | Cross-section   |  Thickness |
| 0 | Volume   |  Cross-section |


## Integration  over the interior

The integrals are always  *volume* integrals. This means that for elements which are of  lower manifold 
dimension than three the "other"  dimension needs to compensate.

For  three-manifold finite elements (tetrahedra and hexahedra) the "other" dimension is always 1.0.
This really means there is no "other" dimension to a volume-like element.

For  finite elements of manifold dimension  less than tthree, the  "other" dimension varies according 
to the model (axially symmetric versus simple  plane 2D) as shown  in the table below.

| Manifold dimension        | Axially symmetric    | Plane 2D |
| ------------- | ------------- | ----- |
| 2 | $2\pi r\times$  |  1.0 |
| 1 | $2\pi r\times$ Thickness  |  Cross-section |
| 0 | $2\pi r  \times$ Cross-section  |  Volume |  

The integral  is approximated with numerical quadrature as

$\int_{\partial \Omega} f dS \approx \sum_q f(\xi_q) J(\xi_q) W_q $ 

Here $f$ is the integrand, $f(\xi_q)$ is the  value of the integrand  at the quadrature point,
$J(\xi_q)$ is the  value of the Jacobian  at the quadrature point.
Importantly, the Jacobian incorporates the "other" dimension,  and therefore it is the  *volume* 
Jacobian. (For the boundary integrals the Jacobian  is computed by the `Jacobianvolume` method.)

## Integration  over the boundary

The integrals are always  *surface* integrals. This means that for elements which are of  lower manifold 
dimension than two the "other"  dimension needs to compensate.

For  two-manifold finite elements (triangles and quadrilaterals) the "other" dimension is always 1.0.
This really means there is no "other" dimension to a surface-like element.

For  finite elements of manifold dimension  less than two, the  "other" dimension varies according 
to the model (axially symmetric versus simple  plane 2D) as shown  in the table below.

| Manifold dimension        | Axially symmetric    | Plane 2D |
| ------------- | ------------- | ----- |
| 1 | $2\pi r$   |  Thickness |
| 0 | $2\pi r  \times$ Thickness  |  Cross-section |

The integral  is approximated with numerical quadrature as

$\int_{\partial \Omega} f dS \approx \sum_q f(\xi_q) J(\xi_q) W_q $ 

Here $f$ is the integrand, $f(\xi_q)$ is the  value of the integrand  at the quadrature point,
$J(\xi_q)$ is the  value of the Jacobian  at the quadrature point.
Importantly, the Jacobian incorporates the "other" dimension,  and therefore it is the  *surface* 
Jacobian. (For the boundary integrals the Jacobian  is computed by the `Jacobiansurface` method.)

### Example: axially symmetric model, line element L2

The surface Jacobian in this case  is  equal to the curve Jacobian times `2*pi*r`.



