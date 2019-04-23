[Table of contents](https://petrkryslucsd.github.io/FinEtools.jl)

# FEM machines

The construction of the matrices and vectors of the discrete form of the weighted residual equation is performed in FEM  machines. (FEM = Finite Element Method.)

As an example consider the weighted-residual form of the heat balance equation

<img src="http://latex.codecogs.com/svg.latex? \int_{V}  \vartheta c_V\frac{\partial T}{\partial t} \; \mathrm{d} V
            +\int_{V}(\mathrm{grad}\vartheta)\; \kappa (\mathrm{grad}T
            )^T\; \mathrm{d} V 
            -\int_{V}  \vartheta Q \; \mathrm{d} V  
            +\int_{S_2} \vartheta\;\overline{q}_{n}\; \mathrm{d} S+ \int_{S_3} \vartheta\;h
            (T-T_a)  \; \mathrm{d} S = 0 " border="0" />

where <img src="http://latex.codecogs.com/svg.latex? \vartheta(x) =0" border="0"/> for <img src="http://latex.codecogs.com/svg.latex? x \in{S_1} " border="0"/>.

The  test function is  taken to be  one  finite element basis function at a time, <img src="http://latex.codecogs.com/svg.latex? \vartheta = N_{\left<j\right>}" border="0"/>, and the trial function is

<img src="http://latex.codecogs.com/svg.latex? T = \sum_{i= 1} ^{N} N_{\left<i\right>} T_i" border="0"/>.

Here by <img src="http://latex.codecogs.com/svg.latex? N_{\left<j\right>}" border="0"/> we mean the basis function constructed on the mesh and associated with the node <img src="http://latex.codecogs.com/svg.latex? \left<j\right>" border="0"/>. We use the notation <img src="http://latex.codecogs.com/svg.latex? \left<j\right>" border="0"/> to mean node number at which the degree of freedom <img src="http://latex.codecogs.com/svg.latex? j" border="0"/> lives.

Now the test function and the trial function is substituted  into the  weighted residual equation.  

## Example:  internal heat generation rate term

For instance,  for the term 

<img src="http://latex.codecogs.com/svg.latex? \int_{V}  \vartheta Q \; \mathrm{d} V  
             " border="0" /> 
             
we obtain

<img src="http://latex.codecogs.com/svg.latex? \int_{V} N_{\left<j\right>} Q \; \mathrm{d} V  
             " border="0" />

This integral evaluates to a number, the heat load  applied to the degree of freedom <img src="http://latex.codecogs.com/svg.latex? j" border="0"/>. When these numbers are evaluated for all  the free degrees of freedom,  they constitute the entries of the global heat load vector. 


Evaluating integrals of this form is so common that there is a module `FEMMBaseModule` with the method `distribloads` that computes and assembles the global vector. For instance to evaluate this heat load vector  on the mesh composed of three-node triangles, for a uniform heat generation rate `Q`, we can write

```julia
fi = ForceIntensity(FFlt[Q]);
F1 = distribloads(FEMMBase(IntegDomain(fes, TriRule(1))), geom, Temp, fi, 3);
```

`IntegDomain(fes, TriRule(1))` constructs integration domain for the  finite elements `fes` using a triangular  integration rule with a single point. `FEMMBase` is the base  FEM  machine,  and all it needs at this point is the integration domain. The method  `distribloads` is defined for the  base FEM machine, the geometry field `geom`, the numbering of the degrees of freedom is taken from the field `Temp`, the internal heat generation rate is defined as the force intensity `fi`, and the integrals  are volume integrals  (3).

## Example: conductivity term

The conductivity term from the weighted residual equation

<img src="http://latex.codecogs.com/svg.latex? \int_{V}(\mathrm{grad}\vartheta)\; \kappa (\mathrm{grad}T
            )^T\; \mathrm{d} V " border="0"/>

is rewritten with the test and trial functions as

<img src="http://latex.codecogs.com/svg.latex? \sum_{i=1}^N \int_{V}(\mathrm{grad}N_{\left<j\right>})\; \kappa (\mathrm{grad}N_{\left<i\right>}
            )^T\; \mathrm{d} V \; T_i" border="0"/>

The sum over the degree of freedom number <img src="http://latex.codecogs.com/svg.latex? i" border="0"/> should be split: some of the  coefficients <img src="http://latex.codecogs.com/svg.latex? T_i" border="0"/> are for free degrees of freedom (<img src="http://latex.codecogs.com/svg.latex? 1 \le i \le  N_{\mathrm{f}}" border="0"/>, with <img src="http://latex.codecogs.com/svg.latex? N_{\mathrm{f}}" border="0"/> being the total number of free degrees of freedom), while some are  fixed (prescribed) for nodes  which are located on the essential boundary condition surface <img src="http://latex.codecogs.com/svg.latex? S_1" border="0"/> (<img src="http://latex.codecogs.com/svg.latex? N_{\mathrm{f}} < i \le N" border="0"/>).

Thus the term splits into two  pieces, 

<img src="http://latex.codecogs.com/svg.latex? \sum_{i=1}^{N_{\mathrm{f}}} \int_{V}(\mathrm{grad}N_{\left<j\right>})\; \kappa (\mathrm{grad}N_{\left<i\right>}
            )^T\; \mathrm{d} V \; T_i" border="0"/>

where the  individual integrals are entries of the conductivity matrix, and

<img src="http://latex.codecogs.com/svg.latex? \sum_{i=N_{\mathrm{f}}+1}^N \int_{V}(\mathrm{grad}N_{\left<j\right>})\; \kappa (\mathrm{grad}N_{\left<i\right>}
            )^T\; \mathrm{d} V \; T_i" border="0"/>

which  will represent heat loads  due to nonzero  prescribed boundary condition.

The FEM machine  for heat conduction  can be created as

```julia
material = MatHeatDiff(thermal_conductivity)
femm = FEMMHeatDiff(IntegDomain(fes, TriRule(1)), material)
```

where we first create a `material` to  deliver the thermal conductivity matrix <img src="http://latex.codecogs.com/svg.latex? \kappa" border="0"/>, and then  we create  the FEM  machine  from the integration domain  for a mesh  consisting of three node triangles, using one-point integration rule, and the material. This  FEM machine  can then be passed to a method

```julia
K = conductivity(femm, geom, Temp)
```

to evaluate the global conductivity matrix `K`, where the geometry comes from the geometry field `geom`, and the temperature field `Temp` provides the  numbering of the degrees of freedom.

The heat load term  due to the  nonzero essential boundary conditions  is evaluated with the method `nzebcloadsconductivity`

```julia
F2 = nzebcloadsconductivity(femm, geom, Temp);
```

where the geometry comes from the geometry field `geom`, and the temperature field `Temp` provides the  numbering of the degrees of freedom and the values of the prescribed (fixed) degrees of freedom. The result is a contribution to the global heat load vector.

## Base FEM machine

The following  operations are provided  by the base FEM machine:

- Integrate  a function expressed in terms of a field. This is typically used  to  evaluate RMS discretization errors.

- Integrate a function of the position. Perhaps the evaluation of the moments of inertia,  or the calculation of the volume.

- Transfer field between meshes of different resolutions.

- Calculate  the distributed-load system vector.

- Construct a field  from integration-point quantities. This is typically used in the postprocessing phase, for instance to construct continuous distribution of stresses in the structure.

## Acoustics FEM machines

There is one for  the interior integrals  and one for  boundary integrals. 
The  machine for the interior integrals can be used to compute:

- Evaluate the acoustic-mass matrix and the acoustic-stiffness matrix.

- Evaluate the load vector corresponding to prescribed pressure  or the prescribed second order  rate of the pressure.

The machine for the boundary integrals can be used to compute:

- Compute  transformation matrix to convert  pressure  to resultant force  or pressure to resultant torque.

- Compute the acoustic  ABC  (absorbing boundary condition) matrix.

##  Heat  conduction FEM machines

There is one for  the interior integrals  and one for  boundary integrals. 
The  machine for the interior integrals can be used to compute:

- Evaluate the conductivity matrix.

- Evaluate the load vector corresponding to prescribed temperature.

The machine for the boundary integrals can be used to compute:

- Compute surface heat transfer  matrix.

- Compute  the heat load vector for surface heat transfer.

- Compute the heat load vector  corresponding to prescribed temperatures on the boundary  with surface heat transfer.

## Linear deformation FEM  machines

For  the base machine for linear deformation, `FEMMDeforLinearBase`, assumes standard isoparametric  finite elements. It evaluates  the interior integrals:

- The stiffness matrix, the mass matrix.

- The load vector corresponding to thermal strains.

Additionally:

- Function to inspect  integration points.

The FEM machine `FEMMDeforLinear` simply stores the data required by the base `FEMMDeforLinearBase`.
 
The machine `FEMMDeforWinkler` is specialized for the boundary integrals for bodies  supported  on continuously distributed springs:

- Compute the stiffness matrix corresponding to the springs.

The  mean-strain FEM machine `FEMMDeforLinearMS` implements advanced hexahedral and tetrahedral elements based on multi-field theory and  energy-sampling  stabilization. It provides functions to compute:

- The stiffness matrix, the mass matrix.

- The load vector corresponding to thermal strains.

Additionally it defines:

- Function to inspect  integration points.









