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

Here <img src="http://latex.codecogs.com/svg.latex? N_{\left<j\right>}" border="0"/> we mean the basis function constructed on the mesh and associated with the node <img src="http://latex.codecogs.com/svg.latex? \left<j\right>" border="0"/>. We use the notation <img src="http://latex.codecogs.com/svg.latex? \left<j\right>" border="0"/> to mean node number at which the degree of freedom <img src="http://latex.codecogs.com/svg.latex? j" border="0"/> lives.

Now the test function and the trial function is substituted  into the  weighted residual equation.  For instance,  for the term <img src="http://latex.codecogs.com/svg.latex? \int_{V}  \vartheta Q \; \mathrm{d} V  
             " border="0" /> we obtain

<img src="http://latex.codecogs.com/svg.latex? \int_{V} N_{\left<j\right>} Q \; \mathrm{d} V  
             " border="0" />

This integral evaluates to a number, the heat load  applied to the degree of freedom <img src="http://latex.codecogs.com/svg.latex? j" border="0"/>. When these numbers are evaluated for all  the free degrees of freedom,  they constitute the entries of the global heat load vector. 


Evaluating integrals of this form is so common that there is a module `FEMMBaseModule` with the method `distribloads` that computes and assembles the global vector. For instance to evaluate this heat load vector  on the mesh composed of three-node triangles, for a uniform heat generation rate `Q`, we can write

```julia
fi = ForceIntensity(FFlt[Q]);
F1 = distribloads(FEMMBase(IntegData(fes, TriRule(1))), geom, Temp, fi, 3);
```
