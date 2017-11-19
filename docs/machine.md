[Table of contents](https://petrkryslucsd.github.io/FinEtools.jl)

# FEM machines

The construction of the matrices and vectors of the discrete form of the weighted residual equation is performed in FEM  machines. (FEM = Finite Element Method.)

As an example consider the weighted-residual form of the heat balance equation

<img src="http://latex.codecogs.com/svg.latex? \int_{V}  \vartheta c_V\frac{\partial T}{\partial t} \; \mathrm{d} V
            +\int_{V}(\mathrm{grad}\vartheta)\; {\kappa} (\mathrm{grad}T
            )^T\; \mathrm{d} V 
            -\int_{V}  \vartheta Q \; \mathrm{d} V  
            +\int_{S_2} \vartheta\;\overline{{q}}_{n}\; \mathrm{d} S+ \int_{S_3} \vartheta\;h
            (T-T_a)  \; \mathrm{d} S = 0 " border="0" />

where <img src="http://latex.codecogs.com/svg.latex? \vartheta(\mbf{x}) =0" border="0"/> for <img src="http://latex.codecogs.com/svg.latex? \mbf{x}\in{S_1} " border="0"/>.