[Table of contents](https://petrkryslucsd.github.io/FinEtools.jl)

# FEM machines

As an example consider the weighted-residual form of the heat balance equation
<img src="http://latex.codecogs.com/svg.latex?    \int_{V}  \vartheta c_V\frac{\partial T}{\partial t} \; d V  +\int_{V}(\mathrm{grad}\vartheta)\;{\kappa} (\mathrm{grad}T)^T\; d V  -\int_{V}  \vartheta Q \; d V +\int_{S_2} \vartheta\;\overline{{q}}_{n}\; d S+\int_{S_3} \vartheta\;h (T-T_a)  \; d S= 0 " border="0"/>
<img src="http://latex.codecogs.com/svg.latex?    \int_{V}  \vartheta c_V\frac{\partial T}{\partial t} \; d V  = 0 " border="0"/>
<img src="http://latex.codecogs.com/svg.latex? \vartheta(\mbf{x}) =0\mbox{ for}\mbf{x}\in{S_1}\;." border="0"/>