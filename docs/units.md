[Table of contents](https://petrkryslucsd.github.io/FinEtools.jl)

# Physical units

The `PhysicalUnitModule` provides a simple function, `phun`, which can help with providing input numbers with the correct conversion between physical units. For instance, it is possible to specify the input data as

```julia
E = 200*phun("GPa");# Young's modulus
nu = 0.3;# Poisson ratio
rho = 8000*phun("KG*M^-3");# mass density
L = 10.0*phun("M"); # side of the square plate
t = 0.05*phun("M"); # thickness of the square plate
```

The usual sets of units are included, `:US`, `:IMPERIAL`, `:CGS`, `:SIMM` (millimeter-based SI units), and `:SI` (meter-based SI units). The resulting  values assigned to the variables are floating-point numbers, for instance

```julia
julia> E = 200*phun("GPa")
2.0e11
```

Numbers output by the simulation can also be converted  to appropriate units for printing as

```julia
julia> E/phun("MPa")
200000.0
```
