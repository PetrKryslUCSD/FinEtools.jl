[Table of contents](https://petrkryslucsd.github.io/FinEtools.jl)

# Algorithms

Solution procedures and other  common operations on FEM models  are expressed  in algorithms. Anything that algorithms can do,  the user of FinEtools  can do manually, but to use an algorithm is convenient.

Algorithms typically (not always) accept a single argument, `modeldata`, a dictionary of data, keyed by Strings. Algorithms  also return `modeldata`,  typically  including additional key/value pairs that represent the data computed by the algorithm.

## Base algorithms

These are not specific to the particular physics at hand. Examples of  algorithms are  Richardson extrapolation,  calculation of the norm of the field, or calculation of the norm  of the difference of two fields. These algorithms are the exceptions, they do not return `modeldata` but rather return directly computed values.

## Acoustics algorithms

At the moment there is one algorithm, for steady-state (harmonic) acoustics.

### Example:  baffled piston

After the mesh  has been generated, the `modeldata` can be set up: Here we begin with  the region.

```julia
material = MatAcoustFluid(bulk, rho)
region1 =  FDataDict("femm"=>FEMMAcoust(IntegDomain(fes, GaussRule(3, 2)), material))
```

We set up a definition of the absorbing boundary condition:

```julia
abc1  =  FDataDict("femm"=>FEMMAcoustSurf(IntegDomain(outer_fes, GaussRule(2, 2)),
          material))
```

The  surface of the piston is associated with a known-flux  boundary condition:

```julia
flux1  =  FDataDict("femm"=>FEMMAcoustSurf(IntegDomain(piston_fes, GaussRule(2, 2)),
          material),  "normal_flux"=> -rho*a_piston+0.0im);
```

And finally we make the model data,

```julia
modeldata =  FDataDict("fens"=>  fens,
                 "omega"=>omega,
                 "regions"=>[region1],
                 "flux_bcs"=>[flux1], "ABCs"=>[abc1])
```

and call  the solver:

```julia
modeldata = FinEtools.AlgoAcoustModule.steadystate(modeldata)
```

When  the algorithm completes, `modeldata["P"]` is the computed pressure field.

## Heat diffusion algorithms

There is an implementation of an algorithm for steady-state heat conduction.

## Linear deformation algorithms

There are algorithms for

- Linear static analysis;
- Export  of the deformed shape for visualization;
- Export  of the nodal and elementwise stress fields for visualization;
- Modal (free-vibration) analysis;
- Export  of modal shapes for visualization;
- Subspace-iteration method implementation.

## Model data

Model data is a dictionary, with string keys, and arbitrary values.
The documentation string for each method of an algorithm lists the required input.
For instance, for the method `linearstatics` of the `AlgoDeforLinearModule`, the
`modeldata` dictionary needs to provide key-value pairs for the finite element node set, and
the regions, the boundary conditions, and so on.

The `modeldata` may be also supplemented with additional key-value pairs inside an algorithm 
and returned for further processing by other algorithms.
