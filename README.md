[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![Build status](https://github.com/PetrKryslUCSD/FinEtools.jl/workflows/CI/badge.svg)](https://github.com/PetrKryslUCSD/FinEtools.jl/actions)
[![Code Coverage](https://codecov.io/gh/PetrKryslUCSD/FinEtools.jl/branch/main/graph/badge.svg)](https://app.codecov.io/gh/PetrKryslUCSD/FinEtools.jl)
[![Documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://petrkryslucsd.github.io/FinEtools.jl/dev)
[![Codebase Graph](https://img.shields.io/badge/Codebase-graph-green.svg)](diagram.svg) <!--(https://github.com/githubocto/repo-visualizer) -->
 
# FinEtools: Finite Element tools in Julia

`FinEtools` is a package for basic operations on finite element meshes: Construction, modification, selection, and evaluation of quantities defined on a mesh. Utilities are provided for maintaining mesh-based data (fields), for defining normals and loads, for working with physical units and coordinate systems, and for integrating over finite element meshes. ![Alt Visualization of acoustic pressure](http://hogwarts.ucsd.edu/~pkrysl/site.images/baffled-piston-aa.png "FinEtools.jl")

The package supports application packages, for instance:

- [Meshing](https://github.com/PetrKryslUCSD/FinEtoolsMeshing.jl);
- [Meshing of biomedical images](https://github.com/PetrKryslUCSD/FinEtoolsVoxelMesher.jl);
- [Linear acoustics](https://github.com/PetrKryslUCSD/FinEtoolsAcoustics.jl);
- [Heat conduction](https://github.com/PetrKryslUCSD/FinEtoolsHeatDiff.jl);
- [Linear stress analysis](https://github.com/PetrKryslUCSD/FinEtoolsDeforLinear.jl);
- [Nonlinear stress analysis](https://github.com/PetrKryslUCSD/FinEtoolsDeforNonlinear.jl);
- [Analysis of flexible beam and shell structures](https://github.com/PetrKryslUCSD/FinEtoolsFlexStructures.jl);
- [Vibration in fluids](https://github.com/PetrKryslUCSD/FinEtoolsVibInFluids.jl).

## News

- 12/31/2023: Update for Julia 1.10.
- 06/19/2023: Introduce DataCache, generic linear and bilinear forms. 
- 04/20/2023: Make all types in the library generic.
- 04/18/2023: Implemented resizing of assembly buffers. Makematrix! resets pointers.
- 04/16/2023: Enabled creation of finite element sets from arbitrary arrays.
- 03/15/2023: Changed strategy when assembling into the COO format.
- 03/10/2023: Genericity of arguments enabled in the assembly module.
- 12/27/2022: Added T4 and T10 cylinder meshing routines.
- 11/03/2022: Added extraction of outer surface of solid.
- 06/29/2022: Reconfigured documentation.
- 06/15/2022: Added point partitioning. Added dual connection matrix implementation.
- 05/03/2022: Updated project to Julia 1.7.2: Julia 1.6 is not supported from version 5.4.0 on.
- 01/26/2022: Added assembler for sparse diagonal matrices (such as the mass
  matrix for explicit finite element methods).
- 01/25/2022: Added export of a ParaView collection of a time sequence of data.


[Past news](oldnews.md)

## Get FinEtools

This package is  registered, and hence one can do just
```julia
] add FinEtools
```
Only version 1.x and the nightly builds of Julia are supported. The best bet is the latest stable Julia.

## Testing

```julia
] test FinEtools
```

## Usage and Documentation

Tutorials in the form
of marked-down Julia source files using the
[Literate](https://github.com/fredrikekre/Literate.jl) workflow are available
and more will  be added in the near future. Each application package has some tutorials. For a complete list refer to [the search](https://github.com/PetrKryslUCSD?tab=repositories&q=FinEtools+Tutorial&type=&language=).

The package has been used to build applications for various purposes. For a complete list refer to [the search](https://github.com/PetrKryslUCSD?tab=repositories&q=FinEtools&type=&language=).

The documentation  is published as [Github pages](https://petrkryslucsd.github.io/FinEtools.jl/latest/).

![Alt Visualization of mechanical stress](http://hogwarts.ucsd.edu/~pkrysl/site.images/ScreenHunter_31%20Feb.%2009%2020.54.jpg "FinEtools.jl")
