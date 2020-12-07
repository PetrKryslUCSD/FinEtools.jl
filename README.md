[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![Build Status](https://img.shields.io/travis/PetrKryslUCSD/FinEtools.jl/master.svg?label=Linux+MacOSX+Windows)](https://travis-ci.com/PetrKryslUCSD/FinEtools.jl) 
[![Coverage Status](https://coveralls.io/repos/github/PetrKryslUCSD/FinEtools.jl/badge.svg?branch=master)](https://coveralls.io/github/PetrKryslUCSD/FinEtools.jl?branch=master) 
[![Latest documentation](https://img.shields.io/badge/docs-latest-blue.svg)](https://petrkryslucsd.github.io/FinEtools.jl/latest)

# FinEtools: Finite Element tools in Julia

`FinEtools` is a package for basic operations on finite element meshes: Construction, modification, selection, and evaluation of quantities defined on a mesh. Utilities are provided for maintaining mesh-based data (fields), for defining normals and loads, for working with physical units and coordinate systems, and for integrating over finite element meshes. ![Alt Visualization of acoustic pressure](http://hogwarts.ucsd.edu/~pkrysl/site.images/baffled-piston-aa.png "FinEtools.jl")

The package supports application packages, for instance:

- [Meshing](https://github.com/PetrKryslUCSD/FinEtoolsMeshing.jl);
- [Meshing of biomedical images](https://github.com/PetrKryslUCSD/FinEtoolsVoxelMesher.jl);
- [Linear acoustics](https://github.com/PetrKryslUCSD/FinEtoolsAcoustics.jl);
- [Heat conduction](https://github.com/PetrKryslUCSD/FinEtoolsHeatDiff.jl);
- [Linear stress analysis](https://github.com/PetrKryslUCSD/FinEtoolsDeforLinear.jl);
- [Nonlinear stress analysis](https://github.com/PetrKryslUCSD/FinEtoolsDeforNonlinear.jl);
- [Nonlinear analysis of flexible beam structures](https://github.com/PetrKryslUCSD/FinEtoolsFlexBeams.jl);
- [Vibration in fluids](https://github.com/PetrKryslUCSD/FinEtoolsVibInFluids.jl).

## News

- 12/07/2020: Export of meshes in the .mesh format with the labels enabled.
- 08/02/2020: Enabled permuting of nodes in the field.
- 02/29/2020: Many new tests added. Code coverage enabled.
- 01/23/2020: Dependencies have been updated to work with Julia 1.3.1.
- 01/02/2020: Matrix multiplication code improved with the help of the `LoopVectorization` package.


[Past news](oldnews.md)

## Get FinEtools

This package is  registered, and hence one can do just
```julia
] add FinEtools
```
Only version 1.x and the nightly builds of Julia are supported.

## Testing

```julia
] test FinEtools
```

## Usage and Documentation

Tutorials in the form
of marked-down Julia source files using the
[Literate](https://github.com/fredrikekre/Literate.jl) workflow are available
and more will  be added in the near future. (Each application package has some tutorials. For instance for the [linear deformation](https://github.com/PetrKryslUCSD/FinEtoolsDeforLinear.jl/tree/master/tutorials).)

The package has been used to build applications for various purposes. In
addition to the list at the top, there are also examples of various  [mesh
generation
tasks](https://github.com/PetrKryslUCSD/FinEtoolsMeshGenerationExamples.git)).

The documentation  is published as [Github pages](https://petrkryslucsd.github.io/FinEtools.jl/latest/).

![Alt Visualization of mechanical stress](http://hogwarts.ucsd.edu/~pkrysl/site.images/ScreenHunter_31%20Feb.%2009%2020.54.jpg "FinEtools.jl")
