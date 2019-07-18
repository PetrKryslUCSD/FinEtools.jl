[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![Build Status](https://travis-ci.org/PetrKryslUCSD/FinEtools.jl.svg?branch=master)](https://travis-ci.org/PetrKryslUCSD/FinEtools.jl) [![codecov.io](http://codecov.io/github/PetrKryslUCSD/FinEtools.jl/coverage.svg?branch=master)](http://codecov.io/github/PetrKryslUCSD/FinEtools.jl?branch=master)
[![Build status](https://ci.appveyor.com/api/projects/status/0qgyw2aa2529fahy?svg=true)](https://ci.appveyor.com/project/PetrKryslUCSD/finetools-jl)  [![Coverage Status](https://coveralls.io/repos/github/PetrKryslUCSD/FinEtools.jl/badge.svg?branch=master)](https://coveralls.io/github/PetrKryslUCSD/FinEtools.jl?branch=master)

[![][docs-latest-img]][docs-latest-url]

[docs-latest-img]: https://img.shields.io/badge/docs-latest-blue.svg
[docs-latest-url]: http://petrkryslucsd.github.io/FinEtools.jl/latest/

# FinEtools: Finite Element tools in Julia

`FinEtools` is a package for basic operations on finite element meshes: Construction, modification, selection, and evaluation of quantities defined on a mesh. Utilities are provided for maintaining mesh-based data (fields), for defining normals and loads, for working with physical units and coordinate systems, and for integrating over finite element meshes. ![Alt Visualization of acoustic pressure](http://hogwarts.ucsd.edu/~pkrysl/site.images/baffled-piston-aa.png "FinEtools.jl") 

The package supports application packages, for instance:

- [Meshing](https://github.com/PetrKryslUCSD/FinEtoolsMeshing.jl);
- [Meshing of biomedical images](https://github.com/PetrKryslUCSD/FinEtoolsVoxelMesher.jl);
- [Linear acoustics](https://github.com/PetrKryslUCSD/FinEtoolsAcoustics.jl);
- [Heat conduction](https://github.com/PetrKryslUCSD/FinEtoolsHeatDiff.jl);
- [Linear stress analysis](https://github.com/PetrKryslUCSD/FinEtoolsDeforLinear.jl);
- [Nonlinear stress analysis](https://github.com/PetrKryslUCSD/FinEtoolsDeforNonlinear.jl).

## News

- 06/11/2019: Applications have been extracted from FinEtools into their own separate packages. This will make the base library lighter and easier to understand.

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

[Tutorials](https://github.com/PetrKryslUCSD/FinEtoolsTutorials.git) in the form
of marked-down Julia source files using the
[Literate](https://github.com/fredrikekre/Literate.jl) workflow are available
and more will  be added in the near future.

The package has been used to build applications for various purposes. In
addition to the list at the top, there are also examples of various  [mesh
generation
tasks](https://github.com/PetrKryslUCSD/FinEtoolsMeshGenerationExamples.git)).

The documentation  is published as [Github pages](https://petrkryslucsd.github.io/FinEtools.jl/latest/).

![Alt Visualization of mechanical stress](http://hogwarts.ucsd.edu/~pkrysl/site.images/ScreenHunter_31%20Feb.%2009%2020.54.jpg "FinEtools.jl")
