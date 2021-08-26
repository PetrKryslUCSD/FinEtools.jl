[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![Build status](https://github.com/PetrKryslUCSD/FinEtools.jl/workflows/CI/badge.svg)](https://github.com/PetrKryslUCSD/FinEtools.jl/actions)
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

- 08/12/2021: Added quadrilateral mesh generation by extrusion of a wire.
- 08/12/2021: Added triangle mesh generation for circular segments.
- 07/14/2021: Implemented the mechanism of the "delegate of"  to facilitate creation of user-defined finite elements.
- 04/08/2021: Implemented export of VTK binary files using WriteVTK.
- 02/20/2021: Implemented extrusion of triangular meshes into tetrahedra.
- 02/07/2021: Tests clean with Julia 1.6 release candidate.
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
