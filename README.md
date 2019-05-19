[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![Build Status](https://travis-ci.org/PetrKryslUCSD/FinEtools.jl.svg?branch=master)](https://travis-ci.org/PetrKryslUCSD/FinEtools.jl) [![codecov.io](http://codecov.io/github/PetrKryslUCSD/FinEtools.jl/coverage.svg?branch=master)](http://codecov.io/github/PetrKryslUCSD/FinEtools.jl?branch=master) 
[![Build status](https://ci.appveyor.com/api/projects/status/0qgyw2aa2529fahy?svg=true)](https://ci.appveyor.com/project/PetrKryslUCSD/finetools-jl)  [![Coverage Status](https://coveralls.io/repos/github/PetrKryslUCSD/FinEtools.jl/badge.svg?branch=master)](https://coveralls.io/github/PetrKryslUCSD/FinEtools.jl?branch=master)

[![][docs-latest-img]][docs-latest-url]

[docs-latest-img]: https://img.shields.io/badge/docs-latest-blue.svg
[docs-latest-url]: http://petrkryslucsd.github.io/FinEtools.jl/latest/

# FinEtools: Finite Element tools in Julia

![Alt Visualization of acoustic pressure](http://hogwarts.ucsd.edu/~pkrysl/site.images/baffled-piston-a.png "FinEtools.jl")

## News
 
- 05/19/2019: Implementation of the computation of the infsup condition for isoparametric, mean-strain, and nodally-integrated solid finite elements.


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

[Tutorials](https://github.com/PetrKryslUCSD/FinEtoolsTutorials.git) in the form of marked-down Julia source files using the
[Literate](https://github.com/fredrikekre/Literate.jl) workflow are available and more will  be added in the near future.

The package comes with examples  of its use 
([heat conduction](https://github.com/PetrKryslUCSD/FinEtoolsHeatConductionExamples.git), 
[linear deformation](https://github.com/PetrKryslUCSD/FinEtoolsLinearDeformationExamples.git), 
[acoustics](https://github.com/PetrKryslUCSD/FinEtoolsAcousticsExamples.git), 
[mesh generation](https://github.com/PetrKryslUCSD/FinEtoolsMeshGenerationExamples.git)). 

The documentation  is published as [Github pages](https://petrkryslucsd.github.io/FinEtools.jl/latest/). 
A use-case package explaining how to integrate FinEtools with  the user's own code is [available here](https://github.com/PetrKryslUCSD/FinEtoolsUseCase).

![Alt Visualization of mechanical stress](http://hogwarts.ucsd.edu/~pkrysl/site.images/ScreenHunter_31%20Feb.%2009%2020.54.jpg "FinEtools.jl")

