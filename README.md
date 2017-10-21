
| Build       | Coverage  |
| ------------- |:-------------:|
| [![Build Status](https://travis-ci.org/PetrKryslUCSD/FinEtools.jl.svg?branch=master)](https://travis-ci.org/PetrKryslUCSD/FinEtools.jl)     |  [![codecov.io](http://codecov.io/github/PetrKryslUCSD/FinEtools.jl/coverage.svg?branch=master)](http://codecov.io/github/PetrKryslUCSD/FinEtools.jl?branch=master) |
|  [![Build status](https://ci.appveyor.com/api/projects/status/0qgyw2aa2529fahy?svg=true)](https://ci.appveyor.com/project/PetrKryslUCSD/finetools-jl)  | [![Coverage Status](https://coveralls.io/repos/github/PetrKryslUCSD/FinEtools.jl/badge.svg?branch=master)](https://coveralls.io/github/PetrKryslUCSD/FinEtools.jl?branch=master) |




# FinEtools: Finite Element tools in Julia


## News

- 10/21/2017: Global RMS error norms can now be evaluated for multi-region  models.

- 10/12/2017: Evaluation of global error norms implemented.

- 08/04/2017: Several Jupyter notebook tutorials are now available. Abaqus mesh import is now also implemented (for continuum elements). Initial implementation  of high-performance mean-strain hexahedra and tetrahedra has been completed.

- 07/23/2017: Export of the finite element model to Abaqus  CAE (from Daussault Systems)
has been implemented. The FinEtools package can of course  handle the finite element calculations that the export enables, but the point is to have an independent verification of the results. Export of continuum-element meshes for  vibration and static stress analysis  is currently included.


## Get FinEtools

Since  this package is not registered, please use cloning:
```julia
Pkg.clone("https://github.com/PetrKryslUCSD/FinEtools.jl")
```
Only version 0.7 (the nightly) builds of Julia  are supported.

## Testing

Pkg.test("FinEtools")

## Usage

The package comes with examples  of its use. Tutorials in notebook form
are available and more will  be added in the near future.

![Alt Sample mesh](http://hogwarts.ucsd.edu/~pkrysl/site.images/ScreenHunter_31%20Feb.%2009%2020.54.jpg "FinEtools.jl")
