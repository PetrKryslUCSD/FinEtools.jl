
| Build       | Coverage  |
| ------------- |:-------------:|
| [![Build Status](https://travis-ci.org/PetrKryslUCSD/FinEtools.jl.svg?branch=master)](https://travis-ci.org/PetrKryslUCSD/FinEtools.jl)     |  [![codecov.io](http://codecov.io/github/PetrKryslUCSD/FinEtools.jl/coverage.svg?branch=master)](http://codecov.io/github/PetrKryslUCSD/FinEtools.jl?branch=master) |
|  [![Build status](https://ci.appveyor.com/api/projects/status/0qgyw2aa2529fahy?svg=true)](https://ci.appveyor.com/project/PetrKryslUCSD/finetools-jl)  | [![Coverage Status](https://coveralls.io/repos/github/PetrKryslUCSD/FinEtools.jl/badge.svg?branch=master)](https://coveralls.io/github/PetrKryslUCSD/FinEtools.jl?branch=master) |




# FinEtools: Finite Element tools in Julia


## News

- 06/08/2018: FinEtools is now a registered package. The release 0.3.0 works with Julia 0.7.0-alpha.

- 03/26/2018: Development of the package for Julia 0.6 is frozen. All development now goes towards 0.7.

- 03/20/2018: FinEtools version to be used with Julia 0.7 and above is at the moment being developed in a separate branch,
path_to_0.7. At this point the package runs and tests cleanly with 0.7.

- 12/10/2017: A use-case package explaining how to integrate FinEtools with  the user's own code is [available here](https://github.com/PetrKryslUCSD/FinEtoolsUseCase).

- 11/18/2017:  The documentation of the design principles  and  organization of the package is published as [Github pages](https://petrkryslucsd.github.io/FinEtools.jl). (At the moment incomplete,  but under rapid development.  Please don't hesitate to contact me if something is unclear.)

- 10/21/2017: Global RMS error norms can now be evaluated for multi-region  models.

- 10/12/2017: Evaluation of global error norms implemented.

- 08/04/2017: Several Jupyter notebook tutorials are now available. Abaqus mesh import is now also implemented (for continuum elements). Initial implementation  of high-performance mean-strain hexahedra and tetrahedra has been completed.

- 07/23/2017: Export of the finite element model to Abaqus  CAE (from Daussault Systems) has been implemented. The FinEtools package can of course  handle the finite element calculations that the export enables, but the point is to have an independent verification of the results. Export of continuum-element meshes for  vibration and static stress analysis  is currently included.


## Get FinEtools

Since  this package is not registered, please use cloning:
```julia
Pkg.clone("https://github.com/PetrKryslUCSD/FinEtools.jl")
```
Only version 0.6 and 0.7 (the nightly) builds of Julia  are supported. 

## Testing

Pkg.test("FinEtools")

## Usage and Documentation

The package comes with examples  of its use. Tutorials in notebook form
are available and more will  be added in the near future.

The documentation  is published as [Github pages](https://petrkryslucsd.github.io/FinEtools.jl). 

![Alt Sample mesh](http://hogwarts.ucsd.edu/~pkrysl/site.images/ScreenHunter_31%20Feb.%2009%2020.54.jpg "FinEtools.jl")
