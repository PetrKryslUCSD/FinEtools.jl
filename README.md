[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![Build Status](https://travis-ci.org/PetrKryslUCSD/FinEtools.jl.svg?branch=master)](https://travis-ci.org/PetrKryslUCSD/FinEtools.jl) [![codecov.io](http://codecov.io/github/PetrKryslUCSD/FinEtools.jl/coverage.svg?branch=master)](http://codecov.io/github/PetrKryslUCSD/FinEtools.jl?branch=master) 
[![Build status](https://ci.appveyor.com/api/projects/status/0qgyw2aa2529fahy?svg=true)](https://ci.appveyor.com/project/PetrKryslUCSD/finetools-jl)  [![Coverage Status](https://coveralls.io/repos/github/PetrKryslUCSD/FinEtools.jl/badge.svg?branch=master)](https://coveralls.io/github/PetrKryslUCSD/FinEtools.jl?branch=master)

# FinEtools: Finite Element tools in Julia

## News

- 08/24/2018: Updated all tutorials for Julia 1.0, using the [Literate](https://github.com/fredrikekre/Literate.jl) workflow.

- 08/23/2018: Added support for NICE (nodally-integrated continuum elements) tetrahedral four-node finite elements.

- 08/09/2018: The toolkit is at this point compatible with Julia 1.0.0 in that
the tests are almost 100% passing. A handful of tests however are failing on Appveyor/Travis because of some LLVM issues, and these are at the moment left out of the test suite.


[Past news](oldnews.md)

## Get FinEtools

This package is  registered, and hence one can do just
```julia
] add FinEtools
```
Only version 0.7 (including the nightly) builds of Julia  are supported. 

## Testing

```julia
] test FinEtools 
```

## Usage and Documentation

The package comes with examples  of its use. Tutorials in notebook form
are available and more will  be added in the near future.

The documentation  is published as [Github pages](https://petrkryslucsd.github.io/FinEtools.jl). 
A use-case package explaining how to integrate FinEtools with  the user's own code is [available here](https://github.com/PetrKryslUCSD/FinEtoolsUseCase).

![Alt Visualization of sample result](http://hogwarts.ucsd.edu/~pkrysl/site.images/ScreenHunter_31%20Feb.%2009%2020.54.jpg "FinEtools.jl")
