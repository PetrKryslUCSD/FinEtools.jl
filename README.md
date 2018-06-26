[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![Build Status](https://travis-ci.org/PetrKryslUCSD/FinEtools.jl.svg?branch=master)](https://travis-ci.org/PetrKryslUCSD/FinEtools.jl) [![codecov.io](http://codecov.io/github/PetrKryslUCSD/FinEtools.jl/coverage.svg?branch=master)](http://codecov.io/github/PetrKryslUCSD/FinEtools.jl?branch=master) 
[![Build status](https://ci.appveyor.com/api/projects/status/0qgyw2aa2529fahy?svg=true)](https://ci.appveyor.com/project/PetrKryslUCSD/finetools-jl)  [![Coverage Status](https://coveralls.io/repos/github/PetrKryslUCSD/FinEtools.jl/badge.svg?branch=master)](https://coveralls.io/github/PetrKryslUCSD/FinEtools.jl?branch=master)

# FinEtools: Finite Element tools in Julia

## News

- 06/19/2018: The testing is at the moment intentionally crippled in order to assist with debugging of a Julia 0.7-alpha issue.
It is always possible to run the tests manually, for instance `include("test/test_miscellaneous.jl")`.

- 06/08/2018: FinEtools is now a registered package. The release 0.3.0 works with Julia 0.7.0-alpha. The notebook examples still need updating for 0.7.

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
