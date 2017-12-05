[Table of contents](https://petrkryslucsd.github.io/FinEtools.jl)

# Make up your own public interface

Here we assume that the FinEtools package is installed. We also assume the user works in his or her own folder, which for simplicity we assume is a package folder in the same tree as the package folder for FinEtools. 

The user may have his or her additions to the FinEtools library, for instance a new material implementation, or a new FEMM (finite element model machine). Additionally, the user writes some code to solve particular problems.

In order to facilitate interactive work at the command line(REPL), it is convenient to have one or two modules so that `using` them allows for the user's code to resolve function names from the FinEtools package and from the user's own code.

Here are two ways in which this can be accomplished.

1. The user exports his or her own additions from the module `add2FinEtools` (the name of this module is not obligatory, it can be anything). In addition, the public interface to the FinEtools package needs to be brought in separately.

```
using FinEtools
using add2FinEtools
```

2. The user may change entirely the public interface to the FinEtools package by selectively including parts of the `FinEtools.jl` file and the code to export his or her own functionality in a single module, let us say `myFinEtools` (this name is arbitrary), so that when the user invokes

```
using myFinEtools
```

all the functions that ACCORDING to the USER constitute the public interface to the FinEtools package and to the user's own additions are exported. 

