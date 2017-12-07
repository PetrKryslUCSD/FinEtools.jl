[Table of contents](https://petrkryslucsd.github.io/FinEtools.jl)

# Tutorials and Examples

## Tutorials

The tutorials are in the Jupyter notebook form. 

## Examples

The examples are in the form of  Julia files with multiple functions, where each function defines one example. Take for instance the example file `Fahy_examples.jl`. This incantation will run all the examples from the example file:

```
include("Fahy_examples.jl"); Fahy_examples.allrun()
```

This will run just a single example from this file:

```
include("Fahy_examples.jl"); Fahy_examples.fahy_H8_example()
```

The example file `Fahy_examples.jl` consists of a module (whose name matches the name of the file), and  the module defines multiple functions,  one for each example, and one to run *all* examples, `allrun`.
