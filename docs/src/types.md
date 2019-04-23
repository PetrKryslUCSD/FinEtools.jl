[Table of contents](https://petrkryslucsd.github.io/FinEtools.jl)

# Arithmetic types

The FinEtools package is fairly strict about typing arguments. The arithmetic types used throughout are `FInt` for integer data,
`FFlt` for floating-point data, and `Complex{FFlt}` for applications that work with complex linear algebra quantities.

The module `FTypesModule` defines these types, and also defines abbreviations for vectors and matrices with entries of these types.

Some algorithms expect input in the form of the data dictionary, `FDataDict`, and also produce output in this form.
