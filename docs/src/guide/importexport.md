[Table of contents](https://petrkryslucsd.github.io/FinEtools.jl)

# Import/export

## Importing

At the moment importing is mostly limited to the mesh data (properties, boundary conditions, analysis of data, etc. are typically not imported).
The following formats of finite element input files can be handled:

- NASTRAN (`.nas` files).
- Abaqus (`.inp` files).

## Exporting

- VTK (`.vtk` so-called legacy files). Export of geometry and fields (nodal and elemental) is supported.
- Abaqus (`.inp` files). Mesh data and selected property, boundary condition, and procedure commands can be handled.
- NASTRAN (`.nas` files). Very basic mesh and select other attributes are handled.
- STL file export of surface data.
- H2Lib triangular-surface export (`.tri` files).
- CSV file export of numerical data is supported.




