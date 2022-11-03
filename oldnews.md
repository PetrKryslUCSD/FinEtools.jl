# Past news

- 12/03/2021: Added a utility to generate a random-looking triangular mesh.

- 11/17/2021: Added utility to distort a block of mesh.

- 09/21/2021: Fixed import problems for Abaqus files.

- 08/29/2021: Mesh export/import in the HDF5 format using DataDrop.

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

- 10/11/2019: Corrected a design blunder in the matrix utilities.

- 06/11/2019: Applications have been extracted from FinEtools into their own separate packages. This will make the base library lighter and easier to understand.

- 05/19/2019: Implementation of the computation of the infsup condition for isoparametric, mean-strain, and nodally-integrated solid finite elements.

- 04/27/2019: The naming and various definitions of abstract types has been unified and streamlined. Because of the ensuing (slight) incompatibilities, the toolkit has been released in the version v2.0.0.

- 04/25/2019: Added support for Documenter-generated guide and user manual.

- 04/13/2019: Added Reverse Cuthill-McKee renumbering.

- 03/07/2019: Meshing functions for circles, spherical surfaces, and cylinders added.

- 11/09/2018: The name IntegData was changed to IntegDomain to better reflect the meaning of this type. Since this is an incompatible change, v1.0.0 tag was released.

- 10/02/2018: The code-coverage computation seems to be broken. The coverage in the FinEtools package hasn't actually changed and it is still at 98%.

- 09/20/2018: Separated out examples and tutorials from the FinEtools library itself.

- 09/02/2018: Updated all examples for Julia 1.0.

- 08/24/2018: Updated all tutorials for Julia 1.0, using the [Literate](https://github.com/fredrikekre/Literate.jl) workflow.

- 08/23/2018: Added support for NICE (nodally-integrated continuum elements) tetrahedral four-node finite elements.

- 08/09/2018: The toolkit is at this point compatible with Julia 1.0.0 in that
the tests are almost 100% passing. A handful of tests however are failing on Appveyor/Travis because of some LLVM issues, and these are at the moment left out of the test suite.

- 06/08/2018: FinEtools is now a registered package. The release 0.3.0 works with Julia 0.7.0-alpha. The notebook examples still need updating for 0.7.

- 03/26/2018: Development of the package for Julia 0.6 is frozen. All development now goes towards 0.7.

- 03/20/2018: FinEtools version to be used with Julia 0.7 and above is at the moment being developed in a separate branch,
path_to_0.7. At this point the package runs and tests cleanly with 0.7.

- 12/10/2017: A use-case package explaining how to integrate FinEtools with  the user's own code is [available here](https://github.com/PetrKryslUCSD/FinEtoolsUseCase).

- 11/18/2017:  The documentation of the design principles  and  organization of the package is published as [Github pages](https://petrkryslucsd.github.io/FinEtools.jl). (At the moment incomplete,  but under rapid development.  Please don't hesitate to contact me if something is unclear.)

- 10/21/2017: Global RMS error norms can now be evaluated for multi-region  models.

- 10/12/2017: Evaluation of global error norms implemented.

- 08/04/2017: Several Jupyter notebook tutorials are now available. Abaqus mesh import is now also implemented (for continuum elements). Initial implementation  of high-performance mean-strain hexahedra and tetrahedra has been completed.

- 07/23/2017: Export of the finite element model to Abaqus  CAE (from Daussault Systems) has been implemented. The FinEtools package can of course  handle the finite element calculations that the export enables, but the point is to have an independent verification of the results. Export of continuum-element meshes for  vibration and static stress analysis  is currently included.
