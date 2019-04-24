using Documenter, FinEtools
using FinEtools.FTypesModule, FinEtools.BoxModule, FinEtools.PhysicalUnitModule, FinEtools.RotationUtilModule, FinEtools.CSysModule, FinEtools.FESetModule, FinEtools.FENodeSetModule, FinEtools.FENodeToFEMapModule, FinEtools.FieldModule, FinEtools.GeneralFieldModule, FinEtools.NodalFieldModule, FinEtools.ElementalFieldModule, FinEtools.MeshUtilModule, FinEtools.MeshSelectionModule, FinEtools.MeshExportModule, FinEtools.MeshModificationModule, FinEtools.MeshImportModule, FinEtools.VectorCacheModule, FinEtools.ForceIntensityModule, FinEtools.SurfaceNormalModule, FinEtools.AssemblyModule, FinEtools.IntegRuleModule, FinEtools.IntegDomainModule, FinEtools.FEMMBaseModule, FinEtools.MeshQuadrilateralModule, FinEtools.MeshLineModule, FinEtools.MeshTriangleModule, FinEtools.MeshHexahedronModule, FinEtools.MeshTetrahedronModule, FinEtools.VoxelBoxModule, FinEtools.VoxelTetMeshingModule, FinEtools.MatHeatDiffModule, FinEtools.MatAcoustFluidModule, FinEtools.FEMMAcoustModule, FinEtools.FEMMAcoustSurfModule, FinEtools.DeforModelRedModule, FinEtools.MatDeforModule, FinEtools.MatDeforElastIsoModule, FinEtools.FEMMDeforLinearBaseModule, FinEtools.FEMMDeforLinearModule, FinEtools.FEMMDeforWinklerModule, FinEtools.FEMMDeforLinearMSModule, FinEtools.FEMMDeforSurfaceDampingModule, FinEtools.FEMMDeforLinearNICEModule, FinEtools.FEMMDeforLinearESNICEModule

makedocs(
   modules = [FinEtools, FinEtools.FTypesModule, FinEtools.BoxModule, FinEtools.PhysicalUnitModule, FinEtools.RotationUtilModule, FinEtools.CSysModule, FinEtools.FESetModule, FinEtools.FENodeSetModule, FinEtools.FENodeToFEMapModule, FinEtools.FieldModule, FinEtools.GeneralFieldModule, FinEtools.NodalFieldModule, FinEtools.ElementalFieldModule, FinEtools.MeshUtilModule, FinEtools.MeshSelectionModule, FinEtools.MeshExportModule, FinEtools.MeshModificationModule, FinEtools.MeshImportModule, FinEtools.VectorCacheModule, FinEtools.ForceIntensityModule, FinEtools.SurfaceNormalModule, FinEtools.AssemblyModule, FinEtools.IntegRuleModule, FinEtools.IntegDomainModule, FinEtools.FEMMBaseModule, FinEtools.MeshQuadrilateralModule, FinEtools.MeshLineModule, FinEtools.MeshTriangleModule, FinEtools.MeshHexahedronModule, FinEtools.MeshTetrahedronModule, FinEtools.VoxelBoxModule, FinEtools.VoxelTetMeshingModule, FinEtools.MatHeatDiffModule, FinEtools.MatAcoustFluidModule, FinEtools.FEMMAcoustModule, FinEtools.FEMMAcoustSurfModule, FinEtools.DeforModelRedModule, FinEtools.MatDeforModule, FinEtools.MatDeforElastIsoModule, FinEtools.FEMMDeforLinearBaseModule, FinEtools.FEMMDeforLinearModule, FinEtools.FEMMDeforWinklerModule, FinEtools.FEMMDeforLinearMSModule, FinEtools.FEMMDeforSurfaceDampingModule, FinEtools.FEMMDeforLinearNICEModule, FinEtools.FEMMDeforLinearESNICEModule],
   doctest = false, clean = true,
 checkdocs = :all,
    format = Documenter.HTML(prettyurls = true),
   authors = "Petr Krysl",
  sitename = "FinEtools.jl",
     pages = Any[
              "Home" => "index.md",
              "Guide" => ["Module structure" => "guide/modules.md",
	              "Arithmetic types" => "guide/types.md",
	              "Physical units" => "guide/units.md",
	              "Mesh entities" => "guide/mesh.md",
	              "Mesh generation" => "guide/meshgen.md",
	              "Selection of mesh entities" => "guide/selection.md",
	              "Field" => "guide/field.md",
	              "Element" => "guide/element.md",
	              "Integration" => "guide/integration.md",
	              "FEM machine" => "guide/machine.md",
	              "Material" => "guide/material.md",
	              "Algorithms" => "guide/algorithms.md",
	              "Querying quadrature-point data" => "guide/query.md",
	              "Postprocessing" => "guide/post.md",
	              "Import/export" => "guide/importexport.md",
	              "Tutorials and examples" => "guide/tutorials.md",],
              "Types and Functions" => Any[
              "man/types.md",
              "man/functions.md"]
             ]
)

# deploydocs(
#     repo = "https://github.com/PetrKryslUCSD/FinEtools.jl.git",
#     target = "build",
#     julia  = "nightly",
#     deps = nothing,
#     make = nothing
# )
