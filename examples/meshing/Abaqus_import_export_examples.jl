module Abaqus_import_export_examples
using FinEtools
using FinEtools.MeshExportModule
using FinEtools.MeshImportModule
using FinEtools.MeshImportModule: import_ABAQUS

function Abaqus_import_export()
    ## Solid cylinder/taper/sphere—-temperature loading; quadratic brick mesh
    
    # The mesh  will be created in a very coarse representation from the
    # key points in the drawing. The first coordinate is radial, the second coordinate is axial.
    rz=[1.     0.;#A
    1.4    0.;#B
    0.995184726672197   0.098017140329561;
    1.393258617341076 0.137223996461385;
    0.980785  0.195090;#
    1.37309939 0.27312645;
    0.956940335732209   0.290284677254462
    1.339716470025092 0.406398548156247
    0.9238795  0.38268;#C
    1.2124  0.7;#D
    0.7071  0.7071;#E
    1.1062  1.045;#F
    0.7071  (0.7071+1.79)/2;#(E+H)/2
    1.      1.39;#G
    0.7071  1.79;#H
    1.      1.79;#I
    ]*phun("M")
    tolerance =1.e-6*phun("M")
    
    # This is the quadrilateral mesh of the cross-section.   It will be modified and
    # refined as  we go.
    fens = FENodeSet(rz);
    fes = FESetQ4([1 2 4 3; 3 4 6 5; 5 6 8 7; 7 8 10 9; 9 10 12 11; 11 12 14 13; 13 14 16 15]);
    
    nref = 1;
    for ref = 1:nref
        fens,fes = Q4refine(fens,fes);
        list = selectnode(fens, distance=1.0+0.1/2^nref, from=[0. 0.], inflate=tolerance);
        fens.xyz[list,:] = FinEtools.MeshUtilModule.ontosphere(fens.xyz[list,:],1.0);
    end
    
    ##
    # The mesh is extruded by sweeping around the axis of symmetry.
    # Only a single layer of elements is generated of internal angle
    # |angslice|.
    nLayers = 7;
    angslice  = 5*pi/16;
    
    ##
    # First the mesh is extruded to a block whose third dimension
    # represents the angular coordinate.
    fens,fes = H8extrudeQ4(fens, fes, nLayers, (rz,k)->[rz[1],rz[2],0.0]-(k)/nLayers*[0.,0.,angslice]);
    # The block is now converted  to the axially symmetric geometry by using the
    # third (angular) coordinate  to sweep out  an axially symmetric domain. The
    # ccoordinates of the nodes at this point are |rza|,  radial distance,
    # Z-coordinate, angle.
    sweep(rza) = [-rza[1]*sin(rza[3]+angslice/2.0), rza[1]*cos(rza[3]+angslice/2.0), rza[2]]
    for j=1:size(fens.xyz,1)
        fens.xyz[j,:] = sweep(fens.xyz[j,:])
    end
    
    AE = AbaqusExporter("LE11NAFEMS_H8");
    HEADING(AE, "LE11NAFEMS: Linear bricks.");
    PART(AE, "part1");
    END_PART(AE);
    ASSEMBLY(AE, "ASSEM1");
    INSTANCE(AE, "INSTNC1", "PART1");
    NODE(AE, fens.xyz);
    ELEMENT(AE, "c3d8rh", "AllElements", 1, connasarray(fes))
    END_INSTANCE(AE);
    END_ASSEMBLY(AE);
    close(AE)
    
    
    output = MeshImportModule.import_ABAQUS("./LE11NAFEMS_H8.inp")
    fens, fes = output["fens"], output["fesets"][1]
    
    File = "LE11NAFEMS_H8.vtk"
    MeshExportModule.vtkexportmesh(File, fens, fes)
    @async run(`"paraview.exe" $File`)
end # Abaqus_import_export


function Abaqus_import_export_LE10_T10()
    output = import_ABAQUS("$(@__DIR__)/" * "NLE10-coarse-T10.inp")
    fens, fes = output["fens"], output["fesets"][1]
    
    File = "LE10NAFEMS_T10.vtk"
    MeshExportModule.vtkexportmesh(File, fens, fes)
    @async run(`"paraview.exe" $File`)
end # Abaqus_import_export_LE10_T10

function allrun()
    println("#####################################################") 
    println("# Abaqus_import_export ")
    Abaqus_import_export()
    println("#####################################################") 
    println("# Abaqus_import_export_LE10_T10 ")
    Abaqus_import_export_LE10_T10()
    return true
end # function allrun

end # module Abaqus_import_export_examples
