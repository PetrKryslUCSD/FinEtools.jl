"""
    MeshExportModule

Module for export of meshes and data defined on meshes.
"""
module MeshExportModule

export vtkexportmesh

using FinEtools.FTypesModule
using FinEtools.FESetModule
using FinEtools.FENodeSetModule

const P1=1
const L2=3
const T3=5
const Q4=9
const T4=10
const H8=12
const L3=21
const T6=22
const Q8=23
const T10=24
const H20=25

VTKtypemap = Dict{DataType, Int}(FESetP1=>P1, FESetL2=>L2, FESetT3=>T3,
    FESetQ4=>Q4, FESetT4=>T4, FESetH8=>H8, FESetQ8=>Q8, FESetL3=>L3, FESetT6=>T6,
    FESetT10=>T10, FESetH20=>H20)

numnodesmap = Dict{Int, Int}(P1=>1, L2=>2, T3=>3,
    Q4=>4, T4=>4, H8=>8, Q8=>8, L3=>3, T6=>6,
    T10=>10, H20=>20)

"""
    vtkexportmesh(theFile::String, fens::FENodeSet, fes::T;
      opts...) where {T<:FESet}

Export mesh to a VTK 1.0 file as an unstructured grid.

opts = keyword argument list, where
scalars = array of tuples, (name, data)
vectors = array of tuples, (name, data)
"""
function vtkexportmesh(theFile::String, fens::FENodeSet, fes::T;
    opts...) where {T<:FESet}
    Cell_type = get(()->error("Unknown VTK type!"), VTKtypemap, typeof(fes));
    return vtkexportmesh(theFile, fes.conn, fens.xyz, Cell_type; opts...)
end

"""
    vtkexportmesh(theFile::String, Connectivity, Points, Cell_type;
      vectors=nothing, vectors_name ="vectors",
      scalars=nothing, scalars_name ="scalars")

Export mesh to a VTK 1.0 file as an unstructured grid.
"""
function vtkexportmesh(theFile::String, Connectivity, Points, Cell_type;
    vectors=nothing, scalars=nothing)
    X = Points
    if size(Points, 2) < 3
        X = zeros(size(Points,1),3)
        for j  = 1:size(Points,1)
            X[j,1:size(Points,2)] =  Points[j,:]
        end
    end
    numnodes = get(()->error("Wrong number of connected nodes!"), numnodesmap, Cell_type);
    @assert numnodes == size(Connectivity,2)

    fid=open(theFile,"w");
    if (fid == -1)
        error(["Could not open " * theFile])
        return nothing
    end
    print(fid,"# vtk DataFile Version 1.0\n");
    print(fid,"FinEtools Export\n");
    print(fid,"ASCII\n");
    print(fid,"\n");
    print(fid,"DATASET UNSTRUCTURED_GRID\n");
    print(fid,"POINTS ", size(X,1), " double\n");
    for i= 1:size(X, 1)
        for j= 1:size(X,2)-1
            print(fid,X[i,j]," ");
        end
        print(fid,X[i,end],"\n");
    end
    print(fid,"\n");
    print(fid,"\n");

    print(fid,"CELLS ",size(Connectivity,1)," ",(size(Connectivity,1)*(size(Connectivity,2)+1)),"\n");
    for i= 1:size(Connectivity, 1)
        print(fid,size(Connectivity,2)," ");
        for j= 1:size(Connectivity,2)-1
            print(fid,Connectivity[i,j]-1," ");
        end
        print(fid,Connectivity[i,end]-1,"\n");
    end
    print(fid,"\n");
    print(fid,"\n");
    print(fid,"CELL_TYPES ",size(Connectivity,1),"\n");
    for i= 1:size(Connectivity,1)
        print(fid,Cell_type,"\n");
    end
    print(fid,"\n");
    print(fid,"\n");


    did_point_data = false

    # First try to write point data
    if (scalars != nothing)
        did_point_data = false
        for  sx = 1:length(scalars)
            name, data = scalars[sx]
            if (size(data, 1) == size(Points, 1)) # point data
                if (!did_point_data)
                    print(fid,"POINT_DATA ",size(data, 1),"\n");
                    did_point_data = true
                end
                print(fid,"SCALARS ", name," double\n");
                print(fid,"LOOKUP_TABLE default\n");
                for j= 1:size(data,  1)
                    print(fid,data[j],"\n");
                end
                print(fid,"\n");
            elseif (size(data, 1) == size(Connectivity, 1)) # cell data
                # Handled below, in the section for cell data
            end
        end
    end

    print(fid,"\n");

    if vectors != nothing
        for  vx = 1:length(vectors)
            name, data = vectors[vx]
            if (!did_point_data)
                print(fid,"POINT_DATA ",size(data, 1),"\n");
                did_point_data = true
            end
            print(fid,"VECTORS ", name," double\n");
            #print(fid,"LOOKUP_TABLE default\n");
            if size(data, 2) < 3
                X = zeros(size(data,1),3)
                for j  = 1:size(data,1)
                    X[j, 1:size(data,2)] =  data[j,:]
                end
            else
                X = data
            end
            for j= 1:size(X,1)
                k=1;
                print(fid,X[j,k]);
                for k = 2:size(X,2)
                    print(fid," ", X[j,k]);
                end
                print(fid,"\n");
            end
            print(fid,"\n");
        end
    end
    print(fid,"\n");

    # Are there any cell data?  If so, write out.
    if (scalars != nothing)
        did_cell_data = false
        for  sx = 1:length(scalars)
            name, data = scalars[sx]
            if (size(data, 1) == size(Points, 1)) # point data
                # Handled above for the case of point data
            elseif (size(data, 1) == size(Connectivity, 1)) # cell data
                if (!did_cell_data)
                    print(fid,"CELL_DATA ",size(data, 1),"\n");
                    did_cell_data = true
                end
                print(fid,"SCALARS ", name," double\n");
                print(fid,"LOOKUP_TABLE default\n");
                for j= 1:size(data,  1)
                    print(fid, data[j],"\n");
                end
                print(fid,"\n");
            end
        end
    end

    fid=close(fid);
    return true
end


end
