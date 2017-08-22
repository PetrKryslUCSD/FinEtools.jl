"""
Module for generating and manipulating voxel 3-D images (boxes of voxels).

This module may be helpful when working with CT scan and other medical images
in the context of finite element analysis carried out with the package
FinEtools.jl .
"""
module VoxelBoxModule

import Base.size

export VoxelBoxVolume, voxeldims, size,
    fillvolume!, fillsolid!,
    intersectionop, unionop, complementop, differenceop,
    solidsphere, solidhalfspace, solidbox, solidcylinder,
    trim,
    vtkexport

mutable struct VoxelBoxVolume{CoordT<:Number,DataT<:Number}
    origin::Array{CoordT,1}
    boxdim::Array{CoordT,1}
    data::Array{DataT,3}
end

function VoxelBoxVolume(::Type{CoordT},::Type{DataT}) where {CoordT<:Number,DataT<:Number}
    V = VoxelBoxVolume(zeros(CoordT, 3), zeros(CoordT, 3), Array{DataT}(0, 0, 0));
    return V
end

function VoxelBoxVolume(::Type{DataT},
    nvox::Array{Int,1}, boxdim::Array{CoordT,1}) where {CoordT<:Number,DataT<:Number}
    origin = zeros(CoordT,3)
    data = zeros(DataT,nvox...)
    V = VoxelBoxVolume(origin, boxdim, data);
    return V
end

function VoxelBoxVolume(data::Array{DataT,3},
    boxdim::Array{CoordT,1})  where {CoordT<:Number,DataT<:Number}
    V = VoxelBoxVolume(CoordT,DataT);
    copy!(V.boxdim, boxdim)
    V.data = deepcopy(data)
    return V
end

function voxeldims(V::VoxelBoxVolume)
    nx, ny, nz = size(V.data)
    voxszx=V.boxdim[1]/nx
    voxszy=V.boxdim[2]/ny
    voxszz=V.boxdim[3]/nz
    return (voxszx, voxszy, voxszz)
end

size(V::VoxelBoxVolume, which) = size(V.data, which)

struct SolidCF{F<:Function}
    f::F
end

function (o::SolidCF)(x::CoordT, y::CoordT, z::CoordT) where {CoordT<:Number}
    return o.f(x, y, z)
end

function intersectionop(o1::SolidCF, o2::SolidCF)
    function f(x, y, z)
        return o1(x, y, z) && o2(x, y, z)
    end
    return SolidCF(f)
end

function unionop(o1::SolidCF, o2::SolidCF)
    function f(x, y, z)
        return o1(x, y, z) || o2(x, y, z)
    end
    return SolidCF(f)
end

function complementop(o1::SolidCF)
    function f(x, y, z)
        return !o1(x, y, z)
    end
    return SolidCF(f)
end

function differenceop(o1::SolidCF, o2::SolidCF)
    function f(x, y, z)
        return o1(x, y, z) && (!o2(x, y, z))
    end
    return SolidCF(f)
end

"""
    solidsphere(center::Tuple{CoordT, CoordT, CoordT}, r::CoordT) where {CoordT<:Number}

Solid sphere.
"""
function solidsphere(center::Tuple{CoordT, CoordT, CoordT}, r::CoordT) where {CoordT<:Number}
    return SolidCF((x, y, z) -> (r*r - ((x-center[1])^2+(y-center[2])^2+(z-center[3])^2) >= 0.0))
end

"""
    solidhalfspace(center::Tuple{CoordT, CoordT, CoordT},
        normal::Tuple{CoordT, CoordT, CoordT}) where {CoordT<:Number}

Solid halfspace.
"""
function solidhalfspace(center::Tuple{CoordT, CoordT, CoordT},
    normal::Tuple{CoordT, CoordT, CoordT}) where {CoordT<:Number}
    nn = sqrt(normal[1]^2 + normal[2]^2 + normal[3]^2)
    n = (normal[1]/nn, normal[2]/nn, normal[3]/nn)
    function g(x, y, z)
        return (x-center[1])*normal[1] + (y-center[2])*normal[2] + (z-center[3])*normal[3]
    end
    return SolidCF((x, y, z) -> (-g(x, y, z) >= 0.0))
end

"""
    solidbox(corner1::Tuple{CoordT, CoordT, CoordT},
        corner2::Tuple{CoordT, CoordT, CoordT}) where {CoordT<:Number}

Solid box  with faces aligned with the global Cartesian axes.
"""
function solidbox(corner1::Tuple{CoordT, CoordT, CoordT},
    corner2::Tuple{CoordT, CoordT, CoordT}) where {CoordT<:Number}
    cmin = (min(corner1[1], corner2[1]), min(corner1[2], corner2[2]), min(corner1[3], corner2[3]))
    cmax = (max(corner1[1], corner2[1]), max(corner1[2], corner2[2]), max(corner1[3], corner2[3]))
    i1 = intersectionop(solidhalfspace(cmin, (-1.0, 0.0, 0.0)), solidhalfspace(cmax, (+1.0, 0.0, 0.0)))
    i2 = intersectionop(solidhalfspace(cmin, (0.0, -1.0, 0.0)), solidhalfspace(cmax, (0.0, +1.0, 0.0)))
    i3 = intersectionop(solidhalfspace(cmin, (0.0, 0.0, -1.0)), solidhalfspace(cmax, (0.0, 0.0, +1.0)))
    return intersectionop(i1, intersectionop(i2, i3))
end

function solidcylinder(center::Tuple{CoordT, CoordT, CoordT},
    axisdirection::Tuple{CoordT, CoordT, CoordT}, radius::CoordT) where {CoordT<:Number}
    nad = sqrt(axisdirection[1]^2 + axisdirection[2]^2 + axisdirection[3]^2)
    ad = (axisdirection[1]/nad, axisdirection[2]/nad, axisdirection[3]/nad)
    function g(x, y, z)
        xcad = (x-center[1])*ad[1] + (y-center[2])*ad[2] + (z-center[3])*ad[3]
        px = x-center[1] - xcad*ad[1]
        py = y-center[2] - xcad*ad[2]
        pz = z-center[3] - xcad*ad[3]
        return ((px^2+py^2+pz^2) - radius*radius)
    end
    return SolidCF((x, y, z) -> (-g(x, y, z) >= 0.0))
end

"""
    fillsolid!(V::VoxelBoxVolume,
        f::SolidCF, fillvalue::DataT) where {DataT<:Number}

Filled a solid using a solid characteristic function.
"""
function fillsolid!(V::VoxelBoxVolume,
    f::SolidCF, fillvalue::DataT) where {DataT<:Number}
    nx, ny, nz = size(V.data)
    voxszx=V.boxdim[1]/nx
    voxszy=V.boxdim[2]/ny
    voxszz=V.boxdim[3]/nz
    for i = 1:nx
        x = V.origin[1] + voxszx*(i+0.5)
        for j = 1:ny
            y = V.origin[2] + voxszy*(j+0.5)
            for k = 1:nz
                z = V.origin[3] + voxszz*(k+0.5)
                if f(x,y,z)
                    V.data[i,j,k] = fillvalue
                end
            end
        end
    end
    return V
end

"""
    fillvolume!(V::VoxelBoxVolume, fillvalue::DataT) where {DataT<:Number}

Fill the volume with a given value.
"""
function fillvolume!(V::VoxelBoxVolume, fillvalue::DataT) where {DataT<:Number}
    fill!(V.data, fillvalue)
    return V
end

function trim(V::VoxelBoxVolume, emptyvalue)
    emptyvalue = convert(eltype(V.data[1]), emptyvalue)
    function sliceisempty(slice, emptyvalue)
        for jx = 1:size(slice, 2)
            for ix = 1:size(slice, 1)
                if slice[ix, jx] != emptyvalue
                    return false
                end
            end
        end
        return true
    end
    xmin = 1
    for i = 1:size(V.data, 1)
        slice = view(V.data, i, :, :)
        if !sliceisempty(slice, emptyvalue)
            break
        end
        xmin = i
    end
    xmax = size(V.data, 1)
    for i = size(V.data, 1):-1:1
        slice = view(V.data, i, :, :)
        if !sliceisempty(slice, emptyvalue)
            break
        end
        xmax = i
    end
    (xmin, xmax)
end


"""
    vtkexport(theFile::String, V::VoxelBoxVolume{CoordT,DataT}) where {CoordT<:Number,DataT<:Number}

Compute.
"""
function vtkexport(theFile::String, V::VoxelBoxVolume{CoordT,DataT}) where {CoordT<:Number,DataT<:Number}
    fid=open(theFile,"w");
    if (fid==-1)
        error(["Could not open " * theFile])
        return nothing
    end
    print(fid,"# vtk DataFile Version 2.0\n");
    print(fid,"Example\n");
    print(fid,"ASCII\n");
    print(fid,"DATASET STRUCTURED_POINTS\n");
    nx, ny, nz = size(V.data)
    print(fid,"DIMENSIONS $(nx) $(ny) $(nz)\n");
    voxszx=V.boxdim[1]/nx
    voxszy=V.boxdim[2]/ny
    voxszz=V.boxdim[3]/nz
    print(fid,"SPACING $(voxszx) $(voxszy) $(voxszz)\n");
    print(fid,"ORIGIN 0 0 0\n");
    print(fid,"POINT_DATA $(nx*ny*nz)\n");
    #     bit , unsigned_char , char , unsigned_short , short , unsigned_int , int ,
    # unsigned_long , long , float , or double
    Types =  Dict{DataType, String}(Int8=>"char" , UInt8=>"unsigned_char",
                                    UInt16=>"unsigned_short", Int16=>"short",
                                    UInt32=>"unsigned_int" , Int32=>"int" ,
                                    UInt64=>"unsigned_long" , Int64=>"long" ,
                                    Float32=>"float", Float64=> "double")
    typed=Types[DataT]
    print(fid,"SCALARS volume_scalars $typed 1\n");
    print(fid,"LOOKUP_TABLE default\n");
    for j=1:(nx*ny*nz)
        print(fid,"$(V.data[j])\n");
    end
     fid=close(fid);
end
export vtkexport


end

# """
#     bar
#
# Compute.
# """
# # Merge the "under" volume into the "over" volume.  The "over" volume voxels have precedence.
# function mergevolume!(Vover::VoxelBoxVolume{CoordT,DataT},
#     Vunder::VoxelBoxVolume{CoordT,DataT},
#     emptyvalue::DataT) where {CoordT<:Number,DataT<:Number}
#     @assert (size(Vover.data) == size(Vunder.data)) "Volumes are incompatible"
#
#     nx, ny, nz = size(Vover.data)
#     voxszx=Vover.dimensions[1]/nx
#     voxszy=Vover.dimensions[2]/ny
#     voxszz=Vover.dimensions[3]/nz
#     for i = 1:nx
#         x = Vover.origin[1] + voxszx*(i+0.5)
#         for j = 1:ny
#             y = Vover.origin[2] + voxszy*(j+0.5)
#             for k = 1:nz
#                 z = Vover.origin[3] + voxszz*(k+0.5)
#                 if (Vunder.data[i,j,k]!=emptyvalue) && (Vover.data[i,j,k]==emptyvalue)
#                     Vover.data[i,j,k] = Vunder.data[i,j,k]
#                 end
#             end
#         end
#     end
#     return Vover
# end
