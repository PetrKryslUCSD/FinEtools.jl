"""
Module for generating and manipulating voxel 3-D images (boxes of voxels).

This module may be helpful when working with CT scan and other medical images
in the context of finite element analysis carried out with the package
FinEtools.jl .
"""
module VoxelBoxModule

import Base.size

mutable struct VoxelBoxVolume{CoordT<:Number,DataT<:Number}
    origin::Array{CoordT,1}
    boxdim::Array{CoordT,1}
    data::Array{DataT,3}
end

function VoxelBoxVolume(::Type{CoordT}, ::Type{DataT}) where {CoordT<:Number,DataT<:Number}
    V = VoxelBoxVolume(zeros(CoordT, 3), zeros(CoordT, 3), fill(zero(DataT), 0, 0, 0));
    return V
end

function VoxelBoxVolume(::Type{DataT}, nvox::Array{Int,1}, boxdim::Array{CoordT,1}) where {CoordT<:Number,DataT<:Number}
    origin = zeros(CoordT,3)
    data = zeros(DataT,nvox...)
    V = VoxelBoxVolume(origin, boxdim, data);
    return V
end

function VoxelBoxVolume(data::Array{DataT,3}, boxdim::Array{CoordT,1})  where {CoordT<:Number,DataT<:Number}
    V = VoxelBoxVolume(CoordT,DataT);
    V.boxdim = deepcopy(boxdim)
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
size(V::VoxelBoxVolume) = size(V.data)

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

"""
    trim(V::VoxelBoxVolume, emptyvalue)

Trim off pieces of the volume that consist only of the empty value.
"""
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
        if !sliceisempty(view(V.data, i, :, :), emptyvalue)
            xmin = i; break
        end
    end
    xmax = size(V.data, 1)
    for i = size(V.data, 1):-1:1
        if !sliceisempty(view(V.data, i, :, :), emptyvalue)
            xmax = i; break
        end
    end
    ymin = 1
    for i = 1:size(V.data, 2)
        if !sliceisempty(view(V.data, :, i, :), emptyvalue)
            ymin = i; break
        end
    end
    ymax = size(V.data, 2)
    for i = size(V.data, 2):-1:1
        if !sliceisempty(view(V.data, :, i, :), emptyvalue)
            ymax = i; break
        end
    end
    zmin = 1
    for i = 1:size(V.data, 3)
        if !sliceisempty(view(V.data, :, :, i), emptyvalue)
            zmin = i; break
        end
    end
    zmax = size(V.data, 3)
    for i = size(V.data, 3):-1:1
        if !sliceisempty(view(V.data, :, :, i), emptyvalue)
            zmax = i; break
        end
    end
    if (xmin == 1) && (xmax == size(V.data, 1)) &&
        (ymin == 1) && (ymax == size(V.data, 2))  &&
        (zmin == 1) && (zmax == size(V.data, 3))
        return V
    end
    origin = V.origin .+ vec([xmin, ymin, zmin]) .* voxeldims(V)
    boxdim = vec([(xmax - xmin) (ymax - ymin) (zmax - zmin)]) .* voxeldims(V)
    return VoxelBoxVolume(origin, boxdim, V.data[xmin:xmax, ymin:ymax, zmin:zmax])
end

"""
    pad(V::VoxelBoxVolume, ipad, jpad, kpad, padvalue)

Pad voxel box with a constant value.
"""
function pad(V::VoxelBoxVolume, ipad, jpad, kpad, padvalue)
    # Adjust the number of pixels
    originalvoxeldims = voxeldims(V)
    data = zeros(eltype(V.data[1]), size(V) .+ (sum(ipad), sum(jpad), sum(kpad)))
    fill!(data, padvalue)
    data[ipad[1]+1:ipad[1]+size(V, 1), jpad[1]+1:jpad[1]+size(V, 2), kpad[1]+1:kpad[1]+size(V, 3)] = V.data
    origin = V.origin .- vec([ipad[1], jpad[1], kpad[1]]) .* voxeldims(V)
    boxdim = [x for x in size(data) .* originalvoxeldims]
    return VoxelBoxVolume(origin, boxdim, data)
end

"""
    threshold(V, threshold_value, voxel_below, voxel_above)

Threshold the data.
"""
function threshold(V, threshold_value, voxel_below, voxel_above)
    for k= 1:size(V, 3)
        for j= 1:size(V, 2)
            for i= 1:size(V, 1)
                if V.data[i,j,k] > threshold_value
                    V.data[i,j,k] = voxel_above
                else
                    V.data[i,j,k] = voxel_below
                end
            end
        end
    end
    return V
end

"""
    vtkexport(theFile::String, V::VoxelBoxVolume{CoordT,DataT}) where {CoordT<:Number,DataT<:Number}

Compute.
"""
function vtkexport(theFile::String, V::VoxelBoxVolume{CoordT,DataT}) where {CoordT<:Number,DataT<:Number}
    open(theFile,"w") do fid
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
    end
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
