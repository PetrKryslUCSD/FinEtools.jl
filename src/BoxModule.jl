"""
    BoxModule

Module for working with bounding boxes.
"""
module BoxModule

export updatebox!, boundingbox

using FinEtools
using FinEtools.FTypesModule

inrange(rangelo::FFlt,rangehi::FFlt,x::FFlt) = ((x>=rangelo) && (x<=rangehi));

"""
    updatebox!(box::AbstractVector, x::AbstractArray)

Update a box with another location, or create a new box.

If the  `box` does not have  the correct dimensions,  it is correctly sized.

`box` = bounding box
    for 1-D `box=[minx,maxx]`, or
    for 2-D `box=[minx,maxx,miny,maxy]`, or
    for 3-D `box=[minx,maxx,miny,maxy,minz,maxz]`
    The `box` is expanded to include the
    supplied location `x`.   The variable `x`  can hold multiple points in rows.
"""
function updatebox!(box::AbstractVector, x::AbstractArray)
    sdim = size(x,2)
    if length(box) < 2*sdim
        box = Array{FFlt}(2*sdim)
        for i = 1:size(x,2)
            box[2*i-1] = Inf;
            box[2*i]   = -Inf;
        end
    end
    for j = 1:size(x,1)
        for i = 1:sdim
            box[2*i-1] = min(box[2*i-1],x[j,i]);
            box[2*i]   = max(box[2*i],x[j,i]);
        end
    end
    return box
end

function updatebox!(box::AbstractVector, x::AbstractVector)
    return updatebox!(box, reshape(x, 1, length(x)))
end

"""
    boundingbox(x::AbstractArray)

Compute the bounding box of the points in `x`.

`x` = holds points, one per row.

Returns `box` = bounding box
    for 1-D `box=[minx,maxx]`, or
    for 2-D `box=[minx,maxx,miny,maxy]`, or
    for 3-D `box=[minx,maxx,miny,maxy,minz,maxz]`
"""
function boundingbox(x::AbstractArray)
    return updatebox!(Array{FFlt}(0), x)
end

end
