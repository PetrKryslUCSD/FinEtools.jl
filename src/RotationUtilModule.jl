"""
    RotationUtilModule

Module to handle two-dimensional and three-dimensional rotations.
"""
module RotationUtilModule

using ..FTypesModule: FInt, FFlt, FCplxFlt, FFltVec, FIntVec, FFltMat, FIntMat, FMat, FVec, FDataDict
import LinearAlgebra: norm

_I3 = [i==j ? one(FFlt) : zero(FFlt) for i=1:3, j=1:3]

"""
    rotmat3!(Rmout::FFltMat, a::FFltVec)

Compute a 3D rotation matrix.
"""
function rotmat3!(Rmout::FFltMat, a::FFltVec)
    m,n = size(Rmout);
    @assert (m == n) && (m == 3) && (length(a) == m)
    na = norm(a);
    if (na == 0.)
        for i = 1:size(Rmout, 1), j = 1:size(Rmout, 2)
            Rmout[i, j] = (i == j ? 1.0 : 0.0)
        end
    else
        thetatilde = zeros(3,3);
        skewmat!(thetatilde, a);
        a = a/na;
        ca = cos(na);
        sa = sin(na);
        aa = a*a';
        copyto!(Rmout, ca * (_I3-aa) + sa/na*thetatilde + aa);
    end
end


"""
    skewmat!(S, theta)

Compute skew-symmetric matrix.
"""
function skewmat!(S, theta)
    @assert  length(theta)== 3 "Input must be a 3-vector"
    S[:,:]=[ 0.0    -theta[3] theta[2];
            theta[3]   0.0   -theta[1];
            -theta[2] theta[1]  0.0];
    return S;
end


# function cross3(theta::AbstractVector{T}, v::AbstractVector{T}) where T
#     @assert (length(theta)== 3) && (length(v)== 3) "Inputs must be 3-vectors"
#     return [ 0.0-theta[3]*v[2]+theta[2]*v[3];
#             theta[3]*v[1]+0.0-theta[1]*v[3];
#             -theta[2]*v[1]+theta[1]*v[2]+0.0];
# end

function cross3!(result::AbstractVector{T}, theta::AbstractVector{T}, v::AbstractVector{T}) where T
    @assert (length(theta)== 3) && (length(v)== 3) "Inputs must be 3-vectors"
    result[1] = -theta[3]*v[2]+theta[2]*v[3];
    result[2] = theta[3]*v[1]-theta[1]*v[3];
    result[3] = -theta[2]*v[1]+theta[1]*v[2];
    return result
end


function cross2(theta::AbstractVector{T}, v::AbstractVector{T}) where T
    @assert (length(theta)== 2) && (length(v)== 2) "Inputs must be 2-vectors"
    return -theta[2]*v[1] +theta[1]*v[2];
end

# """
#     rotmat(theta::FFltMat)
#
# Compute rotation matrix from a rotation vector or from the associated skew
# matrix.
#
# Example:
# ```
# theta=[0.8147;  0.9058;       0.1270]
# julia>        RotationUtilModule.rotmat(theta')
# 3x3 Array{Float64,2}:
#   0.631452  0.22755    0.741276
#   0.422633  0.700499  -0.57505
#  -0.650116  0.676404   0.346161
# ```
# """
# function rotmat(theta::FFltMat)
#     m,n=size(theta);
#     if (m == 3) && (n == 3)         # in this case we are assuming  the input is a skew symmetric matrix
#         # if (norm(theta' + theta) > 1e-9)
#         #     error("theta must be a skew-symmetric matrix");
#         # end
#         a=[ -theta[2,3],  theta[1,3],  -theta[1,2]];
#         # R = expm(thetatilde); # but we are not going to do that, we will use the formula below
#     else                            # we are going to assume that the input is the rotation vector
#         a= theta[:];
#     end
#     na = norm(a);
#     if (na == 0)
#         R = eye(3,3);
#     else
#         thetatilde=zeros(3,3); skewmat!(thetatilde, a);
#         a=a/na;
#         ca = cos(na);
#         sa = sin(na);
#         aa = a*a';
#         R = ca * (eye(3,3)-aa) + sa/na*thetatilde + aa;
#     end
#     return R;
# end


end
