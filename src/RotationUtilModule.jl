"""
    RotationUtilModule

Module to handle two-dimensional and three-dimensional rotations.
"""
module RotationUtilModule

__precompile__(true)

import LinearAlgebra: norm

_I3() where {T}  = [i == j ? one(T) : zero(T) for i in 1:3, j in 1:3]

"""
    rotmat3!(Rmout::Matrix{T}, a::VT) where {T, VT}

Compute a 3D rotation matrix in-place.

`a` = array, vector, or tuple with three floating-point numbers
"""
function rotmat3!(Rmout::Matrix{T}, a::VT) where {T, VT}
    m, n = size(Rmout)
    @assert (m == n) && (m == 3) && (length(a) == m)
    na = norm(a)
    if (na == 0.0)
        for i in axes(Rmout, 1), j in axes(Rmout, 2)
            Rmout[i, j] = (i == j ? 1.0 : 0.0)
        end
    else
        # thetatilde = zeros(3,3);
        # skewmat!(thetatilde, a);
        # a = a/na;
        # ca = cos(na);
        # sa = sin(na);
        # aa = a*a';
        # copyto!(Rmout, ca * (_I3-aa) + sa/na*thetatilde + aa);
        as = (a[1] / na, a[2] / na, a[3] / na)
        ca = cos(na)
        sa = sin(na)
        oca = (1.0 - ca)
        for j in 1:3, i in 1:3
            Rmout[i, j] = oca * as[i] * as[j]
        end
        Rmout[1, 2] += sa * -as[3]
        Rmout[1, 3] += sa * as[2]
        Rmout[2, 1] += sa * as[3]
        Rmout[2, 3] += sa * -as[1]
        Rmout[3, 1] += sa * -as[2]
        Rmout[3, 2] += sa * as[1]
        Rmout[1, 1] += ca
        Rmout[2, 2] += ca
        Rmout[3, 3] += ca
    end
    return Rmout
end


"""
    rotmat3(a::VT) where {VT}

Prepare a rotation matrix from a rotation vector
"""
function rotmat3(a::VT) where {VT}
    Rmout = fill(zero(eltype(a)), 3, 3)
    return rotmat3!(Rmout, a)
end

"""
    skewmat!(S::Matrix{T}, theta::VT) where {T, VT}

Compute skew-symmetric matrix.
"""
function skewmat!(S::Matrix{T}, theta::VT) where {T, VT}
    @assert length(theta) == 3 "Input must be a 3-vector"
    # S[:,:].=[0.0     -theta[3]  theta[2];
    #          theta[3]   0.0    -theta[1];
    #         -theta[2] theta[1]  0.0];
    S[1, 1] = S[2, 2] = S[3, 3] = 0.0
    S[1, 2] = -theta[3]
    S[1, 3] = theta[2]
    S[2, 1] = theta[3]
    S[2, 3] = -theta[1]
    S[3, 1] = -theta[2]
    S[3, 2] = theta[1]
    return S
end

"""
    cross3!(
        result::AbstractVector{T1},
        theta::AbstractVector{T2},
        v::AbstractVector{T3},
    ) where {T1, T2, T3}

Compute the cross product of two vectors in three-space in place.
"""
function cross3!(
    result::AbstractVector{T1},
    theta::AbstractVector{T2},
    v::AbstractVector{T3},
) where {T1, T2, T3}
    @assert (length(theta) == 3) && (length(v) == 3) "Inputs must be 3-vectors"
    result[1] = -theta[3] * v[2] + theta[2] * v[3]
    result[2] = theta[3] * v[1] - theta[1] * v[3]
    result[3] = -theta[2] * v[1] + theta[1] * v[2]
    return result
end


"""
    cross2(theta::AbstractVector{T1}, v::AbstractVector{T2}) where {T1, T2}

Compute the cross product of two vectors in two-space.
"""
function cross2(theta::AbstractVector{T1}, v::AbstractVector{T2}) where {T1, T2}
    @assert (length(theta) == 2) && (length(v) == 2) "Inputs must be 2-vectors"
    return -theta[2] * v[1] + theta[1] * v[2]
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
