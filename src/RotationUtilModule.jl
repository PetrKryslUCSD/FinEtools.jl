"""
    RotationUtilModule

Module to handle two-dimensional and three-dimensional rotations.
"""
module RotationUtilModule

__precompile__(true)

import LinearAlgebra: norm, I

function _I3(one::T = T(1.0), zero::T = T(0.0)) where {T}
    [i == j ? one : zero for i = 1:3, j = 1:3]
end

"""
    rotmat3!(Rmout::Matrix{T}, a::AbstractVector{VT}) where {T,VT}

Compute a 3D rotation matrix in-place.

`a` = array, vector, or tuple with three floating-point numbers
"""
function rotmat3!(Rmout::Matrix{T}, a::AbstractVector{VT}) where {T,VT}
    t = norm(a)
    if t == zero(VT)
        Rmout .= I(3)
    else
        kx, ky, kz = a ./ t
        c = cos(t); s = sin(t)
        Rmout .= [
                c + (1 - c) * kx^2           (1 - c) * kx * ky - s * kz   (1 - c) * kx * kz + s * ky
                (1 - c) * kx * ky + s * kz      c + (1 - c) * ky^2        (1 - c) * ky * kz - s * kx
                (1 - c) * kx * kz - s * ky   (1 - c) * ky * kz + s * kx      c + (1 - c) * kz^2
            ]
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
function skewmat!(S::Matrix{T}, theta::VT) where {T,VT}
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
) where {T1,T2,T3}
    @assert (length(theta) == 3) && (length(v) == 3) "Inputs must be 3-vectors"
    result[1] = -theta[3] * v[2] + theta[2] * v[3]
    result[2] = theta[3] * v[1] - theta[1] * v[3]
    result[3] = -theta[2] * v[1] + theta[1] * v[2]
    return result
end

"""
    cross3!(
        result::AbstractVector{T1},
        theta::Union{AbstractVector{T2}, Tuple{T2, T2, T2}},
        v::Union{AbstractVector{T3}, Tuple{T3, T3, T3}}
    ) where {T1, T2, T3}

Compute the cross product of two vectors in three-space in place.
"""
function cross3!(
    result::AbstractVector{T1},
    theta::Union{AbstractVector{T2},Tuple{T2,T2,T2}},
    v::Union{AbstractVector{T3},Tuple{T3,T3,T3}},
) where {T1,T2,T3}
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
function cross2(theta::AbstractVector{T1}, v::AbstractVector{T2}) where {T1,T2}
    @assert (length(theta) == 2) && (length(v) == 2) "Inputs must be 2-vectors"
    return -theta[2] * v[1] + theta[1] * v[2]
end

end
