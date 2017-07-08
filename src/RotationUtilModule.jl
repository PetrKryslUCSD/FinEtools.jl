"""
    RotationUtilModule

Module to handle two-dimensional and three-dimensional rotations.
"""
module RotationUtilModule

export rotmat,  rotmat3!, skewmat!

using FinEtools.FTypesModule

"""
    rotmat(theta::FFltMat)

Compute rotation matrix from a rotation vector or from the associated skew
matrix.

Example:
```
theta=[0.8147;  0.9058;       0.1270]
julia>        RotationUtilModule.rotmat(theta')
3x3 Array{Float64,2}:
  0.631452  0.22755    0.741276
  0.422633  0.700499  -0.57505
 -0.650116  0.676404   0.346161
```
"""
function rotmat(theta::FFltMat)
    m,n=size(theta);
    if (m == 3) && (n == 3)         # in this case we are assuming  the input is a skew symmetric matrix
        # if (norm(theta' + theta) > 1e-9)
        #     error("theta must be a skew-symmetric matrix");
        # end
        a=[ -theta[2,3],  theta[1,3],  -theta[1,2]];
        # R = expm(thetatilde); # but we are not going to do that, we will use the formula below
    else                            # we are going to assume that the input is the rotation vector
        a= theta[:];
    end
    na = norm(a);
    if (na == 0)
        R = eye(3,3);
    else
        thetatilde=zeros(3,3); skewmat!(thetatilde, a);
        a=a/na;
        ca = cos(na);
        sa = sin(na);
        aa = a*a';
        R = ca * (eye(3,3)-aa) + sa/na*thetatilde + aa;
    end
    return R;
end


"""
    rotmat3!(Rmout::FFltMat, a::FFltVec)

Compute a 3D rotation matrix.
"""
function rotmat3!(Rmout::FFltMat, a::FFltVec)
    m,n = size(Rmout);
    @assert (m == n) && (m == 3) && (length(a) == m)
    na = norm(a);
    if (na == 0.)
        copy!(Rmout, eye(m, n));
    else
        thetatilde = zeros(3,3);
        skewmat!(thetatilde, a);
        a = a/na;
        ca = cos(na);
        sa = sin(na);
        aa = a*a';
        copy!(Rmout, ca * (eye(3,3)-aa) + sa/na*thetatilde + aa);
    end
end


"""
    skewmat!(S, theta)

Compute skew-symmetric matrix.
"""
function skewmat!(S, theta)
    if length(theta)== 3
        S[:,:]=[ 0.0    -theta[3] theta[2];
                theta[3]   0.0   -theta[1];
                -theta[2] theta[1]  0.0];
    elseif  length(theta) == 2
        S[:,:]=[-theta[2] theta[1]];
    else
        error("theta must be a 3-vector");
    end
    return S;
end


function cross(theta,v)
    if length(theta)== 3
        n    =[ 0.0-theta[3]*v[2]+theta[2]*v[3];
                theta[3]*v[1]+0.0-theta[1]*v[3];
                -theta[2]*v[1]+theta[1]*v[2]+0.0];
    elseif  length(theta) == 2
        n    =-theta[2]*v[1] +theta[1]*v[2];
    else
        error("theta must be a 3-vector");
    end
    return n;
end

end
