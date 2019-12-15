"""
    MatrixUtilityModule

Module for general utility matrix product functions.
"""
module MatrixUtilityModule

using ..FTypesModule: FInt, FFlt, FCplxFlt, FFltVec, FIntVec, FFltMat, FIntMat, FMat, FVec, FDataDict

"""
    loc!(loc::FFltMat, ecoords::FFltMat, N::FFltMat)

Compute the location of the quadrature point.

Arguments:
`loc` = matrix of coordinates, overwritten  inside the function
`ecoords` = matrix of the node coordinates for the element.
`N` = matrix of basis function values
"""
function loc!(loc::FFltMat, ecoords::FFltMat, N::FFltMat)
	return mulCAtB!(loc, N, ecoords)
end

"""
    jac!(J::FFltMat, ecoords::FFltMat,gradNparams::FFltMat)

Compute the Jacobian matrix at the quadrature point.

Arguments:
`J` = Jacobian matrix, overwritten  inside the function
`ecoords` = matrix of the node coordinates for the element.
`gradNparams` = matrix of basis function gradients
"""
function jac!(J::FFltMat, ecoords::FFltMat, gradNparams::FFltMat)
    n = size(gradNparams, 1)
    @inbounds for j in 1:size(J, 2)
        @inbounds for i in 1:size(J, 1)
            Ja = 0.0
            @inbounds for k in 1:n
                Ja += ecoords[k, i] * gradNparams[k, j]
            end
            J[i, j] = Ja
        end
    end
    return J
end

"""
    locjac!(loc::FFltMat, J::FFltMat, ecoords::FFltMat, N::FFltMat, gradNparams::FFltMat)

Compute location and Jacobian matrix at the quadrature point.

Arguments:
`loc` = matrix of coordinates, overwritten  inside the function
`J` = Jacobian matrix, overwritten  inside the function
`ecoords` = matrix of the node coordinates for the element.
`N` = matrix of basis function values
`gradNparams` = matrix of basis function gradients
"""
function locjac!(loc::FFltMat, J::FFltMat, ecoords::FFltMat, N::FFltMat, gradNparams::FFltMat)
    return loc!(loc, ecoords, N), jac!(J, ecoords, gradNparams)
end

"""
    add_mggt_ut_only!(Ke::FFltMat, gradN::FFltMat, mult::FFlt)

Add the product `gradN*mult*gradNT*(Jac*w[j])` to the elementwise matrix `Ke`.
The argument `mult` is a scalar.
*Only upper triangle* is computed; the lower triangle is not touched.
(Use `complete_lt!` to complete the lower triangle, if needed.)

The matrix `Ke` is assumed to be suitably initialized.

The matrix `Ke` is modified.  The matrix `gradN` is not modified
inside this function.
"""
function add_mggt_ut_only!(Ke::FFltMat, gradN::FFltMat, mult::FFlt)
    Kedim = size(Ke, 1)
    @assert Kedim == size(Ke, 2) # Square matrix?
    nne, mdim = size(gradN)
    @assert nne == Kedim # compatible matrices?
    @inbounds for nx = 1:Kedim # Do: Ce  =  Ce + gradN*((Jac*w[j]))*gradN' ;
        @inbounds for px = 1:mdim
            a = (mult)*gradN[nx, px]
            @inbounds for mx = 1:nx # only the upper triangle
                Ke[mx, nx] +=  gradN[mx, px] * a
            end
        end
    end
    return true
end


"""
    add_gkgt_ut_only!(Ke::FFltMat, gradN::FFltMat, Jac_w::FFlt,
      kappa_bar::FFltMat, kappa_bargradNT::FFltMat)

Add the product `gradN*kappa_bar*gradNT*(Jac*w[j])` to the elementwise matrix `Ke`.
*Only upper triangle* is computed; the lower triangle is not touched.
(Use `complete_lt!` to complete the lower triangle, if needed.)

The matrix `Ke` is assumed to be suitably initialized.

Upon return,  the matrix `Ke` is updated.  The scratch buffer `kappa_bargradNT`
is overwritten during each call of this function. The matrices `gradN` and
`kappa_bar` are not modified inside this function.
"""
function add_gkgt_ut_only!(Ke::FFltMat, gradN::FFltMat, Jac_w::FFlt,
        kappa_bar::FFltMat, kappa_bargradNT::FFltMat)
    @assert size(Ke, 1) == size(Ke, 2)
    Kedim = size(Ke, 1)
    nne, mdim = size(gradN)
    @assert size(kappa_bar) == (mdim, mdim)
    @assert size(kappa_bargradNT) == (mdim, nne)
    # A_mul_Bt!(kappa_bargradNT, kappa_bar, gradN); # intermediate result
    @inbounds for nx = 1:nne
        @inbounds for mx = 1:mdim
	        accum::FFlt  = 0.0
	        @inbounds @simd for px = 1:mdim
	            accum += kappa_bar[mx, px]*gradN[nx, px]
	        end
	        kappa_bargradNT[mx, nx] = Jac_w*accum
        end
    end
    # Ke = Ke + gradN*(kappa_bar*(Jac*w[j]))*gradN'; only upper triangle
    @inbounds for nx = 1:Kedim
        @inbounds for mx = 1:nx # only the upper triangle
            accum::FFlt  = 0.0
            @inbounds @simd for px = 1:mdim
                accum += gradN[mx, px]*kappa_bargradNT[px, nx]
            end
            Ke[mx, nx] += accum
        end
    end
    return true
end

"""
    complete_lt!(Ke::FFltMat)

Complete the lower triangle of the elementwise matrix `Ke`.

The matrix `Ke` is modified  inside this function. The
upper-triangle  entries  are copied  across the diagonal
to the lower triangle.
"""
function complete_lt!(Ke::FFltMat)
    Kedim = size(Ke, 1)
    @assert Kedim == size(Ke, 2)
    @inbounds for nx = 1:Kedim # complete the lower triangle
        @inbounds for mx = nx+1:Kedim
            Ke[mx, nx] = Ke[nx, mx]
        end
    end
    return true
end

"""
    add_btdb_ut_only!(Ke::FFltMat, B::FFltMat, Jac_w::FFlt,
                  D::FFltMat, DB::FFltMat)

Add the product  `(B'*(D*(Jac*w[j]))*B)`, to the elementwise matrix Ke.
*Only upper triangle* is computed; the lower triangle is not touched.
(Use `complete_lt!` to complete the lower triangle, if needed.)

The matrix Ke is assumed to be suitably initialized.

The matrix Ke is modified.  The matrices B and D are not modified
inside this function. The scratch buffer DB is overwritten
during each call of this function.
"""
function add_btdb_ut_only!(Ke::FFltMat, B::FFltMat, Jac_w::FFlt, D::FFltMat, DB::FFltMat)
    @assert size(Ke, 1) == size(Ke, 2)
    @assert size(B, 1) == size(D, 1)
    nstr, Kedim = size(B)
    @assert size(D) == (nstr, nstr)
    @assert size(DB) == (nstr, Kedim)
    # A_mul_B!(DB, D, B)    # intermediate product
    @inbounds for nx = 1:Kedim
        @inbounds for mx = 1:nstr
            accum::FFlt  = 0.0
            @inbounds for px = 1:nstr
                accum += D[mx, px]*B[px, nx]
            end
            DB[mx, nx] = Jac_w*accum
        end
    end
    #  Ke = Ke + (B'*(D*(Jac*w[j]))*B); only the upper triangle
    @inbounds for nx = 1:Kedim
        @inbounds for mx = 1:nx # only the upper triangle
	        accum::FFlt  = 0.0
	        @inbounds for px = 1:nstr
	            accum += B[px, mx]*DB[px, nx]
	        end
	        Ke[mx, nx] += accum
        end
    end
    return true
end

"""
    add_btsigma!(Fe::FFltVec, B::FFltMat, coefficient::FFlt, sigma::FFltVec)

Add the product  `B'*(sigma*coefficient)`, to the elementwise vector `Fe`.

The vector `Fe` is assumed to be suitably initialized.

The vector `Fe` is modified.  The vector `sigma` is not modified
inside this function. 
"""
function add_btsigma!(Fe::FFltVec, B::FFltMat, coefficient::FFlt, sigma::FFltVec)
    @assert size(B, 1) == length(sigma)
    nstr, Kedim = size(B)
    @inbounds for nx = 1:Kedim
	    accum::FFlt  = 0.0
	    @inbounds for px = 1:nstr
	    	accum += sigma[px]*B[px, nx]
		end
		Fe[nx] += coefficient*accum
	end
    return true
end

"""
    add_nnt_ut_only!(Ke::FFltMat, Nn::FFltMat, Jac_w_coeff::FFlt)

Add the product  `Nn*(Nn'*(coeff*(Jac*w(j)))`, to the elementwise matrix `Ke`.
*Only upper triangle* is computed; the lower triangle is not touched.

The matrix `Ke` is assumed to be suitably initialized.

The matrix `Ke` is modified.  The matrix `Nn` is not modified
inside this function.
"""
function add_nnt_ut_only!(Ke::FMat{T}, N::FFltMat, Jac_w_coeff::T) where {T<:Number}
    Kedim = size(N, 1)
    @assert Kedim == size(Ke, 2)
    @inbounds for nx = 1:Kedim
        a = (Jac_w_coeff)*N[nx]
        @inbounds for mx = 1:nx # only the upper triangle
            Ke[mx,nx] += N[mx] * a
        end
    end
    return true
end

function add_nnt_ut_only!(Ke::FFltMat, N::FFltVec, Jac_w_coeff::FFlt)
    return add_nnt_ut_only!(Ke, reshape(N, length(N), 1), Jac_w_coeff)
end


"""
    add_btv!(elvec::FFltVec, B::FFltMat, sig::FFltMat, Jac_w_coeff::FFlt)

Add product of the strain-displacement matrix transpose times the stress vector.

Note:  the coefficient `Jac_w_coeff` will be typically  NEGATIVE.

The argument `elvec` needs to be suitably  initialized before the first call
(filled with zeros, for instance),  and it is updated upon return.
The arguments `B`, `sig` are not modified.
"""
function add_btv!(elvec::FFltVec, B::FFltMat, sig::FFltVec, Jac_w_coeff::FFlt)
    @assert size(B, 1) == length(sig)
    @assert size(B, 2) == length(elvec)
    Kedim = length(elvec)
    nsig = length(sig)
    @inbounds for mx = 1:Kedim
        @inbounds for px = 1:nsig
            elvec[mx] += B[px, mx]*Jac_w_coeff*sig[px]
        end
    end
    return true
end

"""
    adjugate3!(B, A)

Compute the adjugate matrix of `A`.
"""
function adjugate3!(B, A)
    # % A = rand(3)
    # % det(A)*inv(A)-B
    # B =0*A;
    B[1,1] =(A[5]*A[9]-A[8]*A[6]);
    B[1,2] =-(A[4]*A[9]-A[7]*A[6]);
    B[1,3] =(A[4]*A[8]-A[7]*A[5]);

    B[2,1] =-(A[2]*A[9]-A[8]*A[3]);
    B[2,2] =(A[1]*A[9]-A[7]*A[3]);
    B[2,3] =-(A[1]*A[8]-A[7]*A[2]);

    B[3,1] =(A[2]*A[6]-A[5]*A[3]);
    B[3,2] =-(A[1]*A[6]-A[4]*A[3]);
    B[3,3] =(A[1]*A[5]-A[4]*A[2]);
    return B
end

"""
    symmetrize!(a)

Make the matrix on input symmetric.
"""
function symmetrize!(a)
	@assert size(a, 1) == size(a, 2)
    @inbounds for c in 1:size(a, 2)
    	@inbounds for r in c:size(a, 1)
    		a[c, r] =  a[r, c] += a[c, r]
    	end
    end
    a .= 0.5 * a
end

"""
    mulCAB!(::Val{3}, C, A, B)

Compute the product of matrices `C = A * B`

"""
function mulCAB!(::Val{3}, C, A, B)
	#define MAT3DUTILMULTAB3X3(C, A, B)                                  
	C[1, 1]  = (A[1, 1]*B[1, 1] + A[1, 2]*B[2, 1] + A[1, 3]*B[3, 1]);
	C[1, 2]  = (A[1, 1]*B[1, 2] + A[1, 2]*B[2, 2] + A[1, 3]*B[3, 2]);
	C[1, 3]  = (A[1, 1]*B[1, 3] + A[1, 2]*B[2, 3] + A[1, 3]*B[3, 3]);
	C[2, 1]  = (A[2, 1]*B[1, 1] + A[2, 2]*B[2, 1] + A[2, 3]*B[3, 1]);
	C[2, 2]  = (A[2, 1]*B[1, 2] + A[2, 2]*B[2, 2] + A[2, 3]*B[3, 2]);
	C[2, 3]  = (A[2, 1]*B[1, 3] + A[2, 2]*B[2, 3] + A[2, 3]*B[3, 3]);
	C[3, 1]  = (A[3, 1]*B[1, 1] + A[3, 2]*B[2, 1] + A[3, 3]*B[3, 1]);
	C[3, 2]  = (A[3, 1]*B[1, 2] + A[3, 2]*B[2, 2] + A[3, 3]*B[3, 2]);
	C[3, 3]  = (A[3, 1]*B[1, 3] + A[3, 2]*B[2, 3] + A[3, 3]*B[3, 3]);
	return C
end

"""
    mulCAtB!(::Val{3}, C, A, B)

Compute the product of matrices `C = Transpose(A) * B`

"""
function mulCAtB!(::Val{3}, C, A, B)
	#define MAT3DUTILMULTATB3X3(C, A, B)                                 
	C[1, 1]  = (A[1, 1]*B[1, 1] + A[2, 1]*B[2, 1] + A[3, 1]*B[3, 1]);
	C[1, 2]  = (A[1, 1]*B[1, 2] + A[2, 1]*B[2, 2] + A[3, 1]*B[3, 2]);
	C[1, 3]  = (A[1, 1]*B[1, 3] + A[2, 1]*B[2, 3] + A[3, 1]*B[3, 3]);
	C[2, 1]  = (A[1, 2]*B[1, 1] + A[2, 2]*B[2, 1] + A[3, 2]*B[3, 1]);
	C[2, 2]  = (A[1, 2]*B[1, 2] + A[2, 2]*B[2, 2] + A[3, 2]*B[3, 2]);
	C[2, 3]  = (A[1, 2]*B[1, 3] + A[2, 2]*B[2, 3] + A[3, 2]*B[3, 3]);
	C[3, 1]  = (A[1, 3]*B[1, 1] + A[2, 3]*B[2, 1] + A[3, 3]*B[3, 1]);
	C[3, 2]  = (A[1, 3]*B[1, 2] + A[2, 3]*B[2, 2] + A[3, 3]*B[3, 2]);
	C[3, 3]  = (A[1, 3]*B[1, 3] + A[2, 3]*B[2, 3] + A[3, 3]*B[3, 3]);
	return C
end

"""
    mulCABt!(::Val{3}, C, A, B)

Compute the product of matrices `C = A * Transpose(B)`

"""
function mulCABt!(::Val{3}, C, A, B)
	#define MAT3DUTILMULTABT3X3(C, A, B)                                 
	C[1, 1]  = (A[1, 1]*B[1, 1] + A[1, 2]*B[1, 2] + A[1, 3]*B[1, 3]);
	C[1, 2]  = (A[1, 1]*B[2, 1] + A[1, 2]*B[2, 2] + A[1, 3]*B[2, 3]);
	C[1, 3]  = (A[1, 1]*B[3, 1] + A[1, 2]*B[3, 2] + A[1, 3]*B[3, 3]);
	C[2, 1]  = (A[2, 1]*B[1, 1] + A[2, 2]*B[1, 2] + A[2, 3]*B[1, 3]);
	C[2, 2]  = (A[2, 1]*B[2, 1] + A[2, 2]*B[2, 2] + A[2, 3]*B[2, 3]);
	C[2, 3]  = (A[2, 1]*B[3, 1] + A[2, 2]*B[3, 2] + A[2, 3]*B[3, 3]);
	C[3, 1]  = (A[3, 1]*B[1, 1] + A[3, 2]*B[1, 2] + A[3, 3]*B[1, 3]);
	C[3, 2]  = (A[3, 1]*B[2, 1] + A[3, 2]*B[2, 2] + A[3, 3]*B[2, 3]);
	C[3, 3]  = (A[3, 1]*B[3, 1] + A[3, 2]*B[3, 2] + A[3, 3]*B[3, 3]);
	return C
end

"""
    mulCAtB!(C::FFltMat, A::FFltMat, B::FFltMat)

Compute the matrix `C = A' * B`
"""
function mulCAtB!(C::FFltMat, A::FFltMat, B::FFltMat)
    n = size(B, 1)
    @assert size(C, 1) == size(A, 2)
    @assert size(C, 2) == size(B, 2)
    @assert size(A, 1) == size(B, 1)
    @inbounds for j in 1:size(C, 2)
        @inbounds for i in 1:size(C, 1)
            Ca = 0.0
            @inbounds for k in 1:n
                Ca += A[k, i] * B[k, j]
            end
            C[i, j] = Ca
        end
    end
    return C
end

"""
    mulCAB!(C::FFltMat, A::FFltMat, B::FFltMat)

Compute the matrix `C = A * B`
"""
function mulCAB!(C::FFltMat, A::FFltMat, B::FFltMat)
    n = size(B, 1)
    @assert size(C, 1) == size(A, 1)
    @assert size(C, 2) == size(B, 2)
    @assert size(A, 2) == size(B, 1)
    @inbounds for j in 1:size(C, 2)
        @inbounds for i in 1:size(C, 1)
            Ca = 0.0
            @inbounds for k in 1:n
                Ca += A[i, k] * B[k, j]
            end
            C[i, j] = Ca
        end
    end
    return C
end

"""
    mulCABt!(C::FFltMat, A::FFltMat, B::FFltMat)

Compute the matrix `C = A * B'`
"""
function mulCABt!(C::FFltMat, A::FFltMat, B::FFltMat)
    n = size(B, 2)
    @assert size(C, 1) == size(A, 1)
    @assert size(C, 2) == size(B, 1)
    @assert size(A, 2) == size(B, 2)
    @inbounds for j in 1:size(C, 2)
        @inbounds for i in 1:size(C, 1)
            Ca = 0.0
            @inbounds for k in 1:n
                Ca += A[i, k] * B[j, k]
            end
            C[i, j] = Ca
        end
    end
    return C
end

"""
    detC(::Val{3}, C::FFltMat)

Compute determinant of `C`.
"""
function detC(::Val{3}, C::FFltMat)
	#define MAT3DUTIL_DET_3X3(T3X3) 
	return (  
		(C[1, 1] * C[2, 2] * C[3, 3]) + 
		(C[1, 2] * C[2, 3] * C[3, 1]) +
		(C[1, 3] * C[2, 1] * C[3, 2]) - 
		(C[1, 3] * C[2, 2] * C[3, 1]) - 
		(C[1, 2] * C[2, 1] * C[3, 3]) - 
		(C[1, 1] * C[2, 3] * C[3, 2]))
end

end
