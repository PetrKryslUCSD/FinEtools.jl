"""
    MatrixUtilityModule

Module for general utility matrix product functions.
"""
module MatrixUtilityModule

__precompile__(true)

using DelimitedFiles
using SparseArrays
using LoopVectorization

"""
    loc!(loc::Matrix{T}, ecoords::Matrix{T}, N::Matrix{T}) where {T}

Compute the location of the quadrature point.

Arguments:
`loc` = matrix of coordinates, overwritten  inside the function
`ecoords` = matrix of the node coordinates for the element.
`N` = matrix of basis function values
"""
function loc!(loc::Matrix{T}, ecoords::Matrix{T}, N::Matrix{T}) where {T}
    return mulCAtB!(loc, N, ecoords)
end

"""
    jac!(J::Matrix{T}, ecoords::Matrix{T}, gradNparams::Matrix{T}) where {T}

Compute the Jacobian matrix at the quadrature point.

Arguments:
`J` = Jacobian matrix, overwritten  inside the function
`ecoords` = matrix of the node coordinates for the element.
`gradNparams` = matrix of basis function gradients
"""
function jac!(J::Matrix{T}, ecoords::Matrix{T}, gradNparams::Matrix{T}) where {T}
    return mulCAtB!(J, ecoords, gradNparams)
end

"""
    locjac!(
        loc::Matrix{T},
        J::Matrix{T},
        ecoords::Matrix{T},
        N::Matrix{T},
        gradNparams::Matrix{T},
    ) where {T}

Compute location and Jacobian matrix at the quadrature point.

Arguments:
`loc` = matrix of coordinates, overwritten  inside the function
`J` = Jacobian matrix, overwritten  inside the function
`ecoords` = matrix of the node coordinates for the element.
`N` = matrix of basis function values
`gradNparams` = matrix of basis function gradients
"""
function locjac!(
    loc::Matrix{T},
    J::Matrix{T},
    ecoords::Matrix{T},
    N::Matrix{T},
    gradNparams::Matrix{T},
) where {T}
    return loc!(loc, ecoords, N), jac!(J, ecoords, gradNparams)
end

"""
    add_mggt_ut_only!(Ke::Matrix{T}, gradN::Matrix{T}, mult) where {T}

Add the product `gradN*mult*gradNT` to the matrix `Ke`.

The argument `mult` is a scalar.
*Only upper triangle* is computed; the lower triangle is not touched.
(Use `complete_lt!` to complete the lower triangle, if needed.)

The matrix `Ke` is assumed to be suitably initialized.

The matrix `Ke` is modified.  The matrix `gradN` is not modified
inside this function.
"""
function add_mggt_ut_only!(Ke::Matrix{MT}, gradN::Matrix{T}, mult) where {MT,T}
    Kedim = size(Ke, 1)
    @assert Kedim == size(Ke, 2) # Square matrix?
    nne, mdim = size(gradN)
    @assert nne == Kedim # compatible matrices?
    @inbounds for nx = 1:Kedim # Do: Ce  =  Ce + gradN*((Jac*w[j]))*gradN' ;
        @inbounds for px = 1:mdim
            a = (mult) * gradN[nx, px]
            @inbounds for mx = 1:nx # only the upper triangle
                Ke[mx, nx] += gradN[mx, px] * a
            end
        end
    end
    return true
end

"""
    add_gkgt_ut_only!(
        Ke::Matrix{T},
        gradN::Matrix{T},
        Jac_w::T,
        kappa_bar::Matrix{T},
        kappa_bargradNT::Matrix{T},
    ) where {T}

Add the product `gradN*kappa_bar*gradNT*(Jac*w[j])` to the matrix `Ke`.

*Only upper triangle* is computed; the lower triangle is not touched.
(Use `complete_lt!` to complete the lower triangle, if needed.)

The matrix `Ke` is assumed to be suitably initialized.

Upon return,  the matrix `Ke` is updated.  The scratch buffer `kappa_bargradNT`
is overwritten during each call of this function. The matrices `gradN` and
`kappa_bar` are not modified inside this function.
"""
function add_gkgt_ut_only!(
    Ke::Matrix{T},
    gradN::Matrix{T},
    Jac_w::T,
    kappa_bar::Matrix{T},
    kappa_bargradNT::Matrix{T},
) where {T}
    @assert size(Ke, 1) == size(Ke, 2)
    Kedim = size(Ke, 1)
    nne, mdim = size(gradN)
    @assert size(kappa_bar) == (mdim, mdim)
    @assert size(kappa_bargradNT) == (mdim, nne)
    # A_mul_Bt!(kappa_bargradNT, kappa_bar, gradN); # intermediate result
    @inbounds for nx = 1:nne
        @inbounds for mx = 1:mdim
            accum = T(0.0)
            @inbounds @simd for px = 1:mdim
                accum += kappa_bar[mx, px] * gradN[nx, px]
            end
            kappa_bargradNT[mx, nx] = Jac_w * accum
        end
    end
    # Ke = Ke + gradN*(kappa_bar*(Jac*w[j]))*gradN'; only upper triangle
    @inbounds for nx = 1:Kedim
        @inbounds for mx = 1:nx # only the upper triangle
            accum = T(0.0)
            @inbounds @simd for px = 1:mdim
                accum += gradN[mx, px] * kappa_bargradNT[px, nx]
            end
            Ke[mx, nx] += accum
        end
    end
    return true
end

"""
    complete_lt!(Ke::Matrix{T}) where {T}

Complete the lower triangle of the elementwise matrix `Ke`.

The matrix `Ke` is modified  inside this function. The
upper-triangle  entries  are copied  across the diagonal
to the lower triangle.
"""
function complete_lt!(Ke::Matrix{T}) where {T}
    Kedim = size(Ke, 1)
    @assert Kedim == size(Ke, 2)
    @inbounds for nx = 1:Kedim # complete the lower triangle
        @inbounds for mx = (nx+1):Kedim
            Ke[mx, nx] = Ke[nx, mx]
        end
    end
    return true
end

"""
    add_btdb_ut_only!(Ke::Matrix{T}, B::Matrix{T}, Jac_w::T, D::Matrix{T}, DB::Matrix{T}) where {T}

Add the product  `(B'*(D*(Jac*w[j]))*B)`, to the matrix Ke.

*Only upper triangle* is computed; the lower triangle is not touched.
(Use `complete_lt!` to complete the lower triangle, if needed.)

The matrix Ke is assumed to be suitably initialized.

The matrix Ke is modified.  The matrices B and D are not modified
inside this function. The scratch buffer DB is overwritten
during each call of this function.
"""
function add_btdb_ut_only!(
    Ke::Matrix{T},
    B::Matrix{T},
    Jac_w::T,
    D::Matrix{T},
    DB::Matrix{T},
) where {T}
    @assert size(Ke, 1) == size(Ke, 2)
    @assert size(B, 1) == size(D, 1)
    nstr, Kedim = size(B)
    @assert size(D) == (nstr, nstr)
    @assert size(DB) == (nstr, Kedim)
    # A_mul_B!(DB, D, B)    # intermediate product
    @inbounds for nx = 1:Kedim
        @inbounds for mx = 1:nstr
            accum = T(0.0)
            @inbounds for px = 1:nstr
                accum += D[mx, px] * B[px, nx]
            end
            DB[mx, nx] = Jac_w * accum
        end
    end
    #  Ke = Ke + (B'*(D*(Jac*w[j]))*B); only the upper triangle
    @inbounds for nx = 1:Kedim
        @inbounds for mx = 1:nx # only the upper triangle
            accum = T(0.0)
            @inbounds for px = 1:nstr
                accum += B[px, mx] * DB[px, nx]
            end
            Ke[mx, nx] += accum
        end
    end
    return true
end

"""
    add_b1tdb2!(
        Ke::Matrix{T},
        B1::Matrix{T},
        B2::Matrix{T},
        Jac_w::T,
        D::Matrix{T},
        DB2::Matrix{T},
    ) where {T}

Add the product  `(B1'*(D*(Jac_w))*B2)`, to the matrix `Ke`.

The matrix `Ke` is assumed to be suitably initialized: the results of this
computation are **added**. The matrix `Ke` may be rectangular.

The matrix `D` may be rectangular.

The matrix Ke is modified.  The matrices `B1`, `B2`, and `D` are not modified
inside this function. The scratch buffer `DB` is overwritten
during each call of this function.
"""
function add_b1tdb2!(
    Ke::Matrix{T},
    B1::Matrix{T},
    B2::Matrix{T},
    Jac_w::T,
    D::Matrix{T},
    DB2::Matrix{T},
) where {T}
    Kedim1, Kedim2 = size(Ke)
    Ddim1, Ddim2 = size(D)
    @assert size(B1) == (Ddim1, Kedim1)
    @assert size(B2) == (Ddim2, Kedim2)
    @assert size(DB2) == (Ddim1, Kedim2)
    # A_mul_B!(DB, D, B)    # intermediate product
    @inbounds for nx = 1:Kedim2
        for mx = 1:Ddim1
            accum = T(0.0)
            for px = 1:Ddim2
                accum += D[mx, px] * B2[px, nx]
            end
            DB2[mx, nx] = Jac_w * accum
        end
    end
    #  Ke = Ke + (B1'*(D*(Jac_w))*B2)
    @inbounds for nx = 1:Kedim2
        for mx = 1:Kedim1 # only the upper triangle
            accum = T(0.0)
            for px = 1:Ddim1
                accum += B1[px, mx] * DB2[px, nx]
            end
            Ke[mx, nx] += accum
        end
    end
    return true
end

"""
    add_btsigma!(Fe::Vector{T}, B::Matrix{T}, coefficient::T, sigma::Vector{T}) where {T}

Add the product  `B'*(sigma*coefficient)`, to the elementwise vector `Fe`.

The vector `Fe` is assumed to be suitably initialized.

The vector `Fe` is modified.  The vector `sigma` is not modified
inside this function. 
"""
function add_btsigma!(
    Fe::Vector{T},
    B::Matrix{T},
    coefficient::T,
    sigma::Vector{T},
) where {T}
    @assert size(B, 1) == length(sigma)
    nstr, Kedim = size(B)
    @inbounds for nx = 1:Kedim
        accum = T(0.0)
        @inbounds for px = 1:nstr
            accum += sigma[px] * B[px, nx]
        end
        Fe[nx] += coefficient * accum
    end
    return true
end

"""
    add_nnt_ut_only!(Ke::Matrix{T}, N::Matrix{T}, Jac_w_coeff::T) where {T<:Number}

Add the product  `Nn*(Nn'*(coeff*(Jac*w(j)))`, to the matrix `Ke`.

*Only the upper triangle* is computed; the lower triangle is not touched.

The matrix `Ke` is assumed to be suitably initialized. The 
matrix `Nn` has a single column.

The matrix `Ke` is modified.  The matrix `Nn` is not modified
inside this function.
"""
function add_nnt_ut_only!(Ke::Matrix{T}, N::Matrix{T}, Jac_w_coeff::T) where {T<:Number}
    Kedim = size(N, 1)
    @assert Kedim == size(Ke, 2)
    @inbounds for nx = 1:Kedim
        a = (Jac_w_coeff) * N[nx]
        @inbounds for mx = 1:nx # only the upper triangle
            Ke[mx, nx] += N[mx] * a
        end
    end
    return true
end

"""
    add_n1n2t!(Ke::Matrix{T}, N1::Matrix{T}, N2::Matrix{T}, Jac_w_coeff::T) where {T<:Number}

Add the product  `N1*(N2'*(coeff*(Jac*w(j)))`, to the matrix `Ke`.

The matrix `Ke` is assumed to be suitably initialized. The
matrices `N1` and `N2` have a single column each.

The matrix `Ke` is modified.  The matrix `N1` and `N2` are  not modified inside
this function.
"""
function add_n1n2t!(
    Ke::Matrix{T},
    N1::VecOrMat{T},
    N2::VecOrMat{T},
    Jac_w_coeff::T,
) where {T<:Number}
    @assert length(N1) == size(Ke, 1)
    @assert length(N2) == size(Ke, 2)
    @inbounds for nx = 1:size(Ke, 2)
        a = (Jac_w_coeff) * N2[nx]
        @inbounds for mx = 1:size(Ke, 1) # only the upper triangle
            Ke[mx, nx] += N1[mx] * a
        end
    end
    return true
end

function add_nnt_ut_only!(Ke::Matrix{T}, N::Vector{T}, Jac_w_coeff::T) where {T}
    return add_nnt_ut_only!(Ke, reshape(N, length(N), 1), Jac_w_coeff)
end

"""
    adjugate3!(B, A)

Compute the adjugate matrix of 3x3 matrix `A`.
"""
function adjugate3!(B, A)
    # % A = rand(3)
    # % det(A)*inv(A)-B
    # B =0*A;
    B[1, 1] = (A[5] * A[9] - A[8] * A[6])
    B[1, 2] = -(A[4] * A[9] - A[7] * A[6])
    B[1, 3] = (A[4] * A[8] - A[7] * A[5])

    B[2, 1] = -(A[2] * A[9] - A[8] * A[3])
    B[2, 2] = (A[1] * A[9] - A[7] * A[3])
    B[2, 3] = -(A[1] * A[8] - A[7] * A[2])

    B[3, 1] = (A[2] * A[6] - A[5] * A[3])
    B[3, 2] = -(A[1] * A[6] - A[4] * A[3])
    B[3, 3] = (A[1] * A[5] - A[4] * A[2])
    return B
end

"""
    symmetrize!(a)

Make the matrix on input symmetric.

The operation is in-place.
"""
function symmetrize!(a)
    @assert size(a, 1) == size(a, 2)
    @inbounds for c in axes(a, 2)
        @inbounds for r = c:size(a, 1)
            a[c, r] = a[r, c] += a[c, r]
        end
    end
    a .= 0.5 * a
end

"""
    mulCAB!(::Val{3}, C, A, B)

Compute the product of 3X3 matrices `C = A * B`

"""
function mulCAB!(::Val{3}, C, A, B)
    #define MAT3DUTILMULTAB3X3(C, A, B)                                  
    C[1, 1] = (A[1, 1] * B[1, 1] + A[1, 2] * B[2, 1] + A[1, 3] * B[3, 1])
    C[1, 2] = (A[1, 1] * B[1, 2] + A[1, 2] * B[2, 2] + A[1, 3] * B[3, 2])
    C[1, 3] = (A[1, 1] * B[1, 3] + A[1, 2] * B[2, 3] + A[1, 3] * B[3, 3])
    C[2, 1] = (A[2, 1] * B[1, 1] + A[2, 2] * B[2, 1] + A[2, 3] * B[3, 1])
    C[2, 2] = (A[2, 1] * B[1, 2] + A[2, 2] * B[2, 2] + A[2, 3] * B[3, 2])
    C[2, 3] = (A[2, 1] * B[1, 3] + A[2, 2] * B[2, 3] + A[2, 3] * B[3, 3])
    C[3, 1] = (A[3, 1] * B[1, 1] + A[3, 2] * B[2, 1] + A[3, 3] * B[3, 1])
    C[3, 2] = (A[3, 1] * B[1, 2] + A[3, 2] * B[2, 2] + A[3, 3] * B[3, 2])
    C[3, 3] = (A[3, 1] * B[1, 3] + A[3, 2] * B[2, 3] + A[3, 3] * B[3, 3])
    return C
end

"""
    mulCAtB!(::Val{3}, C, A, B)

Compute the product of 3X3 matrices `C = Transpose(A) * B`

"""
function mulCAtB!(::Val{3}, C, A, B)
    #define MAT3DUTILMULTATB3X3(C, A, B)                                 
    C[1, 1] = (A[1, 1] * B[1, 1] + A[2, 1] * B[2, 1] + A[3, 1] * B[3, 1])
    C[1, 2] = (A[1, 1] * B[1, 2] + A[2, 1] * B[2, 2] + A[3, 1] * B[3, 2])
    C[1, 3] = (A[1, 1] * B[1, 3] + A[2, 1] * B[2, 3] + A[3, 1] * B[3, 3])
    C[2, 1] = (A[1, 2] * B[1, 1] + A[2, 2] * B[2, 1] + A[3, 2] * B[3, 1])
    C[2, 2] = (A[1, 2] * B[1, 2] + A[2, 2] * B[2, 2] + A[3, 2] * B[3, 2])
    C[2, 3] = (A[1, 2] * B[1, 3] + A[2, 2] * B[2, 3] + A[3, 2] * B[3, 3])
    C[3, 1] = (A[1, 3] * B[1, 1] + A[2, 3] * B[2, 1] + A[3, 3] * B[3, 1])
    C[3, 2] = (A[1, 3] * B[1, 2] + A[2, 3] * B[2, 2] + A[3, 3] * B[3, 2])
    C[3, 3] = (A[1, 3] * B[1, 3] + A[2, 3] * B[2, 3] + A[3, 3] * B[3, 3])
    return C
end

"""
    mulCABt!(::Val{3}, C, A, B)

Compute the product of 3X3 matrices `C = A * Transpose(B)`

"""
function mulCABt!(::Val{3}, C, A, B)
    #define MAT3DUTILMULTABT3X3(C, A, B)                                 
    C[1, 1] = (A[1, 1] * B[1, 1] + A[1, 2] * B[1, 2] + A[1, 3] * B[1, 3])
    C[1, 2] = (A[1, 1] * B[2, 1] + A[1, 2] * B[2, 2] + A[1, 3] * B[2, 3])
    C[1, 3] = (A[1, 1] * B[3, 1] + A[1, 2] * B[3, 2] + A[1, 3] * B[3, 3])
    C[2, 1] = (A[2, 1] * B[1, 1] + A[2, 2] * B[1, 2] + A[2, 3] * B[1, 3])
    C[2, 2] = (A[2, 1] * B[2, 1] + A[2, 2] * B[2, 2] + A[2, 3] * B[2, 3])
    C[2, 3] = (A[2, 1] * B[3, 1] + A[2, 2] * B[3, 2] + A[2, 3] * B[3, 3])
    C[3, 1] = (A[3, 1] * B[1, 1] + A[3, 2] * B[1, 2] + A[3, 3] * B[1, 3])
    C[3, 2] = (A[3, 1] * B[2, 1] + A[3, 2] * B[2, 2] + A[3, 3] * B[2, 3])
    C[3, 3] = (A[3, 1] * B[3, 1] + A[3, 2] * B[3, 2] + A[3, 3] * B[3, 3])
    return C
end

"""
    mulCAtB!(C, A, B)

Compute the matrix `C = A' * B`

The use of BLAS is purposefully avoided in order to eliminate contentions of
multi-threaded execution of the library code with the user-level threads.
"""
function mulCAtB!(C, A, B)
    M, N = size(C)
    K = size(B, 1)
    @assert size(C, 1) == size(A, 2)
    @assert size(C, 2) == size(B, 2)
    @assert size(A, 1) == size(B, 1)
    # When the @avx macro is available, this code is faster:
    # z = zero(eltype(C))
    @avx for n in axes(C, 2), m in axes(C, 1)
        Cmn = zero(eltype(C)) # Cmn = z # For some reason this bombed in LoopVectorization 0.4.0
        for k in axes(A, 1)
            Cmn += A[k, m] * B[k, n]
        end
        C[m, n] = Cmn
    end
    # Otherwise, this code works just fine:
    # if mod(M, 2) == 0 && mod(N, 2) == 0
    #     @inbounds for m in 1:2:M
    #     	m1 = m + 1
    #     	@inbounds for n in 1:2:N 
    #     		n1 = n + 1
    # 	    	C11, C21, C12, C22 = 0.0, 0.0, 0.0, 0.0 
    # 	    	@inbounds for k in 1:K
    # 	    		C11 += A[k,m] * B[k,n] 
    # 	    		C21 += A[k,m1] * B[k,n] 
    # 	    		C12 += A[k,m] * B[k,n1] 
    # 	    		C22 += A[k,m1] * B[k,n1]
    # 	    	end
    # 	    	C[m,n] = C11
    # 	    	C[m1,n] = C21
    # 	    	C[m,n1] = C12
    # 	    	C[m1,n1] = C22
    # 	    end
    #     end
    # else
    # 	@inbounds for n in 1:N, m in 1:M 
    #     	Cmn = 0.0
    #     	@inbounds for k in 1:K
    #     		Cmn += A[k,m] * B[k,n]
    #     	end
    #     	C[m,n] = Cmn
    #     end
    # end
    return C
end

"""
    mulCAB!(C, A, B)

Compute the matrix `C = A * B`

The use of BLAS is purposefully avoided in order to eliminate contentions of
multi-threaded execution of the library code with the user-level threads.

Note: See the thread https://discourse.julialang.org/t/ann-loopvectorization/32843/36
"""
function mulCAB!(C, A, B)
    M, N = size(C)
    K = size(B, 1)
    @assert M == size(A, 1)
    @assert N == size(B, 2)
    @assert size(A, 2) == K
    # Surprisingly the code below (nkm) is slower than the one in the order mnk.
    # C .= 0
    # @avx  for n in 1:N, k in 1:K 
    # 	Bkn = B[k,n]
    # 	for m in 1:M
    # 		C[m,n] += A[m,k] * Bkn
    # 	end
    # end
    # z = zero(eltype(C))
    @avx for m = 1:M, n = 1:N
        Cmn = zero(eltype(C))
        for k = 1:K
            Cmn += A[m, k] * B[k, n]
        end
        C[m, n] = Cmn
    end
    return C
end

"""
    mulCABt!(C, A, B)

Compute the matrix `C = A * B'`

The use of BLAS is purposefully avoided in order to eliminate contentions of
multi-threaded execution of the library code with the user-level threads.
"""
function mulCABt!(C, A, B)
    M, N = size(C)
    K = size(B, 2)
    @assert size(C, 1) == size(A, 1)
    @assert size(C, 2) == size(B, 1)
    @assert size(A, 2) == size(B, 2)
    C .= 0
    @avx for n = 1:N, k = 1:K # for k ∈ 1:K, n ∈ 1:N #
        Bnk = B[n, k]
        for m = 1:M
            C[m, n] += A[m, k] * Bnk
        end
    end
    return C
end

"""
    detC(::Val{3}, C::Matrix{T})

Compute determinant of 3X3 `C`.
"""
function detC(::Val{3}, C::Matrix{T}) where {T}
    #define MAT3DUTIL_DET_3X3(T3X3) 
    return (
        (C[1, 1] * C[2, 2] * C[3, 3]) +
        (C[1, 2] * C[2, 3] * C[3, 1]) +
        (C[1, 3] * C[2, 1] * C[3, 2]) - (C[1, 3] * C[2, 2] * C[3, 1]) -
        (C[1, 2] * C[2, 1] * C[3, 3]) - (C[1, 1] * C[2, 3] * C[3, 2])
    )
end

"""
    mulCAB!(C::Vector{T}, A, B::Vector{T})  where {T}

Compute the product `C = A * B`, where `C` and `B` are "vectors".

The use of BLAS is purposefully avoided in order to eliminate contentions of
multi-threaded execution of the library code with the user-level threads.
"""
function mulCAB!(C::Vector{T}, A, B::Vector{T}) where {T}
    M, N = size(A)
    @assert length(C) == M
    @assert length(B) == N
    C .= 0.0
    @avx for n = 1:N
        for m = 1:M
            C[m] += A[m, n] * B[n]
        end
    end
    return C
end

"""
    export_sparse(filnam, M)

Export sparse matrix to a text file.
"""
function export_sparse(filnam, M)
    #readdlm(source, ',')
    open(filnam, "w") do io
        out = [(size(M)..., 0.0)]
        I, J, V = findnz(M)
        out = vcat(out, [(I[i], J[i], V[i]) for i in eachindex(I)])
        writedlm(io, out, ',')
    end
end

"""
    import_sparse(filnam)

Import sparse matrix from a text file.
"""
function import_sparse(filnam)
    B = open(filnam, "r") do io
        rows = readlines(io)
        s = split(replace(rows[1], "," => " "))
        m, n = parse(Int64, s[1]), parse(Int64, s[2])
        I = fill(0, length(rows) - 1)
        J = fill(0, length(rows) - 1)
        V = fill(0.0, length(rows) - 1)
        for i in eachindex(I)
            s = split(replace(rows[i+1], "," => " "))
            I[i] = parse(Int64, s[1])
            J[i] = parse(Int64, s[2])
            V[i] = parse(Float64, s[3])
        end
        sparse(I, J, V, m, n)
    end
    return B
end

"""
    matrix_blocked_ff(A, row_nfreedofs, col_nfreedofs = row_nfreedofs)

Extract the "free-free" partition of a matrix.

The matrix is assumed to be composed of four blocks
```
A = [A_ff A_fd
     A_df A_dd]
```

Here `f` stands for free, and `d` stands for data (i.e. fixed, prescribed, ...).
The size of the `ff` block is `row_nfreedofs, col_nfreedofs`. Neither one of
the blocks is square, unless `row_nfreedofs == col_nfreedofs`.

When `row_nfreedofs == col_nfreedofs`, only the number of rows needs to be given.
"""
function matrix_blocked_ff(A, row_nfreedofs, col_nfreedofs = row_nfreedofs)
    row_nalldofs, col_nalldofs = size(A)
    row_nfreedofs <= row_nalldofs || error("The ff block has too many rows")
    col_nfreedofs <= col_nalldofs || error("The ff block has too many columns")
    row_f_dim = row_nfreedofs
    col_f_dim = col_nfreedofs

    if (row_f_dim > 0 && col_f_dim > 0)
        A_ff = A[1:row_nfreedofs, 1:col_nfreedofs]
    else
        A_ff = spzeros(row_f_dim, col_f_dim)
    end
    return A_ff
end


"""
    matrix_blocked_fd(A, row_nfreedofs, col_nfreedofs = row_nfreedofs)

Extract the "free-data" partition of a matrix.

The matrix is assumed to be composed of four blocks
```
A = [A_ff A_fd
     A_df A_dd]
```

Here `f` stands for free, and `d` stands for data (i.e. fixed, prescribed, ...).
The size of the `ff` block is `row_nfreedofs, col_nfreedofs`. Neither one of
the blocks is square, unless `row_nfreedofs == col_nfreedofs`.

When `row_nfreedofs == col_nfreedofs`, only the number of rows needs to be given.
"""
function matrix_blocked_fd(A, row_nfreedofs, col_nfreedofs = row_nfreedofs)
    row_nalldofs, col_nalldofs = size(A)
    row_nfreedofs <= row_nalldofs || error("The ff block has too many rows")
    col_nfreedofs <= col_nalldofs || error("The ff block has too many columns")
    row_f_dim = row_nfreedofs
    row_d_dim = (row_nfreedofs < row_nalldofs ? row_nalldofs - row_nfreedofs : 0)
    col_f_dim = col_nfreedofs
    col_d_dim = (col_nfreedofs < col_nalldofs ? col_nalldofs - col_nfreedofs : 0)

    if (row_f_dim > 0 && col_d_dim > 0)
        A_fd = A[1:row_nfreedofs, (col_nfreedofs+1):end]
    else
        A_fd = spzeros(row_f_dim, col_d_dim)
    end
    return A_fd
end


"""
    matrix_blocked_df(A, row_nfreedofs, col_nfreedofs = row_nfreedofs)

Extract the "data-free" partition of a matrix.

The matrix is assumed to be composed of four blocks
```
A = [A_ff A_fd
     A_df A_dd]
```

Here `f` stands for free, and `d` stands for data (i.e. fixed, prescribed, ...).
The size of the `ff` block is `row_nfreedofs, col_nfreedofs`. Neither one of
the blocks is square, unless `row_nfreedofs == col_nfreedofs`.

When `row_nfreedofs == col_nfreedofs`, only the number of rows needs to be given.
"""
function matrix_blocked_df(A, row_nfreedofs, col_nfreedofs = row_nfreedofs)
    row_nalldofs, col_nalldofs = size(A)
    row_nfreedofs <= row_nalldofs || error("The ff block has too many rows")
    col_nfreedofs <= col_nalldofs || error("The ff block has too many columns")
    row_f_dim = row_nfreedofs
    row_d_dim = (row_nfreedofs < row_nalldofs ? row_nalldofs - row_nfreedofs : 0)
    col_f_dim = col_nfreedofs
    col_d_dim = (col_nfreedofs < col_nalldofs ? col_nalldofs - col_nfreedofs : 0)

    if (row_d_dim > 0 && col_f_dim > 0)
        A_df = A[(row_nfreedofs+1):end, 1:col_nfreedofs]
    else
        A_df = spzeros(row_d_dim, col_f_dim)
    end
    return A_df
end


"""
    matrix_blocked_dd(A, row_nfreedofs, col_nfreedofs = row_nfreedofs)

Extract the "data-data" partition of a matrix.

The matrix is assumed to be composed of four blocks
```
A = [A_ff A_fd
     A_df A_dd]
```

Here `f` stands for free, and `d` stands for data (i.e. fixed, prescribed, ...).
The size of the `ff` block is `row_nfreedofs, col_nfreedofs`. Neither one of
the blocks is square, unless `row_nfreedofs == col_nfreedofs`.

When `row_nfreedofs == col_nfreedofs`, only the number of rows needs to be given.
"""
function matrix_blocked_dd(A, row_nfreedofs, col_nfreedofs = row_nfreedofs)
    row_nalldofs, col_nalldofs = size(A)
    row_nfreedofs <= row_nalldofs || error("The ff block has too many rows")
    col_nfreedofs <= col_nalldofs || error("The ff block has too many columns")
    row_f_dim = row_nfreedofs
    row_d_dim = (row_nfreedofs < row_nalldofs ? row_nalldofs - row_nfreedofs : 0)
    col_f_dim = col_nfreedofs
    col_d_dim = (col_nfreedofs < col_nalldofs ? col_nalldofs - col_nfreedofs : 0)

    if (row_d_dim > 0 && col_d_dim > 0)
        A_dd = A[(row_nfreedofs+1):end, (col_nfreedofs+1):end]
    else
        A_dd = spzeros(row_d_dim, col_d_dim)
    end
    return A_dd
end

"""
    vector_blocked_f(V, nfreedofs)

Extract the "free" part of a vector.

The vector is composed of two blocks
```
V = [V_f
     V_d]
```
Here `f` stands for free, and `d` stands for data (i.e. fixed, prescribed, ...).
"""
function vector_blocked_f(V, nfreedofs)
    nalldofs = length(V)
    nfreedofs <= nalldofs || error("The f block has too many rows")
    row_f_dim = nfreedofs
    row_d_dim = (nfreedofs < nalldofs ? nalldofs - nfreedofs : 0)
    if (row_f_dim > 0)
        V_f = V[1:nfreedofs]
    else
        V_f = eltype(V)[]
    end
    return V_f
end

"""
    vector_blocked_d(V, nfreedofs)

Extract the "data" part of a vector.

The vector is composed of two blocks
```
V = [V_f
     V_d]
```
Here `f` stands for free, and `d` stands for data (i.e. fixed, prescribed, ...).
"""
function vector_blocked_d(V, nfreedofs)
    nalldofs = length(V)
    nfreedofs <= nalldofs || error("The f block has too many rows")
    row_f_dim = nfreedofs
    row_d_dim = (nfreedofs < nalldofs ? nalldofs - nfreedofs : 0)
    if (row_d_dim > 0)
        V_d = V[(nfreedofs+1):end]
    else
        V_d = eltype(V)[]
    end
    return V_d
end

"""
    setvectorentries!(a, v = zero(eltype(a)))

Set entries of a long vector to a given constant.

The operation runs in parallel on as many threads as are available.
"""
function setvectorentries!(a, v = zero(eltype(a)))
    @inbounds Threads.@threads for k in eachindex(a)
        a[k] = v
    end
    a
end


end
