"""
    MatrixUtilityModule

Module for general utility matrix product functions.
"""
module MatrixUtilityModule

__precompile__(true)

using DelimitedFiles
using SparseArrays
using ..FTypesModule:
    FInt, FFlt, FCplxFlt, FFltVec, FIntVec, FFltMat, FIntMat, FMat, FVec, FDataDict
using LoopVectorization

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
    return mulCAtB!(J, ecoords, gradNparams)
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
function locjac!(
    loc::FFltMat,
    J::FFltMat,
    ecoords::FFltMat,
    N::FFltMat,
    gradNparams::FFltMat,
)
    return loc!(loc, ecoords, N), jac!(J, ecoords, gradNparams)
end

"""
    add_mggt_ut_only!(Ke::FFltMat, gradN::FFltMat, mult::FFlt)

Add the product `gradN*mult*gradNT` to the matrix `Ke`.

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
    @inbounds for nx in 1:Kedim # Do: Ce  =  Ce + gradN*((Jac*w[j]))*gradN' ;
        @inbounds for px in 1:mdim
            a = (mult) * gradN[nx, px]
            @inbounds for mx in 1:nx # only the upper triangle
                Ke[mx, nx] += gradN[mx, px] * a
            end
        end
    end
    return true
end


"""
    add_gkgt_ut_only!(Ke::FFltMat, gradN::FFltMat, Jac_w::FFlt,
      kappa_bar::FFltMat, kappa_bargradNT::FFltMat)

Add the product `gradN*kappa_bar*gradNT*(Jac*w[j])` to the matrix `Ke`.

*Only upper triangle* is computed; the lower triangle is not touched.
(Use `complete_lt!` to complete the lower triangle, if needed.)

The matrix `Ke` is assumed to be suitably initialized.

Upon return,  the matrix `Ke` is updated.  The scratch buffer `kappa_bargradNT`
is overwritten during each call of this function. The matrices `gradN` and
`kappa_bar` are not modified inside this function.
"""
function add_gkgt_ut_only!(
    Ke::FFltMat,
    gradN::FFltMat,
    Jac_w::FFlt,
    kappa_bar::FFltMat,
    kappa_bargradNT::FFltMat,
)
    @assert size(Ke, 1) == size(Ke, 2)
    Kedim = size(Ke, 1)
    nne, mdim = size(gradN)
    @assert size(kappa_bar) == (mdim, mdim)
    @assert size(kappa_bargradNT) == (mdim, nne)
    # A_mul_Bt!(kappa_bargradNT, kappa_bar, gradN); # intermediate result
    @inbounds for nx in 1:nne
        @inbounds for mx in 1:mdim
            accum::FFlt = 0.0
            @inbounds @simd for px in 1:mdim
                accum += kappa_bar[mx, px] * gradN[nx, px]
            end
            kappa_bargradNT[mx, nx] = Jac_w * accum
        end
    end
    # Ke = Ke + gradN*(kappa_bar*(Jac*w[j]))*gradN'; only upper triangle
    @inbounds for nx in 1:Kedim
        @inbounds for mx in 1:nx # only the upper triangle
            accum::FFlt = 0.0
            @inbounds @simd for px in 1:mdim
                accum += gradN[mx, px] * kappa_bargradNT[px, nx]
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
    @inbounds for nx in 1:Kedim # complete the lower triangle
        @inbounds for mx in nx+1:Kedim
            Ke[mx, nx] = Ke[nx, mx]
        end
    end
    return true
end

"""
    add_btdb_ut_only!(Ke::FFltMat, B::FFltMat, Jac_w::FFlt,
                  D::FFltMat, DB::FFltMat)

Add the product  `(B'*(D*(Jac*w[j]))*B)`, to the matrix Ke.

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
    @inbounds for nx in 1:Kedim
        @inbounds for mx in 1:nstr
            accum::FFlt = 0.0
            @inbounds for px in 1:nstr
                accum += D[mx, px] * B[px, nx]
            end
            DB[mx, nx] = Jac_w * accum
        end
    end
    #  Ke = Ke + (B'*(D*(Jac*w[j]))*B); only the upper triangle
    @inbounds for nx in 1:Kedim
        @inbounds for mx in 1:nx # only the upper triangle
            accum::FFlt = 0.0
            @inbounds for px in 1:nstr
                accum += B[px, mx] * DB[px, nx]
            end
            Ke[mx, nx] += accum
        end
    end
    return true
end

"""
    add_b1tdb2!(Ke::FFltMat, B1::FFltMat, B2::FFltMat, Jac_w::FFlt, D::FFltMat, DB2::FFltMat)

Add the product  `(B1'*(D*(Jac_w))*B2)`, to the matrix `Ke`.

The matrix `Ke` is assumed to be suitably initialized: the results of this
computation are **added**. The matrix `Ke` may be rectangular.

The matrix `D` may be rectangular.

The matrix Ke is modified.  The matrices `B1`, `B2`, and `D` are not modified
inside this function. The scratch buffer `DB` is overwritten
during each call of this function.
"""
function add_b1tdb2!(
    Ke::FFltMat,
    B1::FFltMat,
    B2::FFltMat,
    Jac_w::FFlt,
    D::FFltMat,
    DB2::FFltMat,
)
    Kedim1, Kedim2 = size(Ke)
    Ddim1, Ddim2 = size(D)
    @assert size(B1) == (Ddim1, Kedim1)
    @assert size(B2) == (Ddim2, Kedim2)
    @assert size(DB2) == (Ddim1, Kedim2)
    # A_mul_B!(DB, D, B)    # intermediate product
    @inbounds for nx in 1:Kedim2
        for mx in 1:Ddim1
            accum::FFlt = 0.0
            for px in 1:Ddim2
                accum += D[mx, px] * B2[px, nx]
            end
            DB2[mx, nx] = Jac_w * accum
        end
    end
    #  Ke = Ke + (B1'*(D*(Jac_w))*B2)
    @inbounds for nx in 1:Kedim2
        for mx in 1:Kedim1 # only the upper triangle
            accum::FFlt = 0.0
            for px in 1:Ddim1
                accum += B1[px, mx] * DB2[px, nx]
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
    @inbounds for nx in 1:Kedim
        accum::FFlt = 0.0
        @inbounds for px in 1:nstr
            accum += sigma[px] * B[px, nx]
        end
        Fe[nx] += coefficient * accum
    end
    return true
end

"""
    add_nnt_ut_only!(Ke::FFltMat, Nn::FFltMat, Jac_w_coeff::FFlt)

Add the product  `Nn*(Nn'*(coeff*(Jac*w(j)))`, to the matrix `Ke`.

*Only upper triangle* is computed; the lower triangle is not touched.

The matrix `Ke` is assumed to be suitably initialized. The 
matrix `Nn` has a single column.

The matrix `Ke` is modified.  The matrix `Nn` is not modified
inside this function.
"""
function add_nnt_ut_only!(Ke::FMat{T}, N::FFltMat, Jac_w_coeff::T) where {T<:Number}
    Kedim = size(N, 1)
    @assert Kedim == size(Ke, 2)
    @inbounds for nx in 1:Kedim
        a = (Jac_w_coeff) * N[nx]
        @inbounds for mx in 1:nx # only the upper triangle
            Ke[mx, nx] += N[mx] * a
        end
    end
    return true
end

function add_nnt_ut_only!(Ke::FFltMat, N::FFltVec, Jac_w_coeff::FFlt)
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
        @inbounds for r in c:size(a, 1)
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
    @avx for m in 1:M, n in 1:N
        Cmn = zero(eltype(C))
        for k in 1:K
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
    @avx for n in 1:N, k in 1:K # for k ∈ 1:K, n ∈ 1:N #
        Bnk = B[n, k]
        for m in 1:M
            C[m, n] += A[m, k] * Bnk
        end
    end
    return C
end

"""
    detC(::Val{3}, C::FFltMat)

Compute determinant of 3X3 `C`.
"""
function detC(::Val{3}, C::FFltMat)
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
    @avx for n in 1:N
        for m in 1:M
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
        B = sparse(I, J, V, m, n)
    end
    return B
end

end
