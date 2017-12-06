"""
    MatrixUtilityModule

Module for general utility matrix product functions.
"""
module MatrixUtilityModule

using FinEtools.FTypesModule: FInt, FFlt, FCplxFlt, FFltVec, FIntVec, FFltMat, FIntMat, FMat, FVec, FDataDict

"""
    loc!(loc::FFltMat, X::FFltMat, conn::C, N::FFltMat) where {C}

Compute the location of the quadrature point.  

Arguments: 
`loc` = matrix of coordinates, overwritten  inside the function
`X` = matrix of the node coordinates from the entire mesh, index 
    with the connectivity `conn`
`conn` = connectivity, indexes the matrix `X`
`N` = matrix of basis function values
"""
function loc!(loc::FFltMat, X::FFltMat, conn::C, N::FFltMat) where {C}
    n = size(N, 1)
    @inbounds for j = 1:size(loc, 2)
        la = 0.0
        @inbounds for k = 1:n
            la += N[k] * X[conn[k], j]
        end
        loc[j] = la
    end
    return loc
end

"""
    jac!(J::FFltMat, X::FFltMat, conn::C, gradNparams::FFltMat) where {C}

Compute the Jacobian matrix at the quadrature point.  

Arguments: 
`J` = Jacobian matrix, overwritten  inside the function
`X` = matrix of the node coordinates from the entire mesh, index 
    with the connectivity `conn`
`conn` = connectivity, indexes the matrix `X`
`gradNparams` = matrix of basis function gradients
"""
function jac!(J::FFltMat, X::FFltMat, conn::C, gradNparams::FFltMat) where {C}
    n = size(gradNparams, 1)
    @inbounds for j = 1:size(J, 2)
        @inbounds for i = 1:size(J, 1)
            Ja = 0.0
            @inbounds for k = 1:n
                Ja += X[conn[k], i] * gradNparams[k, j]
            end
            J[i, j] = Ja
        end
    end
    return J
end

"""
    locjac!(loc::FFltMat, J::FFltMat, X::FFltMat, conn::C, N::FFltMat, gradNparams::FFltMat) where {C}

Compute location and Jacobian matrix at the quadrature point.  

Arguments: 
`loc` = matrix of coordinates, overwritten  inside the function
`J` = Jacobian matrix, overwritten  inside the function
`X` = matrix of the node coordinates from the entire mesh, index 
    with the connectivity `conn`
`conn` = connectivity, indexes the matrix `X`
`N` = matrix of basis function values
`gradNparams` = matrix of basis function gradients
"""
function locjac!(loc::FFltMat, J::FFltMat, X::FFltMat, conn::C, N::FFltMat, gradNparams::FFltMat) where {C}
    return loc!(loc, X, conn, N), jac!(J, X, conn, gradNparams)
end

"""
    add_mggt_ut_only!(Ke::FFltMat, gradN::FFltMat, mult::FFlt)

Add the product `gradN*kappa_bar*gradNT*(Jac*w[j])` to the elementwise matrix `Ke`.
*Only upper triangle* is computed; the lower triangle is not touched.
(Use `complete_lt!` to complete the lower triangle, if needed.)

The matrix `Ke` is assumed to be suitably initialized.

The matrix `Ke` is modified.  The matrices `gradN` and `kappa_bar` are not modified
inside this function. The scratch buffer `kappa_bargradNT` is overwritten
during each call of this function.
"""
function add_mggt_ut_only!(Ke::FFltMat, gradN::FFltMat, mult::FFlt)
    @assert size(Ke, 1) == size(Ke, 2)
    Kedim = size(Ke, 1)
    nne, mdim = size(gradN)
    mx::FInt = 0
    nx::FInt = 0
    px::FInt = 0
    # A_mul_Bt!(kappa_bargradNT, kappa_bar, gradN); # intermediate result
    @inbounds for nx = 1:Kedim # Do: Ce  =  Ce + gradN*((Jac*w[j]))*gradN' ;
        @inbounds for px = 1:mdim
        @inbounds for mx = 1:nx # only the upper triangle
            Ke[mx, nx] +=  gradN[mx, px]*(mult)*gradN[nx, px]
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
    mx::FInt = 0
    nx::FInt = 0
    px::FInt = 0
    # A_mul_Bt!(kappa_bargradNT, kappa_bar, gradN); # intermediate result
    @inbounds for nx = 1:nne
        @inbounds for mx = 1:mdim
        accum::FFlt  = 0.0
        @inbounds for px = 1:mdim
            accum += kappa_bar[mx, px]*gradN[nx, px]
        end
        kappa_bargradNT[mx, nx] = Jac_w*accum
        end
    end
    # Ke = Ke + gradN*(kappa_bar*(Jac*w[j]))*gradN'; only upper triangle
    @inbounds for nx = 1:Kedim
        @inbounds for mx = 1:nx # only the upper triangle
        accum::FFlt  = 0.0
        @inbounds for px = 1:mdim
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
    @assert size(Ke, 1) == size(Ke, 2)
    Kedim = size(Ke, 1)
    mx::FInt = 0
    nx::FInt = 0
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
    mx::FInt = 0
    nx::FInt = 0
    px::FInt = 0
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
    add_nnt_ut_only!(Ke::FFltMat, N::FFltMat, Jac_w_coeff::FFlt)

Add the product  `Nn*(Nn'*(coeff*(Jac*w(j)))`, to the elementwise matrix `Ke`.
*Only upper triangle* is computed; the lower triangle is not touched.

The matrix `Ke` is assumed to be suitably initialized.

The matrix `Ke` is modified.  The matrices `B` and `D` are not modified
inside this function. The scratch buffer `DB` is overwritten
during each call of this function.
"""
function add_nnt_ut_only!(Ke::FFltMat, N::FFltMat, Jac_w_coeff::FFlt)
    @assert size(Ke, 1) == size(Ke, 2)
    @assert size(N, 1) == size(Ke, 2)
    Kedim = size(N, 1)
    mx::FInt = 0
    nx::FInt = 0
    @inbounds for nx = 1:Kedim
        @inbounds for mx = 1:nx # only the upper triangle
        Ke[mx,nx] = Ke[mx,nx] + N[mx]*(Jac_w_coeff)*N[nx]
        end
    end
    return true
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
    mx::FInt = 0
    px::FInt = 0
    @inbounds for mx = 1:Kedim
        @inbounds for px = 1:nsig
        elvec[mx] = elvec[mx] + B[px, mx]*Jac_w_coeff*sig[px]
        end
    end
    return true
end

end
