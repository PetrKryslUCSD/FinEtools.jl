"""
    FEMMDeforSurfaceDampingModule

Module for operations on the damping associated with absorbing boundary 
conditions (ABC) representation of the effect of infinite extent 
of inviscid fluid next to the surface.
"""
module FEMMDeforSurfaceDampingModule

using FinEtools.FTypesModule: FInt, FFlt, FCplxFlt, FFltVec, FIntVec, FFltMat, FIntMat, FMat, FVec, FDataDict
import FinEtools.FENodeSetModule: FENodeSet
import FinEtools.FESetModule: FESet, nodesperelem, manifdim
import FinEtools.IntegDomainModule: IntegDomain, integrationdata, Jacobiansurface
import FinEtools.FieldModule: ndofs, gatherdofnums!
import FinEtools.NodalFieldModule: NodalField 
import FinEtools.FEMMBaseModule: FEMMAbstractBase
import FinEtools.AssemblyModule: SysvecAssemblerBase, SysmatAssemblerBase, SysmatAssemblerSparseSymm, startassembly!, assemble!, makematrix!
import FinEtools.MatrixUtilityModule: add_nnt_ut_only!, complete_lt!, locjac!
import FinEtools.SurfaceNormalModule: SurfaceNormal, updatenormal!
import LinearAlgebra: norm, cross

"""
    FEMMDeforSurfaceDamping{S<:FESet, F<:Function}

Type for surface damping model.
"""
mutable struct FEMMDeforSurfaceDamping{S<:FESet, F<:Function} <: FEMMAbstractBase
    integdomain::IntegDomain{S, F} # geometry data
end

"""
    dampingABC(self::FEMMDeforSurfaceDamping, assembler::A,
                  geom::NodalField{FFlt}, u::NodalField{T1},
                  impedance::T2, surfacenormal::SurfaceNormal) where {A<:SysmatAssemblerBase, T1<:Number, T2<:Number}

Compute the damping matrix associated with absorbing boundary conditions (ABC) representation of the effect of infinite extent of inviscid fluid next to the surface.
"""
function dampingABC(self::FEMMDeforSurfaceDamping, assembler::A,
    geom::NodalField{FFlt}, u::NodalField{T1},
    impedance::T2, surfacenormal::SurfaceNormal) where {A<:SysmatAssemblerBase, T1<:Number, T2<:Number}
    fes = self.integdomain.fes
    # Constants
    nfes = count(fes); # number of finite elements
    ndn = ndofs(u); # number of degrees of freedom per node
    nne = nodesperelem(fes); # number of nodes per finite element
    sdim = ndofs(geom); # spatial dimension
    mdim = manifdim(fes); # manifold dimension of the finite elements
    Cedim = ndn * nne; # size of damping element matrix
    # Precompute basis f. values + basis f. gradients wrt parametric coor
    npts, Ns, gradNparams, w, pc  =  integrationdata(self.integdomain);
    loc = zeros(FFlt, 1, sdim); # quadrature point coordinate -- used as a buffer
    J = zeros(FFlt, sdim, mdim); # Jacobian matrix -- used as a buffer
    Ce = zeros(T2, Cedim, Cedim); # element damping matrix -- used as a buffer
    Nn = zeros(FFlt, Cedim); # column vector
    dofnums = zeros(FInt, Cedim); # degree of freedom array -- used as a buffer
    # Prepare assembler and temporaries
    startassembly!(assembler, Cedim, Cedim, nfes, u.nfreedofs, u.nfreedofs);
    for i = 1:nfes # loop over finite elements
        fill!(Ce, 0.0); # Initialize element damping matrix
        for j = 1:npts # loop over quadrature points
            locjac!(loc, J, geom.values, fes.conn[i], Ns[j], gradNparams[j])
            Jac = Jacobiansurface(self.integdomain, J, loc, fes.conn[i], Ns[j]);
            n = updatenormal!(surfacenormal, loc, J, self.integdomain.fes.label[i]);
            for k = 1:nne
                Nn[(k-1)*ndn+1:k*ndn] = n * Ns[j][k]
            end
            add_nnt_ut_only!(Ce, Nn, (-1.0)*impedance*Jac*w[j]);
        end # end loop over quadrature points
        complete_lt!(Ce);
        gatherdofnums!(u, dofnums, self.integdomain.fes.conn[i]); # retrieve degrees of freedom
        assemble!(assembler, Ce, dofnums, dofnums); # assemble the element damping matrix
    end # end loop over finite elements
    return makematrix!(assembler);
end

function dampingABC(self::FEMMDeforSurfaceDamping, geom::NodalField{FFlt},
    u::NodalField{T1}, impedance::T2, surfacenormal::SurfaceNormal) where {T1<:Number, T2<:Number}
    assembler = SysmatAssemblerSparseSymm();
    return dampingABC(self, assembler, geom, u, impedance, surfacenormal);
end

end
