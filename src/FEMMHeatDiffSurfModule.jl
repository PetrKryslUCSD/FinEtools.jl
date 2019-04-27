"""
    FEMMHeatDiffSurfModule

Module for operations on boundaries of domains to construct system matrices and
system vectors for linear heat diffusion/conduction.
"""
module FEMMHeatDiffSurfModule

using FinEtools.FTypesModule: FInt, FFlt, FCplxFlt, FFltVec, FIntVec, FFltMat, FIntMat, FMat, FVec, FDataDict
import FinEtools.FENodeSetModule: FENodeSet
import FinEtools.FESetModule: AbstractFESet, nodesperelem, manifdim
import FinEtools.IntegDomainModule: IntegDomain, integrationdata, Jacobiansurface
import FinEtools.FieldModule: ndofs, gatherdofnums!, gatherfixedvalues_asvec!, gathervalues_asvec!
import FinEtools.NodalFieldModule: NodalField 
import FinEtools.AssemblyModule: AbstractSysvecAssembler, AbstractSysmatAssembler, SysmatAssemblerSparseSymm, startassembly!, assemble!, makematrix!, makevector!, SysvecAssembler
import FinEtools.FEMMBaseModule: AbstractFEMM
import FinEtools.MatrixUtilityModule: add_gkgt_ut_only!, add_nnt_ut_only!, complete_lt!, locjac!
import LinearAlgebra: norm, dot, cross

"""
    FEMMHeatDiffSurf{S<:AbstractFESet, F<:Function} <: AbstractFEMM

    Type for heat diffusion finite element modeling machine for boundary integrals.
"""
mutable struct FEMMHeatDiffSurf{S<:AbstractFESet, F<:Function} <: AbstractFEMM
    integdomain::IntegDomain{S, F} # geometry data
    surfacetransfercoeff::FFlt # material object
end

"""
    surfacetransfer(self::FEMMHeatDiffSurf,  assembler::A,
      geom::NodalField{FFlt}, temp::NodalField{FFlt}) where {A<:AbstractSysmatAssembler}

Compute the surface heat transfer matrix.
"""
function surfacetransfer(self::FEMMHeatDiffSurf,  assembler::A, geom::NodalField{FFlt}, temp::NodalField{FFlt}) where {A<:AbstractSysmatAssembler}
    fes = self.integdomain.fes
    # Constants
    nfes = count(fes); # number of finite elements in the set
    ndn = ndofs(temp); # number of degrees of freedom per node
    nne = nodesperelem(fes); # number of nodes for element
    sdim = ndofs(geom);            # number of space dimensions
    mdim = manifdim(fes); # manifold dimension of the element
    Hedim = ndn*nne;             # dimension of the element matrix
    # Precompute basis f. values + basis f. gradients wrt parametric coor
    npts,  Ns,  gradNparams,  w,  pc = integrationdata(self.integdomain);
    # Prepare assembler and temporaries
    He = fill(zero(FFlt), Hedim, Hedim);                # element matrix -- used as a buffer
    dofnums = zeros(FInt, Hedim); # degree of freedom array -- used as a buffer
    loc = fill(zero(FFlt), 1, sdim); # quadrature point location -- used as a buffer
    J = fill(zero(FFlt), sdim, mdim); # Jacobian matrix -- used as a buffer
    startassembly!(assembler, Hedim, Hedim, nfes, temp.nfreedofs, temp.nfreedofs);
    for i = 1:nfes # Loop over elements
        fill!(He,  0.0); # Initialize element matrix
        for j=1:npts # Loop over quadrature points
            locjac!(loc, J, geom.values, fes.conn[i], Ns[j], gradNparams[j]) 
            Jac = Jacobiansurface(self.integdomain, J, loc, fes.conn[i],  Ns[j]);
            add_nnt_ut_only!(He, Ns[j], self.surfacetransfercoeff*Jac*w[j])
        end # Loop over quadrature points
        complete_lt!(He)
        gatherdofnums!(temp, dofnums, fes.conn[i]);# retrieve degrees of freedom
        assemble!(assembler, He, dofnums, dofnums);# assemble symmetric matrix
    end # Loop over elements
    return makematrix!(assembler);
end

function surfacetransfer(self::FEMMHeatDiffSurf{S},
                                   geom::NodalField{FFlt},
                                   temp::NodalField{FFlt}) where {S<:AbstractFESet}
    assembler = SysmatAssemblerSparseSymm()
    return  surfacetransfer(self, assembler, geom, temp);
end

"""
    surfacetransferloads(self::FEMMHeatDiffSurf,  assembler::A,
      geom::NodalField{FFlt}, temp::NodalField{FFlt},
      ambtemp::NodalField{FFlt}) where {A<:AbstractSysvecAssembler}

Compute the load vector corresponding to surface heat transfer.
"""
function surfacetransferloads(self::FEMMHeatDiffSurf,  assembler::A,  geom::NodalField{FFlt}, temp::NodalField{FFlt},  ambtemp::NodalField{FFlt}) where {A<:AbstractSysvecAssembler}
    fes = self.integdomain.fes
    # Constants
    nfes = count(fes); # number of finite elements in the set
    ndn = ndofs(temp); # number of degrees of freedom per node
    nne = nodesperelem(fes); # number of nodes for element
    sdim = ndofs(geom);            # number of space dimensions
    mdim = manifdim(fes); # manifold dimension of the element
    Hedim = ndn*nne;             # dimension of the element matrix
    # Precompute basis f. values + basis f. gradients wrt parametric coor
    npts,  Ns,  gradNparams,  w,  pc = integrationdata(self.integdomain);
    # Prepare assembler and temporaries
    Fe = fill(zero(FFlt), Hedim, 1); # element matrix -- used as a buffer
    dofnums = zeros(FInt, Hedim); # degree of freedom array -- used as a buffer
    loc = fill(zero(FFlt), 1, sdim); # quadrature point location -- used as a buffer
    J = fill(zero(FFlt), sdim, mdim); # Jacobian matrix -- used as a buffer
    pT = fill(zero(FFlt), Hedim);
    startassembly!(assembler,  temp.nfreedofs);
    for i = 1:nfes # Loop over elements
        gathervalues_asvec!(ambtemp, pT, fes.conn[i]);# retrieve ambient temp
        if norm(pT, Inf) != 0.0    # Is the load nonzero?
            fill!(Fe,  0.0); # Initialize element matrix
            for j=1:npts # Loop over quadrature points
                locjac!(loc, J, geom.values, fes.conn[i], Ns[j], gradNparams[j]) 
                Jac = Jacobiansurface(self.integdomain, J, loc, fes.conn[i],  Ns[j]);
                Ta = dot(vec(pT), vec(Ns[j]))
                factor = Ta*self.surfacetransfercoeff*Jac*w[j]
                Fe .+= factor*Ns[j]
            end # Loop over quadrature points
            gatherdofnums!(temp, dofnums, fes.conn[i]); # retrieve degrees of freedom
            assemble!(assembler,  Fe,  dofnums); # assemble element load vector
        end
    end
    F = makevector!(assembler);
    return F
end


function surfacetransferloads(self::FEMMHeatDiffSurf,
                                        geom::NodalField{FFlt},
                                        temp::NodalField{FFlt},
                                        ambtemp::NodalField{FFlt})
    assembler = SysvecAssembler()
    return  surfacetransferloads(self, assembler, geom, temp, ambtemp);
end

"""
    nzebcsurfacetransferloads(self::FEMMHeatDiffSurf, assembler::A,
      geom::NodalField{FFlt}, temp::NodalField{FFlt}) where {A<:AbstractSysvecAssembler}

Compute load vector for nonzero EBC for fixed temperature.
"""
function nzebcsurfacetransferloads(self::FEMMHeatDiffSurf, assembler::A,  geom::NodalField{FFlt}, temp::NodalField{FFlt}) where {A<:AbstractSysvecAssembler}
    fes = self.integdomain.fes
    # Constants
    nfes = count(fes); # number of finite elements in the set
    ndn = ndofs(temp); # number of degrees of freedom per node
    nne = nodesperelem(fes); # number of nodes for element
    sdim = ndofs(geom);            # number of space dimensions
    mdim = manifdim(fes); # manifold dimension of the element
    Hedim = ndn*nne;             # dimension of the element matrix
    # Precompute basis f. values + basis f. gradients wrt parametric coor
    npts,  Ns,  gradNparams,  w,  pc = integrationdata(self.integdomain);
    # Prepare assembler and temporaries
    He = fill(zero(FFlt), Hedim, Hedim);                # element matrix -- used as a buffer
    dofnums = zeros(FInt, Hedim); # degree of freedom array -- used as a buffer
    loc = fill(zero(FFlt), 1, sdim); # quadrature point location -- used as a buffer
    J = fill(zero(FFlt), sdim, mdim); # Jacobian matrix -- used as a buffer
    pT = fill(zero(FFlt), Hedim);
    startassembly!(assembler,  temp.nfreedofs);
    # Now loop over all finite elements in the set
    for i=1:nfes # Loop over elements
        gatherfixedvalues_asvec!(temp, pT, fes.conn[i]);# retrieve element temp
        if norm(pT, Inf) != 0.0    # Is the load nonzero?
            fill!(He,  0.0);
            for j=1:npts # Loop over quadrature points
                locjac!(loc, J, geom.values, fes.conn[i], Ns[j], gradNparams[j]) 
                Jac = Jacobiansurface(self.integdomain, J, loc, fes.conn[i],  Ns[j]);
                add_nnt_ut_only!(He, Ns[j], self.surfacetransfercoeff*Jac*w[j])
            end # Loop over quadrature points
            complete_lt!(He)
            gatherdofnums!(temp, dofnums, fes.conn[i]); # retrieve degrees of freedom
            assemble!(assembler,  vec(He*pT),  dofnums); # assemble element load vector
        end
    end
    F= makevector!(assembler);
    return F
end

function nzebcsurfacetransferloads(self::FEMMHeatDiffSurf,
                                          geom::NodalField{FFlt},
                                          temp::NodalField{FFlt})
    assembler = SysvecAssembler()
    return  nzebcsurfacetransferloads(self, assembler, geom, temp);
end


end
