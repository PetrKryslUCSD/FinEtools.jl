"""
    FEMMAcoustModule

Module for operations on interiors of domains to construct system matrices and
system vectors for linear acoustics.
"""
module FEMMAcoustModule

import Base.Complex

using FinEtools.FTypesModule: FInt, FFlt, FCplxFlt, FFltVec, FIntVec, FFltMat, FIntMat, FMat, FVec, FDataDict
import FinEtools.FENodeSetModule: FENodeSet
import FinEtools.FESetModule: FESet, gradN!, nodesperelem, manifdim
import FinEtools.MatAcoustFluidModule: MatAcoustFluid
import FinEtools.IntegDomainModule: IntegDomain, integrationdata, Jacobianvolume
import FinEtools.FieldModule: ndofs, gatherdofnums!, gatherfixedvalues_asvec!
import FinEtools.NodalFieldModule: NodalField 
import FinEtools.AssemblyModule: SysvecAssemblerBase, SysmatAssemblerBase, SysmatAssemblerSparseSymm, startassembly!, assemble!, makematrix!, SysvecAssembler, makevector!
import FinEtools.FEMMBaseModule: FEMMAbstractBase
import FinEtools.MatrixUtilityModule: add_mggt_ut_only!, add_nnt_ut_only!, complete_lt!, locjac!
import LinearAlgebra: norm

"""
    FEMMAcoust{S<:FESet}

Type for linear acoustics finite element modeling machine.
"""
mutable struct FEMMAcoust{S<:FESet, F<:Function, M} <: FEMMAbstractBase
    integdomain::IntegDomain{S, F} # geometry data 
    material::M # material object
end

function  buffers(self::FEMMAcoust, geom::NodalField{FFlt}, P::NodalField{F}) where {F}
    fes = self.integdomain.fes
    # Constants
    nfes = count(fes); # number of finite elements in the set
    ndn = ndofs(P); # number of degrees of freedom per node
    nne =  nodesperelem(fes); # number of nodes per element
    sdim =  ndofs(geom);            # number of space dimensions
    mdim = manifdim(fes);     # manifold dimension of the element
    Cedim = ndn*nne;          # dimension of the element matrix
    # Prepare assembler and temporaries
    dofnums = fill(zero(FInt), Cedim); # degree of freedom array -- used as a buffer
    loc = fill(zero(FFlt), 1, sdim); # quadrature point location -- used as a buffer
    J = fill(zero(FFlt), sdim, mdim); # Jacobian matrix -- used as a buffer
    gradN = fill(zero(FFlt), nne, mdim); # intermediate result -- used as a buffer
    elmat = fill(zero(FFlt), Cedim, Cedim);       # element matrix -- used as a buffer
    elvec = fill(zero(F), Cedim); # buffer
    elvecfix = fill(zero(F), Cedim); # buffer
    return dofnums, loc, J, gradN, elmat, elvec, elvecfix
end

"""
    acousticmass(self::FEMMAcoust,
      assembler::A, geom::NodalField,
      P::NodalField{T}) where {T<:Number, A<:SysmatAssemblerBase}

Compute the acoustic mass matrix.

Return K as a matrix.
Arguments
self   =  acoustics model
assembler  =  matrix assembler
geom = geometry field
P = acoustic (perturbation) pressure field
"""
function acousticmass(self::FEMMAcoust, assembler::A, geom::NodalField, P::NodalField{T}) where {T<:Number, A<:SysmatAssemblerBase}
    fes = self.integdomain.fes
    dofnums, loc, J, gradN, elmat, elvec, elvecfix =   buffers(self, geom, P)
    # Precompute basis f. values + basis f. gradients wrt parametric coor
    npts, Ns, gradNparams, w, pc  =  integrationdata(self.integdomain);
    Jac = 0.0;
    afactor = 0.0;
    startassembly!(assembler, size(elmat,1), size(elmat,2), count(fes),
        P.nfreedofs, P.nfreedofs);
    for i = 1:count(fes) # Loop over elements
        fill!(elmat, 0.0); # Initialize element matrix
        for j = 1:npts # Loop over quadrature points
            locjac!(loc, J, geom.values, fes.conn[i], Ns[j], gradNparams[j]) 
            Jac = Jacobianvolume(self.integdomain, J, loc, fes.conn[i], Ns[j]);
            # gradient WRT global Cartesian coordinates
            gradN!(fes, gradN, gradNparams[j], J);
            afactor = (Jac*w[j]);
            add_mggt_ut_only!(elmat, gradN, afactor)
        end # Loop over quadrature points
        complete_lt!(elmat)
        gatherdofnums!(P, dofnums, fes.conn[i]);# retrieve degrees of freedom
        assemble!(assembler, elmat, dofnums, dofnums);# assemble symmetric matrix
    end # Loop over elements
    return makematrix!(assembler);
end

function acousticmass(self::FEMMAcoust, geom::NodalField, P::NodalField{T}) where {T<:Number}
    # Make the default assembler object.
    assembler  =  SysmatAssemblerSparseSymm();
    return acousticmass(self, assembler, geom, P);
end

"""
    nzebcloadsacousticmass(self::FEMMAcoust, assembler::A,
      geom::NodalField, P::NodalField{T}) where {T<:Number,
      A<:SysvecAssemblerBase}

Compute load vector for nonzero EBC for fixed pressure..
"""
function nzebcloadsacousticmass(self::FEMMAcoust, assembler::A, geom::NodalField, P::NodalField{T}) where {T<:Number, A<:SysvecAssemblerBase}
    fes = self.integdomain.fes
    dofnums, loc, J, gradN, elmat, elvec, elvecfix = buffers(self, geom, P)
    # Precompute basis f. values + basis f. gradients wrt parametric coor
    npts, Ns, gradNparams, w, pc  =  integrationdata(self.integdomain);
    startassembly!(assembler, P.nfreedofs);
    # Now loop over all finite elements in the set
    for i = 1:count(fes) # Loop over elements
        gatherfixedvalues_asvec!(P, elvecfix, fes.conn[i]);# retrieve element coordinates
        if norm(elvecfix, Inf) !=  0.0     # Is the load nonzero?
            fill!(elmat, 0.0);
            for j = 1:npts # Loop over quadrature points
                locjac!(loc, J, geom.values, fes.conn[i], Ns[j], gradNparams[j]) 
                Jac = Jacobianvolume(self.integdomain, J, loc, fes.conn[i], Ns[j]);
                # gradient WRT global Cartesian coordinates
                gradN!(fes, gradN, gradNparams[j], J);
                afactor = (Jac*w[j]);
                add_mggt_ut_only!(elmat, gradN, afactor)
            end # Loop over quadrature points
            complete_lt!(elmat)
            gatherdofnums!(P, dofnums, fes.conn[i]);# retrieve degrees of freedom
            assemble!(assembler, -elmat*elvecfix, dofnums); # assemble element load vector
        end
    end
    return  makevector!(assembler);
end

function nzebcloadsacousticmass(self::FEMMAcoust, geom::NodalField, P::NodalField{T}) where {T<:Number}
    assembler  =  SysvecAssembler(P.values[1])
    return nzebcloadsacousticmass(self, assembler, geom, P);
end

"""
    acousticstiffness(self::FEMMAcoust, assembler::A,
      geom::NodalField,
      Pddot::NodalField{T}) where {T<:Number,
      A<:SysmatAssemblerBase}

Compute the acoustic stiffness matrix.
"""
function acousticstiffness(self::FEMMAcoust, assembler::A, geom::NodalField, Pddot::NodalField{T}) where {T<:Number, A<:SysmatAssemblerBase}
    fes = self.integdomain.fes
    dofnums, loc, J, gradN, elmat, elvec, elvecfix = buffers(self, geom, Pddot)
    # Precompute basis f. values + basis f. gradients wrt parametric coor
    npts, Ns, gradNparams, w, pc  =  integrationdata(self.integdomain);
    # Material
    bulk_modulus  =   self.material.bulk_modulus;
    mass_density  =   self.material.mass_density;
    c  =  sqrt(bulk_modulus/mass_density); # sound speed
    oc2 = 1.0/c^2;
    startassembly!(assembler, size(elmat,1), size(elmat,2), count(fes),
        Pddot.nfreedofs, Pddot.nfreedofs);
    for i = 1:count(fes) # Loop over elements
        fill!(elmat, 0.0); # Initialize element matrix
        for j = 1:npts # Loop over quadrature points
            locjac!(loc, J, geom.values, fes.conn[i], Ns[j], gradNparams[j]) 
            Jac = Jacobianvolume(self.integdomain, J, loc, fes.conn[i], Ns[j]);
            ffactor = Jac*oc2*w[j]
            add_nnt_ut_only!(elmat, Ns[j], ffactor)
        end # Loop over quadrature points
        complete_lt!(elmat)
        gatherdofnums!(Pddot, dofnums, fes.conn[i]);# retrieve degrees of freedom
        assemble!(assembler, elmat, dofnums, dofnums);# assemble symmetric matrix
    end # Loop over elements
    return makematrix!(assembler);
end

function acousticstiffness(self::FEMMAcoust, geom::NodalField, Pddot::NodalField{T}) where {T<:Number}
    # Make the default assembler object.
    assembler  =  SysmatAssemblerSparseSymm();
    return acousticstiffness(self, assembler, geom, Pddot);
end

"""
    nzebcloadsacousticstiffness(self::FEMMAcoust, assembler::A,
      geom::NodalField,
      Pddot::NodalField{T}) where {T<:Number,
      A<:SysvecAssemblerBase}

Compute load vector for nonzero EBC for fixed second-order pressure rate.
"""
function nzebcloadsacousticstiffness(self::FEMMAcoust, assembler::A, geom::NodalField, Pddot::NodalField{T}) where {T<:Number, A<:SysvecAssemblerBase}
    fes = self.integdomain.fes
    dofnums, loc, J, gradN, elmat, elvec, elvecfix = buffers(self, geom, Pddot)
    # Precompute basis f. values + basis f. gradients wrt parametric coor
    npts, Ns, gradNparams, w, pc  =  integrationdata(self.integdomain);
    # Material
    bulk_modulus  =   self.material.bulk_modulus;
    mass_density  =   self.material.mass_density;
    c  =  sqrt(bulk_modulus/mass_density); # sound speed
    oc2 = 1.0/c^2;
    startassembly!(assembler, Pddot.nfreedofs);
    # Now loop over all finite elements in the set
    for i = 1:count(fes) # Loop over elements
        gatherfixedvalues_asvec!(Pddot, elvecfix, fes.conn[i]);# retrieve element coordinates
        if norm(elvecfix, Inf) !=  0.0  # Is the load nonzero?
            fill!(elmat, 0.0);
            for j = 1:npts # Loop over quadrature points
                locjac!(loc, J, geom.values, fes.conn[i], Ns[j], gradNparams[j]) 
                Jac = Jacobianvolume(self.integdomain, J, loc, fes.conn[i], Ns[j]);
                ffactor = Jac*oc2*w[j]
                add_nnt_ut_only!(elmat, Ns[j], ffactor)
            end # Loop over quadrature points
            complete_lt!(elmat)
            gatherdofnums!(Pddot, dofnums, fes.conn[i]); # retrieve degrees of freedom
            assemble!(assembler, -elmat*elvecfix, dofnums); # assemble element load vector
        end
    end
    return makevector!(assembler);
end

function nzebcloadsacousticstiffness(self::FEMMAcoust, geom::NodalField, Pddot::NodalField{T}) where {T<:Number}
    assembler  =  SysvecAssembler(Pddot.values[1])#T(0.0)
    return nzebcloadsacousticstiffness(self, assembler, geom, Pddot)
end

end
