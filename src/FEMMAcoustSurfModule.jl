"""
    FEMMAcoustSurfModule

Module for operations on boundaries of domains to construct system matrices and
system vectors for linear acoustics.
"""
module FEMMAcoustSurfModule

import Base.Complex

using FinEtools.FTypesModule: FInt, FFlt, FCplxFlt, FFltVec, FIntVec, FFltMat, FIntMat, FMat, FVec, FDataDict
import FinEtools.FENodeSetModule: FENodeSet
import FinEtools.FESetModule: AbstractFESet, gradN!, nodesperelem, manifdim
import FinEtools.MatAcoustFluidModule: MatAcoustFluid, bulkmodulus
import FinEtools.MatModule: massdensity
import FinEtools.IntegDomainModule: IntegDomain, integrationdata, Jacobiansurface
import FinEtools.FieldModule: ndofs, gatherdofnums!
import FinEtools.NodalFieldModule: NodalField 
import FinEtools.GeneralFieldModule: GeneralField 
import FinEtools.AssemblyModule: AbstractSysvecAssembler, AbstractSysmatAssembler, SysmatAssemblerSparseSymm, startassembly!, assemble!, makematrix!, SysmatAssemblerSparse
import FinEtools.FEMMBaseModule: AbstractFEMM
import FinEtools.MatrixUtilityModule: add_mggt_ut_only!, add_nnt_ut_only!, complete_lt!, locjac!
import LinearAlgebra: norm, cross

"""
    FEMMAcoustSurf{S<:AbstractFESet, F<:Function, M, NF<:Function} <: AbstractFEMM

Class for linear acoustics finite element modeling machine.
"""
mutable struct FEMMAcoustSurf{S<:AbstractFESet, F<:Function, M, NF<:Function} <: AbstractFEMM
    integdomain::IntegDomain{S, F} # geometry data
    material::M # material object
    getnormal!::NF # get the  normal to the surface
end


function FEMMAcoustSurf(integdomain::IntegDomain{S, F},  material::M) where {S<:AbstractFESet, F<:Function, M}
    function getnormal!(n::FFltVec, loc::FFltMat, J::FFltMat)
        sdim, mdim = size(J);
        if     mdim == 1 # 1-D fe
        	N = [J[2,1],-J[1,1]];
        elseif     mdim == 2 # 2-D fe
        	N = cross(J[:,1],J[:,2])
        else
        	error("Got an incorrect size of tangents");
        end
        N=N/norm(N)
        copyto!(n, N)
        return n;
    end
    return FEMMAcoustSurf(integdomain, material, getnormal!)
end

"""
    acousticABC(self::FEMMAcoustSurf, assembler::A,
      geom::NodalField,
      Pdot::NodalField{T}) where {T<:Number, A<:AbstractSysmatAssembler}

Compute the acoustic ABC (Absorbing Boundary Condition) matrix.
"""
function acousticABC(self::FEMMAcoustSurf, assembler::A, geom::NodalField, Pdot::NodalField{T}) where {T<:Number, A<:AbstractSysmatAssembler}
    fes = self.integdomain.fes
    # Constants
    nfes = count(fes); # number of finite elements in the set
    ndn = ndofs(Pdot); # number of degrees of freedom per node
    nne =  nodesperelem(fes); # number of nodes per element
    sdim =  ndofs(geom);            # number of space dimensions
    mdim = manifdim(fes);     # manifold dimension of the element
    Dedim = ndn*nne;          # dimension of the element matrix
    # Precompute basis f. values + basis f. gradients wrt parametric coor
    npts, Ns, gradNparams, w, pc  =  integrationdata(self.integdomain);
    # Material
    bulk_modulus  =   bulkmodulus(self.material);
    mass_density  =   massdensity(self.material);
    c  =  sqrt(bulk_modulus/mass_density); # sound speed
    # Prepare assembler and temporaries
    De = fill(zero(FFlt), Dedim, Dedim);                # element matrix -- used as a buffer
    dofnums = fill(zero(FInt), Dedim); # degree of freedom array -- used as a buffer
    loc = fill(zero(FFlt), 1, sdim); # quadrature point location -- used as a buffer
    J = fill(zero(FFlt), sdim, mdim); # Jacobian matrix -- used as a buffer
    startassembly!(assembler, Dedim, Dedim, nfes, Pdot.nfreedofs, Pdot.nfreedofs);
    for i = 1:count(fes) # Loop over elements
        fill!(De, 0.0); # Initialize element matrix
        for j = 1:npts # Loop over quadrature points
            locjac!(loc, J, geom.values, fes.conn[i], Ns[j], gradNparams[j]) 
            Jac = Jacobiansurface(self.integdomain, J, loc, fes.conn[i], Ns[j]);
            ffactor = (Jac/c*w[j])
            add_nnt_ut_only!(De, Ns[j], ffactor)
        end # Loop over quadrature points
        complete_lt!(De)
        gatherdofnums!(Pdot, dofnums, fes.conn[i]);# retrieve degrees of freedom
        assemble!(assembler, De, dofnums, dofnums);# assemble symmetric matrix
    end # Loop over elements
    return makematrix!(assembler);
end

function acousticABC(self::FEMMAcoustSurf, geom::NodalField, Pdot::NodalField{T}) where {T<:Number}
    # Were we supplied assembler object?  If not make a default.
    assembler  =  SysmatAssemblerSparseSymm();
    return acousticABC(self, assembler, geom, Pdot);
end

"""
    pressure2resultantforce(self::FEMMAcoustSurf, assembler::A,
      geom::NodalField,
      P::NodalField{T},
       Force::Field) where {T<:Number, A<:AbstractSysmatAssembler}

Compute the rectangular coupling matrix that transcribes given pressure
on the surface into the resultant force acting on the surface.
"""
function pressure2resultantforce(self::FEMMAcoustSurf, assembler::A, geom::NodalField, P::NodalField{T}, Force::GeneralField) where {T<:Number, A<:AbstractSysmatAssembler}
    fes = self.integdomain.fes
    # Constants
    nfes = count(fes); # number of finite elements in the set
    ndn = ndofs(P); # number of degrees of freedom per node
    nne =  nodesperelem(fes); # number of nodes per element
    sdim =  ndofs(geom);            # number of space dimensions
    mdim = manifdim(fes);     # manifold dimension of the element
    edim = ndn*nne;          # dimension of the element matrix
    # Precompute basis f. values + basis f. gradients wrt parametric coor
    npts, Ns, gradNparams, w, pc  =  integrationdata(self.integdomain);
    transposedNs = [reshape(N, 1, nne) for N in Ns]
    Ge = fill(zero(FFlt), 3, nne); # element coupling matrix -- used as a buffer
    coldofnums = zeros(FInt, 1, edim); # degree of freedom array -- used as a buffer
    rowdofnums = zeros(FInt, 1, 3); # degree of freedom array -- used as a buffer
    loc = fill(zero(FFlt), 1, sdim); # quadrature point location -- used as a buffer
    n = fill(zero(FFlt), 3) # normal vector -- used as a buffer
    J = fill(zero(FFlt), sdim, mdim); # Jacobian matrix -- used as a buffer
    gatherdofnums!(Force, rowdofnums, [1 2 3]);# retrieve degrees of freedom
    startassembly!(assembler, 3, edim, count(fes), 3, P.nfreedofs);
    for i = 1:count(fes) # Loop over elements
        fill!(Ge, 0.0); # Initialize element matrix
        for j = 1:npts # Loop over quadrature points
            locjac!(loc, J, geom.values, fes.conn[i], Ns[j], gradNparams[j]) 
            Jac = Jacobiansurface(self.integdomain, J, loc, fes.conn[i], Ns[j]);
            n = self.getnormal!(n, loc, J);
            ffactor = (Jac*w[j])
            Ge = Ge + (ffactor*n)*transposedNs[j]
        end # Loop over quadrature points
        gatherdofnums!(P, coldofnums, fes.conn[i]);# retrieve degrees of freedom
        assemble!(assembler, Ge, rowdofnums, coldofnums);# assemble unsymmetric matrix
    end # Loop over elements
    return makematrix!(assembler);
end

function pressure2resultantforce(self::FEMMAcoustSurf, geom::NodalField, P::NodalField{T}, Force::GeneralField) where {T<:Number}
    assembler  =  SysmatAssemblerSparse();
    return pressure2resultantforce(self, assembler, geom, P, Force)
end

"""
    pressure2resultanttorque(self::FEMMAcoustSurf, assembler::A,
      geom::NodalField,
      P::NodalField{T},
      Torque::GeneralField, CG::FFltVec) where {T<:Number, A<:AbstractSysmatAssembler}

Compute the rectangular coupling matrix that transcribes given pressure
on the surface into the resultant torque acting on the surface with respect
to the CG.
"""
function pressure2resultanttorque(self::FEMMAcoustSurf, assembler::A, geom::NodalField, P::NodalField{T}, Torque::GeneralField, CG::FFltVec) where {T<:Number,  A<:AbstractSysmatAssembler}
    fes = self.integdomain.fes
    # Constants
    nfes = count(fes); # number of finite elements in the set
    ndn = ndofs(P); # number of degrees of freedom per node
    nne =  nodesperelem(fes); # number of nodes per element
    sdim =  ndofs(geom);            # number of space dimensions
    mdim = manifdim(fes);     # manifold dimension of the element
    edim = ndn*nne;          # dimension of the element matrix
    # Precompute basis f. values + basis f. gradients wrt parametric coor
    npts, Ns, gradNparams, w, pc  =  integrationdata(self.integdomain);
    transposedNs = [reshape(N, 1, nne) for N in Ns]
    Ge = fill(zero(FFlt), 3, nne); # element coupling matrix -- used as a buffer
    coldofnums = zeros(FInt, 1, edim); # degree of freedom array -- used as a buffer
    rowdofnums = zeros(FInt, 1, 3); # degree of freedom array -- used as a buffer
    loc = fill(zero(FFlt), 1, sdim); # quadrature point location -- used as a buffer
    n = fill(zero(FFlt), 3) # normal vector -- used as a buffer
    J = fill(zero(FFlt), sdim, mdim); # Jacobian matrix -- used as a buffer
    gatherdofnums!(Torque, rowdofnums, [1 2 3]);# retrieve degrees of freedom
    startassembly!(assembler, 3, edim, count(fes), 3, P.nfreedofs);
    for i = 1:count(fes) # Loop over elements
        fill!(Ge, 0.0); # Initialize element matrix
        for j = 1:npts # Loop over quadrature points
            locjac!(loc, J, geom.values, fes.conn[i], Ns[j], gradNparams[j]) 
            Jac = Jacobiansurface(self.integdomain, J, loc, fes.conn[i], Ns[j]);
            n = self.getnormal!(n, loc, J);
            ffactor = (Jac*w[j])
            Ge = Ge + (ffactor*cross(vec(vec(loc)-CG), n))*transposedNs[j]
        end # Loop over quadrature points
        gatherdofnums!(P, coldofnums, fes.conn[i]);# retrieve degrees of freedom
        assemble!(assembler, Ge, rowdofnums, coldofnums);# assemble unsymmetric matrix
    end # Loop over elements
    return makematrix!(assembler);
end

function pressure2resultanttorque(self::FEMMAcoustSurf, geom::NodalField, P::NodalField{T}, Torque::GeneralField, CG::FFltVec) where {T<:Number}
    assembler  =  SysmatAssemblerSparse();
    return pressure2resultanttorque(self, assembler, geom, P, Torque, CG)
end

"""
    acousticcouplingpanels(self::FEMMAcoustSurf, geom::NodalField, u::NodalField{T}) where {T}

Compute the acoustic pressure-structure coupling matrix.

The acoustic pressure-nodal force matrix transforms 
the pressure distributed along the surface to forces acting on the nodes
of the finite element model. Its transpose transforms displacements (or velocities, or
accelerations) into the normal component of the displacement (or
velocity, or acceleration) along the surface.

Arguments
`geom`=geometry field
`u` = displacement field 

Notes:
-- `n`=outer normal (pointing into the acoustic medium).
-- The pressures along the surface are assumed constant (uniform) along 
each finite element –- panel. The panel pressures are assumed to
be given the same numbers as the serial numbers of the finite elements in the block.
"""
function acousticcouplingpanels(self::FEMMAcoustSurf, assembler::A, geom::NodalField, u::NodalField{T}) where {A<:AbstractSysmatAssembler, T}
    fes = self.integdomain.fes
    # Constants
    nne =  nodesperelem(fes); # number of nodes per element
    sdim = ndofs(geom);            # number of space dimensions
    mdim = manifdim(fes);     # manifold dimension of the element
    # Precompute basis f. values + basis f. gradients wrt parametric coor
    npts, Ns, gradNparams, w, pc  =  integrationdata(self.integdomain);
    transposedNs = [reshape(N, 1, nne) for N in Ns]
    coldofnums = zeros(FInt, 1, 1); # degree of freedom array -- used as a buffer
    rowdofnums = zeros(FInt, 1, sdim*nne); # degree of freedom array -- used as a buffer
    loc = fill(zero(FFlt), 1, sdim); # quadrature point location -- used as a buffer
    n = fill(zero(FFlt), sdim) # normal vector -- used as a buffer
    J = fill(zero(FFlt), sdim, mdim); # Jacobian matrix -- used as a buffer
    Ge = fill(zero(FFlt), sdim*nne, 1); # Element matrix -- used as a buffer
    startassembly!(assembler, size(Ge, 1), size(Ge, 2), count(fes), u.nfreedofs, count(fes));
    for i = 1:count(fes) # Loop over elements
        fill!(Ge, 0.0); # Initialize element matrix
        for j = 1:npts # Loop over quadrature points
            locjac!(loc, J, geom.values, fes.conn[i], Ns[j], gradNparams[j]) 
            Jac = Jacobiansurface(self.integdomain, J, loc, fes.conn[i], Ns[j]);
            n = self.getnormal!(n, loc, J);
            Ge = Ge + (Jac*w[j])*reshape(reshape(n, sdim, 1)*transposedNs[j], size(Ge, 1), size(Ge, 2))
        end # Loop over quadrature points
        coldofnums[1] = i
        gatherdofnums!(u, rowdofnums, fes.conn[i]);# retrieve degrees of freedom
        assemble!(assembler, Ge, rowdofnums, coldofnums);# assemble unsymmetric matrix
    end # Loop over elements
    return makematrix!(assembler);
end

function acousticcouplingpanels(self::FEMMAcoustSurf, geom::NodalField, u::NodalField{T}) where {T}
    assembler = SysmatAssemblerSparse(); # The matrix is not symmetric
    return acousticcouplingpanels(self, assembler, geom, u)
end

end

