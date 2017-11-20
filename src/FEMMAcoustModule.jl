"""
    FEMMAcoustModule

Module for operations on interiors of domains to construct system matrices and
system vectors for linear acoustics.
"""
module FEMMAcoustModule

export FEMMAcoust
export acousticmass, nzebcloadsacousticmass
export acousticstiffness, nzebcloadsacousticstiffness

import Base.Complex

using FinEtools.FTypesModule
using FinEtools.FESetModule
using FinEtools.IntegDataModule
using FinEtools.FEMMBaseModule
using FinEtools.FieldModule
using FinEtools.NodalFieldModule
using FinEtools.ForceIntensityModule
using FinEtools.MatAcoustFluidModule
using FinEtools.AssemblyModule
using FinEtools.MatrixUtilityModule.add_mggt_ut_only!
using FinEtools.MatrixUtilityModule.add_nnt_ut_only!
using FinEtools.MatrixUtilityModule.complete_lt!

"""
    FEMMAcoust{S<:FESet}

Class for linear acoustics finite element modeling machine.
"""
mutable struct FEMMAcoust{S<:FESet, F<:Function, M} <: FEMMAbstractBase
    integdata::IntegData{S, F} # geometry data finite element modeling machine
    material::M # material object
end

function  buffers(self::FEMMAcoust, geom::NodalField{FFlt}, P::NodalField{F}) where {F}
    integdata = self.integdata
    # Constants
    nfes = count(integdata.fes); # number of finite elements in the set
    ndn = ndofs(P); # number of degrees of freedom per node
    nne =  nodesperelem(integdata.fes); # number of nodes per element
    sdim =  ndofs(geom);            # number of space dimensions
    mdim = manifdim(integdata.fes);     # manifold dimension of the element
    Cedim = ndn*nne;          # dimension of the element matrix
    # Prepare assembler and temporaries
    conn = zeros(FInt, nne, 1); # element nodes -- used as a buffer
    x = zeros(FFlt, nne, sdim); # array of node coordinates -- used as a buffer
    dofnums = zeros(FInt, 1, Cedim); # degree of freedom array -- used as a buffer
    loc = zeros(FFlt, 1, sdim); # quadrature point location -- used as a buffer
    J = eye(FFlt, sdim, mdim); # Jacobian matrix -- used as a buffer
    gradN = zeros(FFlt, nne, mdim); # intermediate result -- used as a buffer
    elmat = zeros(FFlt, Cedim, Cedim);       # element matrix -- used as a buffer
    elvec = zeros(F, Cedim); # buffer
    elvecfix = zeros(F, Cedim); # buffer
    return conn, x, dofnums, loc, J, gradN, elmat, elvec, elvecfix
end

"""
    acousticmass(self::FEMMAcoust,
      assembler::A, geom::NodalFieldModule.NodalField,
      P::NodalFieldModule.NodalField{T}) where {T<:Number, A<:SysmatAssemblerBase}

Compute the acoustic mass matrix.

Return K as a matrix.
Arguments
self   =  acoustics model
assembler  =  matrix assembler
geom = geometry field
P = acoustic (perturbation) pressure field
"""
function acousticmass(self::FEMMAcoust,
    assembler::A, geom::NodalFieldModule.NodalField,
    P::NodalFieldModule.NodalField{T}) where {T<:Number, A<:SysmatAssemblerBase}
    integdata = self.integdata
    conn, x, dofnums, loc, J, gradN, elmat, elvec, elvecfix =
        buffers(self, geom, P)
    # Precompute basis f. values + basis f. gradients wrt parametric coor
    npts, Ns, gradNparams, w, pc  =  integrationdata(integdata);
    Jac = 0.0;
    afactor = 0.0;
    startassembly!(assembler, size(elmat,1), size(elmat,2), count(integdata.fes),
        P.nfreedofs, P.nfreedofs);
    for i = 1:count(integdata.fes) # Loop over elements
        getconn!(integdata.fes, conn, i);# retrieve element node numbers
        gathervalues_asmat!(geom, x, conn);# retrieve element coordinates
        fill!(elmat, 0.0); # Initialize element matrix
        for j = 1:npts # Loop over quadrature points
            At_mul_B!(loc, Ns[j], x);# Quadrature points location
            At_mul_B!(J, x, gradNparams[j]); # calculate the Jacobian matrix
            Jac = Jacobianvolume(integdata, J, loc, conn, Ns[j]);
            # gradient WRT global Cartesian coordinates
            FESetModule.gradN!(integdata.fes, gradN, gradNparams[j], J);
            afactor = (Jac*w[j]);
            add_mggt_ut_only!(elmat, gradN, afactor)
        end # Loop over quadrature points
        complete_lt!(elmat)
        gatherdofnums!(P, dofnums, conn);# retrieve degrees of freedom
        assemble!(assembler, elmat, dofnums, dofnums);# assemble symmetric matrix
    end # Loop over elements
    return makematrix!(assembler);
end

function acousticmass(self::FEMMAcoust,
    geom::NodalFieldModule.NodalField,
    P::NodalFieldModule.NodalField{T}) where {T<:Number}
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
function nzebcloadsacousticmass(self::FEMMAcoust, assembler::A,
    geom::NodalField, P::NodalField{T}) where {T<:Number,
    A<:SysvecAssemblerBase}
    integdata = self.integdata
    conn, x, dofnums, loc, J, gradN, elmat, elvec, elvecfix =
        buffers(self, geom, P)
    # Precompute basis f. values + basis f. gradients wrt parametric coor
    npts, Ns, gradNparams, w, pc  =  integrationdata(integdata);
    startassembly!(assembler, P.nfreedofs);
    # Now loop over all finite elements in the set
    for i = 1:count(integdata.fes) # Loop over elements
        getconn!(integdata.fes, conn, i);# retrieve element node numbers
        gatherfixedvalues_asvec!(P, elvecfix, conn);# retrieve element coordinates
        if norm(elvecfix) !=  0.0     # Is the load nonzero?
            gathervalues_asmat!(geom, x, conn);# retrieve element coordinates
            fill!(elmat, 0.0);
            for j = 1:npts # Loop over quadrature points
                At_mul_B!(loc, Ns[j], x);# Quadrature points location
                At_mul_B!(J, x, gradNparams[j]); # calculate the Jacobian matrix
                Jac = Jacobianvolume(integdata, J, loc, conn, Ns[j]);
                # gradient WRT global Cartesian coordinates
                FESetModule.gradN!(integdata.fes, gradN, gradNparams[j], J);
                afactor = (Jac*w[j]);
                add_mggt_ut_only!(elmat, gradN, afactor)
            end # Loop over quadrature points
            complete_lt!(elmat)
            gatherdofnums!(P, dofnums, conn);# retrieve degrees of freedom
            assemble!(assembler, -elmat*elvecfix, dofnums); # assemble element load vector
        end
    end
    return  makevector!(assembler);
end

function nzebcloadsacousticmass(self::FEMMAcoust,
    geom::NodalField,  P::NodalField{T}) where {T<:Number}
    assembler  =  SysvecAssembler(P.values[1])
    return nzebcloadsacousticmass(self, assembler, geom, P);
end

"""
    acousticstiffness(self::FEMMAcoust, assembler::A,
      geom::NodalFieldModule.NodalField,
      Pddot::NodalFieldModule.NodalField{T}) where {T<:Number,
      A<:SysmatAssemblerBase}

Compute the acoustic stiffness matrix.
"""
function acousticstiffness(self::FEMMAcoust, assembler::A,
    geom::NodalFieldModule.NodalField,
    Pddot::NodalFieldModule.NodalField{T}) where {T<:Number,
    A<:SysmatAssemblerBase}
    integdata = self.integdata
    conn, x, dofnums, loc, J, gradN, elmat, elvec, elvecfix =
        buffers(self, geom, Pddot)
    # Precompute basis f. values + basis f. gradients wrt parametric coor
    npts, Ns, gradNparams, w, pc  =  integrationdata(integdata);
    # Material
    bulk_modulus  =   self.material.bulk_modulus;
    mass_density  =   self.material.mass_density;
    c  =  sqrt(bulk_modulus/mass_density); # sound speed
    oc2 = 1.0/c^2;
    startassembly!(assembler, size(elmat,1), size(elmat,2), count(integdata.fes),
        Pddot.nfreedofs, Pddot.nfreedofs);
    for i = 1:count(integdata.fes) # Loop over elements
        getconn!(integdata.fes, conn, i);# retrieve element node numbers
        gathervalues_asmat!(geom, x, conn);# retrieve element coordinates
        fill!(elmat, 0.0); # Initialize element matrix
        for j = 1:npts # Loop over quadrature points
            At_mul_B!(loc, Ns[j], x);# Quadrature points location
            At_mul_B!(J, x, gradNparams[j]); # calculate the Jacobian matrix
            Jac = Jacobianvolume(integdata, J, loc, conn, Ns[j]);
            ffactor = Jac*oc2*w[j]
            add_nnt_ut_only!(elmat, Ns[j], ffactor)
        end # Loop over quadrature points
        complete_lt!(elmat)
        gatherdofnums!(Pddot, dofnums, conn);# retrieve degrees of freedom
        assemble!(assembler, elmat, dofnums, dofnums);# assemble symmetric matrix
    end # Loop over elements
    return makematrix!(assembler);
end

function acousticstiffness(self::FEMMAcoust,
    geom::NodalFieldModule.NodalField,
    Pddot::NodalFieldModule.NodalField{T}) where {T<:Number}
    # Make the default assembler object.
    assembler  =  SysmatAssemblerSparseSymm();
    return acousticstiffness(self, assembler, geom, Pddot);
end

"""
    nzebcloadsacousticstiffness(self::FEMMAcoust, assembler::A,
      geom::NodalFieldModule.NodalField,
      Pddot::NodalFieldModule.NodalField{T}) where {T<:Number,
      A<:SysvecAssemblerBase}

Compute load vector for nonzero EBC for fixed second-order pressure rate.
"""
function nzebcloadsacousticstiffness(self::FEMMAcoust, assembler::A,
    geom::NodalFieldModule.NodalField,
    Pddot::NodalFieldModule.NodalField{T}) where {T<:Number,
    A<:SysvecAssemblerBase}
    integdata = self.integdata
    conn, x, dofnums, loc, J, gradN, elmat, elvec, elvecfix =
        buffers(self, geom, Pddot)
    # Precompute basis f. values + basis f. gradients wrt parametric coor
    npts, Ns, gradNparams, w, pc  =  integrationdata(integdata);
    # Material
    bulk_modulus  =   self.material.bulk_modulus;
    mass_density  =   self.material.mass_density;
    c  =  sqrt(bulk_modulus/mass_density); # sound speed
    oc2 = 1.0/c^2;
    startassembly!(assembler, Pddot.nfreedofs);
    # Now loop over all finite elements in the set
    for i = 1:count(integdata.fes) # Loop over elements
        getconn!(integdata.fes, conn, i);# retrieve element node numbers
        gatherfixedvalues_asvec!(Pddot, elvecfix, conn);# retrieve element coordinates
        if norm(elvecfix) !=  0.0  # Is the load nonzero?
            gathervalues_asmat!(geom, x, conn);# retrieve element coordinates
            fill!(elmat, 0.0);
            for j = 1:npts # Loop over quadrature points
                At_mul_B!(loc, Ns[j], x);# Quadrature point location
                At_mul_B!(J, x, gradNparams[j]); # calculate the Jacobian matrix
                Jac = Jacobianvolume(integdata, J, loc, conn, Ns[j]);
                ffactor = Jac*oc2*w[j]
                add_nnt_ut_only!(elmat, Ns[j], ffactor)
            end # Loop over quadrature points
            complete_lt!(elmat)
            gatherdofnums!(Pddot, dofnums, conn); # retrieve degrees of freedom
            assemble!(assembler, -elmat*elvecfix, dofnums); # assemble element load vector
        end
    end
    return makevector!(assembler);
end

function nzebcloadsacousticstiffness(self::FEMMAcoust,
    geom::NodalFieldModule.NodalField,
    Pddot::NodalFieldModule.NodalField{T}) where {T<:Number}
    assembler  =  SysvecAssembler(Pddot.values[1])#T(0.0)
    return nzebcloadsacousticstiffness(self, assembler, geom, Pddot)
end

end
