module FEMMAcoustModule

import Base.Complex

using FinEtools.FTypesModule
using FinEtools.FESetModule
using FinEtools.GeoDModule
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
  geod::GeoD{S, F} # geometry data finite element modeling machine
  material::M # material object
end
export FEMMAcoust

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
  geod = self.geod
  # Constants
  const nfes = count(geod.fes); # number of finite elements in the set
  const ndn = ndofs(P); # number of degrees of freedom per node
  const nne =  nodesperelem(geod.fes); # number of nodes per element
  const sdim =  ndofs(geom);            # number of space dimensions
  const mdim = manifdim(geod.fes);     # manifold dimension of the element
  const Cedim = ndn*nne;          # dimension of the element matrix
  # Precompute basis f. values + basis f. gradients wrt parametric coor
  npts, Ns, gradNparams, w, pc  =  integrationdata(geod);
  # Prepare assembler and temporaries
  Ce = zeros(FFlt, Cedim, Cedim);                # element matrix -- used as a buffer
  conn = zeros(FInt, nne, 1); # element nodes -- used as a buffer
  x = zeros(FFlt, nne, sdim); # array of node coordinates -- used as a buffer
  dofnums = zeros(FInt, 1, Cedim); # degree of freedom array -- used as a buffer
  loc = zeros(FFlt, 1, sdim); # quadrature point location -- used as a buffer
  J = eye(FFlt, sdim, mdim); # Jacobian matrix -- used as a buffer
  Jinv = eye(FFlt, sdim, mdim); # Jacobian matrix -- used as a buffer
  gradN = zeros(FFlt, nne, mdim); # intermediate result -- used as a buffer
  Jac = 0.0;
  afactor = 0.0;
  startassembly!(assembler, Cedim, Cedim, nfes, P.nfreedofs, P.nfreedofs);
  for i = 1:nfes # Loop over elements
    getconn!(geod.fes, conn, i);# retrieve element node numbers
    gathervalues_asmat!(geom, x, conn);# retrieve element coordinates
    fill!(Ce, 0.0); # Initialize element matrix
    for j = 1:npts # Loop over quadrature points
      At_mul_B!(loc, Ns[j], x);# Quadrature points location
      At_mul_B!(J, x, gradNparams[j]); # calculate the Jacobian matrix
      Jac = Jacobianvolume(geod, J, loc, conn, Ns[j]);
      # gradient WRT global Cartesian coordinates
      FESetModule.gradN!(geod.fes, gradN, gradNparams[j], J);
      afactor = (Jac*w[j]);
      add_mggt_ut_only!(Ce, gradN, afactor)
    end # Loop over quadrature points
    complete_lt!(Ce)
    gatherdofnums!(P, dofnums, conn);# retrieve degrees of freedom
    assemble!(assembler, Ce, dofnums, dofnums);# assemble symmetric matrix
  end # Loop over elements
  return makematrix!(assembler);
end
export acousticmass

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
  geod = self.geod
  # Constants
  const nfes = count(fes); # number of finite elements in the set
  const ndn = ndofs(Pdot); # number of degrees of freedom per node
  const nne =  nodesperelem(fes); # number of nodes per element
  const sdim =  ndofs(geom);            # number of space dimensions
  const mdim = manifdim(fes);     # manifold dimension of the element
  const Cedim = ndn*nne;          # dimension of the element matrix
  # Precompute basis f. values + basis f. gradients wrt parametric coor
  npts, Ns, gradNparams, w, pc  =  FEMMBaseModule.integrationdata(geod);
  # Prepare assembler and temporaries
  Ce::FFltMat  = zeros(FFlt, Cedim, Cedim);                # element matrix -- used as a buffer
  conn::FIntMat = zeros(FInt, nne, 1); # element nodes -- used as a buffer
  x::FFltMat  = zeros(FFlt, nne, sdim); # array of node coordinates -- used as a buffer
  dofnums::FIntMat = zeros(FInt, 1, Cedim); # degree of freedom array -- used as a buffer
  loc::FFltMat  = zeros(FFlt, 1, sdim); # quadrature point location -- used as a buffer
  J::FFltMat  = eye(FFlt, sdim, mdim); # Jacobian matrix -- used as a buffer
  gradN::FFltMat  = zeros(FFlt, nne, mdim); # intermediate result -- used as a buffer
  pP::FFltMat = zeros(FFlt, Cedim, 1);
  startassembly!(assembler, P.nfreedofs);
  # Now loop over all finite elements in the set
  for i = 1:nfes # Loop over elements
    getconn!(fes, conn, i);# retrieve element node numbers
    gathervaluesasmat!(P, pP, conn);# retrieve element coordinates
    if norm(pP) ! =  0     # Is the load nonzero?
      gathervaluesasmat!(geom, x, conn);# retrieve element coordinates
      fill!(Ce, 0.0);
      for j = 1:npts # Loop over quadrature points
        At_mul_B!(loc, Ns[j], x);# Quadrature points location
        At_mul_B!(J, x, gradNparams[j]); # calculate the Jacobian matrix
        Jac  =  FESetModule.Jacobianvolume(fes, conn, Ns[j], J, x);# Jacobian
        # gradient WRT global Cartesian coordinates
        FESetModule.gradN!(fes, gradN, gradNparams[j], J);#Do: gradN  =  gradNparams[j]/J;
        afactor = (Jac*w[j]);
        add_mggt_ut_only!(Ce, gradN, afactor)
      end # Loop over quadrature points
      complete_lt!(Ce)
      gatherdofnumsasvec!(P, dofnums, conn); # retrieve degrees of freedom
      assemble!(assembler, -Ce*pP, dofnums); # assemble element load vector
    end
  end
  F =  makevector!(assembler);
  return F
end
export nzebcloadsacousticmass

function nzebcloadsacousticmass(self::FEMMAcoust,
  geom::NodalField,  P::NodalField{T}) where {T<:Number}
    assembler  =  SysvecAssembler(P.values[1])
    return nzebcloadsacousticmass(self, assembler, geom, P);
end
export nzebcloadsacousticmass

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
  geod = self.geod
  # Constants
  const nfes = count(geod.fes); # number of finite elements in the set
  const ndn = ndofs(Pddot); # number of degrees of freedom per node
  const nne =  nodesperelem(geod.fes); # number of nodes per element
  const sdim =  ndofs(geom);            # number of space dimensions
  const mdim = manifdim(geod.fes);     # manifold dimension of the element
  const Sedim = ndn*nne;          # dimension of the element matrix
  # Precompute basis f. values + basis f. gradients wrt parametric coor
  npts, Ns, gradNparams, w, pc  =  integrationdata(geod);
  # Material
  bulk_modulus  =   self.material.bulk_modulus;
  mass_density  =   self.material.mass_density;
  c  =  sqrt(bulk_modulus/mass_density); # sound speed
  oc2 = 1.0/c^2;
  # Prepare assembler and temporaries
  Se = zeros(FFlt, Sedim, Sedim);                # element matrix -- used as a buffer
  conn = zeros(FInt, nne, 1); # element nodes -- used as a buffer
  x = zeros(FFlt, nne, sdim); # array of node coordinates -- used as a buffer
  dofnums = zeros(FInt, 1, Sedim); # degree of freedom array -- used as a buffer
  loc = zeros(FFlt, 1, sdim); # quadrature point location -- used as a buffer
  J = eye(FFlt, sdim, mdim); # Jacobian matrix -- used as a buffer
  startassembly!(assembler, Sedim, Sedim, nfes, Pddot.nfreedofs, Pddot.nfreedofs);
  for i = 1:nfes # Loop over elements
    getconn!(geod.fes, conn, i);# retrieve element node numbers
    gathervalues_asmat!(geom, x, conn);# retrieve element coordinates
    fill!(Se, 0.0); # Initialize element matrix
    for j = 1:npts # Loop over quadrature points
      At_mul_B!(loc, Ns[j], x);# Quadrature points location
      At_mul_B!(J, x, gradNparams[j]); # calculate the Jacobian matrix
      Jac = Jacobianvolume(geod, J, loc, conn, Ns[j]);
      ffactor = Jac*oc2*w[j]
      add_nnt_ut_only!(Se, Ns[j], ffactor)
    end # Loop over quadrature points
    complete_lt!(Se)
    gatherdofnums!(Pddot, dofnums, conn);# retrieve degrees of freedom
    assemble!(assembler, Se, dofnums, dofnums);# assemble symmetric matrix
  end # Loop over elements
  return makematrix!(assembler);
end
export acousticstiffness

function acousticstiffness(self::FEMMAcoust,
  geom::NodalFieldModule.NodalField,
  Pddot::NodalFieldModule.NodalField{T}) where {T<:Number}
    # Make the default assembler object.
    assembler  =  SysmatAssemblerSparseSymm();
    return acousticstiffness(self, assembler, geom, Pddot);
end
export acousticstiffness

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
  geod = self.geod
  # Constants
  const nfes = count(geod.fes); # number of finite elements in the set
  const ndn = ndofs(Pdot); # number of degrees of freedom per node
  const nne =  nodesperelem(geod.fes); # number of nodes per element
  const sdim =  ndofs(geom);            # number of space dimensions
  const mdim = manifdim(geod.fes);     # manifold dimension of the element
  const Sedim = ndn*nne;          # dimension of the element matrix
  # Precompute basis f. values + basis f. gradients wrt parametric coor
  npts, Ns, gradNparams, w, pc  =  integrationdata(geod);
  # Material
  bulk_modulus  =   self.material.bulk_modulus;
  mass_density  =   self.material.mass_density;
  c  =  sqrt(bulk_modulus/mass_density); # sound speed
  # Prepare assembler and temporaries
  Se::FFltMat  = zeros(FFlt, Sedim, Sedim);                # element matrix -- used as a buffer
  conn::FIntMat = zeros(FInt, nne, 1); # element nodes -- used as a buffer
  x::FFltMat  = zeros(FFlt, nne, sdim); # array of node coordinates -- used as a buffer
  dofnums::FIntMat = zeros(FInt, 1, Sedim); # degree of freedom array -- used as a buffer
  loc::FFltMat  = zeros(FFlt, 1, sdim); # quadrature point location -- used as a buffer
  J::FFltMat  = eye(FFlt, sdim, mdim); # Jacobian matrix -- used as a buffer
  pPddot::FMat{T} = zeros(T, Sedim, 1);
  startassembly!(assembler, Pddot.nfreedofs);
  # Now loop over all finite elements in the set
  for i = 1:nfes # Loop over elements
    getconn!(fes, conn, i);# retrieve element node numbers
    gathervalues_asmat!(Pddot, pPddot, conn);# retrieve element coordinates
    if norm(pPddot) ! =  0     # Is the load nonzero?
      gathervalues_asmat!(geom, x, conn);# retrieve element coordinates
      fill!(Se, 0.0);
      for j = 1:npts # Loop over quadrature points
        At_mul_B!(loc, Ns[j], x);# Quadrature points location
        At_mul_B!(J, x, gradNparams[j]); # calculate the Jacobian matrix
        Jac = Jacobianvolume(geod, J, loc, conn, Ns[j]);
        ffactor = Jac*oc2*w[j]
        add_nnt_ut_only!(Se, Ns[j], ffactor)
      end # Loop over quadrature points
      complete_lt!(Se)
      gatherdofnums!(Pddot, dofnums, conn); # retrieve degrees of freedom
      assemble!(assembler, -Se*pPddot, dofnums); # assemble element load vector
    end
  end
  F =  makevector!(assembler);
  return F
end
export nzebcloadsacousticstiffness

function nzebcloadsacousticstiffness(self::FEMMAcoust,
  geom::NodalFieldModule.NodalField,
  Pddot::NodalFieldModule.NodalField{T}) where {T<:Number}
    assembler  =  SysvecAssembler(Pddot.values[1])#T(0.0)
    return nzebcloadsacousticstiffness(self, assembler, geom, Pddot)
end
export nzebcloadsacousticstiffness

end
