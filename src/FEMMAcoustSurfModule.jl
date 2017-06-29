module FEMMAcoustSurfModule

import Base.Complex

using FinEtools.FTypesModule
using FinEtools.FESetModule
using FinEtools.GeoDModule
using FinEtools.FEMMBaseModule
using FinEtools.FieldModule
using FinEtools.NodalFieldModule
using FinEtools.GeneralFieldModule
using FinEtools.ForceIntensityModule
using FinEtools.RotationUtilModule
using FinEtools.AssemblyModule
using FinEtools.MatrixUtilityModule.add_mggt_ut_only!
using FinEtools.MatrixUtilityModule.add_nnt_ut_only!
using FinEtools.MatrixUtilityModule.complete_lt!

"""
    FEMMAcoustSurf{S<:FESet, F<:Function, M} <: FEMMAbstractBase

Class for linear acoustics finite element modeling machine.
"""
type FEMMAcoustSurf{S<:FESet, F<:Function, M, NF<:Function} <: FEMMAbstractBase
  geod::GeoD{S, F} # geometry data finite element modeling machine
  material::M # material object
  getnormal!::NF # get the  normal to the surface
end
export FEMMAcoustSurf

function FEMMAcoustSurf(geod::GeoD{S, F},
  material::M) where {S<:FESet, F<:Function, M}
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
    copy!(n, N)
    return n;
  end
  return FEMMAcoustSurf(geod, material, getnormal!)
end

"""
    acousticABC(self::FEMMAcoustSurf, assembler::A,
      geom::NodalFieldModule.NodalField,
      Pdot::NodalFieldModule.NodalField{T}) where {T<:Number, A<:SysmatAssemblerBase}

Compute the acoustic ABC (Absorbing Boundary Condition) matrix.
"""
function acousticABC(self::FEMMAcoustSurf, assembler::A,
  geom::NodalFieldModule.NodalField,
  Pdot::NodalFieldModule.NodalField{T}) where {T<:Number, A<:SysmatAssemblerBase}
  geod = self.geod
  # Constants
  const nfes = count(geod.fes); # number of finite elements in the set
  const ndn = ndofs(Pdot); # number of degrees of freedom per node
  const nne =  nodesperelem(geod.fes); # number of nodes per element
  const sdim =  ndofs(geom);            # number of space dimensions
  const mdim = manifdim(geod.fes);     # manifold dimension of the element
  const Dedim = ndn*nne;          # dimension of the element matrix
  # Precompute basis f. values + basis f. gradients wrt parametric coor
  npts, Ns, gradNparams, w, pc  =  integrationdata(geod);
  # Material
  bulk_modulus  =   self.material.bulk_modulus;
  mass_density  =   self.material.mass_density;
  c  =  sqrt(bulk_modulus/mass_density); # sound speed
  # Prepare assembler and temporaries
  De = zeros(FFlt, Dedim, Dedim);                # element matrix -- used as a buffer
  conn = zeros(FInt, nne, 1); # element nodes -- used as a buffer
  x = zeros(FFlt, nne, sdim); # array of node coordinates -- used as a buffer
  dofnums = zeros(FInt, 1, Dedim); # degree of freedom array -- used as a buffer
  loc = zeros(FFlt, 1, sdim); # quadrature point location -- used as a buffer
  J = eye(FFlt, sdim, mdim); # Jacobian matrix -- used as a buffer
  startassembly!(assembler, Dedim, Dedim, nfes, Pdot.nfreedofs, Pdot.nfreedofs);
  for i = 1:count(geod.fes) # Loop over elements
    getconn!(geod.fes, conn, i);# retrieve element node numbers
    gathervalues_asmat!(geom, x, conn);# retrieve element coordinates
    fill!(De, 0.0); # Initialize element matrix
    for j = 1:npts # Loop over quadrature points
      At_mul_B!(loc, Ns[j], x);# Quadrature points location
      At_mul_B!(J, x, gradNparams[j]); # calculate the Jacobian matrix
      Jac = Jacobiansurface(geod, J, loc, conn, Ns[j]);
      ffactor = (Jac/c*w[j])
      add_nnt_ut_only!(De, Ns[j], ffactor)
    end # Loop over quadrature points
    complete_lt!(De)
    gatherdofnums!(Pdot, dofnums, conn);# retrieve degrees of freedom
    assemble!(assembler, De, dofnums, dofnums);# assemble symmetric matrix
  end # Loop over elements
  return makematrix!(assembler);
end
export acousticABC

function acousticABC(self::FEMMAcoustSurf,
  geom::NodalFieldModule.NodalField,
  Pdot::NodalFieldModule.NodalField{T}) where {T<:Number}
  # Were we supplied assembler object?  If not make a default.
  assembler  =  SysmatAssemblerSparseSymm();
  return acousticABC(self, assembler, geom, Pdot);
end
export acousticABC

"""
    pressure2resultantforce(self::FEMMAcoustSurf, assembler::A,
      geom::NodalFieldModule.NodalField,
      P::NodalFieldModule.NodalField{T},
       Force::Field) where {T<:Number, A<:SysmatAssemblerBase}

Compute.

Compute the rectangular coupling matrix that transcribes given pressure
on the surface into the resultant force acting on the surface.
"""
function pressure2resultantforce(self::FEMMAcoustSurf, assembler::A,
  geom::NodalFieldModule.NodalField,
  P::NodalFieldModule.NodalField{T},
   Force::GeneralField) where {T<:Number, A<:SysmatAssemblerBase}
  geod = self.geod
  # Constants
  const nfes = count(geod.fes); # number of finite elements in the set
  const ndn = ndofs(P); # number of degrees of freedom per node
  const nne =  nodesperelem(geod.fes); # number of nodes per element
  const sdim =  ndofs(geom);            # number of space dimensions
  const mdim = manifdim(geod.fes);     # manifold dimension of the element
  const edim = ndn*nne;          # dimension of the element matrix
  # Precompute basis f. values + basis f. gradients wrt parametric coor
  npts, Ns, gradNparams, w, pc  =  integrationdata(geod);
  Ge = zeros(FFlt, 3, nne); # element coupling matrix -- used as a buffer
  conn = zeros(FInt, nne, 1); # element nodes -- used as a buffer
  x = zeros(FFlt, nne, sdim); # array of node coordinates -- used as a buffer
  coldofnums = zeros(FInt, 1, edim); # degree of freedom array -- used as a buffer
  rowdofnums = zeros(FInt, 1, 3); # degree of freedom array -- used as a buffer
  loc = zeros(FFlt, 1, sdim); # quadrature point location -- used as a buffer
  n = zeros(FFlt, 3) # normal vector -- used as a buffer
  J = eye(FFlt, sdim, mdim); # Jacobian matrix -- used as a buffer
  gatherdofnums!(Force, rowdofnums, [1 2 3]);# retrieve degrees of freedom
  startassembly!(assembler, 3, edim, count(geod.fes), 3, P.nfreedofs);
  for i = 1:count(geod.fes) # Loop over elements
    getconn!(geod.fes, conn, i);# retrieve element node numbers
    gathervalues_asmat!(geom, x, conn);# retrieve element coordinates
    fill!(Ge, 0.0); # Initialize element matrix
    for j = 1:npts # Loop over quadrature points
      At_mul_B!(loc, Ns[j], x);# Quadrature points location
      At_mul_B!(J, x, gradNparams[j]); # calculate the Jacobian matrix
      Jac = Jacobiansurface(geod, J, loc, conn, Ns[j]);
      n = self.getnormal!(n, loc, J);
      ffactor = (Jac*w[j])
      Ge = Ge + (ffactor*n)*transpose(Ns[j])
    end # Loop over quadrature points
    gatherdofnums!(P, coldofnums, conn);# retrieve degrees of freedom
    assemble!(assembler, Ge, rowdofnums, coldofnums);# assemble unsymmetric matrix
  end # Loop over elements
  return makematrix!(assembler);
end
export pressure2resultantforce

function pressure2resultantforce(self::FEMMAcoustSurf,
  geom::NodalFieldModule.NodalField,
  P::NodalFieldModule.NodalField{T},
  Force::GeneralField) where {T<:Number}
  assembler  =  SysmatAssemblerSparse();
  return pressure2resultantforce(self, assembler, geom, P, Force)
end

"""
    pressure2resultanttorque(self::FEMMAcoustSurf, assembler::A,
      geom::NodalFieldModule.NodalField,
      P::NodalFieldModule.NodalField{T},
      Torque::GeneralField, CG::FFltVec) where {T<:Number, A<:SysmatAssemblerBase}

Compute the rectangular coupling matrix that transcribes given pressure
on the surface into the resultant torque acting on the surface with respect
to the CG.
"""
function pressure2resultanttorque(self::FEMMAcoustSurf, assembler::A,
  geom::NodalFieldModule.NodalField,
  P::NodalFieldModule.NodalField{T},
  Torque::GeneralField, CG::FFltVec) where {T<:Number, A<:SysmatAssemblerBase}
  geod = self.geod
  # Constants
  const nfes = count(geod.fes); # number of finite elements in the set
  const ndn = ndofs(P); # number of degrees of freedom per node
  const nne =  nodesperelem(geod.fes); # number of nodes per element
  const sdim =  ndofs(geom);            # number of space dimensions
  const mdim = manifdim(geod.fes);     # manifold dimension of the element
  const edim = ndn*nne;          # dimension of the element matrix
  # Precompute basis f. values + basis f. gradients wrt parametric coor
  npts, Ns, gradNparams, w, pc  =  integrationdata(geod);
  Ge = zeros(FFlt, 3, nne); # element coupling matrix -- used as a buffer
  conn = zeros(FInt, nne, 1); # element nodes -- used as a buffer
  x = zeros(FFlt, nne, sdim); # array of node coordinates -- used as a buffer
  coldofnums = zeros(FInt, 1, edim); # degree of freedom array -- used as a buffer
  rowdofnums = zeros(FInt, 1, 3); # degree of freedom array -- used as a buffer
  loc = zeros(FFlt, 1, sdim); # quadrature point location -- used as a buffer
  n = zeros(FFlt, 3) # normal vector -- used as a buffer
  J = eye(FFlt, sdim, mdim); # Jacobian matrix -- used as a buffer
  gatherdofnums!(Torque, rowdofnums, [1 2 3]);# retrieve degrees of freedom
  startassembly!(assembler, 3, edim, count(geod.fes), 3, P.nfreedofs);
  for i = 1:count(geod.fes) # Loop over elements
    getconn!(geod.fes, conn, i);# retrieve element node numbers
    gathervalues_asmat!(geom, x, conn);# retrieve element coordinates
    fill!(Ge, 0.0); # Initialize element matrix
    for j = 1:npts # Loop over quadrature points
      At_mul_B!(loc, Ns[j], x);# Quadrature points location
      At_mul_B!(J, x, gradNparams[j]); # calculate the Jacobian matrix
      Jac = Jacobiansurface(geod, J, loc, conn, Ns[j]);
      n = self.getnormal!(n, loc, J);
      ffactor = (Jac*w[j])
      Ge = Ge + (ffactor*cross(vec(vec(loc)-CG), n))*transpose(Ns[j])
    end # Loop over quadrature points
    gatherdofnums!(P, coldofnums, conn);# retrieve degrees of freedom
    assemble!(assembler, Ge, rowdofnums, coldofnums);# assemble unsymmetric matrix
  end # Loop over elements
  return makematrix!(assembler);
end
export pressure2resultanttorque

function pressure2resultanttorque(self::FEMMAcoustSurf,
  geom::NodalFieldModule.NodalField,
  P::NodalFieldModule.NodalField{T},
  Torque::GeneralField, CG::FFltVec) where {T<:Number}
  assembler  =  SysmatAssemblerSparse();
  return pressure2resultanttorque(self, assembler, geom, P, Torque, CG)
end

end
