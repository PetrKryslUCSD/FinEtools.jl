module FEMMHeatDiffModule

using FinEtools.FTypesModule
using FinEtools.FESetModule
using FinEtools.FESetModule.gradN!
using FinEtools.CSysModule
using FinEtools.GeoDModule
using FinEtools.FEMMBaseModule
using FinEtools.FieldModule
using FinEtools.NodalFieldModule
using FinEtools.ForceIntensityModule
using FinEtools.MatHeatDiffModule
using FinEtools.AssemblyModule
using FinEtools.MatrixUtilityModule.add_gkgt_ut_only!
using FinEtools.MatrixUtilityModule.complete_lt!
using FinEtools.MatrixUtilityModule.mv_product!

# Class for heat diffusion finite element modeling machine.
type FEMMHeatDiff{S<:FESet, F<:Function, M<:MatHeatDiff} <: FEMMAbstractBase
    geod::GeoD{S, F} # geometry data finite element modeling machine
    material::M # material object
end
export FEMMHeatDiff

function  buffers(self::FEMMHeatDiff, geom::NodalField{FFlt},
  temp::NodalField{FFlt})
  # Constants
  const nfes = count(self.geod.fes); # number of finite elements in the set
  const ndn = ndofs(temp); # number of degrees of freedom per node
  const nne = nodesperelem(self.geod.fes); # number of nodes for element
  const sdim = ndofs(geom);   # number of space dimensions
  const mdim = manifdim(self.geod.fes); # manifold dimension of the element
  const Kedim = ndn*nne;      # dimension of the element matrix
  elmat = zeros(FFlt, Kedim, Kedim); # buffer
  elvec = zeros(FFlt, Kedim); # buffer
  elvecfix = zeros(FFlt, Kedim); # buffer
  conn = zeros(FInt, nne); # buffer
  x = zeros(FFlt, nne, sdim); # buffer
  dofnums = zeros(FInt, 1, Kedim); # buffer
  loc = zeros(FFlt, 1, sdim); # buffer
  J = eye(FFlt, sdim, mdim); # buffer
  RmTJ = zeros(FFlt, mdim, mdim); # buffer
  gradN = zeros(FFlt, nne, mdim); # buffer
  kappa_bargradNT = zeros(FFlt, mdim, nne); # buffer
  return conn, x, dofnums, loc, J, RmTJ, gradN, kappa_bargradNT, elmat, elvec, elvecfix
end

"""
    conductivity(self::FEMMHeatDiff,
      assembler::A, geom::NodalField{FFlt},
      temp::NodalField{FFlt}) where {A<:SysmatAssemblerBase}

Compute the conductivity matrix.
"""
function conductivity(self::FEMMHeatDiff,
  assembler::A, geom::NodalField{FFlt},
  temp::NodalField{FFlt}) where {A<:SysmatAssemblerBase}
  geod = self.geod
  npts,  Ns,  gradNparams,  w,  pc = integrationdata(geod);
  # Thermal conductivity matrix is in local  material coordinates.
  kappa_bar =  self.material.thermal_conductivity;
  # Prepare assembler and buffers
  conn, x, dofnums, loc, J, RmTJ, gradN, kappa_bargradNT, elmat =
            buffers(self, geom, temp)
  startassembly!(assembler, size(elmat,1), size(elmat,2), count(geod.fes),
   temp.nfreedofs, temp.nfreedofs);
  for i = 1:count(geod.fes) # Loop over elements
    getconn!(geod.fes, conn, i);
    gathervalues_asmat!(geom, x, conn);# retrieve element coordinates
    fill!(elmat,  0.0); # Initialize element matrix
    for j=1:npts # Loop over quadrature points
      At_mul_B!(loc, Ns[j], x);# Quadrature point location
      At_mul_B!(J, x, gradNparams[j]); # Jacobian matrix
      Jac = Jacobianvolume(geod, J, loc, conn, Ns[j]);
      updatecsmat!(geod.mcsys, loc, J, geod.fes.label[i]);
      At_mul_B!(RmTJ,  geod.mcsys.csmat,  J); # local Jacobian matrix
      gradN!(geod.fes, gradN, gradNparams[j], RmTJ);
      # Add the product gradN*kappa_bar*gradNT*(Jac*w[j])
      factor::FFlt = (Jac*w[j])
      add_gkgt_ut_only!(elmat, gradN, factor, kappa_bar, kappa_bargradNT)
    end # Loop over quadrature points
    complete_lt!(elmat)
    gatherdofnums!(temp, dofnums, conn);# retrieve degrees of freedom
    assemble!(assembler, elmat, dofnums, dofnums);# assemble symmetric matrix
  end # Loop over elements
  return makematrix!(assembler);
end
export conductivity

function conductivity(self::FEMMHeatDiff,
  geom::NodalField{FFlt},  temp::NodalField{FFlt})
    assembler = SysmatAssemblerSparseSymm();
    return conductivity(self, assembler, geom, temp);
end
export conductivity

"""
    nzebcloadsconductivity(self::FEMMHeatDiff,
      assembler::A,  geom::NodalField{FFlt},
      temp::NodalField{FFlt}) where {A<:SysvecAssemblerBase}

Compute load vector for nonzero EBC of prescribed temperature.
"""
function nzebcloadsconductivity(self::FEMMHeatDiff,
  assembler::A,  geom::NodalField{FFlt},
  temp::NodalField{FFlt}) where {A<:SysvecAssemblerBase}
  geod = self.geod
  npts,  Ns,  gradNparams,  w,  pc = integrationdata(geod);
  # Thermal conductivity matrix is in local  material coordinates.
  kappa_bar = self.material.thermal_conductivity;
  # Prepare assembler and buffers
  conn, x, dofnums, loc, J, RmTJ, gradN, kappa_bargradNT, elmat, elvec, elvecfix =
            buffers(self, geom, temp)
  startassembly!(assembler,  temp.nfreedofs);
  # Now loop over all finite elements in the set
  for i = 1:count(geod.fes) # Loop over elements
    getconn!(geod.fes, conn, i);
    gathervalues_asvec!(temp, elvecfix, conn);# retrieve element coordinates
    if norm(elvecfix) != 0.     # Is the load nonzero?
      gathervalues_asmat!(geom, x, conn);# retrieve element coordinates
      fill!(elmat,  0.0);
      for j=1:npts # Loop over quadrature points
        At_mul_B!(loc, Ns[j], x);# Quadrature point location
        At_mul_B!(J, x, gradNparams[j]); # Jacobian matrix
        Jac = Jacobianvolume(geod, J, loc, conn, Ns[j]);
        updatecsmat!(geod.mcsys, loc, J, geod.fes.label[i]);
        At_mul_B!(RmTJ,  geod.mcsys.csmat,  J); # local Jacobian matrix
        gradN!(geod.fes, gradN, gradNparams[j], RmTJ);
        # Add the product gradN*kappa_bar*gradNT*(Jac*w[j])
        factor::FFlt = (Jac*w[j])
        add_gkgt_ut_only!(elmat, gradN, factor, kappa_bar, kappa_bargradNT)
      end # Loop over quadrature points
      complete_lt!(elmat)
      mv_product!(elvec, elmat, elvecfix) # compute  the load vector
      gatherdofnums!(temp, dofnums, conn); # retrieve degrees of freedom
      assemble!(assembler,  -elvec,  dofnums); # assemble element load vector
    end
  end
  return makevector!(assembler);
end
export nzebcloadsconductivity

function nzebcloadsconductivity(self::FEMMHeatDiff,
  geom::NodalField{FFlt},   temp::NodalField{FFlt})
    assembler = SysvecAssembler()
    return  nzebcloadsconductivity(self, assembler, geom, temp);
end
export nzebcloadsconductivity


end
