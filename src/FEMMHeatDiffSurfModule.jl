"""
    FEMMHeatDiffSurfModule

Module for operations on boundaries of domains to construct system matrices and
system vectors for linear heat diffusion/conduction.
"""
module FEMMHeatDiffSurfModule

export  FEMMHeatDiffSurf
export surfacetransfer,  surfacetransferloads,  nzebcsurfacetransferloads

using FinEtools.FTypesModule
using FinEtools.FESetModule
using FinEtools.CSysModule
using FinEtools.GeoDModule
using FinEtools.FEMMBaseModule
using FinEtools.FieldModule
using FinEtools.NodalFieldModule
using FinEtools.ForceIntensityModule
using FinEtools.MatHeatDiffModule
using FinEtools.AssemblyModule
using FinEtools.MatrixUtilityModule.add_nnt_ut_only!
using FinEtools.MatrixUtilityModule.complete_lt!

# Class for heat diffusion finite element modeling machine.
mutable struct FEMMHeatDiffSurf{S<:FESet, F<:Function} <: FEMMAbstractBase
  geod::GeoD{S, F} # geometry data finite element modeling machine
  surfacetransfercoeff::FFlt # material object
end

"""
    surfacetransfer(self::FEMMHeatDiffSurf,  assembler::A,
      geom::NodalField{FFlt}, temp::NodalField{FFlt}) where {A<:SysmatAssemblerBase}

Compute the surface heat transfer matrix.
"""
function surfacetransfer(self::FEMMHeatDiffSurf,  assembler::A,
  geom::NodalField{FFlt}, temp::NodalField{FFlt}) where {A<:SysmatAssemblerBase}
  geod = self.geod
  # Constants
  const nfes = count(geod.fes); # number of finite elements in the set
  const ndn = ndofs(temp); # number of degrees of freedom per node
  const nne = nodesperelem(geod.fes); # number of nodes for element
  const sdim = ndofs(geom);            # number of space dimensions
  const mdim = manifdim(geod.fes); # manifold dimension of the element
  const Hedim = ndn*nne;             # dimension of the element matrix
  # Precompute basis f. values + basis f. gradients wrt parametric coor
  npts,  Ns,  gradNparams,  w,  pc = integrationdata(geod);
  # Prepare assembler and temporaries
  He = zeros(FFlt, Hedim, Hedim);                # element matrix -- used as a buffer
  conn = zeros(FInt, nne, 1); # element nodes -- used as a buffer
  x = zeros(FFlt, nne, sdim); # array of node coordinates -- used as a buffer
  dofnums = zeros(FInt, 1, Hedim); # degree of freedom array -- used as a buffer
  loc = zeros(FFlt, 1, sdim); # quadrature point location -- used as a buffer
  J = eye(FFlt, sdim, mdim); # Jacobian matrix -- used as a buffer
  startassembly!(assembler, Hedim, Hedim, nfes, temp.nfreedofs, temp.nfreedofs);
  for i = 1:nfes # Loop over elements
    getconn!(geod.fes, conn, i);
    gathervalues_asmat!(geom, x, conn);# retrieve element coordinates
    fill!(He,  0.0); # Initialize element matrix
    for j=1:npts # Loop over quadrature points
      At_mul_B!(loc, Ns[j], x);# Quadrature points location
      At_mul_B!(J,  x,  gradNparams[j]); # calculate the Jacobian matrix
      Jac = Jacobiansurface(geod, J, loc, conn,  Ns[j]);
      add_nnt_ut_only!(He, Ns[j], self.surfacetransfercoeff*Jac*w[j])
    end # Loop over quadrature points
    complete_lt!(He)
    gatherdofnums!(temp, dofnums, conn);# retrieve degrees of freedom
    assemble!(assembler, He, dofnums, dofnums);# assemble symmetric matrix
  end # Loop over elements
  return makematrix!(assembler);
end

function surfacetransfer(self::FEMMHeatDiffSurf{S},
                                   geom::NodalField{FFlt},
                                   temp::NodalField{FFlt}) where {S<:FESet}
    assembler = SysmatAssemblerSparseSymm()
    return  surfacetransfer(self, assembler, geom, temp);
end

"""
    surfacetransferloads(self::FEMMHeatDiffSurf,  assembler::A,
      geom::NodalField{FFlt}, temp::NodalField{FFlt},
      ambtemp::NodalField{FFlt}) where {A<:SysvecAssemblerBase}

Compute the load vector corresponding to surface heat transfer.
"""
function surfacetransferloads(self::FEMMHeatDiffSurf,  assembler::A,
  geom::NodalField{FFlt}, temp::NodalField{FFlt},
  ambtemp::NodalField{FFlt}) where {A<:SysvecAssemblerBase}
  geod = self.geod
  # Constants
  const nfes = count(geod.fes); # number of finite elements in the set
  const ndn = ndofs(temp); # number of degrees of freedom per node
  const nne = nodesperelem(geod.fes); # number of nodes for element
  const sdim = ndofs(geom);            # number of space dimensions
  const mdim = manifdim(geod.fes); # manifold dimension of the element
  const Hedim = ndn*nne;             # dimension of the element matrix
  # Precompute basis f. values + basis f. gradients wrt parametric coor
  npts,  Ns,  gradNparams,  w,  pc = integrationdata(geod);
  # Prepare assembler and temporaries
  Fe = zeros(FFlt, Hedim, 1); # element matrix -- used as a buffer
  conn = zeros(FInt, nne, 1); # element nodes -- used as a buffer
  x = zeros(FFlt, nne, sdim); # array of node coordinates -- used as a buffer
  dofnums = zeros(FInt, 1, Hedim); # degree of freedom array -- used as a buffer
  loc = zeros(FFlt, 1, sdim); # quadrature point location -- used as a buffer
  J = eye(FFlt, sdim, mdim); # Jacobian matrix -- used as a buffer
  pT = zeros(FFlt, Hedim, 1);
  startassembly!(assembler,  temp.nfreedofs);
  for i = 1:nfes # Loop over elements
    getconn!(geod.fes, conn, i);
    gathervalues_asvec!(ambtemp, pT, conn);# retrieve element coordinates
    if norm(pT) != 0.0    # Is the load nonzero?
      gathervalues_asmat!(geom, x, conn);# retrieve element coordinates
      fill!(Fe,  0.0); # Initialize element matrix
      for j=1:npts # Loop over quadrature points
        At_mul_B!(loc, Ns[j], x);# Quadrature points location
        At_mul_B!(J,  x,  gradNparams[j]); # calculate the Jacobian matrix
        Jac = Jacobiansurface(geod, J, loc, conn,  Ns[j]);
        Ta = dot(pT, Ns[j])
        factor = Ta*self.surfacetransfercoeff*Jac*w[j]
        @inbounds for  mx = 1:Hedim
          Fe[mx] += factor*Ns[j][mx]
        end
      end # Loop over quadrature points
      gatherdofnums!(temp, dofnums, conn); # retrieve degrees of freedom
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
      geom::NodalField{FFlt}, temp::NodalField{FFlt}) where {A<:SysvecAssemblerBase}

Compute load vector for nonzero EBC for fixed temperature.
"""
function nzebcsurfacetransferloads(self::FEMMHeatDiffSurf, assembler::A,
  geom::NodalField{FFlt}, temp::NodalField{FFlt}) where {A<:SysvecAssemblerBase}
  geod = self.geod
  # Constants
  const nfes = count(geod.fes); # number of finite elements in the set
  const ndn = ndofs(temp); # number of degrees of freedom per node
  const nne = nodesperelem(geod.fes); # number of nodes for element
  const sdim = ndofs(geom);            # number of space dimensions
  const mdim = manifdim(geod.fes); # manifold dimension of the element
  const Hedim = ndn*nne;             # dimension of the element matrix
  # Precompute basis f. values + basis f. gradients wrt parametric coor
  npts,  Ns,  gradNparams,  w,  pc = integrationdata(geod);
  # Prepare assembler and temporaries
  He = zeros(FFlt, Hedim, Hedim);                # element matrix -- used as a buffer
  conn = zeros(FInt, nne, 1); # element nodes -- used as a buffer
  x = zeros(FFlt, nne, sdim); # array of node coordinates -- used as a buffer
  dofnums = zeros(FInt, 1, Hedim); # degree of freedom array -- used as a buffer
  loc = zeros(FFlt, 1, sdim); # quadrature point location -- used as a buffer
  J = eye(FFlt, sdim, mdim); # Jacobian matrix -- used as a buffer
  pT = zeros(FFlt, Hedim, 1);
  startassembly!(assembler,  temp.nfreedofs);
  # Now loop over all finite elements in the set
  for i=1:nfes # Loop over elements
    getconn!(geod.fes, conn, i);
    gathervalues_asvec!(temp, pT, conn);# retrieve element coordinates
    if norm(pT) != 0.0    # Is the load nonzero?
      gathervalues_asmat!(geom, x, conn);# retrieve element coordinates
      fill!(He,  0.0);
      for j=1:npts # Loop over quadrature points
        At_mul_B!(loc, Ns[j], x);# Quadrature points location
        At_mul_B!(J,  x,  gradNparams[j]); # calculate the Jacobian matrix
        Jac = Jacobiansurface(geod, J, loc, conn,  Ns[j]);
        add_nnt_ut_only!(He, Ns[j], self.surfacetransfercoeff*Jac*w[j])
      end # Loop over quadrature points
      complete_lt!(He)
      gatherdofnums!(temp, dofnums, conn); # retrieve degrees of freedom
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
