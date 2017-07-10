"""
    FEMMDeforLinearBaseModule

Base module for operations on interiors of domains to construct system matrices and
system vectors for linear deformation models.
"""
module FEMMDeforLinearBaseModule

export FEMMDeforLinearAbstract
export stiffness, nzebcloadsstiffness, mass, thermalstrainloads,
       inspectintegpoints

using FinEtools.FTypesModule
using FinEtools.FESetModule
using FinEtools.CSysModule
using FinEtools.GeoDModule
using FinEtools.FEMMBaseModule
using FinEtools.FieldModule
using FinEtools.NodalFieldModule
using FinEtools.ElementalFieldModule
using FinEtools.ForceIntensityModule
using FinEtools.AssemblyModule
using FinEtools.DeforModelRedModule
using FinEtools.MatDeforModule

abstract type FEMMDeforLinearAbstract <: FEMMAbstractBase end

function buffers(self::FEMMDeforLinearAbstract, geom::NodalField, u::NodalField)
  const ndn = ndofs(u); # number of degrees of freedom per node
  const nne = nodesperelem(self.geod.fes); # number of nodes for element
  const sdim = ndofs(geom);            # number of space dimensions
  const mdim = manifdim(self.geod.fes); # manifold dimension of the element
  const nstrs = nstsstn(self.mr);  # number of stresses
  const elmatdim = ndn*nne;             # dimension of the element matrix
  # Prepare buffers
  elmat = zeros(FFlt, elmatdim, elmatdim);      # element matrix -- buffer
  conn = zeros(FInt, nne, 1); # element nodes -- buffer
  x = zeros(FFlt, nne, sdim); # array of node coordinates -- buffer
  dofnums = zeros(FInt, 1, elmatdim); # degree of freedom array -- buffer
  loc = zeros(FFlt, 1, sdim); # quadrature point location -- buffer
  J = eye(FFlt, sdim, mdim); # Jacobian matrix -- buffer
  csmatTJ = zeros(FFlt, mdim, mdim); # intermediate result -- buffer
  gradN = zeros(FFlt, nne, mdim); # intermediate result -- buffer
  D = zeros(FFlt, nstrs, nstrs); # material stiffness matrix -- buffer
  B = zeros(FFlt, nstrs, elmatdim); # strain-displacement matrix -- buffer
  DB = zeros(FFlt, nstrs, elmatdim); # strain-displacement matrix -- buffer
  elvecfix = zeros(FFlt, elmatdim, 1); # vector of prescribed displ. -- buffer
  elvec = zeros(FFlt, elmatdim); # element vector -- buffer
  return conn, x, dofnums, loc, J, csmatTJ, gradN, D, B, DB, elmat, elvec, elvecfix
end

"""
    stiffness(self::FEMMDeforLinear, assembler::A,
          geom::NodalField{FFlt},
          u::NodalField{T}) where {A<:SysmatAssemblerBase, T<:Number}

Compute and assemble  stiffness matrix.
"""
function stiffness(self::FEMMDeforLinearAbstract, assembler::A,
      geom::NodalField{FFlt},
      u::NodalField{T}) where {A<:SysmatAssemblerBase, T<:Number}
  return spzeros(u.nfreedofs, u.nfreedofs);
end

function stiffness(self::FEMMDeforLinearAbstract,
            geom::NodalField{FFlt},  u::NodalField{T}) where {T<:Number}
  assembler = SysmatAssemblerSparseSymm();
  return stiffness(self, assembler, geom, u);
end

"""
    nzebcloadsstiffness(self::FEMMDeforLinear,  assembler::A,
      geom::NodalField{FFlt},
      u::NodalField{T}) where {A<:SysvecAssemblerBase, T<:Number}

Compute load vector for nonzero EBC for fixed displacement.
"""
function nzebcloadsstiffness(self::FEMMDeforLinearAbstract,  assembler::A,
  geom::NodalField{FFlt},
  u::NodalField{T}) where {A<:SysvecAssemblerBase, T<:Number}
  return zeros(FFlt, u.nfreedofs);
end

function nzebcloadsstiffness(self::FEMMDeforLinearAbstract,
  geom::NodalField{FFlt},
  u::NodalField{T}) where {T<:Number}
  assembler = SysvecAssembler()
  return  nzebcloadsstiffness(self, assembler, geom, u);
end

"""
    mass(self::FEMMDeforLinearAbstract,  assembler::A,
      geom::NodalField{FFlt},
      u::NodalField{T}) where {A<:SysmatAssemblerBase, T<:Number}

Compute the consistent mass matrix

This is a general routine for the abstract linear-deformation  FEMM.
"""
function mass(self::FEMMDeforLinearAbstract,  assembler::A,
  geom::NodalField{FFlt},
  u::NodalField{T}) where {A<:SysmatAssemblerBase, T<:Number}
  geod = self.geod
  npts,  Ns,  gradNparams,  w,  pc = integrationdata(geod);
  conn, x, dofnums, loc, J, csmatTJ, gradN, D, B, DB, elmat =
                buffers(self, geom, u)  # Prepare buffers
  rho::FFlt = self.material.mass_density; # mass density
  NexpTNexp = Array{FFltMat}(1, npts);# basis f. matrix -- buffer
  ndn = ndofs(u)
  for j = 1:npts # This quantity is the same for all quadrature points
    Nexp = zeros(FFlt, ndn, size(elmat,1))
    for l1 = 1:length(conn)
      Nexp[1:ndn, (l1-1)*ndn+1:(l1)*ndn] = eye(ndn)*Ns[j][l1];
    end
    NexpTNexp[j] = Nexp'*Nexp;
  end
  startassembly!(assembler,  size(elmat,1),  size(elmat,2),  count(geod.fes),
  u.nfreedofs,  u.nfreedofs);
  for i = 1:count(geod.fes) # Loop over elements
    getconn!(geod.fes, conn, i);
    gathervalues_asmat!(geom, x, conn);# retrieve element coordinates
    fill!(elmat,  0.0); # Initialize element matrix
    for j = 1:npts # Loop over quadrature points
      At_mul_B!(loc, Ns[j], x);# Quadrature point location
      At_mul_B!(J, x, gradNparams[j]); # calculate the Jacobian matrix
      Jac = Jacobianvolume(geod, J, loc, conn, Ns[j]);
      thefactor::FFlt =(rho*Jac*w[j]);
      @inbounds for nx = 1:size(elmat,1) # Do: Me = Me + (Nexp'*Nexp) * (rho * Jac * w(j));
        @inbounds for mx = 1:size(elmat,2)
          elmat[mx, nx] = elmat[mx, nx] + NexpTNexp[j][mx, nx]*thefactor
        end
      end
    end # Loop over quadrature points
    gatherdofnums!(u,  dofnums,  conn);# retrieve degrees of freedom
    assemble!(assembler,  elmat,  dofnums,  dofnums);# assemble symmetric matrix
  end # Loop over elements
  return makematrix!(assembler);
end

function mass(self::FEMMDeforLinearAbstract,
  geom::NodalField{FFlt},  u::NodalField{T}) where {T<:Number}
  assembler = SysmatAssemblerSparseSymm();
  return mass(self, assembler, geom, u);
end

"""
    thermalstrainloads(self::FEMMDeforLinearAbstract, assembler::A,
        geom::NodalField{FFlt}, u::NodalField{T},
        dT::NodalField{FFlt}) where {A<:SysvecAssemblerBase, T<:Number}

Compute the thermal-strain load vector.
"""
function  thermalstrainloads(self::FEMMDeforLinearAbstract, assembler::A,
    geom::NodalField{FFlt}, u::NodalField{T},
    dT::NodalField{FFlt}) where {A<:SysvecAssemblerBase, T<:Number}
  return zeros(FFlt, u.nfreedofs);
end

function thermalstrainloads(self::FEMMDeforLinearAbstract,
              geom::NodalField{FFlt},
              u::NodalField{T}, dT::NodalField{FFlt}) where {T<:Number}
    assembler = SysvecAssembler()
    return  thermalstrainloads(self, assembler, geom, u, dT);
end

"""
    inspectintegpoints(self::FEMMDeforLinearAbstract,
      geom::NodalField{FFlt},  u::NodalField{T},
      dT::NodalField{FFlt},
      felist::FIntVec,
      inspector::F,  idat, quantity=:Cauchy;
      context...) where {T<:Number, F<:Function}

Inspect integration point quantities.

`geom` - reference geometry field
`u` - displacement field
`dT` - temperature difference field
`felist` - indexes of the finite elements that are to be inspected:
     The fes to be included are: `fes[felist]`.
`context`    - structure: see the update!() method of the material.
`inspector` - functionwith the signature
        idat = inspector(idat, j, conn, x, out, loc);
   where
    `idat` - a structure or an array that the inspector may
           use to maintain some state,  for instance minimum or maximum of
           stress, `j` is the element number, `conn` is the element connectivity,
           `out` is the output of the update!() method,  `loc` is the location
           of the integration point in the *reference* configuration.
### Return
The updated inspector data is returned.
"""
function inspectintegpoints(self::FEMMDeforLinearAbstract,
  geom::NodalField{FFlt},  u::NodalField{T},
  dT::NodalField{FFlt},
  felist::FIntVec,
  inspector::F,  idat, quantity=:Cauchy;
  context...) where {T<:Number, F<:Function}
  return idat; # return the updated inspector data
end

function inspectintegpoints(self::FEMMDeforLinearAbstract,
  geom::NodalField{FFlt},  u::NodalField{T},
  felist::FIntVec,
  inspector::F, idat, quantity=:Cauchy;
  context...) where {T<:Number, F<:Function}
  dT = NodalField(zeros(FFlt, nnodes(geom), 1)) # zero difference in temperature
  return inspectintegpoints(self, geom, u, dT, felist,
            inspector, idat, quantity; context...);
end

end
