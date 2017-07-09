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
    mass(self::FEMMDeforLinear,  assembler::A,
      geom::NodalField{FFlt},  u::NodalField{T}) where {A<:SysmatAssemblerBase, T<:Number}

Compute the mass matrix.
"""
function mass(self::FEMMDeforLinearAbstract,  assembler::A,
  geom::NodalField{FFlt},
  u::NodalField{T}) where {A<:SysmatAssemblerBase, T<:Number}
  return spzeros(u.nfreedofs, u.nfreedofs);
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
