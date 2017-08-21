"""
    FEMMBaseModule

Module for comments/base operations on interiors and boundaries of domains.
"""
module FEMMBaseModule

export FEMMAbstractBase, FEMMBase
export associategeometry!, integratefieldfunction, integratefunction, distribloads
export connectionmatrix, fieldfromintegpoints, elemfieldfromintegpoints

using FinEtools
using FinEtools.FTypesModule
using FinEtools.FESetModule
using FinEtools.IntegRuleModule
using FinEtools.NodalFieldModule
using FinEtools.ForceIntensityModule
using FinEtools.CSysModule
using FinEtools.GeoDModule
using FinEtools.AssemblyModule

"""
    FEMMAbstractBase

Abstract base type for all finite element models.
"""
abstract type FEMMAbstractBase; end

"""
    FEMMBase

Class for base finite element modeling machine.
"""
mutable struct FEMMBase{S<:FESet, F<:Function} <: FEMMAbstractBase
  geod::GeoD{S, F} # geometry data finite element modeling machine
end

"""
    associategeometry!(self::FEMMAbstractBase,  geom::NodalField{FFlt})

Associate geometry field with the FEMM.

There may be operations that could benefit from pre-computations
that involve a geometry field. If so, associating the geometry
field gives the FEMM a chance to save on repeated computations.

Geometry field is normally passed into any routine that evaluates some forms
(integrals) over the mesh.  Whenever the geometry passed into a routine is not
consistent with the one for which `associategeometry!()` was called before,
`associategeometry!()` needs to be called with the new geometry field.
"""
function associategeometry!(self::FEMMAbstractBase,  geom::NodalField{FFlt})
  return self # default is no-op
end

"""
    integratefieldfunction(self::FEMMAbstractBase,
      geom::NodalField{FFlt},  afield::NodalField{T},
      fh::F,  initial::R;  m::FInt=-1)  where {T<:Number,R, F<:Function}

Integrate a nodal-field function over the discrete manifold.
"""
function integratefieldfunction(self::FEMMAbstractBase,
  geom::NodalField{FFlt},  afield::NodalField{T},
  fh::F,  initial::R;  m::FInt=-1)  where {T<:Number, R, F<:Function}
  geod = self.geod                # finite elements
  # Constants
  nfes = count(geod.fes); # number of finite elements in the set
  ndn = ndofs(afield); # number of degrees of freedom per node
  nne = nodesperelem(geod.fes); # number of nodes per element
  sdim = ndofs(geom);            # number of space dimensions
  mdim = manifdim(geod.fes);     # manifold dimension of the element
  # Precompute basis f. values + basis f. gradients wrt parametric coor
  npts, Ns, gradNparams, w, pc = integrationdata(geod);
  conn = zeros(FInt,nne,1); # element nodes -- used as a buffer
  x = zeros(FFlt,nne,sdim); # array of node coordinates -- used as a buffer
  a = zeros(FFlt,nne,ndn); # array of field DOFS-- used as a buffer
  loc = zeros(FFlt,1,sdim); # quadrature point location -- used as a buffer
  val = zeros(FFlt,1,ndn); # field value at the point -- used as a buffer
  J = eye(FFlt,sdim,mdim); # Jacobian matrix -- used as a buffer
  if m >= 0
    # Either the manifold dimension was supplied
  else
    m=mdim;# ...Or it is implied
  end
  result = initial;           # initial value for the result
  for i=1:count(geod.fes) #Now loop over all fes in the block
    getconn!(geod.fes,conn,i);
    gathervalues_asmat!(geom,x,conn);# retrieve element coordinates
    gathervalues_asmat!(afield,a,conn);# retrieve element dofs
    for j = 1:npts #Loop over all integration points
      At_mul_B!(loc, Ns[j], x);# Quadrature point location
      At_mul_B!(val, Ns[j], a);# Field value at the quadrature point
      At_mul_B!(J, x, gradNparams[j]); # calculate the Jacobian matrix
      Jac = Jacobianmdim(geod, J, loc, conn,  Ns[j], m);
      result = result + fh(loc,val)*Jac*w[j];
    end
  end
  return result
end

"""
    integratefunction(self::FEMMAbstractBase,
        geom::NodalField{FFlt}, fh::F, m::FInt = -1) where {F<:Function}

Integrate a function over the discrete manifold.

Integrate some scalar function over the geometric cells. The function takes a
single argument, the position vector.

When the scalar
function returns just +1 (such as `(x) ->  1.0`), the result
measures the volume (number of points, length, area, 3-D volume,
according to the manifold dimension). When the function returns the mass
density, the method measures the mass, when the function returns the
x-coordinate equal measure the static moment with respect to the y- axis,
and so on.

### Example:
Compute the volume of the mesh and then its center of gravity:
```
V = integratefunction(femm, geom, (x) ->  1.0)
Sx = integratefunction(femm, geom, (x) ->  x[1])
Sy = integratefunction(femm, geom, (x) ->  x[2])
Sz = integratefunction(femm, geom, (x) ->  x[3])
CG = vec([Sx Sy Sz]/V)
```
Compute the moments of inertia:
```
% Now compute the moments of inertia
Ixx = =measure(femm,geom,@(x)(norm(x-CG)^2*eye(3)-(x-CG)'*(x-CG)))
mass=V*rhos;
Inertia=I*rhos;
```
"""
function integratefunction(self::FEMMAbstractBase,
  geom::NodalField{FFlt}, fh::F, m::FInt = -1) where {F<:Function}
  geod =self.geod;
  if m < 0
    m = manifdim(geod.fes);  # native  manifold dimension
  end
  # Constants
  nfes = count(geod.fes); # number of finite elements in the set
  nne = nodesperelem(geod.fes); # number of nodes per element
  sdim = ndofs(geom);            # number of space dimensions
  mdim = manifdim(geod.fes);     # manifold dimension of the element
  # Precompute basis f. values + basis f. gradients wrt parametric coor
  npts, Ns, gradNparams, w, pc = integrationdata(geod);
  conn = zeros(FInt, nne, 1); # element nodes -- used as a buffer
  x = zeros(FFlt, nne, sdim); # array of node coordinates -- used as a buffer
  loc = zeros(FFlt, 1, sdim); # quadrature point location -- used as a buffer
  J = eye(FFlt, sdim, mdim); # Jacobian matrix -- used as a buffer
  result = 0.0;# Initialize the result
  for i = 1:count(geod.fes)  # Now loop over all fes in the set
    getconn!(geod.fes, conn, i);
    gathervalues_asmat!(geom, x, conn);# retrieve element coordinates
    for j=1:npts #Loop over all integration points
      At_mul_B!(loc, Ns[j], x);# Quadrature points location
      At_mul_B!(J, x, gradNparams[j]); # calculate the Jacobian matrix
      Jac = Jacobianmdim(geod, J, loc, conn,  Ns[j], m);
      result = result + fh(vec(loc))*Jac*w[j];
    end
  end
  return result
end

"""
    distribloads(self::FEMM, assembler::A,
      geom::NodalField{FFlt},
      P::NodalField{T},
      fi::ForceIntensity,
      m::FInt) where {FEMM<:FEMMAbstractBase, T<:Number, A<:SysvecAssemblerBase}

Compute the distributed-load vector.

fi=force intensity object
m= manifold dimension, 2= surface, 3= volume
"""
function distribloads(self::FEMM, assembler::A,
  geom::NodalField{FFlt},
  P::NodalField{T},
  fi::ForceIntensity,
  m::FInt) where {FEMM<:FEMMAbstractBase, T<:Number, A<:SysvecAssemblerBase}
  geod = self.geod;         # the geometry data
  # Constants
  nfes = count(geod.fes); # number of finite elements in the set
  ndn = ndofs(P); # number of degrees of freedom per node
  nne = nodesperelem(geod.fes); # number of nodes per element
  sdim = ndofs(geom);            # number of space dimensions
  mdim = manifdim(geod.fes);     # manifold dimension of the element
  Cedim = ndn*nne;       # dimension of the element matrix/vector
  # Precompute basis f. values + basis f. gradients wrt parametric coor
  npts, Ns, gradNparams, w, pc = integrationdata(geod);
  # Prepare some buffers:
  conn = zeros(FInt,nne,1); # element nodes -- used as a buffer
  x = zeros(FFlt,nne,sdim); # array of node coordinates -- used as a buffer
  dofnums = zeros(FInt,1,Cedim); # degree of freedom array -- used as a buffer
  loc = zeros(FFlt,1,sdim); # quadrature point location -- used as a buffer
  J = eye(FFlt,sdim,mdim); # Jacobian matrix -- used as a buffer
  Fe = zeros(T,Cedim);
  startassembly!(assembler, P.nfreedofs);
  for i = 1:nfes # Loop over elements
    getconn!(geod.fes, conn, i);
    gathervalues_asmat!(geom, x, conn);# retrieve element coordinates
    fill!(Fe, 0.0);
    for j = 1:npts
      At_mul_B!(loc, Ns[j], x);# Quadrature point location
      At_mul_B!(J, x, gradNparams[j]); # calculate the Jacobian matrix
      Jac = Jacobianmdim(geod, J, loc, conn,  Ns[j], m);
      updateforce!(fi, loc, J, geod.fes.label[i]); # retrieve the applied load
      Factor::FFlt = (Jac * w[j]);
      NkxF::FFlt = 0.0
      rx::FInt=1;
      for kx = 1:nne # all the nodes
        NkxF = Ns[j][kx] * Factor
        for mx = 1:ndn   # all the degrees of freedom
          Fe[rx] = Fe[rx] + NkxF * fi.force[mx];
          rx = rx+1;    # next component of the vector
        end
      end
    end
    gatherdofnums!(P, dofnums, conn);
    assemble!(assembler, Fe, dofnums);
  end
  F = makevector!(assembler);
  return F
end

function distribloads(self::FEMM,
  geom::NodalField{FFlt},
  P::NodalField{T},
  fi::ForceIntensity,
  m::FInt) where {FEMM<:FEMMAbstractBase, T<:Number}
  assembler = SysvecAssembler(0.0*P.values[1])#T(0.0))
  return distribloads(self, assembler, geom, P, fi, m)
end

"""
    connectionmatrix(self::FEMM, nnodes::FInt) where {FEMM<:FEMMAbstractBase}

Compute the connection matrix.

The matrix has a nonzero in all the rows and columns which correspond to nodes
connected by some finite element.
"""
function connectionmatrix(self::FEMM, nnodes::FInt) where {FEMM<:FEMMAbstractBase}
  nfes = size(self.fes.conn,1)
  nconns = size(self.fes.conn,2)
  N = nfes*nconns*nconns
  rb = FInt[]; sizehint!(rb, N)
  cb = FInt[]; sizehint!(cb, N)
  vb = ones(FInt, N);
  @inbounds for  j = 1:nfes
    @inbounds for  k = 1:nconns
      append!(rb, self.fes.conn[j, :])
      @inbounds for  m = 1:nconns
        push!(cb, self.fes.conn[j, k])
      end
    end
  end
  return sparse(rb, cb, vb, nnodes, nnodes)
end


struct InverseDistanceInspectorData
  d::FFltVec
  sum_inv_dist::FFltVec
  sum_quant_inv_dist::FFltMat
end

"""
    fieldfromintegpoints(self::FEMM,
      geom::NodalField{FFlt},  u::NodalField{T},
      dT::NodalField{FFlt},  quantity::Symbol,  component::FInt;
      context...) where {FEMM<:FEMMAbstractBase, T<:Number}

Construct nodal field from integration points.

Input arguments
`geom`     - reference geometry field
`u`        - displacement field
`dT`       - temperature difference field
`quantity`   - this is what you would assign to the 'quantity' argument
           of the material update!() method.
`component`- component of the 'quantity' array: see the material update()
           method.
Output argument
 - the new field that can be used to map values to colors and so on
"""
function fieldfromintegpoints(self::FEMM,
  geom::NodalField{FFlt},  u::NodalField{T},
  dT::NodalField{FFlt},  quantity::Symbol,  component::FInt;
  context...) where {FEMM<:FEMMAbstractBase, T<:Number}
  geod = self.geod
  # Constants
  nne = nodesperelem(geod.fes); # number of nodes for element
  sdim = ndofs(geom);            # number of space dimensions
  # Container of intermediate results
  idat = InverseDistanceInspectorData(
  zeros(FFlt, nne),
  zeros(FFlt, nnodes(geom)),
  zeros(FFlt, nnodes(geom), length(component))
  );
  # This is an inverse-distance interpolation inspector.
  function idi_inspector(idat, elnum, conn, xe,  out,  xq)
    # xe = coordinates of the nodes of the element
    # xq = coordinate of the quadrature point
    mindn = Inf
    for jjj = 1:length(idat.d)
      idat.d[jjj]  = sum((vec(xe[jjj, :])-vec(xq)).^2);
      if (idat.d[jjj] > 0.0)
        mindn = min(mindn, idat.d[jjj])
      end
    end
    mindn = mindn/1.0e9;
    for jjj = 1:length(idat.d)
      invdjjj  = 1.0/(idat.d[jjj]+mindn);
      quant =  out[component]
      for kkk = 1:length(quant)
        idat.sum_quant_inv_dist[conn[jjj], kkk] += invdjjj*quant[kkk];
      end
      idat.sum_inv_dist[conn[jjj]] += invdjjj;
    end
    return idat
  end
  # Loop over cells to interpolate to nodes
  idat = inspectintegpoints(self,  geom,  u,  dT,
            collect(1:count(geod.fes)),  idi_inspector,  idat, quantity;
            context...);
  # The data for the field to be constructed is initialized
  nvals = zeros(FFlt, nnodes(geom), length(component));
  # compute the data array
  for kkk = 1:size(nvals, 1)
    for j = 1:length(component)
      if (idat.sum_inv_dist[kkk] > 0.0)
        nvals[kkk, j] = idat.sum_quant_inv_dist[kkk, j]/idat.sum_inv_dist[kkk];
      end
    end
  end
  # Make the field
  return NodalField(nvals);
end

function fieldfromintegpoints(self::FEMM,
  geom::NodalField{FFlt},  u::NodalField{T},
  quantity::Symbol,  component::FInt;
  context...) where {FEMM<:FEMMAbstractBase, T<:Number}
  dT = NodalField(zeros(FFlt, nnodes(geom), 1)) # zero difference in temperature
  return fieldfromintegpoints(self, geom, u, dT, quantity, component; context...)
end


struct MeanValueInspectorData
  n_quant::FIntVec
  sum_quant_value::FFltMat
end

"""
    elemfieldfromintegpoints(self::FEMM,
      geom::NodalField{FFlt},  u::NodalField{T},
      dT::NodalField{FFlt},  quantity::Symbol,  component::FInt;
      context...) where {FEMM<:FEMMAbstractBase, T<:Number}

Construct nodal field from integration points.

Input arguments
`geom`     - reference geometry field
`u`        - displacement field
`dT`       - temperature difference field
`quantity`   - this is what you would assign to the 'quantity' argument
           of the material update!() method.
`component`- component of the 'quantity' array: see the material update()
           method.
Output argument
 - the new field that can be used to map values to colors and so on
"""
function elemfieldfromintegpoints(self::FEMM,
  geom::NodalField{FFlt},  u::NodalField{T},
  dT::NodalField{FFlt},  quantity::Symbol,  component::FInt;
  context...) where {FEMM<:FEMMAbstractBase, T<:Number}
  geod = self.geod
  # Constants
  nne = nodesperelem(geod.fes); # number of nodes for element
  sdim = ndofs(geom);            # number of space dimensions
  # Container of intermediate results
  idat = MeanValueInspectorData(
      zeros(FInt, count(geod.fes)),
      zeros(FFlt, count(geod.fes), length(component))
  );
  # This is an mean-value interpolation inspector. The mean of the
  # quadrature-point quantities is reported per element.
  function mv_inspector(idat, elnum, conn, xe,  out,  xq)
    # xe = coordinates of the nodes of the element
    # xq = coordinate of the quadrature point
    idat.n_quant[elnum] += 1
    quant = out[component]
    for kkk = 1:length(quant)
      idat.sum_quant_value[elnum, kkk] += quant[kkk];
    end
    return idat
  end
  # Loop over cells to interpolate to nodes
  idat = inspectintegpoints(self,  geom,  u,  dT,
            collect(1:count(geod.fes)), mv_inspector,  idat, quantity;
            context...);
  # The data for the field to be constructed is initialized
  evals = zeros(FFlt, count(geod.fes), length(component));
  # compute the data array
  for j = 1:size(evals, 1)
    for kkk = 1:size(evals, 2)
      evals[j, kkk] = idat.sum_quant_value[j, kkk]/idat.n_quant[j];
    end
  end
  # Make the field
  return ElementalField(evals);
end

function elemfieldfromintegpoints(self::FEMM,
  geom::NodalField{FFlt},  u::NodalField{T},
  quantity::Symbol,  component::FInt; context...) where {FEMM<:FEMMAbstractBase,
  T<:Number}
  dT = NodalField(zeros(FFlt, nnodes(geom), 1)) # zero difference in temperature
  return elemfieldfromintegpoints(self, geom, u, dT, quantity, component; context...)
end

end
