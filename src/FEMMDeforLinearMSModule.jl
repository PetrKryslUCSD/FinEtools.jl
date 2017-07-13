"""
    FEMMDeforLinearMSModule

Module for operations on interiors of domains to construct system matrices and
system vectors for linear deformation models:  mean-strain  formulation.
"""
module FEMMDeforLinearMSModule

export FEMMDeforLinearMSH8, FEMMDeforLinearMST10
export stiffness, nzebcloadsstiffness, thermalstrainloads,
       inspectintegpoints

using FinEtools.FTypesModule
using FinEtools.FESetModule
using FinEtools.FESetModule.gradN!
using FinEtools.CSysModule
using FinEtools.GeoDModule
using FinEtools.FEMMBaseModule
using FinEtools.FEMMDeforLinearBaseModule
using FinEtools.FieldModule
using FinEtools.NodalFieldModule
using FinEtools.ElementalFieldModule
using FinEtools.ForceIntensityModule
using FinEtools.AssemblyModule
using FinEtools.DeforModelRedModule
using FinEtools.MatDeforModule
using FinEtools.MatDeforElastIsoModule
using FinEtools.FENodeToFEMapModule
using FinEtools.MatrixUtilityModule.add_btdb_ut_only!
using FinEtools.MatrixUtilityModule.complete_lt!
using FinEtools.MatrixUtilityModule.mv_product!
using FinEtools.MatrixUtilityModule.add_btv!
import FinEtools.FEMMDeforLinearBaseModule.stiffness
import FinEtools.FEMMDeforLinearBaseModule.nzebcloadsstiffness
import FinEtools.FEMMDeforLinearBaseModule.mass
import FinEtools.FEMMDeforLinearBaseModule.thermalstrainloads
import FinEtools.FEMMDeforLinearBaseModule.inspectintegpoints
import FinEtools.FEMMBaseModule.associategeometry!

abstract type FEMMDeforLinearAbstractMS <: FEMMDeforLinearAbstract end

mutable struct FEMMDeforLinearMSH8{MR<:DeforModelRed,
  S<:FESetH8, F<:Function, M<:MatDefor} <: FEMMDeforLinearAbstractMS
  mr::Type{MR}
  geod::GeoD{S, F} # geometry data finite element modeling machine
  material::M # material object
  stabilization_material::MatDeforElastIso
  phis::FFltVec
  function FEMMDeforLinearMSH8(mr::Type{MR}, geod::GeoD{S, F},
    material::M) where {MR<:DeforModelRed,
    S<:FESetH8, F<:Function, M<:MatDefor}
    @assert mr === material.mr "Model reduction is mismatched"
    @assert (mr === DeforModelRed3D) "3D model required"
    stabilization_material = make_stabilization_material(material)
    return new{MR, S, F, M}(mr, geod,  material, stabilization_material,
      zeros(FFlt, 1))
  end
end

function make_stabilization_material(material::M) where {M}
  ns = fieldnames(typeof(material))
  E = 0.0; nu = 0.0
  if :E in ns
    E = material.E
    if material.nu < 0.3
      nu = material.nu
    else
      nu = 0.3 + (material.nu - 0.3) / 2.0
    end
  else
    if :E1 in ns
      E = min(material.E1, material.E2, material.E3)
      nu = min(material.nu12, material.nu13, material.nu23)
    else
      error("No clues on how to construct the stabilization material")
    end
  end
  return  MatDeforElastIso(material.mr, 0.0, E, nu, 0.0)
end

function buffers1(self::FEMMDeforLinearAbstractMS, geom::NodalField, npts::FInt)
  const nne = nodesperelem(self.geod.fes); # number of nodes for element
  const sdim = ndofs(geom);            # number of space dimensions
  const mdim = manifdim(self.geod.fes); # manifold dimension of the element
    # Prepare buffers
  conn = zeros(FInt, nne, 1); # element nodes -- buffer
  x = zeros(FFlt, nne, sdim); # array of node coordinates -- buffer
  loc = zeros(FFlt, 1, sdim); # quadrature point location -- buffer
  J = eye(FFlt, sdim, mdim); # Jacobian matrix -- buffer
  csmatTJ = zeros(FFlt, mdim, mdim); # intermediate result -- buffer
  gradN = zeros(FFlt, nne, mdim);
  return conn, x, loc, J, csmatTJ, gradN
end

function buffers2(self::FEMMDeforLinearAbstractMS, geom::NodalField, u::NodalField, npts::FInt)
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
  AllgradN = Array{FFltMat}(npts);
  for ixxxx = 1:npts
    AllgradN[ixxxx] = zeros(FFlt, nne, mdim);
  end
  Jac = zeros(FFlt, npts);
  MeangradN = zeros(FFlt, nne, mdim); # intermediate result -- buffer
  D = zeros(FFlt, nstrs, nstrs); # material stiffness matrix -- buffer
  Dstab = zeros(FFlt, nstrs, nstrs); # material stiffness matrix -- buffer
  B = zeros(FFlt, nstrs, elmatdim); # strain-displacement matrix -- buffer
  DB = zeros(FFlt, nstrs, elmatdim); # strain-displacement matrix -- buffer
  Bbar = zeros(FFlt, nstrs, elmatdim); # strain-displacement matrix -- buffer
  elvecfix = zeros(FFlt, elmatdim, 1); # vector of prescribed displ. -- buffer
  elvec = zeros(FFlt, elmatdim); # element vector -- buffer
  return conn, x, dofnums, loc, J, csmatTJ, AllgradN, MeangradN, Jac,
            D, Dstab, B, DB, Bbar, elmat, elvec, elvecfix
end

"""
    associategeometry!(self::FEMMAbstractBase,  geom::NodalField{FFlt})

Associate geometry field with the FEMM.

Compute the  correction factors to account for  the shape of the  elements.
"""
function associategeometry!(self::F,  geom::NodalField{FFlt}) where {F<:FEMMDeforLinearMSH8}
  geod = self.geod
  npts,  Ns,  gradNparams,  w,  pc = integrationdata(geod);
  conn, x, loc, J, csmatTJ, gradN = buffers1(self, geom, npts)
  self.phis = zeros(FFlt, count(geod.fes))
  for i = 1:count(geod.fes) # Loop over elements
    getconn!(geod.fes, conn, i);
    gathervalues_asmat!(geom, x, conn);# retrieve element coordinates
    for j = 1:npts # Loop over quadrature points
      At_mul_B!(J, x, gradNparams[j]); # calculate the Jacobian matrix
      At_mul_B!(csmatTJ, geod.mcsys.csmat, J); # local Jacobian matrix
      gradN!(geod.fes, gradN, gradNparams[j], csmatTJ);
      h2 = diag(transpose(csmatTJ)*csmatTJ)
      cap_phi = (2 * (1 + self.stabilization_material.nu) * (minimum(h2) / maximum(h2)))  # Plane stress
      phi = cap_phi / (1 + cap_phi)
      self.phis[i] = max(self.phis[i], phi)
    end # Loop over quadrature points
  end # Loop over elements
  return self
end

"""
stiffness(self::FEMMDeforLinearAbstractMS, assembler::A,
      geom::NodalField{FFlt},
      u::NodalField{T}) where {A<:SysmatAssemblerBase, T<:Number}

Compute and assemble  stiffness matrix.
"""
function stiffness(self::FEMMDeforLinearAbstractMS, assembler::A,
      geom::NodalField{FFlt},
      u::NodalField{T}) where {A<:SysmatAssemblerBase, T<:Number}
  geod = self.geod
  npts,  Ns,  gradNparams,  w,  pc = integrationdata(geod);
  conn, x, dofnums, loc, J, csmatTJ, AllgradN, MeangradN, Jac,
    D, Dstab, B, DB, Bbar, elmat, elvec, elvecfix = buffers2(self, geom, u, npts)
  self.material.tangentmoduli!(self.material, D, 0.0, 0.0, loc, 0)
  self.stabilization_material.tangentmoduli!(self.stabilization_material,
    Dstab, 0.0, 0.0, loc, 0)
  startassembly!(assembler, size(elmat, 1), size(elmat, 2), count(geod.fes),
    u.nfreedofs, u.nfreedofs);
  for i = 1:count(geod.fes) # Loop over elements
    getconn!(geod.fes, conn, i);
    gathervalues_asmat!(geom, x, conn);# retrieve element coordinates
    # NOTE: the coordinate system should be evaluated at a single point within the
    # element in order for the derivatives to be consistent at all quadrature points
    loc .= mean(x, 1)
    updatecsmat!(geod.mcsys, loc, J, geod.fes.label[i]);
    vol = 0.0; # volume of the element
    fill!(MeangradN, 0.0) # mean basis function gradients
    for j = 1:npts # Loop over quadrature points
      At_mul_B!(J, x, gradNparams[j]); # calculate the Jacobian matrix
      Jac[j] = Jacobianvolume(geod, J, loc, conn, Ns[j]);
      At_mul_B!(csmatTJ, geod.mcsys.csmat, J); # local Jacobian matrix
      gradN!(geod.fes, AllgradN[j], gradNparams[j], csmatTJ);
      dvol = Jac[j]*w[j]
      MeangradN .= MeangradN .+ AllgradN[j]*dvol
      vol = vol + dvol
    end # Loop over quadrature points
    MeangradN .= MeangradN/vol
    Blmat!(self.mr, Bbar, Ns[1], MeangradN, loc, geod.mcsys.csmat);
    fill!(elmat,  0.0); # Initialize element matrix
    add_btdb_ut_only!(elmat, Bbar, vol, D, DB)
    add_btdb_ut_only!(elmat, Bbar, -self.phis[i]*vol, Dstab, DB)
    for j = 1:npts # Loop over quadrature points
      Blmat!(self.mr, B, Ns[j], AllgradN[j], loc, geod.mcsys.csmat);
      add_btdb_ut_only!(elmat, B, self.phis[i]*Jac[j]*w[j], Dstab, DB)
    end # Loop over quadrature points
    complete_lt!(elmat)
    gatherdofnums!(u, dofnums, conn); # retrieve degrees of freedom
    assemble!(assembler, elmat, dofnums, dofnums); # assemble symmetric matrix
  end # Loop over elements
  return makematrix!(assembler);
end

# function stiffness(self::FEMMDeforLinearAbstractMS,
#             geom::NodalField{FFlt},  u::NodalField{T}) where {T<:Number}
#   assembler = SysmatAssemblerSparseSymm();
#   return stiffness(self, assembler, geom, u);
# end

# """
#     nzebcloadsstiffness(self::FEMMDeforLinear,  assembler::A,
#       geom::NodalField{FFlt},
#       u::NodalField{T}) where {A<:SysvecAssemblerBase, T<:Number}
#
# Compute load vector for nonzero EBC for fixed displacement.
# """
# function nzebcloadsstiffness(self::FEMMDeforLinear,  assembler::A,
#   geom::NodalField{FFlt},
#   u::NodalField{T}) where {A<:SysvecAssemblerBase, T<:Number}
#   geod = self.geod
#   npts,  Ns,  gradNparams,  w,  pc = integrationdata(geod);
#   conn, x, dofnums, loc, J, csmatTJ, gradN, D, B, DB, elmat, elvec, elvecfix =
#                 buffers(self, geom, u)  # Prepare buffers
#   self.material.tangentmoduli!(self.material, D, 0.0, 0.0, loc, 0)
#   startassembly!(assembler,  u.nfreedofs);
#   for i = 1:count(geod.fes) # Loop over elements
#     getconn!(geod.fes, conn, i);
#     gathervalues_asvec!(u, elvecfix, conn);# retrieve element displacement vector
#     if norm(elvecfix) != 0     # Is the load nonzero?
#       gathervalues_asmat!(geom, x, conn);# retrieve element coordinates
#       fill!(elmat,  0.0); # Initialize element matrix
#       for j = 1:npts # Loop over quadrature points
#         At_mul_B!(loc, Ns[j], x);# Quadrature points location
#         At_mul_B!(J, x, gradNparams[j]); # calculate the Jacobian matrix
#         Jac = Jacobianvolume(geod, J, loc, conn, Ns[j]);
#         updatecsmat!(geod.mcsys, loc, J, geod.fes.label[i]);
#         At_mul_B!(csmatTJ, geod.mcsys.csmat, J); # local Jacobian matrix
#         gradN!(geod.fes, gradN, gradNparams[j], csmatTJ);
#         Blmat!(self.mr, B, Ns[j], gradN, loc, geod.mcsys.csmat);
#         add_btdb_ut_only!(elmat, B, Jac*w[j], D, DB)
#       end # Loop over quadrature points
#       complete_lt!(elmat)
#       gatherdofnums!(u, dofnums, conn); # retrieve degrees of freedom
#       mv_product!(elvec, elmat, elvecfix)
#       assemble!(assembler,  -elvec,  dofnums); # assemble element load vector
#     end
#   end # Loop over elements
#   return makevector!(assembler);
# end
#
# function nzebcloadsstiffness(self::FEMMDeforLinear,
#               geom::NodalField{FFlt},
#               u::NodalField{T}) where {T<:Number}
#     assembler = SysvecAssembler()
#     return  nzebcloadsstiffness(self, assembler, geom, u);
# end
#

# """
#     thermalstrainloads(self::FEMMDeforLinear, assembler::A,
#         geom::NodalField{FFlt}, u::NodalField{T},
#         dT::NodalField{FFlt}) where {A<:SysvecAssemblerBase, T<:Number}
#
# Compute the thermal-strain load vector.
# """
# function  thermalstrainloads(self::FEMMDeforLinear, assembler::A,
#     geom::NodalField{FFlt}, u::NodalField{T},
#     dT::NodalField{FFlt}) where {A<:SysvecAssemblerBase, T<:Number}
#   geod = self.geod
#   npts,  Ns,  gradNparams,  w,  pc = integrationdata(geod);
#   conn, x, dofnums, loc, J, csmatTJ, gradN, D, B, DB, elmat, elvec, elvecfix =
#                 buffers(self, geom, u)# Prepare buffers
#   t= 0.0
#   dt = 0.0
#   DeltaT = zeros(FFlt, length(conn))
#   strain = zeros(FFlt, nstsstn(self.mr)); # total strain -- buffer
#   thstrain = zeros(FFlt, nthstn(self.mr)); # thermal strain -- buffer
#   thstress = zeros(FFlt, nstsstn(self.mr)); # thermal stress -- buffer
#   startassembly!(assembler,  u.nfreedofs);
#   for i = 1:count(geod.fes) # Loop over elements
#     getconn!(geod.fes, conn, i);
#     gathervalues_asvec!(dT, DeltaT, conn);# retrieve element temperatures
#     if norm(DeltaT) != 0     # Is the thermal increment nonzero?
#       gathervalues_asmat!(geom, x, conn);# retrieve element coordinates
#       fill!(elvec,  0.0); # Initialize element matrix
#       for j = 1:npts # Loop over quadrature points
#         At_mul_B!(loc, Ns[j], x);# Quadrature points location
#         At_mul_B!(J, x, gradNparams[j]); # calculate the Jacobian matrix
#         Jac = Jacobianvolume(geod, J, loc, conn, Ns[j]);
#         updatecsmat!(geod.mcsys, loc, J, geod.fes.label[i]);
#         At_mul_B!(csmatTJ,  geod.mcsys.csmat,  J); # local Jacobian matrix
#         gradN!(geod.fes, gradN, gradNparams[j], csmatTJ);#Do: gradN = gradNparams[j]/csmatTJ;
#         Blmat!(self.mr, B, Ns[j], gradN, loc, geod.mcsys.csmat);# strains in mater cs, displ in global cs
#         self.material.thermalstrain!(self.material, thstrain,
#           dot(vec(Ns[j]), DeltaT))
#         thstress = self.material.update!(self.material, thstress,
#           thstress, strain, thstrain, t, dt, loc, geod.fes.label[i], :nothing)
#         add_btv!(elvec, B, thstress, (-1)*(Jac*w[j]))
#       end
#       gatherdofnums!(u, dofnums, conn); # retrieve degrees of freedom
#       assemble!(assembler,  elvec,  dofnums); # assemble element load vector
#     end
#   end # Loop over elements
#   return makevector!(assembler);
# end
#
# function thermalstrainloads(self::FEMMDeforLinear,
#               geom::NodalField{FFlt},
#               u::NodalField{T}, dT::NodalField{FFlt}) where {T<:Number}
#     assembler = SysvecAssembler()
#     return  thermalstrainloads(self, assembler, geom, u, dT);
# end
#
"""
    inspectintegpoints(self::FEMMDeforLinear,
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
function inspectintegpoints(self::FEMMDeforLinearAbstractMS,
  geom::NodalField{FFlt},  u::NodalField{T},
  dT::NodalField{FFlt},
  felist::FIntVec,
  inspector::F,  idat, quantity=:Cauchy;
  context...) where {T<:Number, F<:Function}
  geod = self.geod
  npts,  Ns,  gradNparams,  w,  pc = integrationdata(geod);
  conn, x, dofnums, loc, J, csmatTJ, AllgradN, MeangradN, Jac,
  D, Dstab, B, DB, Bbar, elmat, elvec, elvecfix = buffers2(self, geom, u, npts)
  MeanN = deepcopy(Ns[1])
  # Sort out  the output requirements
  outputcsys = geod.mcsys; # default: report the stresses in the material coord system
  for arg in context
    sy,  val = arg
    if sy == :outputcsys
      outputcsys = val
    end
  end
  t= 0.0
  dt = 0.0
  dTe = zeros(FFlt, length(conn), 1) # nodal temperatures -- buffer
  ue = zeros(FFlt, size(elmat, 1), 1); # array of node displacements -- buffer
  qpdT = 0.0; # node temperature increment
  qpstrain = zeros(FFlt, nstsstn(self.mr), 1); # total strain -- buffer
  qpthstrain = zeros(FFlt, nthstn(self.mr)); # thermal strain -- buffer
  qpstress = zeros(FFlt, nstsstn(self.mr)); # stress -- buffer
  out1 = zeros(FFlt, nstsstn(self.mr)); # stress -- buffer
  out =  zeros(FFlt, nstsstn(self.mr));# output -- buffer
  # Loop over  all the elements and all the quadrature points within them
  for ilist = 1:length(felist) # Loop over elements
    i = felist[ilist];
    getconn!(geod.fes, conn, i);
    gathervalues_asmat!(geom, x, conn);# retrieve element coordinates
    gathervalues_asvec!(u, ue, conn);# retrieve element displacements
    gathervalues_asvec!(dT, dTe, conn);# retrieve element temperature increments
    # NOTE: the coordinate system should be evaluated at a single point within the
    # element in order for the derivatives to be consistent at all quadrature points
    loc .= mean(x, 1)
    updatecsmat!(geod.mcsys, loc, J, geod.fes.label[i]);
    vol = 0.0; # volume of the element
    fill!(MeangradN, 0.0) # mean basis function gradients
    fill!(MeanN, 0.0) # mean basis function gradients
    for j = 1:npts # Loop over quadrature points
      At_mul_B!(J, x, gradNparams[j]); # calculate the Jacobian matrix
      Jac[j] = Jacobianvolume(geod, J, loc, conn, Ns[j]);
      At_mul_B!(csmatTJ, geod.mcsys.csmat, J); # local Jacobian matrix
      gradN!(geod.fes, AllgradN[j], gradNparams[j], csmatTJ);
      dvol = Jac[j]*w[j]
      MeangradN .= MeangradN .+ AllgradN[j]*dvol
      MeanN .= MeanN .+ Ns[j]*dvol
      vol = vol + dvol
    end # Loop over quadrature points
    MeangradN .= MeangradN/vol
    Blmat!(self.mr, Bbar, MeanN, MeangradN, loc, geod.mcsys.csmat);
    MeanN .= MeanN/vol
    qpdT = dot(vec(dTe), vec(MeanN));# Quadrature point temperature increment
    # Quadrature point quantities
    A_mul_B!(qpstrain, Bbar, ue); # strain in material coordinates
    self.material.thermalstrain!(self.material, qpthstrain, qpdT)
    # Material updates the state and returns the output
    out = self.material.update!(self.material, qpstress, out,
    vec(qpstrain), qpthstrain, t, dt, loc, geod.fes.label[i], quantity)
    if (quantity == :Cauchy)   # Transform stress tensor,  if that is "out"
      (length(out1) >= length(out)) || (out1 = zeros(length(out)))
      rotstressvec(self.mr, out1, out, geod.mcsys.csmat')# To global coord sys
      rotstressvec(self.mr, out, out1, outputcsys.csmat)# To output coord sys
    end
    # Call the inspector
    idat = inspector(idat, i, conn, x, out, loc);
  end # Loop over elements
  return idat; # return the updated inspector data
end

function inspectintegpoints(self::FEMMDeforLinearAbstractMS,
  geom::NodalField{FFlt},  u::NodalField{T},
  felist::FIntVec,
  inspector::F, idat, quantity=:Cauchy;
  context...) where {T<:Number, F<:Function}
  dT = NodalField(zeros(FFlt, nnodes(geom), 1)) # zero difference in temperature
  return inspectintegpoints(self, geom, u, dT, felist,
            inspector, idat, quantity; context...);
end


end
