module MatDeforElastOrthoModule

using FinEtools.FTypesModule
using FinEtools.DeforModelRedModule
using FinEtools.MatDeforModule

"""
    MatDeforElastOrtho

Linear orthotropic elasticity  material.


```tangentmoduli!(self::MatDeforElastOrtho,
  ms::MatDeforElastOrthoMS, D::FFltMat,
  t::FFlt, dt::FFlt, loc::FFltMat, label::FInt)
```


"""
struct  MatDeforElastOrtho{MR<:DeforModelRed,
  MTAN<:Function, MUPD<:Function, MTHS<:Function} <: MatDefor
  mr::Type{MR}
  mass_density::FFlt # mass density
  E1::FFlt # Young's modulus for material direction 1
  E2::FFlt # Young's modulus for material direction 2
  E3::FFlt # Young's modulus for material direction 3
  nu12::FFlt
  nu13::FFlt
  nu23::FFlt
  G12::FFlt
  G13::FFlt
  G23::FFlt
  CTE1::FFlt # three thermal expansion coefficients
  CTE2::FFlt # three thermal expansion coefficients
  CTE3::FFlt # three thermal expansion coefficients
  D::FFltMat # cached matrix of tangent moduli
  tangentmoduli!::MTAN
  update!::MUPD
  thermalstrain!::MTHS
end
export MatDeforElastOrtho

function MatDeforElastOrtho(mr::Type{DeforModelRed3D},
  mass_density::FFlt, E1::FFlt, E2::FFlt, E3::FFlt,
  nu12::FFlt, nu13::FFlt, nu23::FFlt,
  G12::FFlt, G13::FFlt, G23::FFlt,
  CTE1::FFlt, CTE2::FFlt, CTE3::FFlt)
  C = [  1.0/E1      -nu12/E1    -nu13/E1  0.0   0.0   0.0;
       -nu12/E1      1.0/E2      -nu23/E2  0.0   0.0   0.0;
       -nu13/E1    -nu23/E2       1.0/E3   0.0   0.0   0.0;
       0.0           0.0           0.0   1/G12   0.0   0.0;
       0.0           0.0           0.0     0.0 1/G13   0.0;
       0.0           0.0           0.0     0.0   0.0  1/G23];
  D = inv(C)
  if (rank(D)<6)
    error("Non-positive definite D!");
  end
  ev = eigvals(D)
  for  e in ev
    #println("$e")
    if (e < 0.0)
      error("Indefinite D!");
    end
  end
  return MatDeforElastOrtho(mr, mass_density, E1, E2, E3,
  nu12, nu13, nu23, G12, G13, G23, CTE1, CTE2, CTE3, D,
  tangentmoduli3d!, update3d!, thermalstrain3d!)
end

################################################################################
# 3-D solid model
################################################################################

"""
    tangentmoduli3d!(self::MatDeforElastOrtho,
      ms::MatDeforElastOrthoMS, D::FFltMat,
      t::FFlt, dt::FFlt, loc::FFltMat, label::FInt)

Calculate the material stiffness matrix.

`D` = matrix of tangent moduli, 6 x 6, supplied as a buffer and overwritten.
"""
function tangentmoduli3d!(self::MatDeforElastOrtho,
  D::FFltMat,
  t::FFlt, dt::FFlt, loc::FFltMat, label::FInt)
  copy!(D,self.D);
  return D
end

"""
    update3d!(self::MatDeforElastOrtho,
      output::FFltMat, stress::FFltVec,
      strain::FFltVec, thstrain::FFltVec=zeros(6,1), t::FFlt= 0.0, dt::FFlt= 0.0,
      loc::FFltMat=zeros(3,1), label::FInt=0, quantity=:nothing)

Update material state.

`strain` = strain vector,
`thstrain` = thermal strain vector,
`t` = current time,
`dt` = current time step,
`loc` = location of the quadrature point in global Cartesian coordinates,
`label` = label of the finite element in which the quadrature point is found.

These quantities get updated or defined:

`stress` = Cauchy stress, defined upon return
`output` = array for outputs, needs to be pre-allocated, defined upon return
"""
function update3d!(self::MatDeforElastOrtho,
  output::FFltMat, stress::FFltVec,
  strain::FFltVec, thstrain::FFltVec=zeros(6,1), t::FFlt= 0.0, dt::FFlt= 0.0,
  loc::FFltMat=zeros(3,1), label::FInt=0, quantity=:nothing)
  A_mul_B!(stress, self.D, strain-thstrain);
  if quantity == :nothing
    #Nothing to be copied to the output array
  elseif quantity == :pressure || quantity == :Pressure
    if length(output) < 1
      output = zeros(1,1)
    end
    copy!(output,  [-sum(stress[1:3])/3.]);
  elseif quantity == :princCauchy || quantity == :princcauchy
    t=zeros(FFlt,3,3)
    t = stress6vto3x3t!(stress,t);
    ep=eig(t);
    if length(output) < 3
      output = zeros(3,1)
    end
    copy!(output,  sort(ep[1]));
  elseif quantity==:vonMises || quantity==:vonmises || quantity==:von_mises || quantity==:vm
    s1=stress[1]; s2=stress[2]; s3=stress[3];
    s4=stress[4]; s5=stress[5]; s6=stress[6];
    if length(output) < 1
      output = zeros(1,1)
    end
    copy!(output,  [sqrt(1.0/2*((s1-s2)^2+(s1-s3)^2+(s2-s3)^2+6*(s4^2+s5^2+s6^2)))])
  end
  return output
end

"""
    thermalstrain3d!(self::MatDeforElastOrtho, thstrain::FFltMat, dT= 0.0)

Compute thermal strain from the supplied temperature increment.
"""
function thermalstrain3d!(self::MatDeforElastOrtho, thstrain::FFltVec, dT= 0.0)
  thstrain[1] = self.CTE1*dT
  thstrain[2] = self.CTE2*dT
  thstrain[3] = self.CTE3*dT
  thstrain[4] = 0.0
  thstrain[5] = 0.0
  thstrain[6] = 0.0
  return thstrain
end

end

#
# ################################################################################
# # Plane strain model
#
# function tangentmoduli!{P<:PropertyDeformationLinear}(::Type{DeforModelRed2DStrain},
#                         self::MaterialDeformationLinear{P},
#                         D::FFltMat; context...)
#     # # Calculate the material stiffness matrix.
#     # #
#     # # Arguments
#     # #     m=material
#     # #     context=structure with mandatory and optional fields that are required
#     # # by the specific material implementation.
#     # #
#     # # the output arguments are
#     # #     D=matrix 6x6 in the local material orientation Cartesian basis
#     # #
#
#     D3d=zeros(FFlt,6,6)
#     tangentmoduli3d!(self.property, D3d; context...);
#     const ix=[1, 2, 4];
#     for i=1:length(ix)
#         for j=1:length(ix)
#             D[j,i]=D3d[ix[j],ix[i]];
#         end
#     end
#     return D
# end
# export tangentmoduli!
#
#
# function update!{MR<:DeforModelRed2DStrain}(::Type{MR},
#                         self::MaterialDeformationLinear, ms; context...)
#         # Update material state.
#         #
#         # function [out, newms] = update (self, ms, context)
#         #
#         # Update material state.  Return the updated material state, and the
#         # requested quantity (default is the stress).
#         #   Call as:
#         #     [out,newms] = update(m, ms, context)
#         #  where
#         #     m=material
#         #     ms = material state
#         #     context=structure
#         #        with mandatory fields
#         #           strain=strain vector  in the local material
#         #               directions (which may be the same as the global coordinate
#         #               directions)
#         #        and optional fields
#         #           output=type of quantity to output, and interpreted by the
#         #               particular material; [] is returned when the material does not
#         #               recognize the requested quantity to indicate uninitialized
#         #               value.  It can be tested with isempty ().
#         #                  output ='Cauchy' - Cauchy stress; this is the default
#         #                      when output type is not specified.
#         #                  output ='princCauchy' - principal Cauchy stress;
#         #                  output ='pressure' - pressure;
#         #                  output ='vol_strain' - volumetric strain;
#         #           outputRm=optional orientation matrix in which output should
#         #               supplied
#         #
#         #   It is assumed that stress is output in m-component vector
#         #           form. m=3 for plane stress, m=4 for plane strain or axially
#         #           symmetric
#         #   The output arguments are
#         #     out=requested quantity
#         #     newms=new material state; don't forget that if the update is final
#         #           the material state newms must be assigned and stored.  Otherwise
#         #           the material update is lost!
#     #
#     Ev=FFlt[]                  # it is empty: we need to get it from context
#     output=:Cauchy
#     dT= 0.0
#     for arg in context
#         sy, val = arg
#         if sy==:strain
#             Ev=val
#         elseif sy==:output
#             output=val
#         elseif sy==:dT
#             dT=val
#         end
#     end
#     D=zeros(FFlt,3,3)
#     tangentmoduli!(MR,self,D; context...)
#     tSigma = thermalstress(MR,self; context...);
#     stress = D * Ev + tSigma;
#     # Need to output the along-the-thickness  stress
#     CTE=self.property.CTE
#     D3d=zeros(FFlt,6,6)
#     tangentmoduli3d!(self.property, D3d; context...);
#     sz=D3d[3,1:2]*Ev[1:2]-dT[1]*D3d[3,1:2]*CTE[1:2]-dT[1]*D3d[3,3]*CTE[3];
#     # sigmax, sigmay, tauxy, sigmaz
#     stress = [vec(stress[1:3]), vec(sz)];
#     if output==:Cauchy
#         out = stress;
#     elseif output==:pressure
#         out = -sum(stress[[1,2,4]])/3;
#     elseif output==:volstrain
#         out = sum(Ev[1:2]);
#     elseif output==:princCauchy
#         t=zeros(FFlt,3,3)
#         t = stress4vto3x3t!(stress,t);
#         ep=eig(t);
#         out =sort(ep[1]);
#     else
#         out = stress;
#     end
#     newms = ms;
#     return out,newms;
# end
#
# function thermalstress{MR<:DeforModelRed2DStrain}(::Type{MR},
#                         self::MaterialDeformationLinear; context...)
# # Calculate vector of thermal stress components.
# #
# # function v = thermal_stress(self,context)
# #
# #   Call as:
# #     v = thermal_stress(m,context)
# #  where
# #     m=material
# #     context=structure; see the update() method
#     #
#     for arg in context
#         sy, val = arg
#         if sy==:dT
#             dT=val
#             D=zeros(FFlt,3,3)
#             tangentmoduli!(MR, self, D; context...);# local material stiffness
#             # there is no sigmaz in this vector!
#             CTE=self.property.CTE
#             v = -dT[1]*[D[1:2,1:2]*CTE[1:2]+CTE[3]*D[1:2,3]; 0];
#             return v
#         end
#     end
#     return  zeros(FFlt, 3, 1);
# end
# export thermalstress
#
# ################################################################################
# # Plane stress model
#
# function tangentmoduli!{P<:PropertyDeformationLinear,
#     MR<:DeforModelRed2DStress}(::Type{MR},
#                         self::MaterialDeformationLinear{P},
#                         D::FFltMat; context...)
#     # # Calculate the material stiffness matrix.
#     # #
#     # # Arguments
#     # #     m=material
#     # #     context=structure with mandatory and optional fields that are required
#     # # by the specific material implementation.
#     # #
#     # # the output arguments are
#     # #     D=matrix 6x6 in the local material orientation Cartesian basis
#     # #
#
#     D3d=zeros(FFlt,6,6)
#     tangentmoduli3d!(self.property, D3d; context...);
#     Dt = D3d[1:2, 1:2]- reshape(D3d[1:2,3], (2,1))*reshape(D3d[3,1:2], (1,2))/D3d[3,3];
#     const ix=[1, 2, 4];
#     for i=1:2
#         for j=1:2
#             D[j,i]= Dt[j,i];
#         end
#     end
#     for i=1:3
#         D[3,i]= D[i,3]= D3d[4,ix[i]];
#     end
#     return D
# end
# export tangentmoduli!
#
# function update!{MR<:DeforModelRed2DStress}(::Type{MR},
#                         self::MaterialDeformationLinear, ms; context...)
#         # Update material state.
#         #
#         # function [out, newms] = update (self, ms, context)
#         #
#         # Update material state.  Return the updated material state, and the
#         # requested quantity (default is the stress).
#         #   Call as:
#         #     [out,newms] = update(m, ms, context)
#         #  where
#         #     m=material
#         #     ms = material state
#         #     context=structure
#         #        with mandatory fields
#         #           strain=strain vector  in the local material
#         #               directions (which may be the same as the global coordinate
#         #               directions)
#         #        and optional fields
#         #           output=type of quantity to output, and interpreted by the
#         #               particular material; [] is returned when the material does not
#         #               recognize the requested quantity to indicate uninitialized
#         #               value.  It can be tested with isempty ().
#         #                  output ='Cauchy' - Cauchy stress; this is the default
#         #                      when output type is not specified.
#         #                  output ='princCauchy' - principal Cauchy stress;
#         #                  output ='pressure' - pressure;
#         #                  output ='vol_strain' - volumetric strain;
#         #           outputRm=optional orientation matrix in which output should
#         #               supplied
#         #
#         #   It is assumed that stress is output in m-component vector
#         #           form. m=3 for plane stress, m=4 for plane strain or axially
#         #           symmetric
#         #   The output arguments are
#         #     out=requested quantity
#         #     newms=new material state; don't forget that if the update is final
#         #           the material state newms must be assigned and stored.  Otherwise
#         #           the material update is lost!
#     #
#     Ev=FFlt[]                  # it is empty: we need to get it from context
#     output=:Cauchy
#     for arg in context
#         sy, val = arg
#         if sy==:strain
#             Ev=val
#         elseif sy==:output
#             output=val
#         end
#     end
#     D=zeros(FFlt,3,3)
#     tangentmoduli!(MR,self,D; context...)
#     tSigma = thermalstress(MR,self; context...);
#     stress = D * Ev + tSigma;
#     if output==:Cauchy
#         out = stress;
#     elseif output==:pressure
#         out = -sum(stress(1:2))/3;
#     elseif output==:volstrain
#         out = sum(Ev(1:2));     # actually this is probably incorrect
#         #for plane stress: the transverse strain is not accounted for
#     elseif output==:princCauchy
#         t=zeros(FFlt,3,3)
#         t = stress3vto3x3t!(stress,t);
#         ep=eig(t);
#         out =sort(ep[1]);
#     else
#         out = stress;
#     end
#     newms = ms;
#     return out,newms;
# end
#
# function thermalstress{MR<:DeforModelRed2DStress}(::Type{MR},
#                         self::MaterialDeformationLinear; context...)
# # Calculate vector of thermal stress components.
# #
# # function v = thermal_stress(self,context)
# #
# #   Call as:
# #     v = thermal_stress(m,context)
# #  where
# #     m=material
# #     context=structure; see the update() method
#     #
#     for arg in context
#         sy, val = arg
#         if sy==:dT
#             dT=val
#             D=zeros(FFlt,3,3)
#             tangentmoduli!(MR, self, D; context...);# local material stiffness
#             v = -D*(dT[1]*[self.property.CTE[1:2],0.0]);
#             return v
#             # switch self.reduction
#             #     case 'axisymm'
#             #         D = tangent_moduli(self, context);
#             #         v = -D*context.dT*[alphas(1:3).*ones(3, 1); 0];
#             #     case 'strain'
#             #         D=  self.property.tangent_moduli(context);% need 3-D
#             #         v = -context.dT*[D(1:2, 1:2)*(alphas(1:2).*ones(2,1))+...
#             #             alphas(3)*D(1:2,3); 0];
#
#         end
#     end
#     return  zeros(FFlt, 3, 1);
# end
# export thermalstress
#
#
#
# ################################################################################
# # Axially symmetric model
#
# function tangentmoduli!(::Type{DeforModelRed2DAxisymm},
#                         self::MaterialDeformationLinear,
#                         D::FFltMat; context...)
#     # # Calculate the material stiffness matrix.
#     # #
#     # # Arguments
#     # #     m=material
#     # #     context=structure with mandatory and optional fields that are required
#     # # by the specific material implementation.
#     # #
#     # # the output arguments are
#     # #     D=matrix 6x6 in the local material orientation Cartesian basis
#     # #
#
#     D3d=zeros(FFlt,6,6)
#     tangentmoduli3d!(self.property, D3d; context...);
#     for i=1:4
#         for j=1:4
#             D[j,i]= D3d[j,i];
#         end
#     end
#     return D
# end
# export tangentmoduli!
#
# function update!{MR<:DeforModelRed2DAxisymm}(::Type{MR},
#                         self::MaterialDeformationLinear, ms; context...)
#         # Update material state.
#         #
#         # function [out, newms] = update (self, ms, context)
#         #
#         # Update material state.  Return the updated material state, and the
#         # requested quantity (default is the stress).
#         #   Call as:
#         #     [out,newms] = update(m, ms, context)
#         #  where
#         #     m=material
#         #     ms = material state
#         #     context=structure
#         #        with mandatory fields
#         #           strain=strain vector  in the local material
#         #               directions (which may be the same as the global coordinate
#         #               directions)
#         #        and optional fields
#         #           output=type of quantity to output, and interpreted by the
#         #               particular material; [] is returned when the material does not
#         #               recognize the requested quantity to indicate uninitialized
#         #               value.  It can be tested with isempty ().
#         #                  output ='Cauchy' - Cauchy stress; this is the default
#         #                      when output type is not specified.
#         #                  output ='princCauchy' - principal Cauchy stress;
#         #                  output ='pressure' - pressure;
#         #                  output ='vol_strain' - volumetric strain;
#         #           outputRm=optional orientation matrix in which output should
#         #               supplied
#         #
#         #   It is assumed that stress is output in m-component vector
#         #           form. m=3 for plane stress, m=4 for plane strain or axially
#         #           symmetric
#         #   The output arguments are
#         #     out=requested quantity
#         #     newms=new material state; don't forget that if the update is final
#         #           the material state newms must be assigned and stored.  Otherwise
#         #           the material update is lost!
#     #
#     Ev=FFlt[]                  # it is empty: we need to get it from context
#     output=:Cauchy
#     for arg in context
#         sy, val = arg
#         if sy==:strain
#             Ev=val
#         elseif sy==:output
#             output=val
#         end
#     end
#     D=zeros(FFlt,4,4)
#     tangentmoduli!(MR,self,D; context...)
#     tSigma = thermalstress(MR,self; context...);
#     stress = D * Ev + tSigma;
#     if output==:Cauchy
#         out = stress;
#     elseif output==:pressure
#         out = -sum(stress[[1,2,4]])/3;
#     elseif output==:volstrain
#         out = sum(Ev[[1,2,4]]);
#     elseif output==:princCauchy
#         t=zeros(FFlt,3,3)
#         t = stress4vto3x3t!(stress[[1,2,4,3]],t);
#         ep=eig(t);
#         out =sort(ep[1]);
#     else
#         out = stress;
#     end
#     newms = ms;
#     return out,newms;
# end
#
# function thermalstress{MR<:DeforModelRed2DAxisymm}(::Type{MR},
#                         self::MaterialDeformationLinear; context...)
# # Calculate vector of thermal stress components.
# #
# # function v = thermal_stress(self,context)
# #
# #   Call as:
# #     v = thermal_stress(m,context)
# #  where
# #     m=material
# #     context=structure; see the update() method
#     #
#     for arg in context
#         sy, val = arg
#         if sy==:dT
#             dT=val
#             D=zeros(FFlt,4,4)
#             tangentmoduli!(MR, self, D; context...);# local material stiffness
#             v = -D*[self.property.CTE; 0.0]*dT;
#             return v
#             #     case 'strain'
#             #         D=  self.property.tangent_moduli(context);% need 3-D
#             #         v = -context.dT*[D(1:2, 1:2)*(alphas(1:2).*ones(2,1))+...
#             #             alphas(3)*D(1:2,3); 0];
#
#         end
#     end
#     return  zeros(FFlt, 4, 1);
# end
# export thermalstress
#
# ################################################################################
# # 1D model
#
# function tangentmoduli!{P<:PropertyDeformationLinear}(::Type{DeforModelRed1D},
#                         self::MaterialDeformationLinear{P},
#                         D::FFltMat; context...)
#     # # Calculate the material stiffness matrix.
#     # #
#     # # Arguments
#     # #     m=material
#     # #     context=structure with mandatory and optional fields that are required
#     # # by the specific material implementation.
#     # #
#     # # the output arguments are
#     # #     D=matrix 6x6 in the local material orientation Cartesian basis
#     # #
#
#     D3d=zeros(FFlt,6,6)
#     tangentmoduli3d!(self.property, D3d; context...);
#     D[1,1] = D3d[1, 1]- D3d[1,2:3]*D3d[2:3,2:3]\D3d[2:3,1];
#     return D
# end
# export tangentmoduli!
