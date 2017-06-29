"""
    IntegRuleModule

Module for numerical integration rules.
"""
module IntegRuleModule

using FinEtools.FTypesModule

abstract type IntegRule end
export IntegRule

"""
    TriRule

Class of triangular quadrature rule.
Used for integration on the standard triangle.
"""
type TriRule <: IntegRule
  npts::FInt
  param_coords::FFltMat
  weights::FFltMat
end
export TriRule

"""
    TriRule(npts=1)

Type for triangular quadrature rule.  Used for integration of the standard
triangle, which is between 0 and 1 in both parametric coordinates.
`npts` = number of points (1-- one-point rule, 3 -- three-point rule,
     6 -- six point rule, 10 -- Strang 10 point, order 13, degree of precision 7, rule).
"""
function TriRule(npts=1)
  if npts==1 # integrates exactly linear polynomials
    param_coords = [1./3. 1./3.];
    weights = reshape([1.0]/2.0,1,1);
  elseif npts==3 # integrates exactly quadratic polynomials
    param_coords = [ 2/3 1/6; 1/6 2/3; 1/6 1/6 ];
    weights = [1/3 1/3 1/3]/2;
  elseif npts==6 # integrates exactly quartic polynomials
    param_coords = [ 0.816847572980459 0.091576213509771;
    0.091576213509771 0.816847572980459;
    0.091576213509771 0.091576213509771;
    0.108103018168070 0.445948490915965;
    0.445948490915965 0.108103018168070;
    0.445948490915965 0.445948490915965];
    weights = [0.109951743655322*[1, 1, 1] 0.223381589678011*[1, 1, 1] ]/2;
  elseif npts==10 # integrates exactly quartic polynomials
    param_coords = [0.333333333333333  0.333333333333333
    0.479308067841923  0.260345966079038
    0.260345966079038  0.479308067841923
    0.260345966079038  0.260345966079038
    0.869739794195568  0.065130102902216
    0.065130102902216  0.869739794195568
    0.065130102902216  0.065130102902216
    0.638444188569809  0.312865496004875
    0.638444188569809  0.048690315425316
    0.312865496004875  0.638444188569809
    0.312865496004875  0.048690315425316
    0.048690315425316  0.638444188569809
    0.048690315425316  0.312865496004875
    ];
    weights = [ -0.149570044467670
    0.175615257433204
    0.175615257433204
    0.175615257433204
    0.053347235608839
    0.053347235608839
    0.053347235608839
    0.077113760890257
    0.077113760890257
    0.077113760890257
    0.077113760890257
    0.077113760890257
    0.077113760890257
    ]'/2;
  else
    #nothing doing: this input is wrong
    error( "Unknown number of integration points" )
  end
  rule  = TriRule(npts, param_coords, weights)
  return  rule
end


"""
    GaussRule

Class of the Gauss rule.
The rule is applicable for a tensor product of  intervals -1 <=x<= +1.
"""
type GaussRule <: IntegRule
    order::FInt
    dim::FInt
    npts::FInt
    param_coords::FFltMat
    weights::FFltMat
end
export GaussRule

"""
  GaussRule(dim=1, order=1)

Gauss rule.
"""
function GaussRule(dim=1, order=1)
  @assert 1 <= dim <= 3
  if     (order==1)
    param_coords = [ 0.0 ];
    weights  = [ 2.0 ];
  elseif (order==2)
    param_coords = [ -0.577350269189626 0.577350269189626 ];
    weights  = [ 1.0 1.0 ];
  elseif (order==3)
    param_coords = [ -0.774596669241483  0.0  0.774596669241483 ];
    weights  = [ 0.5555555555555556 0.8888888888888889 0.5555555555555556 ];
  elseif (order==4)
    param_coords = [ -0.86113631159405  -0.33998104358486   0.33998104358486   0.86113631159405];
    weights  = [ 0.34785484513745   0.65214515486255   0.65214515486255   0.34785484513745];
  else
    error("Gauss rule of order $(order) not available")
  end
  if (dim==1)
    param_coords = transpose(param_coords);
    weights = weights;
    npts=order
    param_coords= reshape(param_coords',npts,1)
    weights= reshape(weights,npts,1)
  elseif (dim==2)
    pc=param_coords;
    w=weights;
    npts=order^2
    param_coords=zeros(FFlt,npts,2);
    weights=zeros(FFlt,npts,1);
    r=1
    for i=1:order
      for j=1:order
        param_coords[r,:]=[pc[i] pc[j]];
        weights[r]=w[i]*w[j];
        r=r+1
      end
    end
  elseif (dim==3)
    pc=param_coords;
    w=weights;
    npts=order^3
    param_coords=zeros(FFlt,npts,3);
    weights=zeros(FFlt,npts,1);
    r=1
    for i=1:order
      for j=1:order
        for k=1:order
          param_coords[r,:]=[pc[i] pc[j] pc[k]];
          weights[r]=w[i]*w[j]*w[k];
          r=r+1
        end
      end
    end
  else
    error("Gauss rule of dimension $(dim) not available")
  end
  return GaussRule(order,dim,npts,param_coords,weights)
end


"""
    TetRule

Class of tetrahedral quadrature rule.
Used for integration on the standard tetrahedron.
"""
type TetRule <: IntegRule
  npts::FInt
  param_coords::FFltMat
  weights::FFltMat
end
export TetRule

"""
    TetRule(npts=1)

Tetrahedral integration rule.

npts=number of points (1-- one-point rule, 4 -- four-point rule, 5 -- five point
rule).
"""
function TetRule(npts=1)
  if npts == 1 # integrates exactly linear polynomials
    param_coords = reshape([0.25,0.25,0.25],1,3);
    weights = reshape([1.0]/6.0,1,1);
  elseif npts == 4 # integrates exactly quadratic polynomials
    param_coords = [[0.13819660 0.13819660 0.13819660];
    [0.58541020 0.13819660 0.13819660];
    [0.13819660 0.58541020 0.13819660];
    [0.13819660 0.13819660 0.58541020]];;
    weights = [ 0.041666666666666666667   0.041666666666666666667   0.041666666666666666667   0.041666666666666666667];
  elseif npts == 5 #  Zienkiewicz #3.
    a =   1.0 / 6.0;
    b =   0.25;
    c =   0.5;
    d = - 0.8;
    e =   0.45;
    param_coords = [[b b b];
    [c a a];
    [a c a];
    [a a c];
    [a a a]];
    weights = [d  e  e  e  e]/6;
  else
    #nothing doing: this input is wrong
    error( "Unknown number of integration points" )
  end

  return  TetRule(npts, param_coords, reshape(weights,length(weights),1))
end

end
