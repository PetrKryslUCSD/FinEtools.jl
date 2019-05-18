"""
	IntegRuleModule

Module for numerical integration rules.
"""
module IntegRuleModule

using FinEtools.FTypesModule: FInt, FFlt, FCplxFlt, FFltVec, FIntVec, FFltMat, FIntMat, FMat, FVec, FDataDict

"""
    AbstractIntegRule

Abstract type for integration rule.
"""
abstract type AbstractIntegRule end

"""
	TriRule <: AbstractIntegRule

Triangular quadrature rule for integration on the standard triangle.
"""
struct TriRule <: AbstractIntegRule
    npts::FInt
    param_coords::Array{Float64, 2}
    weights::Array{Float64, 2}
end

"""
    TriRule(npts=1)

Type for triangular quadrature rule.  Used for integration of the standard
triangle, which is between 0 and 1 in both parametric coordinates.
`npts` = number of points (1-- one-point rule, 3 -- three-point rule, 6 -- six
 point rule, 9 --nine point rule, 10 -- Strang 10 point, order 13, degree of
 precision 7, rule), 12 and 13--twelve- and thirteen-point rule.
"""
function TriRule(npts=1)
    if npts == 1 # integrates exactly linear polynomials
        param_coords = [1.0/3. 1.0/3.];
        weights = reshape([1.0]/2.0,1,1);
    elseif npts == 3 # integrates exactly quadratic polynomials
        param_coords = [ 2.0/3 1.0/6; 1.0/6 2.0/3; 1.0/6 1.0/6 ];
        weights = [1.0/3 1.0/3 1.0/3]/2;
    elseif npts == 4 # integrates exactly quadratic polynomials
        param_coords = [0.333333333333333   0.333333333333333
                        0.200000000000000   0.200000000000000
                        0.600000000000000   0.200000000000000
                        0.200000000000000   0.600000000000000];
        weights = [-0.281250000000000
                    0.260416666666667
                    0.260416666666667
                    0.260416666666667];
    elseif npts == 6 # integrates exactly quartic polynomials
        param_coords = [ 0.816847572980459 0.091576213509771;
                        0.091576213509771 0.816847572980459;
                        0.091576213509771 0.091576213509771;
                        0.108103018168070 0.445948490915965;
                        0.445948490915965 0.108103018168070;
                        0.445948490915965 0.445948490915965];
        weights = [0.109951743655322*[1, 1, 1] 0.223381589678011*[1, 1, 1] ]/2;
    elseif npts == 7 # integrates exactly ? polynomials
        param_coords = [0.101286507323456   0.101286507323456
                        0.797426958353087   0.101286507323456
                        0.101286507323456   0.797426958353087
                        0.470142064105115   0.470142064105115
                        0.059715871789770   0.470142064105115
                        0.470142064105115   0.059715871789770
                        0.333333333333333   0.333333333333333];
        weights = [0.062969590272414
                    0.062969590272414
                    0.062969590272414
                    0.066197076394253
                    0.066197076394253
                    0.066197076394253
                    0.112500000000000];
   	elseif npts == 9 # integrates exactly ? polynomials
        param_coords = [   0.437525248383384   0.437525248383384
            0.124949503233232   0.437525248383384
            0.437525248383384   0.124949503233232
            0.165409927389841   0.037477420750088
            0.037477420750088   0.165409927389841
            0.797112651860071   0.165409927389841
            0.165409927389841   0.797112651860071
            0.037477420750088   0.797112651860071
            0.797112651860071   0.037477420750088];
        weights = [   0.205950504760887
            0.205950504760887
            0.205950504760887
            0.063691414286223
            0.063691414286223
            0.063691414286223
            0.063691414286223
            0.063691414286223
            0.063691414286223] ./ 2;
 	elseif npts == 12 # integrates exactly ? polynomials
        param_coords = [   0.063089014491502   0.063089014491502
            0.873821971016996   0.063089014491502
            0.063089014491502   0.873821971016996
            0.249286745170910   0.249286745170910
            0.501426509658179   0.249286745170910
            0.249286745170910   0.501426509658179
            0.310352451033785   0.053145049844816
            0.053145049844816   0.310352451033785
            0.636502499121399   0.310352451033785
            0.310352451033785   0.636502499121399
            0.053145049844816   0.636502499121399
            0.636502499121399   0.053145049844816];
        weights = [0.050844906370207
            0.050844906370207
            0.050844906370207
            0.116786275726379
            0.116786275726379
            0.116786275726379
            0.082851075618374
            0.082851075618374
            0.082851075618374
            0.082851075618374
            0.082851075618374
            0.082851075618374] ./ 2;
    elseif npts == 13 # integrates exactly ? polynomials
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
    rule  = TriRule(FInt(npts), reshape(param_coords, size(param_coords, 1), 2), reshape(weights, length(weights), 1))
    return  rule
end


"""
    GaussRule <: AbstractIntegRule

The Gauss rul, applicable for a tensor product of  intervals -1 <=x<= +1.
"""
struct GaussRule <: AbstractIntegRule
    dim::FInt
    order::FInt
    npts::FInt
    param_coords::Array{Float64, 2}
    weights::Array{Float64, 2}
end

"""
    GaussRule(dim=1, order=1)

Gauss rule.
"""
function GaussRule(dim=1, order=1)
    @assert 1 <= dim <= 3 "Gauss rule of dimension $(dim) not available"
    if     (order==1)
        param_coords = vec([ 0.0 ]);
        weights  = vec([ 2.0 ]);
    elseif (order==2)
        param_coords = vec([ -0.577350269189626 0.577350269189626 ]);
        weights  = vec([ 1.0 1.0 ]);
    elseif (order==3)
        param_coords = vec([ -0.774596669241483  0.0  0.774596669241483 ]);
        weights  = vec([ 0.5555555555555556 0.8888888888888889 0.5555555555555556 ]);
    elseif (order==4)
        param_coords = vec([ -0.86113631159405  -0.33998104358486   0.33998104358486   0.86113631159405]);
        weights  = vec([ 0.34785484513745   0.65214515486255   0.65214515486255   0.34785484513745]);
    elseif (order==5)
        param_coords = [ -0.906179845938664        -0.538469310105683         0.000000000000000         0.538469310105683         0.906179845938664];
        weights  = [0.236926885056189        0.478628670499367        0.568888888888889        0.478628670499367        0.236926885056189];
    elseif (order==6)
        param_coords = vec([-0.932469514203152
            -0.661209386466264
            -0.238619186083197
            0.238619186083197
            0.661209386466264
            0.932469514203152]);
        weights  = vec([0.171324492379171
            0.360761573048139
            0.467913934572691
            0.467913934572692
            0.360761573048139
            0.171324492379170]);
    elseif (order==7)
        param_coords = vec([-0.949107912342758
            -0.741531185599394
            -0.405845151377397
            -0.000000000000000
            0.405845151377397
            0.741531185599395
            0.949107912342758]);
        weights  = vec([0.129484966168870
            0.279705391489277
            0.381830050505119
            0.417959183673469
            0.381830050505119
            0.279705391489276
            0.129484966168870]);
    elseif (order==8)
        param_coords = vec([-0.960289856497536
            -0.796666477413627
            -0.525532409916329
            -0.183434642495650
            0.183434642495650
            0.525532409916329
            0.796666477413627
            0.960289856497536]);
        weights  = vec([0.101228536290376
            0.222381034453374
            0.313706645877887
            0.362683783378362
            0.362683783378362
            0.313706645877887
            0.222381034453374
            0.101228536290376]);
    elseif (order==9)
        param_coords = vec([-0.968160239507626
            -0.836031107326636
            -0.613371432700590
            -0.324253423403809
            0.000000000000000
            0.324253423403809
            0.613371432700591
            0.836031107326636
            0.968160239507626]);
        weights  = vec([0.081274388361574
            0.180648160694858
            0.260610696402935
            0.312347077040002
            0.330239355001259
            0.312347077040002
            0.260610696402936
            0.180648160694857
            0.081274388361574]);
    elseif (order==10)
        param_coords = vec([-0.973906528517171
            -0.865063366688985
            -0.679409568299025
            -0.433395394129247
            -0.148874338981631
            0.148874338981631
            0.433395394129247
            0.679409568299024
            0.865063366688984
            0.973906528517172]);
        weights  = vec([0.066671344308688
            0.149451349150581
            0.219086362515981
            0.269266719309996
            0.295524224714752
            0.295524224714753
            0.269266719309996
            0.219086362515982
            0.149451349150581
            0.066671344308688]);
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
    else # (dim==3)
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
    end
    return GaussRule(order, dim, FInt(npts), reshape(param_coords, size(param_coords, 1), dim), reshape(weights, length(weights), 1))
end


"""
    TetRule <: AbstractIntegRule

Tetrahedral quadrature rule, used for integration on the standard tetrahedron.
"""
struct TetRule <: AbstractIntegRule
    npts::FInt
    param_coords::Array{Float64, 2}
    weights::Array{Float64, 2}
end

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


"""
    PointRule <: AbstractIntegRule

Point quadrature rule, used for integration on the standard "point" shape.
"""
struct PointRule <: AbstractIntegRule
    npts::FInt
    param_coords::Array{Float64, 2}
    weights::Array{Float64, 2}
end

"""
    PointRule()

POINT integration rule.
"""
function PointRule()
    return PointRule(1, reshape([1.0], 1, 1), reshape([1.0], 1, 1))
end

"""
    SimplexRule <: AbstractIntegRule

Simplex quadrature rule.

Used for integration on the standard triangle or the standard tetrahedron.
"""
struct SimplexRule <: AbstractIntegRule
    npts::FInt
    param_coords::Array{Float64, 2}
    weights::Array{Float64, 2}
end

"""
    SimplexRule(dim=1, npts=1)

Return simplex rule, appropriate for the manifold dimension `dim`.
"""
function SimplexRule(dim=1, npts=1)
    @assert 0 <= dim <= 3 "Simplex rule of dimension $(dim) not available"

    if dim == 0
        r = PointRule()
        return SimplexRule(r.npts, r.param_coords, r.weights)
    elseif dim == 1
        r = GaussRule(dim, npts)
        return SimplexRule(r.npts, r.param_coords, r.weights)
    elseif dim == 2
        r = TriRule(npts)
        return SimplexRule(r.npts, r.param_coords, r.weights)
    else
        r = TetRule(npts)
        return SimplexRule(r.npts, r.param_coords, r.weights)
    end
end


"""
    TrapezoidalRule <: AbstractIntegRule

The trapezoidal rule.

The rule is applicable for a tensor product of  intervals -1 <=x<= +1.
"""
struct TrapezoidalRule <: AbstractIntegRule
    dim::FInt
    npts::FInt
    param_coords::Array{Float64, 2}
    weights::Array{Float64, 2}
end

"""
    TrapezoidalRule(dim=1)

Trapezoidal rule.
"""
function TrapezoidalRule(dim=1)
    @assert 1 <= dim <= 3 "Trapezoidal rule of dimension $(dim) not available"
    order = 2 # order of the trapezoidal rule: fixed
    param_coords = [ -1.0 1.0 ];
    weights  = [ 1.0 1.0 ];
    if (dim==1)
        param_coords = transpose(param_coords);
        weights = weights;
        npts=order^1
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
    else # (dim==3)
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
    end
    return TrapezoidalRule(dim,npts,param_coords,weights)
end

"""
    NodalSimplexRule <: AbstractIntegRule

The nodal-quadrature simplex rule.

The rule is applicable for line segments, triangles, tetrahedra.

!!! note

The quadrature points for a nodal quadrature rule must be listed in the order 
in which the nodes are used in the definition of the element!
"""
struct NodalSimplexRule <: AbstractIntegRule
    dim::FInt
    npts::FInt
    param_coords::Array{Float64, 2}
    weights::Array{Float64, 2}
end

"""
    NodalSimplexRule(dim=1)

Nodal-quadrature simplex rule.
"""
function NodalSimplexRule(dim=1)
    @assert 1 <= dim <= 3 "Nodal-quadrature simplex rule of dimension $(dim) not available"
    order = 1 # order of the rule
    if (dim==1)
        param_coords = reshape([-1.0 1.0]',2,1)
        weights = reshape([1.0 1.0],2,1)
    elseif (dim==2)
        param_coords = FFlt[0 0; 1 0; 0 1];
        weights = FFlt[1/3 1/3 1/3]'/2;
    else # (dim==3)
        param_coords = FFlt[0 0 0; 1 0 0; 0 1 0; 0 0 1];
        weights = FFlt[1, 1, 1, 1] /6.0 / 4;
    end
    npts = length(weights)
    return NodalSimplexRule(dim, npts, reshape(param_coords, npts, dim), reshape(weights, npts, 1))
end

"""
    NodalTensorProductRule <: AbstractIntegRule

The tensor-product nodal-quadrature rule.

The rule is applicable for line segments, quadrilaterals, hexahedra.

!!! note

The quadrature points for a nodal quadrature rule must be listed in the order 
in which the nodes are used in the definition of the element!
"""
struct NodalTensorProductRule <: AbstractIntegRule
    dim::FInt
    npts::FInt
    param_coords::Array{Float64, 2}
    weights::Array{Float64, 2}
end

"""
    NodalTensorProductRule(dim=1)

Nodal-quadrature tensor-product rule.
"""
function NodalTensorProductRule(dim=1)
    @assert 1 <= dim <= 3 "Nodal-quadrature tensor-product rule of dimension $(dim) not available"
    order = 1 # order of the rule
    if (dim==1)
        param_coords = FFlt[-1.0 -1.0];
            	weights = fill(1.0, 2);
    elseif (dim==2)
    	param_coords = FFlt[-1.0 -1.0; 1.0 -1.0; 1.0  1.0;-1.0  1.0];
    	weights = fill(1.0, 4);
    else # (dim==3)
        param_coords = FFlt[-1.0 -1.0 -1.0; 1.0 -1.0 -1.0; 1.0  1.0 -1.0;-1.0  1.0 -1.0;-1.0 -1.0  1.0; 1.0 -1.0  1.0; 1.0  1.0  1.0;-1.0  1.0  1.0];
        weights = fill(1.0, 8);
    end
    npts = length(weights)
    return NodalTensorProductRule(dim, npts, reshape(param_coords, npts, dim), reshape(weights, npts, 1))
end

end