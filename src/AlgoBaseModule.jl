"""
    AlgoBaseModule

Module for base  algorithms.
"""
module AlgoBaseModule

import ..FEMMBaseModule: integratefieldfunction, transferfield!
import LinearAlgebra: norm, dot
import Statistics: mean

function _keymatch(key::String, allowed_keys::Array{String})
    matched_key = nothing
    for j in eachindex(allowed_keys)
        m = match(Regex("^$key"), allowed_keys[j])
        if (m != nothing)
            matched_key = allowed_keys[j]
            break
        end
    end
    return matched_key
end

function dcheck!(d::Dict{String,Any}, recognized_keys::Array{String})
    notmatched = fill("", 0)
    for key in keys(d)
        matched_key = _keymatch(key, recognized_keys)
        if matched_key == nothing
            push!(notmatched, "Key \"$key\" not matched")
        else
            if key != matched_key
                push!(
                    notmatched,
                    "Key \"$key\" not fully matched (partial match \"$matched_key\")",
                )
            end
        end
    end
    return notmatched
end

"""
    richextrapol(solns::T, params::T; lower_conv_rate = 0.001, upper_conv_rate = 10.0) where {T<:AbstractArray{Tn} where {Tn}}

Richardson extrapolation.

# Arguments
- `solns` =  array of three solution values
- `params` = array of values of three parameters for which the `solns` have been obtained. 

The assumption is that the error of the solution is expanded in a Taylor series,
and only the first term in the Taylor series is kept.
    ```
        qex - qapprox ~ C param^beta
    ```
Here `qex` is the true solution, `qapprox` is an approximate solution, `param`
is the element size, or the relative element size, in other words the parameter
of the extrapolation, and `beta` is the convergence rate. The constant `C` is 
the third unknown quantity in this expansion. If we obtain three successive 
approximations, we can solve for the three unknown quantities, `qex`, `beta`, 
and `C`.

It is assumed that the first solution is obtained for the largest value of the
extrapolation parameter, while the last solution in the list is obtained for
the smallest value of the extrapolation parameter: 
    ```
        params[1] > params[2] > params[3]
    ```

# Output
- `solnestim`= estimate of the asymptotic solution from the data points in the
  `solns` array
- `beta`= convergence rate
- `c` = constant in the estimate `error=c*h^beta`
- `maxresidual` = maximum residual after equations from which the above quantities were
  solved (this is a measure of how accurately was the system solved).
"""
function richextrapol(
    solns::T,
    params::T;
    lower_conv_rate = 0.001,
    upper_conv_rate = 10.0,
) where {T<:AbstractArray{Tn} where {Tn}}
    # These two constants may needs to be tweaked in special cases. They are the
    # lower and upper bound on the convergence rate.
    lower, upper = lower_conv_rate, upper_conv_rate
    if !((params[1] > params[2]) && (params[2] > params[3]))
        error("Extrapolation parameter in the wrong order: must be largest to smallest")
    end
    solnn = maximum(abs.(solns))
    nsolns = solns ./ solnn # Normalize data for robust calculation
    nhs1, nhs2, nhs3 = params ./ maximum(params) # Normalize the parameter values
    napproxerror1, napproxerror2 = diff(nsolns) # Normalized approximate errors
    napperr1, napperr2 =
        (napproxerror1, napproxerror2) ./ max(abs(napproxerror1), abs(napproxerror2))
    maxfun = -Inf
    minfun = Inf
    for y in lower:lower:upper
        a = napperr1 * (nhs2^y - nhs3^y) - napperr2 * (nhs1^y - nhs2^y)
        maxfun = max(maxfun, a)
        minfun = min(minfun, a)
    end
    if minfun * maxfun > 0.0
        error("Convergence rate function does not cross zero value: no solution")
    end
    napperr1, napperr2 = (napperr1, napperr2) ./ (maxfun - minfun)  # Normalize 
    fun = y -> napperr1 * (nhs2^y - nhs3^y) - napperr2 * (nhs1^y - nhs2^y)
    # x = collect(lower:lower:upper)
    # y = fun.(x)
    # p = plot(x=x, y = y, Geom.line);
    # draw(PDF("File.pdf", 16cm, 9cm), p)
    tolx, tolf = 1.0e-6 * lower, 1.0e-12  # Note: the function value is normalized to 1.0
    fl = fun(lower)
    fu = fun(upper)
    if fl * fu > 0.0
        error("Convergence rate couldn't be solved, no bracket found")
    end
    beta = bisect(fun, lower, upper, fl, fu, tolx, tolf)
    beta = (beta[1] + beta[2]) / 2.0
    c = napproxerror1 / (params[1]^beta - params[2]^beta)
    nestimtrueerror = c * params[3]^beta
    solnestim = (nsolns[3] + nestimtrueerror) * solnn
    c = c * solnn # adjust for not-normalized data
    # just to check things, calculate the residual
    maxresidual = 0.0
    for I in 1:3
        maxresidual = max(maxresidual, abs((solnestim - solns[I]) - c * params[I]^beta)) # this should be close to zero
    end
    return solnestim, beta, c, maxresidual
end

"""
    richextrapoluniform(solns::T, params::T) where {T<:AbstractArray{Tn} where {Tn}}

Richardson extrapolation.

# Argument
- `solns` =  array of solution values
- `params` = array of values of parameters for which the `solns` have been
  obtained. This function is applicable only to fixed (uniform) ratio between
  the mesh sizes, `params[1]/params[2) = params[2)/params[3)`.
  
# Output 
- `solnestim`= estimate of the asymptotic solution from the data points in the
  `solns` array
- `beta`= convergence rate
- `c` = constant in the estimate `error=c*h^beta`
- `residual` = residual after equations from which the above quantities were
  solved (this is a measure of how accurately was the system solved).
"""
function richextrapoluniform(solns::T, params::T) where {T<:AbstractArray{Tn} where {Tn}}
    @assert abs(params[1] / params[2] - params[2] / params[3]) < 1e-6 "Parameter pair ratio not fixed"
    nsolns = solns ./ solns[1]
    solnestim =
        ((
            -(nsolns[1] * nsolns[3] - nsolns[2]^2) /
            (2 * nsolns[2] - nsolns[1] - nsolns[3])
        )) * solns[1]
    if (solnestim - solns[1]) <= 0
        beta =
            log((solnestim - solns[2]) / (solnestim - solns[3])) /
            log(params[2] / params[3])
    else
        beta =
            log((solnestim - solns[1]) / (solnestim - solns[3])) /
            log(params[1] / params[3])
    end
    # just to check things, calculate the residual
    c = (solnestim - solns[1]) / params[1]^beta
    residual = fill(zero(solns[1]), 3)
    for I in 1:3
        residual[I] = (solnestim - solns[I]) - c * params[I]^beta # this should be close to zero
    end
    return solnestim, beta, c, residual
end

"""
    bisect(fun, xl, xu, fl, fu, tolx, tolf)

Implementation of the bisection method.

Tolerance both on `x` and on `f(x)` is used.
- `fun` = function,
- `xl`,`xu`= lower and upper value of the bracket,
- `fl`,`fu`= function value at the lower and upper limit of the bracket.
The true values must have opposite signs (that is they must constitute a bracket). Otherwise  this algorithm will fail.
- `tolx`= tolerance on the location of the root,
- `tolf`= tolerance on the function value
"""
function bisect(fun, xl, xu, fl, fu, tolx, tolf)
    iter = 0
    if (xl > xu)
        temp = xl
        xl = xu
        xu = temp
    end
    # fl = fun(xl);
    # fu = fun(xu);
    # @assert fl*fu < 0.0 "Need to get a bracket"
    while true
        xr = (xu + xl) / 2.0 # bisect interval
        fr = fun(xr) # value at the midpoint
        if (fr * fl < 0.0)
            xu = xr
            fu = fr# upper --> midpoint
        elseif (fr == 0.0)
            xl = xr
            xu = xr# exactly at the root
        else
            xl = xr
            fl = fr# lower --> midpoint
        end
        if (abs(xu - xl) < tolx) || (abs(fr) < tolf)
            break # We are done
        end
        iter = iter + 1
    end
    return xl, xu
end

"""
    bisect(fun, xl, xu, tolx, tolf)

Implementation of the bisection method.

Tolerance both on `x` and on `f(x)` is used.
- `fun` = function,
- `xl`= lower value of the bracket,
- `xu`= upper Value of the bracket,
- `tolx`= tolerance on the location of the root,
- `tolf`= tolerance on the function value
"""
function bisect(fun, xl, xu, tolx, tolf)
    iter = 0
    if (xl > xu)
        temp = xl
        xl = xu
        xu = temp
    end
    fl = fun(xl)
    fu = fun(xu)
    @assert fl * fu < 0.0 "Need to get a bracket"
    while true
        xr = (xu + xl) / 2.0 # bisect interval
        fr = fun(xr) # value at the midpoint
        if (fr * fl < 0.0)
            xu = xr
            fu = fr# upper --> midpoint
        elseif (fr == 0.0)
            xl = xr
            xu = xr# exactly at the root
        else
            xl = xr
            fl = fr# lower --> midpoint
        end
        if (abs(xu - xl) < tolx) || (abs(fr) < tolf)
            break # We are done
        end
        iter = iter + 1
    end
    return xl, xu
end

"""
    fieldnorm(modeldata)

Compute norm of the target field.  

# Argument
- `modeldata` = data dictionary, mandatory keys:
    + `fens` = finite element node set
    + `regions` = array of regions
    + `targetfields` = array of fields, one for each region
    + `geom` = geometry field
    + `elementsize` = representative element size,

# Output
- Norm of the field as floating-point scalar.
"""
function fieldnorm(modeldata)
    fens = modeldata["fens"]
    regions = modeldata["regions"]
    targetfields = modeldata["targetfields"]
    geom = modeldata["geom"]

    @assert length(regions) == length(targetfields)

    fnorm = 0.0 # Initialize the norm of the difference
    for i in eachindex(regions)
        # Compute the addition to the norm of the field
        fnorm += integratefieldfunction(
            regions[i]["femm"],
            geom,
            targetfields[i],
            (x, v) -> norm(v)^2,
            0.0,
        )
    end

    return sqrt(fnorm)
end

"""
    fielddiffnorm(modeldatacoarse, modeldatafine)

Compute norm of the difference of the fields.  

# Arguments
- `modeldatacoarse`, `modeldatafine` = data dictionaries.

For both the "coarse"- and "fine"-mesh `modeldata` the data dictionaries need to contain the mandatory keys:
- `"fens"` = finite element node set
- `"regions"` = array of regions
- `"targetfields"` = array of fields, one for each region
- `"geom"` = geometry field
- `"elementsize"` = representative element size,
- `"geometricaltolerance"` = geometrical tolerance (used in field transfer;
  refer to the documentation of `transferfield!`)
- `"parametrictolerance"` = parametric tolerance (used in field transfer; refer
  to the documentation of `transferfield!`)

# Output
- Norm of the field as floating-point scalar.
"""
function fielddiffnorm(modeldatacoarse, modeldatafine)
    # Load coarse-mesh data
    fenscoarse = modeldatacoarse["fens"]
    regionscoarse = modeldatacoarse["regions"]
    targetfieldscoarse = modeldatacoarse["targetfields"]

    # Load fine-mesh data
    fensfine = modeldatafine["fens"]
    regionsfine = modeldatafine["regions"]
    targetfieldsfine = modeldatafine["targetfields"]

    # Tolerances
    geometricaltolerance = modeldatacoarse["geometricaltolerance"]
    parametrictolerance = get(modeldatacoarse, "parametrictolerance", 0.01)

    geom = modeldatafine["geom"]

    @assert length(regionscoarse) == length(regionsfine)
    @assert length(regionscoarse) == length(targetfieldscoarse)
    @assert length(regionsfine) == length(targetfieldsfine)

    diffnorm = 0.0 # Initialize the norm of the difference
    for i in eachindex(regionscoarse)
        # Transfer the result from the coarse mesh to the fine-mesh
        ffine = targetfieldsfine[i]
        fcoarse = targetfieldscoarse[i]
        fcoarsetransferred = deepcopy(ffine)
        fcoarsetransferred = transferfield!(
            fcoarsetransferred,
            fensfine,
            regionsfine[i]["femm"].integdomain.fes,
            fcoarse,
            fenscoarse,
            regionscoarse[i]["femm"].integdomain.fes,
            geometricaltolerance;
            parametrictolerance = parametrictolerance,
        )
        # Form the difference  field
        diffff = deepcopy(fcoarsetransferred)
        diffff.values[:] = ffine.values - fcoarsetransferred.values
        # Compute the addition to the norm of the difference
        diffnorm += integratefieldfunction(
            regionsfine[i]["femm"],
            geom,
            diffff,
            (x, v) -> norm(v)^2,
            0.0,
        )
    end

    return sqrt(diffnorm)
end

"""
    evalconvergencestudy(modeldatasequence)

Evaluate a convergence study from a model-data sequence.  

- `modeldatasequence` = array of `modeldata` dictionaries. At least two must be included.

Refer to methods `fieldnorm` and `fielddiffnorm` for details 
on the required keys in the dictionaries.

# Output
- `elementsizes` = element size array, 
- `errornorms` = norms of the error, 
- `convergencerate` = rate of convergence
"""
function evalconvergencestudy(modeldatasequence)
    # Find the element sizes
    elementsizes = [md["elementsize"] for md in modeldatasequence]

    # The  finest solution will provide the norm used in the normalization of the errors
    finestsolnorm = fieldnorm(modeldatasequence[end])
    # Compute the approximate errors  as the differences of successive solutions
    diffnorms = Float64[]
    for i in 1:(length(modeldatasequence)-1)
        push!(diffnorms, fielddiffnorm(modeldatasequence[i], modeldatasequence[i+1]))
    end
    # Normalize the errors
    errornorms = diffnorms ./ finestsolnorm
    # Compute the convergence rate
    f = log.(vec(errornorms))
    A = hcat(log.(vec(elementsizes[1:end-1])), ones(size(f)))
    p = A \ f
    convergencerate = p[1]

    return elementsizes, errornorms, convergencerate
end

"""
    conjugategradient(A::MT, b::Vector{T}, x0::Vector{T}, maxiter) where {MT, T<:Number}

Compute one or more iterations of the conjugate gradient process.  
"""
function conjugategradient(A::MT, b::Vector{T}, x0::Vector{T}, maxiter) where {MT, T<:Number}
    x = deepcopy(x0)
    gt = deepcopy(x0)
    d = deepcopy(x0)
    gt .= A * x - b# transpose of gradient: g = x'*A-b';
    @. d = gt
    for iter in 1:maxiter
        Ad = A * d
        rho = dot(d, Ad)
        alpha = -dot(gt, d) / rho # alpha =(-g*d)/(d'*A*d);
        @. x += alpha * d
        gt .= A * x - b# negative transpose of gradient
        beta = dot(gt, Ad) / rho # beta =(g*A*d)/(d'*A*d);
        @. d = beta * d - gt
    end
    return x
end

"""
    qtrap(ps::VecOrMat{T}, xs::VecOrMat{T}) where {T<:Number}

Compute the area under the curve given by a set of parameters along 
an interval and the values of the 'function' at those parameter values.  
The parameter values need not be uniformly distributed.

Trapezoidal rule is used to evaluate the integral. The 'function' is 
assumed to vary linearly inbetween the given points.
"""
function qtrap(ps::VecOrMat{T}, xs::VecOrMat{T}) where {T<:Number}
    @assert length(ps) == length(xs)
    num = T(0.0)
    for i in 1:length(ps)-1
        num += (ps[i+1] - ps[i]) * (xs[i+1] + xs[i]) / 2.0
    end
    return num
end

"""
    qcovariance(ps::VecOrMat{T}, xs::VecOrMat{T}, ys::VecOrMat{T}; ws = nothing) where {T<:Number}

Compute the covariance for two 'functions' given by the arrays `xs` and `ys` 
at the values of the parameter `ps`. `ws` is the optional weights vector;  
if it is not supplied, uniformly distributed default weights are assumed.  

Notes: 

– The mean is subtracted from both functions. 
– This function is not particularly efficient: it computes the mean of both
  functions and it allocates arrays instead of overwriting the contents of the
  arguments.
"""
function qcovariance(ps::VecOrMat{T}, xs::VecOrMat{T}, ys::VecOrMat{T}; ws = nothing) where {T<:Number}
    @assert length(ps) == length(xs) == length(ys)
    if (ws == nothing)
        ws = ones(T, length(ps))
    end
    @assert length(ws) == length(xs)
    ws[:] = ws ./ qtrap(ps, ws) # Make sure the integral of the weight is equal to 1.0
    xmean = qtrap(ps, xs .* ws)
    ymean = qtrap(ps, ys .* ws)
    zxs = deepcopy(xs) .- xmean
    zys = deepcopy(ys) .- ymean
    return qtrap(ps, zxs .* zys .* ws)
end

"""
    qvariance(ps, xs; ws = nothing)

Compute the variance of a function given by the array `xs` at 
the values of the parameter `ps`. `ws` is the optional weights vector  
with unit default weights.    
"""
qvariance(ps, xs; ws = nothing) = qcovariance(ps, xs, xs; ws = ws)

"""
    penaltyebc!(K, F, dofnums, prescribedvalues, penfact)

Apply penalty essential boundary conditions.

# Arguments
- `K` = stiffness matrix
- `F` = global load vector 
- `dofnums`, `prescribedvalues` = arrays computed by `prescribeddofs()`
- `penfact` = penalty multiplier, in relative terms: how many times the
  maximum absolute value of the diagonal elements should the penalty term be?

# Output
- Updated matrix `K` and vector `F`.
"""
function penaltyebc!(K, F, dofnums, prescribedvalues, penfact)
    maxdiagK = 0.0
    for j in eachindex(dofnums)
        gj = dofnums[j]
        maxdiagK = max(abs(K[gj, gj]), maxdiagK)
    end
    penalty = penfact * maxdiagK
    for j in eachindex(dofnums)
        gj = dofnums[j]
        K[gj, gj] = K[gj, gj] + penalty
        F[gj] = F[gj] + penalty * prescribedvalues[j]
    end
    return K, F
end

end
