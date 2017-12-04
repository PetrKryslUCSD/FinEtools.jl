"""
    AlgoBaseModule

Module for base  algorithms.
"""
module AlgoBaseModule

using FinEtools
using FinEtools.FEMMBaseModule

function _keymatch(key::String, allowed_keys::Array{String})
    matched_key = nothing
    for  j = 1:length(allowed_keys)
        m = match(Regex("^$key"), allowed_keys[j])
        if (m != nothing)
            matched_key = allowed_keys[j]
            break
        end
    end
    return matched_key
end

function dcheck!(d::FDataDict, recognized_keys::Array{String})
    notmatched = Array{String}(0)
    for  key in keys(d)
        matched_key = _keymatch(key, recognized_keys)
        if matched_key == nothing
            push!(notmatched, "Key \"$key\" not matched")
        else
            if key !=  matched_key
                push!(notmatched, "Key \"$key\" not fully matched (partial match \"$matched_key\")")
            end
        end
    end
    return notmatched
end

"""
    richextrapol(solns::AbstractArray{Real}, params::AbstractArray{Real})

Richardson extrapolation.

`solns` =  array of solution values
`params` = array of values of parameters for which the solns have been obtained

This function is applicable only to fixed ratio between the mesh sizes,
  `params[1]/params[2) = params[2)/params[3)`.
Output:
solnestim= estimate of the asymptotic solution from the data points in the solns array
beta= convergence rate
c = constant in the estimate "error=c*h^beta"
residual = residual after equations from which the above quantities were
     solved (this is a measure of how accurately was the system solved).
"""
function richextrapol(solns::FFltVec, params::FFltVec)
    @assert abs(params[1]/params[2] - params[2]/params[3]) < 1e-6 "Parameter pair ratio not fixed"
    nsolns = solns./solns[1];
    solnestim = ((-(nsolns[1]*nsolns[3]-nsolns[2]^2)/(2*nsolns[2]-nsolns[1]-nsolns[3])))*solns[1];
    if (solnestim-solns[1]) <= 0
        beta = log((solnestim-solns[2])/(solnestim-solns[3]))/log(params[2]/params[3]);
    else
        beta = log((solnestim-solns[1])/(solnestim-solns[3]))/log(params[1]/params[3]);
    end
    # just to check things, calculate the residual
    c = (solnestim-solns[1])/params[1]^beta;
    residual = zeros(FFlt,3)
    for I =1:3
        residual[I] = (solnestim-solns[I])-c*params[I]^beta; # this should be close to zero
    end
    return solnestim, beta, c, residual
end


"""
    fieldnorm(modeldata)

Compute norm of the target field.  

`modeldata` = data dictionary, mandatory keys:
    - "fens" = finite element node set
    - "regions" = array of regions
    - "targetfields" = array of fields, one for each region
    - "geom" = geometry field
    - "elementsize" = representative element size,
"""
function fieldnorm(modeldata)
    fens = modeldata["fens"]
    regions = modeldata["regions"]
    targetfields = modeldata["targetfields"]
    geom = modeldata["geom"]

    @assert length(regions) == length(targetfields)

    fnorm = 0.0 # Initialize the norm of the difference
    for i = 1:length(regions)
        # Compute the addition to the norm of the field
        fnorm += integratefieldfunction(regions[i]["femm"], geom, targetfields[i], (x, v) -> norm(v)^2, 0.0)
    end

    return sqrt(fnorm)
end

"""
    fielddiffnorm(modeldatacoarse, modeldatafine)

Compute norm of the difference of the target fields.  

For both the "coarse"- and "fine"-mesh `modeldata` the data dictionaries need to contain the mandatory keys:
- "fens" = finite element node set
- "regions" = array of regions
- "targetfields" = array of fields, one for each region
- "geom" = geometry field
- "elementsize" = representative element size,
- "geometricaltolerance" = geometrical tolerance (used in field transfer; refer to the documentation of `transferfield!`)
- "parametrictolerance" = parametric tolerance (used in field transfer; refer to the documentation of `transferfield!`)
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
    parametrictolerance = get(modeldatacoarse, "parametrictolerance", 0.01);
    
    geom = modeldatafine["geom"]

    @assert length(regionscoarse) == length(regionsfine)
    @assert length(regionscoarse) == length(targetfieldscoarse)
    @assert length(regionsfine) == length(targetfieldsfine)

    diffnorm = 0.0 # Initialize the norm of the difference
    for i = 1:length(regionscoarse)
        # Transfer the result from the coarse mesh to the fine-mesh
        ffine = targetfieldsfine[i]
        fcoarse = targetfieldscoarse[i]
        fcoarsetransferred = deepcopy(ffine)
        fcoarsetransferred = transferfield!(fcoarsetransferred,
            fensfine, regionsfine[i]["femm"].integdata.fes, fcoarse,
            fenscoarse, regionscoarse[i]["femm"].integdata.fes, geometricaltolerance; parametrictolerance = parametrictolerance)
        # Form the difference  field
        diffff = deepcopy(fcoarsetransferred)
        diffff.values[:] = ffine.values - fcoarsetransferred.values
        # Compute the addition to the norm of the difference
        diffnorm += integratefieldfunction(regionsfine[i]["femm"], geom, diffff, (x, v) -> norm(v)^2, 0.0)
    end

    return sqrt(diffnorm)
end

"""
    evalconvergencestudy(modeldatasequence, File)

Evaluate a convergence study from a model-data sequence.  

`modeldatasequence` = array of `modeldata` dictionaries.

Refer to methods `fieldnorm` and `fielddiffnorm` for details on the required keys in the dictionaries.
"""
function evalconvergencestudy(modeldatasequence)
    # Find the element sizes
    elementsizes = [md["elementsize"] for md in modeldatasequence]
    
    # The  finest solution will provide the norm used in the normalization of the errors
    finestsolnorm = fieldnorm(modeldatasequence[end])
    # Compute the approximate errors  as the differences of successive solutions
    diffnorms = FFlt[]
    for i = 1:(length(modeldatasequence) - 1)
        push!(diffnorms, fielddiffnorm(modeldatasequence[i], modeldatasequence[i + 1]))
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

end
