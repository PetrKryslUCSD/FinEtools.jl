"""
    AlgoBaseModule

Module for base  algorithms.
"""
module AlgoBaseModule

using FinEtools.FTypesModule

export FDataDict
export dcheck!, richextrapol

# Type for the model-data packaging system (used by all FinEtools algorithms)
const FDataDict = Dict{String, Any}

function keymatch(key::String, allowed_keys::Array{String})
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
        matched_key = keymatch(key, recognized_keys)
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


end
