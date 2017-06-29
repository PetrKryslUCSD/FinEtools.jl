module AlgoBaseModule

# Type for the model-data packaging system (used by all FinEtools algorithms)
const FDataDict = Dict{String, Any}
export FDataDict

function dmatch(key::String, allowed_keys::Array{String})
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
  for  key in keys(d)
    matched_key = dmatch(key, recognized_keys)
    if matched_key == nothing
      # warn("No key matches \"$key\"")
    else
      if key !=  matched_key
        warn("Key \"$key\" not fully matched (partial match \"$matched_key\")")
        # d[matched_key] = d[key]
      end
    end
  end
end
export dcheck!

end
