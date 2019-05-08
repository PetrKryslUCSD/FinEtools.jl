"""
    AlgoHeatDiffModule

Module for algorithms in linear heat conduction/diffusion  models.
"""
module AlgoHeatDiffModule

using FinEtools.FTypesModule: FInt, FFlt, FCplxFlt, FFltVec, FIntVec, FFltMat, FIntMat, FMat, FVec, FDataDict
import FinEtools.AlgoBaseModule: dcheck!
import FinEtools.FieldModule: ndofs, setebc!, numberdofs!, applyebc!, scattersysvec!
import FinEtools.NodalFieldModule: NodalField, nnodes
import FinEtools.FEMMBaseModule: associategeometry!, distribloads
import FinEtools.FEMMHeatDiffModule: conductivity, nzebcloadsconductivity
import FinEtools.FEMMHeatDiffSurfModule: surfacetransfer, surfacetransferloads, nzebcsurfacetransferloads
import FinEtools.ForceIntensityModule: ForceIntensity
import FinEtools.MeshSelectionModule: connectednodes
import SparseArrays: spzeros
import LinearAlgebra: cholesky

"""
    steadystate(modeldata::FDataDict)

Steady-state heat conduction solver.

# Argument
`modeldata` = dictionary with items

- `"fens"`  = finite element node set
- `"regions"`  = array of region dictionaries
- `"essential_bcs"` = array of essential boundary condition dictionaries
- `"convection_bcs"` = array of convection boundary condition dictionaries
- `"flux_bcs"` = array of flux boundary condition dictionaries

For each region (connected piece of the domain made of a particular material),
mandatory, the  region dictionary  contains items:
- `"femm"` = finite element mmodel machine (mandatory);
- `"Q"` = material internal heat generation rate (optional; default  0.0)

For essential boundary conditions (optional) each dictionary
would hold
  + `"temperature"` = fixed (prescribed) temperature (scalar),  or
            a function with signature
                function T = f(x)
            If not given, zero temperatures assumed.
  + `"node_list"` = list of nodes on the boundary to which the condition applies
            (mandatory)

For convection boundary conditions (optional) each dictionary
may hold
  + `"femm"` = finite element mmodel machine (mandatory);
  + `"ambient_temperature"` = fixed (prescribed) ambient temperature (scalar)
        If not given, zero temperatures assumed.

For flux boundary conditions (optional) each dictionary
would hold
  + `"femm"` = finite element mmodel machine (mandatory);
  + `"normal_flux"` = normal component of the flux through the boundary (scalar)
        Positive  when outgoing.

# Output
`modeldata`= the dictionary on input is augmented with
- `"geom"` = the nodal field that is the geometry
- `"temp"` = the nodal field that is the computed temperature
"""
function steadystate(modeldata::FDataDict)
    # Lists of recognized keys for the data dictionaries:
    modeldata_recognized_keys = ["fens", "regions", "essential_bcs", "convection_bcs", "flux_bcs"]
    essential_bcs_recognized_keys = ["temperature", "node_list"]
    convection_bcs_recognized_keys = ["femm", "ambient_temperature"]
    flux_bcs_recognized_keys = ["femm", "normal_flux"]
    regions_recognized_keys = ["femm", "Q"]

    # Check that keys defined  in the dictionary are matched
    dcheck!(modeldata, modeldata_recognized_keys)

    # Extract the nodes
    fens = get(()->error("Must get finite element nodes (fens)!"), modeldata, "fens")

    # Construct the geometry field
    geom = NodalField(fens.xyz)

    # Construct the temperature field
    temp = NodalField(zeros(size(fens.xyz,1),1))

    # Apply the essential boundary conditions on the temperature field
    essential_bcs = get(modeldata, "essential_bcs", nothing);
    if (essential_bcs != nothing)
        for j = 1:length(essential_bcs)
            ebc = essential_bcs[j]
            dcheck!(ebc, essential_bcs_recognized_keys)
            fenids = get(()->error("Must get node list!"), ebc, "node_list");
            temperature = get(ebc, "temperature", nothing);
            T_fixed = zeros(FFlt,length(fenids)); # default is  zero temperature
            if (temperature != nothing) # if it is nonzero,
                if (typeof(temperature) <: Function) # it could be a function
                    for k = 1:length(fenids)
                        T_fixed[k] = temperature(geom.values[fenids[k],:])[1];
                    end
                else # or it could be a constant
                    fill!(T_fixed, temperature);
                end
            end
            setebc!(temp, fenids[:], true, 1, T_fixed);
            applyebc!(temp);
        end
    end

    # Number the equations
    numberdofs!(temp)           #,Renumbering_options); # NOT DONE <<<<<<<<<<<<<<<<<

    # Initialize the heat loads vector
    F = zeros(FFlt,temp.nfreedofs);

    # Construct the system conductivity matrix
    K = spzeros(temp.nfreedofs,temp.nfreedofs); # (all zeros, for the moment)
    regions = get(()->error("Must get regions!"), modeldata, "regions")
    for i = 1:length(regions)
        region = regions[i]
        dcheck!(region, regions_recognized_keys)
        femm = region["femm"];
        # Add up all the conductivity matrices for all the regions
        K = K + conductivity(femm, geom, temp);
        Q = get(region, "Q", [0.0]);
        if (typeof(Q) <: Function)
            fi = ForceIntensity(FFlt, 1, Q);
        else
            fi = ForceIntensity(Q);
        end
        F = F + distribloads(femm, geom, temp, fi, 3);
        # Loads due to the essential boundary conditions on the temperature field
        essential_bcs = get(modeldata, "essential_bcs", nothing);
        if (essential_bcs != nothing)
        F = F + nzebcloadsconductivity(femm, geom, temp);
        end
    end

    # Process the convection boundary condition
    convection_bcs = get(modeldata, "convection_bcs", nothing);
    if (convection_bcs != nothing)
        amb = deepcopy(temp); # create the ambient temperature field
        for i = 1:length(convection_bcs)
            convbc = convection_bcs[i]
            dcheck!(convbc, convection_bcs_recognized_keys)
            femm = get(()->error("Must get femm!"), convbc, "femm");
            # Apply the prescribed ambient temperature
            fenids = connectednodes(femm.integdomain.fes);
            fixed = ones(length(fenids));
            T_fixed = zeros(FFlt, length(fenids)); # default is zero
            ambient_temperature = get(convbc, "ambient_temperature", nothing);
            if ambient_temperature != nothing  # if given as nonzero
                if (typeof(ambient_temperature) <: Function) # given by function
                    for k = 1:length(fenids)
                        T_fixed[k] = ambient_temperature(geom.values[fenids[k],:])[1];
                    end
                else # it could be a constant
                    fill!(T_fixed, ambient_temperature);
                end
            end
            setebc!(amb, fenids[:], true, 1, T_fixed);
            applyebc!(amb);
            femm = convbc["femm"];
            K = K + surfacetransfer(femm, geom, temp);
            F = F + surfacetransferloads(femm, geom, temp, amb);
            # Note that EBC will contribute through the surface heat transfer matrix
            essential_bcs = get(modeldata, "essential_bcs", nothing);
            # If any essential boundary condition defined, the convection BC could
            # contribute a load term
            if (essential_bcs != nothing)
                F = F + nzebcsurfacetransferloads(femm, geom, temp);
            end
        end
    end

    # # Process the flux boundary condition
    flux_bcs = get(modeldata, "flux_bcs", nothing);
    if (flux_bcs != nothing)
        for j = 1:length(flux_bcs)
            fluxbc = flux_bcs[j]
            dcheck!(fluxbc, flux_bcs_recognized_keys)
            normal_flux = fluxbc["normal_flux"];
            if (typeof(normal_flux) <: Function)
                fi = ForceIntensity(FFlt, 1, normal_flux);
            else
                if typeof(normal_flux) <: AbstractArray
                else
                    normal_flux = FFlt[normal_flux]
                end
                fi = ForceIntensity(normal_flux);
            end
            femm = fluxbc["femm"]
            # Note the sign  which reflects the formula (negative sign
            # in front of the integral)
            F = F - distribloads(femm, geom, temp, fi, 2);
        end
    end

    # Solve for the temperatures
    K = cholesky(K);
    U = K\F;
    scattersysvec!(temp, U[:])

    # Update the model data
    setindex!(modeldata, geom, "geom");
    setindex!(modeldata, temp, "temp");
    return modeldata            # ... And return the updated model data
end

end
