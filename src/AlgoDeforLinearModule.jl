"""
    AlgoDeforLinearModule

Module for algorithms used in linear deformation models.
"""
module AlgoDeforLinearModule

using FinEtools.FTypesModule: FInt, FFlt, FCplxFlt, FFltVec, FIntVec, FFltMat, FIntMat, FMat, FVec, FDataDict
import FinEtools.AlgoBaseModule: dcheck!
import Arpack: eigs
import SparseArrays: spzeros
import LinearAlgebra: mul!
my_A_mul_B!(C, A, B) = mul!(C, A, B)
import FinEtools.FieldModule: AbstractField, ndofs, setebc!, numberdofs!, applyebc!, scattersysvec!
import FinEtools.NodalFieldModule: NodalField, nnodes
import FinEtools.FEMMBaseModule: associategeometry!, distribloads, fieldfromintegpoints, elemfieldfromintegpoints
import FinEtools.FEMMDeforLinearBaseModule: stiffness, mass, nzebcloadsstiffness, thermalstrainloads, inspectintegpoints
import FinEtools.FEMMDeforLinearMSModule: stiffness, mass, nzebcloadsstiffness, thermalstrainloads, inspectintegpoints
import FinEtools.DeforModelRedModule: stresscomponentmap
import FinEtools.ForceIntensityModule: ForceIntensity
import FinEtools.MeshModificationModule: meshboundary
import FinEtools.MeshExportModule: vtkexportmesh
import LinearAlgebra: eigen, qr, dot, cholesky

"""
    AlgoDeforLinearModule.linearstatics(modeldata::FDataDict)

Algorithm for static linear deformation (stress) analysis.

`modeldata` = dictionary with values for keys

* "fens"  = finite element node set
* "regions"  = array of region dictionaries
* "essential_bcs" = array of essential boundary condition dictionaries
* "traction_bcs" = array of traction boundary condition dictionaries
* "temperature_change" = dictionary of data for temperature change

For each region (connected piece of the domain made of a particular material),
mandatory, the  region dictionary  contains values for keys:
* "femm" = finite element model machine (mandatory);

For essential boundary conditions (optional) each dictionary
would hold
  + "displacement" = fixed (prescribed) displacement (scalar),  or
            a function with signature
                function w = f(x)
            If not given, zero displacement assumed.
  + "component" = which component is prescribed  (1, 2, 3)?
  + "node_list" = list of nodes on the boundary to which the condition applies
            (mandatory)

For traction boundary conditions (optional) each dictionary
would hold
  + "femm" = finite element model machine (mandatory);
  + "traction_vector" = traction vector,  either  a constant or  a function
        Positive  when outgoing.

Output:
modeldata= the dictionary on input is augmented with
* "geom" = the nodal field that is the geometry
* "u" = the nodal field that is the computed displacement
* "dT" = the nodal field that is the temperature increment
"""
function linearstatics(modeldata::FDataDict)

    # For traction boundary conditions (optional):
    # model_data.boundary_conditions.traction = cell array of struct,
    #           each piece of surface with traction boundary condition gets one
    #           element of the array with a struct with the attributes
    #     traction=traction (vector), supply a zero for component in which
    #           the boundary condition is inactive
    #     fes = finite element set on the boundary to which
    #                       the condition applies
    #     integration_rule= integration rule
    #
    # For body loads (optional):
    # model_data.body_load = cell array of struct,
    #          each piece of the domain can have each its own body load
    #     force  = force density vector
    #     fes = finite element set to which the load applies
    #     integration_rule= integration rule
    #
    # For multi point constraints (MPC) (optional):
    # model_data.mpc= cell array of structs, each for one MPC.
    #      node_list = list of node numbers involved in the MPC,
    #      dof_list= numbers of degrees of freedom for the nodes above,
    #      umultipliers=multipliers for the nodes above,
    #      penfact=the penalty factor to multiply  the constraint matrix,
    #          The MPC looks like this: sum_i m_i u_{dof(i),node(i)} =0
    #          where m_i is the multiplier.

    # Lists of recognized keys for the data dictionaries:
    modeldata_recognized_keys = ["fens", "regions",
    "essential_bcs",  "traction_bcs", "temperature_change",
    "factorize"]
    essential_bcs_recognized_keys = ["displacement", "node_list", "component"]
    traction_bcs_recognized_keys = ["femm", "traction_vector"]
    regions_recognized_keys = ["femm", "body_load"]
    temperature_change_recognized_keys = ["temperature"]

    # Extract the nodes
    fens=get(()->error("Must get fens!"), modeldata, "fens")

    # Construct the geometry field
    geom = NodalField(fens.xyz)

    # Construct the displacement field
    u = NodalField(zeros(nnodes(geom), ndofs(geom)))

    # Construct the temperature field
    temp = NodalField(zeros(nnodes(geom), 1))

    modeldata["timing"] = FDataDict()

    tstart = time()
    # Apply the essential boundary conditions on the displacement field
    essential_bcs = get(modeldata, "essential_bcs", nothing);
    if (essential_bcs != nothing)
        for j = 1:length(essential_bcs)
            ebc = essential_bcs[j]
            dcheck!(ebc, essential_bcs_recognized_keys)
            fenids = get(()->error("Must get node list!"), ebc, "node_list");
            displacement = get(ebc, "displacement", nothing);
            u_fixed = zeros(FFlt, length(fenids)); # default is  zero displacement
            if (displacement != nothing) # if it is nonzero,
                if (typeof(displacement) <: Function) # it could be a function
                    for k = 1:length(fenids)
                        u_fixed[k] = displacement(geom.values[fenids[k],:])[1];
                    end
                else # or it could be a constant
                    fill!(u_fixed, displacement);
                end
            end
            component = get(ebc, "component", 0); # which component?
            setebc!(u, fenids[:], true, component, u_fixed);
        end
        applyebc!(u);
    end

    # Number the equations
    numberdofs!(u)           #,Renumbering_options); # NOT DONE <<<<<<<<<<<<<<<<<
    modeldata["timing"]["essential_bcs"] = time() - tstart

    tstart = time()
    # Initialize the heat loads vector
    F = zeros(FFlt,u.nfreedofs);
    # Construct the system stiffness matrix
    K = spzeros(u.nfreedofs,u.nfreedofs); # (all zeros, for the moment)
    regions = get(()->error("Must get region list!"), modeldata, "regions")
    for i = 1:length(regions)
        region = regions[i]
        dcheck!(region, regions_recognized_keys)
        femm = region["femm"];
        # Give the  FEMM a chance  to precompute  geometry-related quantities
        femm = associategeometry!(femm, geom);
        # Add up all the conductivity matrices for all the regions
        K = K + stiffness(femm, geom, u);
        # Loads due to the essential boundary conditions on the displacement field
        essential_bcs = get(modeldata, "essential_bcs", nothing);
        if (essential_bcs != nothing) # there was at least one EBC applied
            F = F + nzebcloadsstiffness(femm, geom, u);
        end
    end
    modeldata["timing"]["stiffness"] = time() - tstart

    # # Process the body load
    # body_load = get(modeldata, "body_load", nothing);
    # if (body_load  !=nothing)
    #     for j=1:length(model_data.body_load)
    #         body_load =model_data.body_load{j};
    #         femm = femm_deformation_linear (struct ('material',[],...
    #             'fes',body_load.fes,...
    #             'integration_rule',body_load.integration_rule));
    #         fi= force_intensity(struct('magn',body_load.force));
    #         F = F + distrib_loads(femm, sysvec_assembler, geom, u, fi, 3);
    #     end
    #     clear body_load fi  femm
    # end

    tstart = time()
    # Process the traction boundary condition
    traction_bcs = get(modeldata, "traction_bcs", nothing);
    if (traction_bcs != nothing)
        for j=1:length(traction_bcs)
            tractionbc = traction_bcs[j]
            dcheck!(tractionbc, traction_bcs_recognized_keys)
            traction_vector = tractionbc["traction_vector"];
            if (typeof(traction_vector) <: Function)
                fi = ForceIntensity(FFlt, ndofs(geom), traction_vector);
            else
                fi = ForceIntensity(traction_vector);
            end
            femm = tractionbc["femm"]
            F = F + distribloads(femm, geom, u, fi, 2);
        end
    end
    modeldata["timing"]["traction_bcs"] = time() - tstart

    tstart = time()
    # Process the thermal strain  loading
    temperature_change = get(modeldata, "temperature_change", nothing);
    if (temperature_change != nothing)
        dcheck!(temperature_change, temperature_change_recognized_keys)
        # zero temperature change is a reasonable default
        temp = NodalField(zeros(size(fens.xyz,1),1))
        temperature = get(temperature_change, "temperature", nothing);
        if (temperature != nothing) # if it is nonzero,
            if (typeof(temperature) <: Function) # it could be a function
                for k = 1:count(fens)
                    temp.values[k] = temperature(geom.values[k,:])[1];
                end
            else # or it could be a constant
                fill!(temp.values, temperature);
            end
        end
        for i = 1:length(regions)
            region = regions[i]
            femm = region["femm"];
            F = F + thermalstrainloads(femm, geom, u, temp)
        end
    end
    modeldata["timing"]["temperature_change"] = time() - tstart

    # # Process the nodal force boundary condition
    # if (isfield(model_data.boundary_conditions, 'nodal_force' ))
    #     for j=1:length(model_data.boundary_conditions.nodal_force)
    #         nodal_force =model_data.boundary_conditions.nodal_force{j};
    #         femm = femm_deformation_linear (struct ('material',[],...
    #             'fes',fe_set_P1(struct('conn',reshape(nodal_force.node_list,[],1))),...
    #             'integration_rule',point_rule));
    #         fi= force_intensity(struct('magn',nodal_force.force));
    #         F = F + distrib_loads(femm, sysvec_assembler, geom, u, fi, 0);
    #     end
    #     clear nodal_force fi femm
    # end

    # # Apply multi point constraints
    # if isfield(model_data,'mpc')
    #     for i=1:length(model_data.mpc)
    #         mpc =model_data.mpc{i};
    #         dofnums=0*mpc.umultipliers;# Construct an array of the degree of freedom numbers
    #         for kx=1:length(mpc.node_list)
    #             dofnums(kx)=u.dofnums(mpc.node_list(kx),mpc.dof_list(kx));
    #         end
    #         # Now call the utility function to calculate the constraint matrix
    #         [Kmpc,Fmpc]=apply_penalty_mpc(u.nfreedofs,dofnums,mpc.umultipliers,0.0,mpc.penfact);
    #         K = K + Kmpc;
    #         F = F + Fmpc;
    #     end
    #     clear Kmpc Fmpc
    # end

    tstart = time()
    # Solve the system of linear algebraic equations
    K = cholesky(K);
    U = K\F;
    scattersysvec!(u, U[:])
    modeldata["timing"]["solution"] = time() - tstart


    # Update the model data
    setindex!(modeldata, geom, "geom");
    setindex!(modeldata, u, "u");
    setindex!(modeldata, temp, "temp");
    setindex!(modeldata, temp, "dT");
    setindex!(modeldata, dot(F,U)/2, "work");
    return modeldata            # ... And return the updated model data
end

"""
    AlgoDeforLinearModule.exportdeformation(modeldata::FDataDict)

Algorithm for exporting of the deformation for visualization in Paraview.

`modeldata` = dictionary with values for keys

* "fens"  = finite element node set
* "regions"  = array of region dictionaries
* "geom" = geometry field
* "u" = displacement field, or
* "us" = array of  tuples (name, displacement field)
* "postprocessing" = dictionary  with values for keys
  + "boundary_only" = should only the boundary of the  regions be rendered?
                      Default is render the interior.
  + "file" = name of the  postprocessing file

For each region (connected piece of the domain made of a particular material),
mandatory, the  region dictionary  contains values for keys:
* "femm" = finite element mmodel machine (mandatory);

Output: modeldata updated with
* modeldata["postprocessing"]["exported"] = array of data dictionaries, one for
        each exported file. The data is stored with the keys:
    - "file" - names of exported file
    - "field" - nodal or elemental field
"""
function exportdeformation(modeldata::FDataDict)
    modeldata_recognized_keys = ["fens", "regions", "geom", "u", "postprocessing"]
    postprocessing_recognized_keys = ["boundary_only", "file"]
    # Defaults
    boundary_only = false;
    ffile = "deformation"
    dcheck!(modeldata, modeldata_recognized_keys)

    # Let's have a look at what's been specified
    postprocessing = get(modeldata, "postprocessing", nothing);
    if (postprocessing != nothing)
        dcheck!(postprocessing, postprocessing_recognized_keys)
        boundary_only =  get(postprocessing, "boundary_only", boundary_only);
        ffile =  get(postprocessing, "file", ffile);
    end

    fens = get(()->error("Must get fens!"), modeldata, "fens")
    geom = get(()->error("Must get geometry field!"), modeldata, "geom");
    u = get(modeldata, "u", nothing);
    us = get(modeldata, "us", nothing);
    if us == nothing
        us = [("u", u)]
    end

    # Export one file for each region
    modeldata["postprocessing"]["exported"] = Array{FDataDict, 1}()
    regions = get(()->error("Must get region!"), modeldata, "regions")
    for i = 1:length(regions)
        region = regions[i]
        femm = region["femm"]
        rfile = ffile * "$i" * ".vtk";
        vectors = Tuple{String, FFltMat}[]
        for ixxxx = 1:length(us)
            push!(vectors, (us[ixxxx][1], us[ixxxx][2].values))
        end
        if boundary_only
            bfes = meshboundary(femm.integdomain.fes);
            vtkexportmesh(rfile, fens, bfes;  vectors=vectors)
        else
            vtkexportmesh(rfile, fens, femm.integdomain.fes; vectors=vectors)
        end
        ed = FDataDict("file"=>rfile, "field"=>u, "region"=>i,
            "type"=>"displacement")
        push!(modeldata["postprocessing"]["exported"], ed)
    end

    return modeldata
end

"""
    AlgoDeforLinearModule.exportstress(modeldata::FDataDict)

Algorithm for exporting of the stress for visualization in Paraview.

`modeldata` = dictionary with values for keys

* "fens"  = finite element node set
* "regions"  = array of region dictionaries
* "geom" = geometry field
* "u" = displacement field
* "postprocessing" = dictionary  with values for keys
  + "boundary_only" = should only the boundary of the  regions be rendered?
                      Default is render the interior.
  + "file" = name of the  postprocessing file
  + "quantity" = quantity to be exported (default :Cauchy)
  + "component" = which component of the quantity?
  + "outputcsys" = output coordinate system
  + "inspectormeth" = inspector method to pass to `inspectintegpoints()`
  + "extrap" = method for extrapolating from the quadrature points to the nodes
    within one element

For each region (connected piece of the domain made of a particular material),
mandatory, the  region dictionary  contains values for keys:
* "femm" = finite element mmodel machine (mandatory);

Output: modeldata updated with
* modeldata["postprocessing"]["exported"] = array of data dictionaries, one for
        each exported file. The data is stored with the keys:
    - "file" - name of exported file
    - "field" - nodal field
"""
function exportstress(modeldata::FDataDict)
    modeldata_recognized_keys = ["fens", "regions", "geom", "u",
    "dT", "postprocessing"]
    postprocessing_recognized_keys = ["boundary_only", "file", "quantity",
    "component", "outputcsys", "nodevalmethod", "reportat"]
    # Defaults
    boundary_only = false;
    ffile = "stress"
    dcheck!(modeldata, modeldata_recognized_keys)
    quantity = :Cauchy
    component = 1
    outputcsys = nothing
    reportat = :default
    nodevalmethod = :invdistance
    # Let's have a look at what's been specified
    postprocessing = get(modeldata, "postprocessing", nothing);
    if (postprocessing != nothing)
        dcheck!(postprocessing, postprocessing_recognized_keys)
        boundary_only = get(postprocessing, "boundary_only", boundary_only);
        ffile = get(postprocessing, "file", ffile);
        quantity = get(postprocessing, "quantity", quantity);
        component = get(postprocessing, "component", component);
        outputcsys = get(postprocessing, "outputcsys", outputcsys);
        nodevalmethod = get(postprocessing, "nodevalmethod", nodevalmethod);
        reportat = get(postprocessing, "reportat", reportat);
    end

    fens = get(()->error("Must get fens!"), modeldata, "fens")
    geom = get(()->error("Must get geometry field!"), modeldata, "geom");
    u = get(()->error("Must get displacement field!"), modeldata, "u");
    dT = get(modeldata, "dT", nothing);

    context = []
    if (outputcsys != nothing)
        push!(context, (:outputcsys, outputcsys))
    end
    if (nodevalmethod != nothing)
        push!(context, (:nodevalmethod, nodevalmethod))
    end
    if (reportat != nothing)
        push!(context, (:reportat, reportat))
    end

    # Export a file for each region
    modeldata["postprocessing"]["exported"] = Array{FDataDict, 1}()
    regions = get(()->error("Must get region!"), modeldata, "regions")
    for i = 1:length(regions)
        region = regions[i]
        femm = region["femm"]
        rfile = ffile * "-" * string(quantity) * string(component) * "-region $i" * ".vtk";
        if (typeof(component) == Symbol)
            componentnum = stresscomponentmap(femm.mr)[component]
        else
            componentnum = component
        end
        componentname = length(componentnum) > 1 ? "" : "$(componentnum)"
        # Note that we are creating a field  separately for each region.  This is
        # important  for the following reason: if the regions were of different
        # materials, or if they were of the same material but with different material
        # axes orientation, averaging across the material interface  would not make
        # sense.
        if (dT != nothing)
            fld = fieldfromintegpoints(femm, geom, u, dT, quantity, componentnum;
            context...)
        else
            fld = fieldfromintegpoints(femm, geom, u, quantity, componentnum;
            context...)
        end
        if boundary_only
            bfes = meshboundary(femm.integdomain.fes);
            vtkexportmesh(rfile, fens, bfes;
            scalars=[(string(quantity)*componentname, fld.values)],
            vectors=[("u", u.values)])
        else
            vtkexportmesh(rfile, fens, femm.integdomain.fes;
            scalars=[(string(quantity)*componentname, fld.values)],
            vectors=[("u", u.values)])
        end
        ed = FDataDict("file"=>rfile, "field"=>fld, "region"=>i,
            "type"=>"nodal stress",
            "quantity"=>quantity, "component"=>component, "outputcsys"=>outputcsys,
            "nodevalmethod"=>nodevalmethod, "reportat"=>reportat)
        push!(modeldata["postprocessing"]["exported"], ed)
    end

    return modeldata
end

"""
    AlgoDeforLinearModule.exportstresselementwise(modeldata::FDataDict)


Algorithm for exporting of the elementwise stress for visualization in Paraview.

`modeldata` = dictionary with values for keys

* "fens"  = finite element node set
* "regions"  = array of region dictionaries
* "geom" = geometry field
* "u" = displacement field
* "postprocessing" = dictionary  with values for keys
  + "boundary_only" = should only the boundary of the  regions be rendered?
                      Default is render the interior.
  + "file" = name of the  postprocessing file
  + "quantity" = quantity to be exported (default :Cauchy)
  + "component" = which component of the quantity?
  + "outputcsys" = output coordinate system

For each region (connected piece of the domain made of a particular material),
mandatory, the  region dictionary  contains values for keys:
* "femm" = finite element mmodel machine (mandatory);

Output: modeldata updated with
* modeldata["postprocessing"]["exported"] = array of data dictionaries, one for
        each exported file. The data is stored with the keys:
    - "file" - name of exported file
    - "field" - elemental field
"""
function exportstresselementwise(modeldata::FDataDict)
    modeldata_recognized_keys = ["fens", "regions", "geom", "u",
    "dT", "postprocessing"]
    postprocessing_recognized_keys = ["boundary_only", "file", "quantity",
    "component", "outputcsys" ]
    # Defaults
    boundary_only = false;
    ffile = "stress"
    dcheck!(modeldata, modeldata_recognized_keys)
    quantity = :Cauchy
    component = 1
    outputcsys = nothing
    # Let's have a look at what's been specified
    postprocessing = get(modeldata, "postprocessing", nothing);
    if (postprocessing != nothing)
        dcheck!(postprocessing, postprocessing_recognized_keys)
        boundary_only = get(postprocessing, "boundary_only", boundary_only);
        ffile = get(postprocessing, "file", ffile);
        quantity = get(postprocessing, "quantity", quantity);
        component = get(postprocessing, "component", component);
        outputcsys = get(postprocessing, "outputcsys", outputcsys);
    end

    fens = get(()->error("Must get fens!"), modeldata, "fens")
    geom = get(()->error("Must get geometry field!"), modeldata, "geom");
    u = get(()->error("Must get displacement field!"), modeldata, "u");
    dT = get(modeldata, "dT", nothing);

    context = []
    if (outputcsys != nothing)
        push!(context, (:outputcsys, outputcsys))
    end

    # Export a file for each region
    modeldata["postprocessing"]["exported"] = Array{FDataDict, 1}()
    regions = get(()->error("Must get region!"), modeldata, "regions")
    for i = 1:length(regions)
        region = regions[i]
        femm = region["femm"]
        rfile = ffile * "-" * string(quantity) * string(component) * "-region $i" * ".vtk";
        if (typeof(component) == Symbol)
            componentnum = stresscomponentmap(femm.mr)[component]
        else
            componentnum = component
        end
        componentname = length(componentnum) > 1 ? "" : "$(componentnum)"
        # Note that we are creating a field  separately for each region.  This is
        # important  for the following reason: if the regions were of different
        # materials, or if they were of the same material but with different material
        # axes orientation, averaging across the material interface  would not make
        # sense.
        if (dT != nothing)
            fld = elemfieldfromintegpoints(femm, geom, u, dT, quantity, componentnum;
            context...)
        else
            fld = elemfieldfromintegpoints(femm, geom, u, quantity, componentnum;
            context...)
        end
        if boundary_only
            bfes = meshboundary(femm.integdomain.fes);
            vtkexportmesh(rfile, fens, bfes;
            scalars=[(string(quantity)*componentname, fld.values)],
            vectors=[("u", u.values)])
        else
            vtkexportmesh(rfile, fens, femm.integdomain.fes;
            scalars=[(string(quantity)*componentname, fld.values)],
            vectors=[("u", u.values)])
        end
        ed = FDataDict("file"=>rfile, "field"=>fld, "region"=>i,
            "type"=>"elemental stress",
            "quantity"=>quantity, "component"=>component, "outputcsys"=>outputcsys)
        push!(modeldata["postprocessing"]["exported"], ed)
    end

    return modeldata
end

"""
    AlgoDeforLinearModule.modal(modeldata::FDataDict)

Modal (free-vibration) analysis solver.


`modeldata` = dictionary with values for keys

* "fens"  = finite element node set
* "regions"  = array of region dictionaries
* "essential_bcs" = array of essential boundary condition dictionaries

For each region (connected piece of the domain made of a particular material),
mandatory, the  region dictionary  contains values for keys:
* "femm" = finite element mmodel machine (mandatory);

For essential boundary conditions (optional) each dictionary
would hold
  + "displacement" = fixed (prescribed) displacement (scalar): only zero
        displacement is  allowed for modal analysis.
  + "component" = which component is prescribed  (1, 2, 3)?
  + "node_list" = list of nodes on the boundary to which the condition applies
            (mandatory)

Control parameters:
* "neigvs" = number of eigenvalues/eigenvectors to compute
* "omega_shift"= angular frequency shift for mass shifting
* "use_lumped_mass" = true or false?  (Default is false: consistent mass)


Output:
modeldata= the dictionary on input is augmented with
* "geom" = the nodal field that is the geometry
* "u" = the nodal field that is the computed displacement
* "neigvs" = Number of computed eigenvectors
* "W" = Computed eigenvectors, neigvs columns
* "omega" =  Computed angular frequencies, array of length neigvs
* "raw_eigenvalues" = Raw computed eigenvalues
"""
function modal(modeldata::FDataDict)

    # For multi point constraints (MPC) (optional):
    # model_data.mpc= cell array of structs, each for one MPC.
    #      mpc.node_list = list of node numbers involved in the MPC,
    #      mpc.dof_list= numbers of degrees of freedom for the nodes above,
    #      mpc.umultipliers=multipliers for the nodes above,
    #      mpc.penfact=the penalty factor to multiply  the constraint matrix,
    #          The MPC looks like this: sum_i m_i u_{dof(i),node(i)} =0
    #          where m_i is the multiplier.


    # Lists of recognized keys for the data dictionaries:
    modeldata_recognized_keys = ["fens", "regions",
    "essential_bcs", "neigvs", "omega_shift", "use_lumped_mass"]
    essential_bcs_recognized_keys = ["displacement", "node_list", "component"]
    regions_recognized_keys = ["femm", "femm_stiffness", "femm_mass", "body_load"]

    neigvs = get(modeldata, "neigvs", 7); # Number of eigenvalues

    omega_shift = get(modeldata, "omega_shift", 0.0); # Mass shifting

    use_factorization = get(modeldata, "use_factorization", false); # Factorization?

    use_lumped_mass = get(modeldata, "use_lumped_mass", false); # Lumped mass?

    # Extract the nodes
    fens = get(()->error("Must get fens!"), modeldata, "fens")

    # Construct the geometry field
    geom = NodalField(fens.xyz)

    # Construct the displacement field
    u = NodalField(zeros(nnodes(geom),ndofs(geom)))

    # Apply the essential boundary conditions on the displacement field
    essential_bcs = get(modeldata, "essential_bcs", nothing);
    if (essential_bcs != nothing)
        for j = 1:length(essential_bcs)
            ebc = essential_bcs[j]
            dcheck!(ebc, essential_bcs_recognized_keys)
            fenids = get(()->error("Must get node list!"), ebc, "node_list");
            displacement = get(ebc, "displacement", nothing);
            u_fixed = zeros(FFlt, length(fenids)); # only zero displacement accepted
            component = get(ebc, "component", 0); # which component?
            setebc!(u, fenids[:], true, component, u_fixed);
        end
        applyebc!(u);
    end

    # Number the equations
    numberdofs!(u)           #,Renumbering_options); # NOT DONE <<<<<<<<<<<<<<<<<

    # Construct the system stiffness matrix
    K = spzeros(u.nfreedofs,u.nfreedofs); # (all zeros, for the moment)
    regions = get(()->error("Must get region list!"), modeldata, "regions")
    for i = 1:length(regions)
        region = regions[i]
        dcheck!(region, regions_recognized_keys)
        if "femm_stiffness"  in keys(region)
            femm = region["femm_stiffness"];
        else
            femm = get(()->error("Must get femm or femm_stiffness!"), region, "femm")
        end
        # Give the  FEMM a chance  to precompute  geometry-related quantities
        femm = associategeometry!(femm, geom);
        # Add up all the stiffness matrices for all the regions
        K = K + stiffness(femm, geom, u);
    end

    # Construct the system mass matrix
    M = spzeros(u.nfreedofs,u.nfreedofs); # (all zeros, for the moment)
    regions = get(()->error("Must get region list!"), modeldata, "regions")
    for i = 1:length(regions)
        region = regions[i]
        dcheck!(region, regions_recognized_keys)
        if "femm_mass"  in keys(region)
            femm = region["femm_mass"];
        else
            femm = get(()->error("Must get femm or femm_mass!"), region, "femm")
        end
        # Give the  FEMM a chance  to precompute  geometry-related quantities
        femm = associategeometry!(femm, geom);
        # Add up all the mass matrices for all the regions
        M = M + mass(femm, geom, u);
    end

    # Options for the eigenproblem solution

    # Solve
    # if (~ use_factorization )
    #     # This is one way of solving the eigenvalue problem, just pass the matrices
    #     [W,Omega]= eigs(K+omega_shift*M, M, neigvs, 'SM', evopts);
    # else
    # This form uses the factorized matrix and has the potential of being much faster
    # Factorize the left-hand side matrix for efficiency (Choleski)
    # [mA,status] = chol(K+omega_shift*M,'lower');#,'vector',prm
    # if ( status ~= 0 ) error('Choleski factorization failed'), end
    # clear K; # Not needed anymore
    # mAt= mA';
    # [W,Omega]= eigs(@(bv)mAt\(mA\bv), u.nfreedofs, M, neigvs, 'SM', evopts);
    #          [W,Omega]= eigen(full(K+omega_shift*M), full(M));

    d,v,nev,nconv = eigs(K+omega_shift*M, M; nev=neigvs, which=:SM)
    broadcast!(+, d, d, -omega_shift);

    modeldata["raw_eigenvalues"] = d;
    #    Subtract the mass-shifting Angular frequency
    if any(imag(d) .!= 0.0)
        d=real(d);
    end
    if any(real(d) .< 0.0)
        d = abs.(d);
    end
    #    Sort  the angular frequencies by magnitude.  Make sure all
    #    imaginary parts of the eigenvalues are removed.
    ix =sortperm(d);

    # Update the model data: store geometry
    modeldata["geom"] = geom;
    # Store the displacement field
    modeldata["u"] = u;
    # Number of computed eigenvectors
    modeldata["neigvs"] = nev;
    #  Computed eigenvectors: we are ignoring the imaginary part here
    #  because the modal analysis is presumed to have been performed for
    #  an undamped structure
    modeldata["W"] = real(v[:,ix]);
    #  Computed angular frequencies
    modeldata["omega"] = sqrt.(d[ix]);
    return modeldata
end

"""
    AlgoDeforLinearModule.exportmode(modeldata::FDataDict)

Algorithm for exporting of the mmode shape for visualization in Paraview.

`modeldata` = dictionary with values for keys

* "fens"  = finite element node set
* "regions"  = array of region dictionaries
* "geom" = geometry field
* "u" = displacement field
* "W" = Computed free-vibration eigenvectors, neigvs columns
* "omega" =  Computed free-vibration angular frequencies, array of length neigvs
* "postprocessing" = dictionary  with values for keys
  + "boundary_only" = should only the boundary of the  regions be rendered?
                      Default is render the interior.
  + "file" = name of the  postprocessing file
  + "mode" = which mode should be visualized?
  + "component" = which component of the quantity?
  + "outputcsys" = output coordinate system

For each region (connected piece of the domain made of a particular material),
mandatory, the  region dictionary  contains values for keys:
* "femm" = finite element mmodel machine (mandatory);

Output: modeldata updated with
* modeldata["postprocessing"]["exported"] = see `exportdeformation()`
"""
function exportmode(modeldata::FDataDict)
    modeldata_recognized_keys = ["fens", "regions", "geom", "u",
    "omega", "W",
    "postprocessing"]
    postprocessing_recognized_keys = ["boundary_only", "file", "mode"]
    mode = 1;
    dcheck!(modeldata, modeldata_recognized_keys)

    # Let's have a look at what's been specified
    postprocessing = get(modeldata, "postprocessing", nothing);
    if (postprocessing != nothing)
        dcheck!(postprocessing, postprocessing_recognized_keys)
        mode =  get(postprocessing, "mode", mode);
    end

    omega = modeldata["omega"]

    # Scatter the desired mode
    W = modeldata["W"]
    if typeof(mode)<:Int
        @assert 0 < mode <= length(omega) "Invalid mode number $mode"
        scattersysvec!(modeldata["u"], W[:,mode])
    else
        us = Tuple{String, AbstractField}[]
        u = modeldata["u"]
        for ixxxx in mode
            @assert 0 < ixxxx <= length(omega) "Invalid mode number $ixxxx"
            scattersysvec!(u, W[:,ixxxx])
            push!(us, ("mode_$(ixxxx)", deepcopy(u)))
        end
        modeldata["us"] = us
    end

    return exportdeformation(modeldata)
end

"""
    ssit(K, M; nev::Int=6, evshift::FFlt = 0.0,
        v0::FFltMat = Array{FFlt}(0, 0),
        tol::FFlt = 1.0e-3, maxiter::Int = 300, verbose::Bool=false)

Subspace  Iteration (block inverse power) method.

Block inverse power method for k smallest eigenvalues of the generalized
eigenvalue problem
           `K*v= lambda*M*v`

# Arguments
* `K` =  square symmetric stiffness matrix (if necessary mass-shifted),
* `M` =  square symmetric mass matrix,

# Keyword arguments
* `v0` =  initial guess of the eigenvectors (for instance random),
* `nev` = the number of eigenvalues sought
* `tol` = relative tolerance on the eigenvalue, expressed in terms of norms of the
      change of the eigenvalue estimates from iteration to iteration.
* `maxiter` =  maximum number of allowed iterations
* `withrr` = with Rayleigh-Ritz problem solved to improve the subspace?  (default
    is false)
* `verbose` = verbose? (default is false)

#  Return
* `labm` = computed eigenvalues,
* `v` = computed eigenvectors,
* `nconv` = number of converged eigenvalues
* `niter` = number of iterations taken
* `nmult` = ignore this output
* `lamberr` = eigenvalue errors, defined as  normalized  differences  of
    successive  estimates of the eigenvalues
"""
function ssit(K, M; nev::Int=6, evshift::FFlt = 0.0,
    v0::FFltMat = Array{FFlt}(0, 0),
    tol::FFlt = 1.0e-3, maxiter::Int = 300, withrr::Bool=false,
    verbose::Bool=false)
    @assert nev >= 1
    v = deepcopy(v0)
    if isempty(v0)
        v = rand(size(K, 1), nev)
    end
    @assert nev <= size(v, 2)
    nvecs = size(v, 2)  # How many eigenvalues are iterated?
    plamb = zeros(nvecs)  # previous eigenvalue
    lamb = zeros(nvecs)
    lamberr = zeros(nvecs)
    converged = falses(nvecs)  # not yet
    niter = 0
    nconv = 0
    nmult = 0
    factor = cholesky(K+evshift*M)
    Kv = zeros(size(K, 1), size(v, 2))
    Mv = zeros(size(M, 1), size(v, 2))
    for i = 1:maxiter
        u = factor\(M*v)
        factorization = qr(u)  # ; full=falseeconomy factorization
        v = Array(factorization.Q)
        my_A_mul_B!(Kv, K, v)
        my_A_mul_B!(Mv, M, v)
        for j = 1:nvecs
            lamb[j] = dot(v[:, j], Kv[:, j]) / dot(v[:, j], Mv[:, j])
            lamberr[j] = abs(lamb[j] - plamb[j])/abs(lamb[j])
            converged[j] = lamberr[j] <= tol
        end
        nconv = length(findall(converged[1:nev]))
        verbose && println("nconv = $(nconv)")
        if nconv >= nev # converged on all requested eigenvalues
            break
        end
        if withrr
            decomp = eigen(transpose(v)*Kv, transpose(v)*Mv)
            ix = sortperm(abs.(decomp.values))
            rrd = decomp.values[ix]
            rrv = decomp.vectors[:, ix]
            v = v*rrv
        end
        plamb, lamb = lamb, plamb # swap the eigenvalue arrays
        niter = niter + 1
    end
    return lamb, v, nconv, niter, nmult, lamberr
end
# (d,[v,],nconv,niter,nmult,resid)
# eigs returns the nev requested eigenvalues in d, the corresponding Ritz vectors
# v (only if ritzvec=true), the number of converged eigenvalues nconv, the number
# of iterations niter and the number of matrix vector multiplications nmult, as
# well as the final residual vector resid.

end
