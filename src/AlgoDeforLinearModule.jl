module AlgoDeforLinearModule

using FinEtools
using FinEtools.AlgoBaseModule.dcheck!

"""
    AlgoDeforLinearModule.linearstatics(modeldata::FDataDict)

Algorithm for static linear deformation (stress) analysis.

`modeldata` = dictionary with values for keys

- "fens"  = finite element node set
- "regions"  = array of region dictionaries
- "essential_bcs" = array of essential boundary condition dictionaries
- "traction_bcs" = array of traction boundary condition dictionaries
- "temperature_change" = dictionary of data for temperature change

For each region (connected piece of the domain made of a particular material),
mandatory, the  region dictionary  contains values for keys:
- "femm" = finite element mmodel machine (mandatory);

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
  + "femm" = finite element mmodel machine (mandatory);
  + "traction_vector" = traction vector,  either  a constant or  a function
        Positive  when outgoing.

Output:
modeldata= the dictionary on input is augmented with
- "geom" = the nodal field that is the geometry
- "u" = the nodal field that is the computed displacement
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

  # Initialize the heat loads vector
  F = zeros(FFlt,u.nfreedofs);

  # Construct the system stiffness matrix
  K = spzeros(u.nfreedofs,u.nfreedofs); # (all zeros, for the moment)
  regions = get(()->error("Must get region list!"), modeldata, "regions")
  for i = 1:length(regions)
    region = regions[i]
    dcheck!(region, regions_recognized_keys)
    femm = region["femm"];
    # # Give the  FEMM a chance  to precompute  geometry-related quantities
    # region.femm = associate_geometry(region.femm,geom);
    # Add up all the conductivity matrices for all the regions
    K = K + stiffness(femm, geom, u);
    # Loads due to the essential boundary conditions on the displacement field
    essential_bcs = get(modeldata, "essential_bcs", nothing);
    if (essential_bcs != nothing) # there was at least one EBC applied
      F = F + nzebcloadsstiffness(femm, geom, u);
    end
  end

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

  # Process the traction boundary condition
  traction_bcs = get(modeldata, "traction_bcs", nothing);
  if (traction_bcs != nothing)
    for j=1:length(traction_bcs)
      tractionbc = traction_bcs[j]
      dcheck!(tractionbc, traction_bcs_recognized_keys)
      traction_vector = tractionbc["traction_vector"];
      if (typeof(traction_vector) <: Function)
        fi = ForceIntensity(FFlt, ndofn(geom), traction_vector);
      else
        fi = ForceIntensity(traction_vector);
      end
      femm = tractionbc["femm"]
      F = F + distribloads(femm, geom, u, fi, 2);
    end
  end

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

  # Solve the system of linear algebraic equations
  K = cholfact(K);
  U = K\F;
  scattersysvec!(u, U[:])


  # Update the model data
  setindex!(modeldata, geom, "geom");
  setindex!(modeldata, u, "u");
  setindex!(modeldata, dot(F,U)/2, "work");
  return modeldata            # ... And return the updated model data
end

"""
    AlgoDeforLinearModule.exportdeformation(modeldata::FDataDict)

Export the deformation for visualization in Paraview.
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
  u = get(()->error("Must get displacement field!"), modeldata, "u");

  # Plot the surface for each region
  regions = get(()->error("Must get region!"), modeldata, "regions")
  for i = 1:length(regions)
    region = regions[i]
    femm = region["femm"]
    rfile = ffile * "$i" * ".vtk";
    if boundary_only
      bfes = meshboundary(femm.geod.fes);
      vtkexportmesh(rfile, fens, bfes;  vectors=[("u", u.values)])
    else
      vtkexportmesh(rfile, fens, femm.geod.fes; vectors=[("u", u.values)])
    end
  end

  return modeldata
end

"""
    AlgoDeforLinearModule.exportstress(modeldata::FDataDict)

Export the stress/deformation for visualization in Paraview.
"""
function exportstress(modeldata::FDataDict)
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

  # Plot the surface for each region
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
      bfes = meshboundary(femm.geod.fes);
      vtkexportmesh(rfile, fens, bfes;
        scalars=[(string(quantity)*string(component), fld.values)],
        vectors=[("u", u.values)])
    else
      vtkexportmesh(rfile, fens, femm.geod.fes;
        scalars=[(string(quantity)*string(component), fld.values)],
        vectors=[("u", u.values)])
    end
  end

  return modeldata
end

"""
    AlgoDeforLinearModule.exportstresselementwise(modeldata::FDataDict)

Export the elementwise stress/deformation for visualization in Paraview.
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

  # Plot the surface for each region
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
      bfes = meshboundary(femm.geod.fes);
      vtkexportmesh(rfile, fens, bfes;
        scalars=[(string(quantity)*string(component), fld.values)],
        vectors=[("u", u.values)])
    else
      vtkexportmesh(rfile, fens, femm.geod.fes;
        scalars=[(string(quantity)*string(component), fld.values)],
        vectors=[("u", u.values)])
    end
  end

  return modeldata
end

function modal(modeldata::FDataDict)
  # Modal (free-vibration) analysis solver.
  #
  # function model_data = deformation_linear_modal_analysis(model_data)
  #
  # Arguments
  # model_data = struct  with fields as follows.
  #
  # model_data.fens = finite element node set (mandatory)
  #
  # For each region (connected piece of the domain made of a particular material),
  # mandatory:
  # model_data.region= cell array of struct with the attributes, each region
  #           gets a struct with attributes
  #     fes= finite element set that covers the region
  #     integration_rule =integration rule; alternatively, one could
  #           specify separate integration rules for the stiffness and
  #           for the mass matrix.
  #     integration_rule_stiffness =integration rule for the stiffness
  #     integration_rule_mass =integration rule for the mass
  #     property = material property hint (optional): 'isotropic' (default),
  #           'orthotropic',...
  #        For isotropic property (default):
  #     E = material Young's modulus
  #     nu = material Poisson's ratio
  #        For orthotropic property:
  #     E1, E2, E3 = material Young's modulus
  #           in the three material directions
  #     G12, G13, G23 = material shear modulus
  #           between the three material directions
  #     nu12, nu13, nu23 = material Poisson's ratio
  #     rho = mass density (optional for statics, mandatory for dynamics)
  #        If the region is for a two-dimensional model (plane strain, plane
  #        stress, or axially symmetric) the attribute reduction  needs to be
  #        specified.
  #     reduction = 'strain' or 'stress' or 'axisymm'
  #
  #        If the material orientation matrix is not the identity and needs
  #        to be supplied,  include the attribute
  #     Rm= constant orientation matrix or a handle  to a function to compute
  #           the orientation matrix (see the class femm_base).
  #
  # Example:
  #     clear region
  #     region.E =E;
  #     region.nu=nu;
  #     region.fes= fes;
  #     region.integration_rule = gauss_rule (struct('dim', 3, 'order', 2));
  #     model_data.region{1} =region;
  #
  # Example:
  #     clear region
  #     region.property = 'orthotropic';
  #     region.rho =rho;
  #     region.E1 =E1;     region.E2 =E2;     region.E3 =E3;
  #     region.G12=G12;     region.G13=G13;     region.G23=G23;
  #     region.nu12=nu12;     region.nu13=nu13;     region.nu23=nu23;
  #     region.fes= fes;# set of finite elements for the interior of the domain
  #     region.integration_rule = gauss_rule (struct('dim', 3, 'order', 2));
  #     region.Rm =@LayerRm;
  #     model_data.region{1} =region;
  #
  # For essential boundary conditions (optional):
  # model_data.boundary_conditions.essential = cell array of struct,
  #           each piece of surface with essential boundary condition gets one
  #           element of the array with a struct with the attributes
  #           as recognized by the set_ebc() method of nodal_field
  #     component=which component is affected (if not supplied,
  #           or if supplied as empty default is all components)
  #     fixed_value=is always  0.0 (only that makes sense in modal analysis)
  #     fes = finite element set on the boundary to which
  #                       the condition applies
  #               or alternatively
  #     node_list = list of nodes on the boundary to which
  #                       the condition applies
  #           Only one of fes and node_list needs to be given.
  #
  # Example:
  #     clear essential
  #     essential.component= [1,3];
  #     essential.fixed_value= 0;
  #     essential.node_list = [[fenode_select(fens, struct('box', [0,a,0,0,-Inf,Inf],...
  #         'inflate',tolerance))],[fenode_select(fens, struct('box', [0,a,b,b,-Inf,Inf],...
  #         'inflate',tolerance))]];;
  #     model_data.boundary_conditions.essential{2} = essential;
  #
  # Example:
  #     clear essential
  #     essential.component= [1];
  #     essential.fixed_value= 0;
  #     essential.fes = mesh_boundary(fes);
  #     model_data.boundary_conditions.essential{1} = essential;
  #
  # For multi point constraints (MPC) (optional):
  # model_data.mpc= cell array of structs, each for one MPC.
  #      mpc.node_list = list of node numbers involved in the MPC,
  #      mpc.dof_list= numbers of degrees of freedom for the nodes above,
  #      mpc.umultipliers=multipliers for the nodes above,
  #      mpc.penfact=the penalty factor to multiply  the constraint matrix,
  #          The MPC looks like this: sum_i m_i u_{dof(i),node(i)} =0
  #          where m_i is the multiplier.
  #
  # Control parameters:
  # model_data.neigvs = number of eigenvalues/eigenvectors to compute
  # model_data.omega_shift= angular frequency shift for mass shifting
  # model_data.renumber = true or false flag (default is true)
  # model_data.renumbering_method = optionally choose the renumbering
  #       method  (symrcm or symamd)
  # model_data.use_lumped_mass = true or false?  (Default is false: consistent
  #       mass)
  #
  # Output
  # model_data = updated structure that was given on input
  # The struct model_data on output incorporates in addition to the input the fields
  #     model_data.geom = geometry field
  #     model_data.u = displacement field
  #     model_data.neigvs=Number of computed eigenvectors
  #     model_data.W = Computed eigenvectors, neigvs columns
  #     model_data.Omega=  Computed angular frequencies, array of length neigvs

  # What is the model reduction scheme?  Default is 3-D.
  modelreduction = get(modeldata, "modelreduction", DeformationModelReduction3D);

  neigvs = get(modeldata, "neigvs", 7); # Number of eigenvalues

  omega_shift = get(modeldata, "omega_shift", 0.0); # Mass shifting

  use_factorization = get(modeldata, "use_factorization", false); # Factorization?

  use_lumped_mass = get(modeldata, "use_lumped_mass", false); # Lumped mass?

  # # Should we renumber the nodes to minimize the cost of the solution of
  # # the coupled linear algebraic equations?
  # if (renumber) && (use_factorization)
  #     renumbering_method = 'symamd'; # default choice
  #     if ( isfield(model_data,'renumbering_method'))
  #         renumbering_method  =model_data.renumbering_method;;
  #     end
  #     # Run the renumbering algorithm
  #     model_data =renumber_mesh(model_data, renumbering_method);;
  #     # Save  the renumbering  (permutation of the nodes)
  #     clear  Renumbering_options; Renumbering_options.node_perm  =model_data.node_perm;
  # end

  # Extract the nodes
  fens=get(()->error("Must get fens!"), modeldata, "fens")

  # Construct the geometry field
  geom = NodalField(name ="geom",data =fens.xyz)

  # Construct the displacement field
  u = NodalField(name ="u",data =zeros(nnodes(geom),ndofn(geom)))

  # Apply the essential boundary conditions on the displacement field
  boundary_conditions = get(()->error("Must get boundary conditions!"), modeldata, "boundary_conditions");
  essential = get(boundary_conditions, "essential", nothing);
  if (essential!= nothing)
    for j=1:length(essential)
      fes = get(essential[j], "fes", nothing);
      if (fes!= nothing)
        fenids= connectednodes(fes);
      else
        fenids= get(()->error("Must get node list!"), essential[j], "node_list");
      end
      component=get(essential[j], "component", 1:ndofn(u)); # which components?  Default is all
      u_fixed=0.0;  # nonzero displacement does not make sense in modal analysis
      for k=1:length(component)
        for n=1:length(fenids)
          setebc!(u,fenids[n],true,k,u_fixed);
        end
      end
      applyebc!(u);
    end
  end


  # Number the equations
  numberdofs!(u)           #,Renumbering_options); # NOT DONE <<<<<<<<<<<<<<<<<

  # Initialize the heat loads vector
  F =zeros(FFlt,u.nfreedofs);

  # Create the finite element models for the regions (if not supplied as input)
  region=get(()->error("Must get region!"), modeldata, "region")
  for i=1:length(region)
    femm=get(region[i], "femm", nothing);
    if (femm==nothing) # see if you got the specialized FEMMs
      femm_stiffness=get(region[i], "femm_stiffness", nothing);
      femm_mass=get(region[i], "femm_mass", nothing);
      if (femm_stiffness==nothing) || (femm_mass==nothing)
        #need to construct the FEMM
        # Construct the property  and material  objects
        property=get(()->error("Must get property!"), region[i], "property");
        mater=MaterialDeformationLinear(property);
        Rm=get(region[i], "Rm", MaterialOrientation());
        # This is the model object for the current region: note that we supply
        # integration  rule and the  material orientation matrix
        fes=get(()->error("Must get finite elements!"), region[i], "fes");
        integration_rule_stiffness=get(region[i], "integration_rule_stiffness", nothing);
        integration_rule_mass=get(region[i], "integration_rule_mass", nothing);
        if (integration_rule_stiffness==nothing) || (integration_rule_mass==nothing)
          integration_rule=get(()->error("Must get integration rule!"), region[i], "integration_rule");
          integration_rule_stiffness=integration_rule
          integration_rule_mass=integration_rule
        end
        femm_stiffness=FEMMDeformationLinear(FEMMBase(fes,integration_rule_stiffness,Rm),mater);
        femm_mass=FEMMDeformationLinear(FEMMBase(fes,integration_rule_mass,Rm),mater);
      end
    else   # assume both should be the same
      femm_stiffness=femm;
      femm_mass=femm
    end
    setindex!(region[i], femm_stiffness, "femm_stiffness");
    setindex!(region[i], femm_mass, "femm_mass");
  end
  setindex!(modeldata, region, "region"); # put the generated data back


  # Construct the system stiffness matrix
  K=  spzeros(u.nfreedofs,u.nfreedofs); # (all zeros, for the moment)
  M=  spzeros(u.nfreedofs,u.nfreedofs); # (all zeros, for the moment)
  region=get(()->error("Must get region!"), modeldata, "region")
  for i=1:length(region)
    femm_stiffness=get(region[i], "femm_stiffness", nothing);
    femm_mass=get(region[i], "femm_mass", nothing);
    # # Give the  FEMM a chance  to precompute  geometry-related quantities
    # region.femm = associate_geometry(region.femm,geom);
    K = K + stiffness(modelreduction, femm_stiffness, geom, u);
    M = M + mass(modelreduction, femm_mass, geom, u);
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
  #          [W,Omega]= eig(full(K+omega_shift*M), full(M));
  #  end                         #

  d,v,nev,nconv =eigs(K+omega_shift*M, M; nev=neigvs, which=:SM)
  d = d - omega_shift;

  #    Subtract the mass-shifting Angular frequency
  if any(imag(d) .!= 0.0)
    warn("Some complex angular frequencies detected");
    warn(" $(d)");
    d=real(d);
  end
  if any(real(d) .< 0.0)
    warn("Some negative angular frequencies detected");
    warn(" $(d)");
    d=abs(d);
  end
  #    Sort  the angular frequencies by magnitude.  Make sure all
  #    imaginary parts of the eigenvalues are removed.
  ix =sortperm(d);

  # Update the model data: store geometry
  modeldata["geom"] = geom;
  # Store the displacement field
  modeldata["u"] = u;
  # Number of computed eigenvectors
  modeldata["neigvs"]=nev;
  #  Computed eigenvectors: we are ignoring the imaginary part here
  #  because the modal analysis is presumed to have been performed for
  #  an undamped structure
  modeldata["W"] = real(v[:,ix]);
  #  Computed angular frequencies
  modeldata["omega"]=sqrt(d[ix]);
  return modeldata
end

function exportmode(modeldata::FDataDict)
  # Export the modal deformation for visualization in Paraview.
  #
  # Arguments
  # model_data= model data as produced by deformation_linear_statics()
  # model_data.postprocessing= optional struct with optional fields
  #      gv = graphic viewer; if not supplied, a graphic viewer is created
  #           and returned in options.gv
  #      u_scale = deflection scale, default 1.0;
  #      modelist= default is 1:model_data.neigvs.
  #      save_frame= should we save images for the modes displayed? default false;
  #      frame_name= name for the mode images; default 'frame_name';
  #      camera  = camera, default is [] which means use the default orientation
  #           of the view;
  #      cmap= colormap (default: jet)
  #
  # Output
  # model_data = structure on input updated with
  # model_data.postprocessing.gv=graphic viewer used to display the data


  # Defaults
  mode=1;
  boundary_only= false;
  ffile= "deformation"
  # Let's have a look at what's been specified
  postprocessing = get(modeldata, "postprocessing", nothing);
  if (postprocessing!=nothing)
    mode=  get(postprocessing, "mode", mode);
    boundary_only=  get(postprocessing, "boundary_only", boundary_only);
    ffile=  get(postprocessing, "file", ffile);
  end

  omega=modeldata["omega"]
  if (length(omega)<mode) || (mode<0)
    error("Invalid node number $mode")
  end


  fens=get(()->error("Must get fens!"), modeldata, "fens")
  geom = get(()->error("Must get geometry field!"), modeldata, "geom");
  u = get(()->error("Must get displacement field!"), modeldata, "u");

  # Scatter the desired mode
  W=modeldata["W"]
  scattersysvec!(u,W[:,mode])

  # Plot the surface for each region
  region=get(()->error("Must get region!"), modeldata, "region")
  for i=1:length(region)
    femm=get(region[i], "femm", nothing);
    if (femm==nothing)
      femm=get(region[i], "femm_stiffness", nothing);
    end
    if (femm==nothing)
      error("No FEMM")
    end
    rfile = ffile * "$i" * ".vtk";
    if boundary_only
      bfes= meshboundary(femm.femmbase.fes);
      vtkexportmesh(rfile, fens, bfes;  vectors=u.values, vectors_name ="EV$mode")
    else
      vtkexportmesh(rfile, fens, femm.femmbase.fes;  vectors=u.values, vectors_name ="EV$mode")
    end
  end

  return modeldata
end

end
