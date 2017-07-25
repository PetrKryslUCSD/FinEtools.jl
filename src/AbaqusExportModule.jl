"""
AbaqusExportModule

Module to export an abaqus INP file.

## Example:
Very simple  definition of a model.
```
AE = AbaqusExporter();
open("f" * ".inp");
HEADING(["My first example ");
NODE(model_data.geom.values);
ELEMENT(ElType,"All",1,fes.conn);
ELEMENT("SFM3D4","Traction",count(fes)+1,model_data.boundary_conditions.traction{1}.fes.conn);
NSET_NSET("Corner",nc);
ORIENTATION("Orientation", axis, aangle);
MATERIAL("Material");
ELASTIC(E1,E2,E3,nu12,nu13,nu23,G12,G13,G23);
SOLID_SECTION("Material","Orientation","All");
SURFACE_SECTION("Traction");
STEP_PERTURBATION("Linear solution");
BOUNDARY_direct(model_data.u.is_fixed,model_data.u.fixed_values);
DLOAD("Traction",model_data.boundary_conditions.traction{1}.traction);
NODE_PRINT("Corner");
END_STEP();
close();
```
## Example:
In this example we define a part, and an instance within an assembly.
```

AE=Abaqus_exporter;
open(["f-" num2str(ref) ".inp"]);
HEADING([mfilename " " "ElType=" ElType " " "ref=" num2str(ref)]);
PART("part1");
END_PART();
ASSEMBLY("ASSEM1");
INSTANCE("INSTNC1","PART1");
NODE(model_data.geom.values);
ELEMENT(ElType,"All",1,fes.conn);
ELEMENT(SurfElType,"Traction",count(fes)+1,model_data.boundary_conditions.traction{1}.fes.conn);
NSET_NSET("Corner",nc);
ORIENTATION("Orientation", Rm(:,1), Rm(:,2));
SOLID_SECTION("Material","Orientation","All");
SURFACE_SECTION("Traction");
NSET_NSET("xfix",find(model_data.u.is_fixed(:,1)));
NSET_NSET("yfix",find(model_data.u.is_fixed(:,2)));
NSET_NSET("zfix",find(model_data.u.is_fixed(:,3)));
END_INSTANCE();
END_ASSEMBLY();
MATERIAL("Material");
ELASTIC(E1,E2,E3,nu12,nu13,nu23,G12,G13,G23);
STEP_PERTURBATION("Linear solution");
DLOAD("ASSEM1.INSTNC1.Traction",model_data.boundary_conditions.traction{1}.traction);
BOUNDARY("ASSEM1.INSTNC1.xfix",1);
BOUNDARY("ASSEM1.INSTNC1.yfix",2);
BOUNDARY("ASSEM1.INSTNC1.zfix",3);
NODE_PRINT("ASSEM1.INSTNC1.Corner");
ENERGY_PRINT();
END_STEP();
close();
```

"""
module AbaqusExportModule

using FinEtools.FTypesModule
using FinEtools.FESetModule
using FinEtools.FENodeSetModule
import Base.open
import Base.close

export AbaqusExporter
export close
export HEADING, COMMENT
export PART, END_PART, ASSEMBLY, END_ASSEMBLY, INSTANCE, END_INSTANCE
export NODE, ELEMENT, NSET_NSET, ELSET_ELSET, ORIENTATION, MATERIAL,
  ELASTIC, EXPANSION, DENSITY, SECTION_CONTROLS, SOLID_SECTION,
  SURFACE_SECTION
export STEP_PERTURBATION_STATIC, STEP_FREQUENCY, BOUNDARY, DLOAD, CLOAD,
  TEMPERATURE,
  END_STEP,  NODE_PRINT, EL_PRINT,  ENERGY_PRINT

mutable struct AbaqusExporter
  filename::AbstractString
  ios::IO
  element_range::Tuple{Int64, Int64}
  function AbaqusExporter(filename::AbstractString)
    if match(r".*\.inp$", filename) == nothing
      filename = filename * ".inp"
    end
    ios = open(filename,"w+")
    return new(deepcopy(filename), ios, (typemax(Int64), 0))
  end
end

"""
    COMMENT(self::AbaqusExporter, Text::AbstractString)

Write out the `**` comment option.
"""
function COMMENT(self::AbaqusExporter, Text::AbstractString)
  println(self.ios, "**" * Text)
end

"""
    HEADING(self::AbaqusExporter, Text::AbstractString)

Write out the `*HEADING` option.
"""
function HEADING(self::AbaqusExporter, Text::AbstractString)
  println(self.ios, "*HEADING");
  println(self.ios, Text)
end

"""
    PART(self::AbaqusExporter, NAME::AbstractString)

Write out the `*PART` option.
"""
function PART(self::AbaqusExporter, NAME::AbstractString)
  println(self.ios, "*PART, NAME=" * NAME)
end

"""
    END_PART(self::AbaqusExporter)

Write out the `*END PART` option.
"""
function END_PART(self::AbaqusExporter)
  println(self.ios, "*END PART")
end

"""
    ASSEMBLY(self::AbaqusExporter, NAME::AbstractString)

Write out the `*ASSEMBLY` option.
"""
function ASSEMBLY(self::AbaqusExporter, NAME::AbstractString)
  println(self.ios, "*ASSEMBLY, NAME=" * NAME)
end

"""
    END_ASSEMBLY(self::AbaqusExporter)

Write out the `*END ASSEMBLY` option.
"""
function END_ASSEMBLY(self::AbaqusExporter)
  println(self.ios, "*END ASSEMBLY")
end

"""
    INSTANCE(self::AbaqusExporter, NAME::AbstractString, PART::AbstractString)

Write out the `*INSTANCE` option.
"""
function INSTANCE(self::AbaqusExporter, NAME::AbstractString, PART::AbstractString)
  println(self.ios, "*INSTANCE, NAME=" * NAME * ", PART=" * PART)
end

"""
    END_INSTANCE(self::AbaqusExporter)

Write out the `*END INSTANCE` option.
"""
function END_INSTANCE(self::AbaqusExporter)
  println(self.ios, "*END INSTANCE")
end

"""
    NODE(self::AbaqusExporter, xyz::AbstractArray{T, 2}) where {T}

Write out the `*NODE` option.

`xyz`=array of node coordinates
"""
function NODE(self::AbaqusExporter, xyz::AbstractArray{T, 2}) where {T}
  println(self.ios, "*NODE")
  for j = 1:size(xyz,1)
    print(self.ios, "$j,$(xyz[j, 1])")
    for ixxxx = 2:1:size(xyz,2)
      print(self.ios, ",$(xyz[j, ixxxx])")
    end
    print(self.ios, "\n")
  end
end

"""
    ELEMENT(self::AbaqusExporter, TYPE::AbstractString, ELSET::AbstractString,
      start::Integer, conn::AbstractArray{T, 2}) where {T<:Integer}

Write out the `*ELEMENT` option.

`TYPE`= element type code,
`ELSET`= element set to which the elements belong,
`start`= start the element numbering at this integer,
`conn`= connectivity array for the elements, one row per element
"""
function ELEMENT(self::AbaqusExporter, TYPE::AbstractString, ELSET::AbstractString,
  start::Integer, conn::AbstractArray{T, 2}) where {T<:Integer}
  # Check that start is valid
  @assert start > 0 "The  starting element number must be > 0"
  # Check that the current element number range is disjoint from the
  # range of the elements on input: they must not overlap
  @assert (self.element_range[2] < start) ||
          (self.element_range[1] > start+size(conn, 1)-1) "Elements  must be given unique numbers"
  # Update the element range
  self.element_range = (min(self.element_range[1], start),
                        max(self.element_range[2], start+size(conn, 1)-1))
  println(self.ios, "*ELEMENT, TYPE =" * TYPE * ", ELSET=" * ELSET)
  for j=1:size(conn,1)
    print(self.ios, "$(j+start-1),")
    for ixxxx = 1:size(conn,2)-1
      if ixxxx>15
        print(self.ios, "\n")
      end
      print(self.ios, "$(conn[j,ixxxx]),")
    end
    print(self.ios, "$(conn[j,size(conn,2)])\n")
  end
end

"""
    NSET_NSET(self::AbaqusExporter, NSET::AbstractString,
      n::AbstractVector{T}) where {T<:Integer}

Write out the `*NSET` option.

`NSET` = name of the set,
`n` = array of the node numbers
"""
function NSET_NSET(self::AbaqusExporter, NSET::AbstractString, n::AbstractVector{T}) where {T<:Integer}
  println(self.ios, "*NSET, NSET=" * NSET)
  for j=1:length(n)-1
    println(self.ios, "$(n[j]),")
  end
  println(self.ios, "$(n[length(n)])")
end

"""
    ELSET_ELSET(self::AbaqusExporter, ELSET::AbstractString, n::AbstractArray{T, 1}) where {T<:Integer}

Write out the `*ELSET` option.

`ELSET` = name of the set,
`n` = array of the node numbers
"""
function ELSET_ELSET(self::AbaqusExporter, ELSET::AbstractString, n::AbstractArray{T, 1}) where {T<:Integer}
  println(self.ios, "*ELSET, ELSET=" * ELSET)
  for j=1:length(n)-1
    println(self.ios, "$(n[j]),")
  end
  println(self.ios, "$(n[length(n)])")
end

"""
    ORIENTATION(self::AbaqusExporter, ORIENTATION::AbstractString,
      a::AbstractArray{T,1}, b::AbstractArray{T,1})

Write out the `*ORIENTATION` option.

Invoke at level: Part,  Part instance,  Assembly
"""
function ORIENTATION(self::AbaqusExporter, ORIENTATION::AbstractString,
  a::AbstractArray{T,1}, b::AbstractArray{T,1}) where {T<:Real}
  println(self.ios, "*ORIENTATION,NAME=" * ORIENTATION)
  println(self.ios, "$(a[1]),$(a[2]),$(a[3]),$(b[1]),$(b[2]),$(b[3])")
  println(self.ios, "1,0.0")
end

"""
    MATERIAL(self::AbaqusExporter, MATERIAL::AbstractString)

Write out the `*MATERIAL` option.
"""
function MATERIAL(self::AbaqusExporter, MATERIAL::AbstractString)
  println(self.ios, "*MATERIAL,NAME=" * MATERIAL);
end

"""
    ELASTIC(self::AbaqusExporter, E1::F, E2::F, E3::F, nu12::F, nu13::F, nu23::F,
      G12::F, G13::F, G23::F) where {F}

Write out the `*ELASTIC,TYPE=ENGINEERING CONSTANTS` option.
"""
function ELASTIC(self::AbaqusExporter, E1::F, E2::F, E3::F, nu12::F, nu13::F, nu23::F,
  G12::F, G13::F, G23::F) where {F}
  println(self.ios, "*ELASTIC,TYPE=ENGINEERING CONSTANTS");
  println(self.ios, "$E1,$E2,$E3,$nu12,$nu13,$nu23,$G12,$G13,");
  println(self.ios, "$G23");
end

"""
    ELASTIC(self::AbaqusExporter, E::F, nu::F) where {F}

Write out the `*ELASTIC,TYPE=ISOTROPIC` option.
"""
function ELASTIC(self::AbaqusExporter, E::F, nu::F) where {F}
  println(self.ios, "*ELASTIC,TYPE=ISOTROPIC ");
  println(self.ios, "$E,$nu");
end

"""
    EXPANSION(self::AbaqusExporter, CTE::F) where {F}

Write out the `*EXPANSION` option.
"""
function EXPANSION(self::AbaqusExporter, CTE::F) where {F}
  println(self.ios, "*EXPANSION");
  println(self.ios, "$CTE,");
end

"""
    DENSITY(self::AbaqusExporter, rho)

Write out the `*DENSITY` option.
"""
function DENSITY(self::AbaqusExporter, rho)
  println(self.ios, "*DENSITY ");
  println(self.ios, "$rho");
end

"""
    SECTION_CONTROLS(self::AbaqusExporter, NAME::AbstractString,
      OPTIONAL::AbstractString)

Write out the `*SECTION CONTROLS` option.

`OPTIONAL` = string, for instance
                     HOURGLASS=ENHANCED
"""
function SECTION_CONTROLS(self::AbaqusExporter, NAME::AbstractString,
  OPTIONAL::AbstractString)
  println(self.ios, "*SECTION CONTROLS, NAME=" * NAME * "," * OPTIONAL);
end

"""
    SOLID_SECTION(self::AbaqusExporter, MATERIAL::AbstractString,
      ORIENTATION::AbstractString, ELSET::AbstractString,
      CONTROLS::AbstractString)

Write out the `*SOLID SECTION` option.

Level: Part,  Part instance

"""
function SOLID_SECTION(self::AbaqusExporter, MATERIAL::AbstractString,
  ORIENTATION::AbstractString, ELSET::AbstractString,
  CONTROLS::AbstractString)
  if "CONTROLS" == ""
    println(self.ios, "*SOLID SECTION,MATERIAL=" * MATERIAL *
      ",ORIENTATION =" * ORIENTATION * ",ELSET=" * ELSET);
  else
    println(self.ios, "*SOLID SECTION,MATERIAL=" * MATERIAL *
      ",ORIENTATION =" * ORIENTATION * ",ELSET=" * ELSET *
      ", CONTROLS  =" * CONTROLS);
  end
end

"""
    SOLID_SECTION(self::AbaqusExporter, MATERIAL::AbstractString,
      ORIENTATION::AbstractString, ELSET::AbstractString)

Write out the `*SOLID SECTION` option.

Level: Part,  Part instance

"""
function SOLID_SECTION(self::AbaqusExporter, MATERIAL::AbstractString,
  ORIENTATION::AbstractString, ELSET::AbstractString)
  println(self.ios, "*SOLID SECTION,MATERIAL=" * MATERIAL *
      ",ORIENTATION =" * ORIENTATION * ",ELSET=" * ELSET);
end

"""
    HOURGLASS(self::AbaqusExporter, KIND::AbstractString, VALUE::F) where {F}

Write out the `*HOURGLASS` option.

```
Example:
*SOLID SECTION,ELSET=SOLID3,MATERIAL=MAT,CONTROL=A
*HOURGLASS STIFFNESS
5.E8
```
"""
function HOURGLASS(self::AbaqusExporter, KIND::AbstractString, VALUE::F) where {F}
  println(self.ios, "*HOURGLASS " * KIND);
  println(self.ios, "$( VALUE )");
end

"""
    SURFACE_SECTION(self::AbaqusExporter, ELSET::AbstractString)

Write out the `*SURFACE SECTION` option.
"""
function SURFACE_SECTION(self::AbaqusExporter, ELSET::AbstractString)
  println(self.ios, "*SURFACE SECTION, ELSET=" * ELSET);
end

"""
    STEP_PERTURBATION_STATIC(self::AbaqusExporter)

Write out the `*STEP,PERTURBATION` option for linear static analysis.
"""
function STEP_PERTURBATION_STATIC(self::AbaqusExporter)
  println(self.ios, "*STEP,PERTURBATION");
  println(self.ios, "*STATIC");
end

"""
    STEP_PERTURBATION_BUCKLE(self::AbaqusExporter, neigv::Integer)

Write out the `*STEP,PERTURBATION` option for linear buckling analysis.
"""
function STEP_PERTURBATION_BUCKLE(self::AbaqusExporter, neigv::Integer)
  println(self.ios, "*STEP, name=Buckling, nlgeom=NO, perturbation");
  println(self.ios, "*BUCKLE");
  println(self.ios, "$(neigv), , , , ");
end

"""
    STEP_FREQUENCY(self::AbaqusExporter, Nmodes::Integer)

Write out the `*STEP,FREQUENCY` option.
"""
function STEP_FREQUENCY(self::AbaqusExporter, Nmodes::Integer)
  println(self.ios, "*STEP");
  println(self.ios, "*FREQUENCY, EIGENSOLVER=LANCZOS");
  println(self.ios, "$(Nmodes), , ,-1.E6 \n");
end

"""
    BOUNDARY(self::AbaqusExporter, is_fixed::AbstractArray{B,2},
      fixed_value::AbstractArray{F,2}) where {B, F}

Write out the `*BOUNDARY` option.

`is_fixed`= array of Boolean flags (true means fixed, or prescribed),  one row per node,
`fixed_value`=array of displacements to which the corresponding displacement components is fixed
"""
function BOUNDARY(self::AbaqusExporter, is_fixed::AbstractArray{B,2},
  fixed_value::AbstractArray{F,2}) where {B, F}
  println(self.ios, "*BOUNDARY");
  for j=1:size(is_fixed,1)
    for k=1:size(is_fixed,2)
      #<node number>, <first dof>, <last dof>, <magnitude of displacement>
      if is_fixed[j,k]
        println(self.ios, "$j,$k,$k,$fixed_value[j,k]");
      end
    end
  end
end

"""
    BOUNDARY(self::AbaqusExporter, NSET::AbstractString, dof::Integer)

Write out the `*BOUNDARY` option to fix displacements at zero for a node set.

Invoke at Level: Model,  Step

`NSET`= node set,
`dof`=Degree of freedom, 1, 2, 3
"""
function BOUNDARY(self::AbaqusExporter, NSET::AbstractString, dof::Integer)
  println(self.ios, "*BOUNDARY");
  println(self.ios, NSET * ",$dof");
end

"""
    BOUNDARY(self::AbaqusExporter, NSET::AbstractString, dof::Integer,
      value::F) where {F}

Write out the `*BOUNDARY` option to fix displacements at nonzero value for a
node set.

`NSET`= node set,
`dof`=Degree of freedom, 1, 2, 3
"""
function BOUNDARY(self::AbaqusExporter, NSET::AbstractString, dof::Integer, value::F) where {F}
  println(self.ios, "*BOUNDARY,TYPE=DISPLACEMENT");
  println(self.ios, NSET * ",$dof,$dof,$value");
end

"""
    DLOAD(self::AbaqusExporter, ELSET::AbstractString,
      traction::AbstractVector{F}) where {F}

Write out the `*DLOAD` option.
"""
function DLOAD(self::AbaqusExporter, ELSET::AbstractString,
  traction::AbstractVector{F}) where {F}
  println(self.ios, "*DLOAD, follower=NO");
  nt = norm(traction)
  print(self.ios, ELSET * ",TRVEC," * "$(nt),")
  println(self.ios, "$(traction[1]/nt),$(traction[2]/nt),$(traction[3]/nt)");
end

"""
    CLOAD(self::AbaqusExporter, NSET::AbstractString, dof::Integer,
      magnitude::F) where {F}

Write out the `*CLOAD` option.

NSET=Number of node
dof= 1, 2, 3,
magnitude= signed multiplier
"""
function CLOAD(self::AbaqusExporter, NSET::AbstractString, dof::Integer,
  magnitude::F) where {F}
  println(self.ios, "*CLOAD");
  println(self.ios, NSET * "$dof,$magnitude");
end

"""
    TEMPERATURE(self::AbaqusExporter, nlist::AbstractArray{I, 1},
      tlist::AbstractArray{F, 1}) where {I, F}

Write out the `*TEMPERATURE` option.
"""
function TEMPERATURE(self::AbaqusExporter, Classifier::AbstractString,
  nlist::AbstractArray{I, 1},
  tlist::AbstractArray{F, 1}) where {I, F}
  @assert length(nlist) == length(tlist)
  println(self.ios, "*TEMPERATURE");
  for ixxxx = 1:length(nlist)
    println(self.ios, Classifier * "$(nlist[ixxxx]),$(tlist[ixxxx])");
  end
end

"""
    NODE_PRINT(self::AbaqusExporter, NSET::AbstractString)

Write out the `*NODE PRINT` option.
"""
function NODE_PRINT(self::AbaqusExporter, NSET::AbstractString; KEYS = "U")
  println(self.ios, "*NODE PRINT, NSET=" * NSET);
  println(self.ios, KEYS);
end

"""
    EL_PRINT(self::AbaqusExporter, ELSET::AbstractString, KEYS::AbstractString)

Write out the `*EL PRINT` option.
"""
function EL_PRINT(self::AbaqusExporter, ELSET::AbstractString, KEYS::AbstractString)
  println(self.ios, "*EL PRINT, ELSET=" * ELSET * ", POSITION=INTEGRATION POINTS, SUMMARY= YES");
  println(self.ios, "");
  println(self.ios, KEYS);
end

"""
    ENERGY_PRINT(self::AbaqusExporter)

Write out the `*ENERGY PRINT` option.
"""
function ENERGY_PRINT(self::AbaqusExporter)
  println(self.ios, "*ENERGY PRINT");
end

"""
    END_STEP(self::AbaqusExporter)

Write out the `*END STEP` option.
"""
function END_STEP(self::AbaqusExporter)
  println(self.ios, "*END STEP");
end

"""
    close(self::AbaqusExporter)

Close  the strain opened by the exporter.
"""
function close(self::AbaqusExporter)
  close(self.ios)
end

end # module
