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
export open, close
export HEADING
export PART, END_PART, ASSEMBLY, END_ASSEMBLY, INSTANCE, END_INSTANCE
export NODE, ELEMENT, NSET_NSET, ELSET_ELSET, ORIENTATION, MATERIAL,
  ELASTIC, ELASTIC_ISOTROPIC, DENSITY, SECTION_CONTROLS

"""
    AbaqusExporter(filename::AbstractString)

Open an output stream for the Abaqus export.
"""
AbaqusExporter(filename::AbstractString) = open(filename * ".inp","w+")

"""
    HEADING(self::IOStream, Text::AbstractString)

Write out the `*HEADING` option.
"""
function HEADING(self::IOStream, Text::AbstractString)
  println(self, "*HEADING");
  println(self, Text)
end

"""
    PART(self::IOStream, NAME::AbstractString)

Write out the `*PART` option.
"""
function PART(self::IOStream, NAME::AbstractString)
  println(self, "*PART, NAME=" * NAME)
end

"""
    END_PART(self::IOStream)

Write out the `*END PART` option.
"""
function END_PART(self::IOStream)
  println(self, "*END PART")
end

"""
    ASSEMBLY(self::IOStream, NAME::AbstractString)

Write out the `*ASSEMBLY` option.
"""
function ASSEMBLY(self::IOStream, NAME::AbstractString)
  println(self, "*ASSEMBLY, NAME=" * NAME)
end

"""
    END_ASSEMBLY(self::IOStream)

Write out the `*END ASSEMBLY` option.
"""
function END_ASSEMBLY(self::IOStream)
  println(self, "*END ASSEMBLY")
end

"""
    INSTANCE(self::IOStream, NAME::AbstractString, PART::AbstractString)

Write out the `*INSTANCE` option.
"""
function INSTANCE(self::IOStream, NAME::AbstractString, PART::AbstractString)
  println(self, "*INSTANCE, NAME=" * NAME * ", PART=" * PART)
end

"""
    END_INSTANCE(self::IOStream)

Write out the `*END INSTANCE` option.
"""
function END_INSTANCE(self::IOStream)
  println(self, "*END INSTANCE")
end

"""
    NODE(self::IOStream, xyz::AbstractArray{T, 2}) where {T}

Write out the `*NODE` option.

`xyz`=array of node coordinates
"""
function NODE(self::IOStream, xyz::AbstractArray{T, 2}) where {T}
  println(self, "*NODE")
  for j = 1:size(xyz,1)
    println(self, "$j,$(xyz[j, 1]),$(xyz[j, 2]),$(xyz[j, 3])")
  end
end

"""
    ELEMENT(self::IOStream, TYPE::AbstractString, ELSET::AbstractString,
      start::Integer, conn::AbstractArray{T, 2}) where {T<:Integer}

Write out the `*ELEMENT` option.

`TYPE`= element type code,
`ELSET`= element set to which the elements belong,
`start`= start the element numbering at this integer,
`conn`= connectivity array for the elements, one row per element
"""
function ELEMENT(self::IOStream, TYPE::AbstractString, ELSET::AbstractString,
  start::Integer, conn::AbstractArray{T, 2}) where {T<:Integer}
  println(self, "*ELEMENT, TYPE =" * TYPE * ", ELSET=" * ELSET)
  nn9=size(conn,2); nn8=min(nn9,15);
  for j=1:size(conn,1)
    print(self, "$(j+start-1)")
    for ixxxx = 1:min(15, size(conn,2))
      print(self, ",$(conn[j,ixxxx])")
    end
    if (nn9>15)
      print(self, ",\n")
      for ixxxx = 15+1:size(conn,2)-1
        print(self, "$(conn[j,ixxxx]),")
      end
      for ixxxx = size(conn,2)
        print(self, "$(conn[j,ixxxx])")
      end
    end
    print(self, "\n")
  end
end

"""
    NSET_NSET(self::IOStream, NSET::AbstractString,
      n::AbstractVector{T}) where {T<:Integer}

Write out the `*NSET` option.

`NSET` = name of the set,
`n` = array of the node numbers
"""
function NSET_NSET(self::IOStream, NSET::AbstractString, n::AbstractVector{T}) where {T<:Integer}
  println(self, "*NSET, NSET=" * NSET)
  for j=1:length(n)-1
    println(self, "$(n[j]),")
  end
  println(self, "$(n[length(n)])")
end

"""
    ELSET_ELSET(self::IOStream, ELSET::AbstractString, n::AbstractArray{T, 1}) where {T<:Integer}

Write out the `*ELSET` option.

`ELSET` = name of the set,
`n` = array of the node numbers
"""
function ELSET_ELSET(self::IOStream, ELSET::AbstractString, n::AbstractArray{T, 1}) where {T<:Integer}
  println(self, "*ELSET, ELSET=" * ELSET)
  for j=1:length(n)-1
    println(self, "$(n[j]),")
  end
  println(self, "$(n[length(n)])")
end

"""
    ORIENTATION(self::IOStream, ORIENTATION::AbstractString,
      a::AbstractArray{T,1}, b::AbstractArray{T,1})

Write out the `*ORIENTATION` option.
"""
function ORIENTATION(self::IOStream, ORIENTATION::AbstractString,
  a::AbstractArray{T,1}, b::AbstractArray{T,1}) where {T<:Real}
  println(self, "*ORIENTATION,NAME=" * ORIENTATION)
  println(self, "$(a[1]),$(a[2]),$(a[3]),$(b[1]),$(b[2]),$(b[3])")
  println(self, "1.0, 0.0")
end

"""
    MATERIAL(self::IOStream, MATERIAL::AbstractString)

Write out the `*MATERIAL` option.
"""
function MATERIAL(self::IOStream, MATERIAL::AbstractString)
  println(self, "*MATERIAL,NAME=" * MATERIAL);
end

"""
    ELASTIC(self::IOStream, E1::F, E2::F, E3::F, nu12::F, nu13::F, nu23::F,
      G12::F, G13::F, G23::F) where {F}

Write out the `*ELASTIC,TYPE=ENGINEERING CONSTANTS` option.
"""
function ELASTIC(self::IOStream, E1::F, E2::F, E3::F, nu12::F, nu13::F, nu23::F,
  G12::F, G13::F, G23::F) where {F}
  println(self, "*ELASTIC,TYPE=ENGINEERING CONSTANTS");
  println(self, "$E1,$E2,$E3,$nu12,$nu13,$nu23,$G12,$G13,");
  println(self, "$G23");
end

"""
    ELASTIC(self::IOStream, E::F, nu::F) where {F}
    
Write out the `*ELASTIC,TYPE=ISOTROPIC` option.
"""
function ELASTIC(self::IOStream, E::F, nu::F) where {F}
  println(self, "*ELASTIC,TYPE=ISOTROPIC ");
  println(self, "$E,$nu");
end

"""
    DENSITY(self::IOStream, rho)

Write out the `*DENSITY` option.
"""
function DENSITY(self::IOStream, rho)
  println(self, "*DENSITY ");
  println(self, "$rho");
end

"""
    SECTION_CONTROLS(self::IOStream, NAME::AbstractString,
      OPTIONAL::AbstractString)

Write out the `*SECTION CONTROLS` option.

`OPTIONAL` = string, for instance
                     HOURGLASS=ENHANCED
"""
function SECTION_CONTROLS(self::IOStream, NAME::AbstractString,
  OPTIONAL::AbstractString)
  println(self, "*SECTION CONTROLS, NAME=" * NAME * "," * OPTIONAL);
end

"""
    SOLID_SECTION(self::IOStream, MATERIAL::AbstractString,
      ORIENTATION::AbstractString, ELSET::AbstractString,
      CONTROLS::AbstractString)

Write out the `*SOLID SECTION` option.
"""
function SOLID_SECTION(self::IOStream, MATERIAL::AbstractString,
  ORIENTATION::AbstractString, ELSET::AbstractString,
  CONTROLS::AbstractString)
  if "CONTROLS" == ""
    println(self, "*SOLID SECTION,MATERIAL=" * MATERIAL *
      ",ORIENTATION =" * ORIENTATION * ",ELSET=" * ELSET);
  else
    println(self, "*SOLID SECTION,MATERIAL=" * MATERIAL *
      ",ORIENTATION =" * ORIENTATION * ",ELSET=" * ELSET *
      ", CONTROLS  =" * CONTROLS);
  end
end

"""
    HOURGLASS(self::IOStream, KIND::AbstractString, VALUE::F) where {F}

Write out the `*HOURGLASS` option.

```
Example:
*SOLID SECTION,ELSET=SOLID3,MATERIAL=MAT,CONTROL=A
*HOURGLASS STIFFNESS
5.E8
```
"""
function HOURGLASS(self::IOStream, KIND::AbstractString, VALUE::F) where {F}
  println(self, "*HOURGLASS " * KIND);
  println(self, "$( VALUE )");
end

"""
    SURFACE_SECTION(self::IOStream, ELSET::AbstractString)

Write out the `*SURFACE SECTION` option.
"""
function SURFACE_SECTION(self::IOStream, ELSET::AbstractString)
  println(self, "*SURFACE SECTION, ELSET=" * ELSET);
end

"""
    STEP_PERTURBATION_STATIC(self::IOStream)

Write out the `*STEP,PERTURBATION` option for linear static analysis.
"""
function STEP_PERTURBATION_STATIC(self::IOStream)
  println(self, "*STEP,PERTURBATION");
  println(self, "*STATIC");
end

"""
    STEP_PERTURBATION_BUCKLE(self::IOStream, neigv::Integer)

Write out the `*STEP,PERTURBATION` option for linear buckling analysis.
"""
function STEP_PERTURBATION_BUCKLE(self::IOStream, neigv::Integer)
  println(self, "*STEP, name=Buckling, nlgeom=NO, perturbation");
  println(self, "*BUCKLE");
  println(self, "$(neigv), , , , ");
end

"""
    STEP_FREQUENCY(self::IOStream, Nmodes::Integer)

Write out the `*STEP,FREQUENCY` option.
"""
function STEP_FREQUENCY(self::IOStream, Nmodes::Integer)
  println(self, "*STEP");
  println(self, "*FREQUENCY, EIGENSOLVER=LANCZOS");
  println(self, "$(Nmodes), , ,-1.E6 \n");
end

"""
    BOUNDARY(self::IOStream, is_fixed::AbstractArray{B,2},
      fixed_value::AbstractArray{F,2}) where {B, F}

Write out the `*BOUNDARY` option.

`is_fixed`= array of Boolean flags (true means fixed, or prescribed),  one row per node,
`fixed_value`=array of displacements to which the corresponding displacement components is fixed
"""
function BOUNDARY(self::IOStream, is_fixed::AbstractArray{B,2},
  fixed_value::AbstractArray{F,2}) where {B, F}
  println(self, "*BOUNDARY");
  for j=1:size(is_fixed,1)
    for k=1:size(is_fixed,2)
      #<node number>, <first dof>, <last dof>, <magnitude of displacement>
      if is_fixed[j,k]
        println(self, "$j,$k,$k,$fixed_value[j,k]");
      end
    end
  end
end

"""
    BOUNDARY(self::IOStream, NSET::AbstractString, dof::Integer)

Write out the `*BOUNDARY` option to fix displacements at zero for a node set.

`NSET`= node set,
`dof`=Degree of freedom, 1, 2, 3
"""
function BOUNDARY(self::IOStream, NSET::AbstractString, dof::Integer)
  println(self, "*BOUNDARY");
  println(self, NSET * ",$dof");
end

"""
    BOUNDARY(self::IOStream, NSET::AbstractString, dof::Integer,
      value::F) where {F}

Write out the `*BOUNDARY` option to fix displacements at nonzero value for a
node set.

`NSET`= node set,
`dof`=Degree of freedom, 1, 2, 3
"""
function BOUNDARY(self::IOStream, NSET::AbstractString, dof::Integer, value::F) where {F}
  println(self, "*BOUNDARY,TYPE=DISPLACEMENT");
  println(self, NSET * ",$dof,$dof,$value");
end

"""
    DLOAD(self::IOStream, ELSET::AbstractString,
      traction::AbstractVector{F}) where {F}

Write out the `*DLOAD` option.
"""
function DLOAD(self::IOStream, ELSET::AbstractString,
  traction::AbstractVector{F}) where {F}
  println(self, "*DLOAD, follower=NO");
  nt = norm(traction)
  print(self, ELSET * ",TRVEC," * "$(nt),")
  println(self, "$(traction[1]/nt),$(traction[2]/nt),$(traction[3]/nt)");
end

"""
    CLOAD(self::IOStream, NSET::AbstractString, dof::Integer,
      magnitude::F) where {F}

Write out the `*CLOAD` option.

NSET=Number of node
dof= 1, 2, 3,
magnitude= signed multiplier
"""
function CLOAD(self::IOStream, NSET::AbstractString, dof::Integer,
  magnitude::F) where {F}
  println(self, "*CLOAD");
  println(self, NSET * "$dof,$magnitude");
end

"""
    NODE_PRINT(self::IOStream, NSET::AbstractString)

Write out the `*NODE PRINT` option.
"""
function NODE_PRINT(self::IOStream, NSET::AbstractString; KEYS = "U")
  println(self, "*NODE PRINT, NSET=" * NSET);
  println(self, KEYS);
end

"""
    EL_PRINT(self::IOStream, ELSET::AbstractString, KEYS::AbstractString)

Write out the `*EL PRINT` option.
"""
function EL_PRINT(self::IOStream, ELSET::AbstractString, KEYS::AbstractString)
  println(self, "*EL PRINT, ELSET=" * ELSET * ", POSITION=INTEGRATION POINTS, SUMMARY= YES");
  println(self, "");
  println(self, KEYS);
end

"""
    ENERGY_PRINT(self::IOStream)

Write out the `*ENERGY PRINT` option.
"""
function ENERGY_PRINT(self::IOStream)
  println(self, "*ENERGY PRINT");
end

"""
    END_STEP(self::IOStream)

Write out the `*END STEP` option.
"""
function END_STEP(self::IOStream)
  println(self, "*END STEP");
end


end # module

# module mmmmmmmmmm
# using FinEtools
# using Base.Test
# using AbaqusExportModule
# function test()
#   AE = AbaqusExporter("SimpleExample");
#   HEADING(AE, "My first example ");
#   PART(AE, "part1");
#   END_PART(AE);
#   ASSEMBLY(AE, "ASSEM1");
#   INSTANCE(AE, "INSTNC1","PART1");
#   NODE(AE, rand(4,3));
#   ELEMENT(AE, "c3d8", "Allh8", 1, rand(UInt8, 4, 8)+1)
#   ELEMENT(AE, "c3d20", "Allh20", 1, rand(UInt8, 4, 20)+1)
#   NSET_NSET(AE, "Corner", vec([1 2 3 4]));
#   END_INSTANCE(AE);
#   END_ASSEMBLY(AE);
#   close(AE)
# end
# end
# using mmmmmmmmmm
# mmmmmmmmmm.test()
