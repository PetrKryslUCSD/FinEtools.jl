"""
    MeshExportModule

Module for export of meshes and data defined on meshes.
"""
module MeshExportModule

__precompile__(true)

function _match_type(d, t)
    for k in keys(d)
        if t <: k
            return d[k]
        end
    end
    return nothing
end

module VTK
################################################################################
# VTK export
################################################################################
import ...FESetModule:
    AbstractFESet,
    FESetP1,
    FESetL2,
    FESetT3,
    FESetQ4,
    FESetT4,
    FESetH8,
    FESetQ8,
    FESetQ9,
    FESetL3,
    FESetT6,
    FESetT10,
    FESetH20,
    connasarray
import ...FENodeSetModule: FENodeSet
import Base.close
using ..MeshExportModule: _match_type

using Printf
import LinearAlgebra: norm, cross

const P1 = 1
const L2 = 3
const T3 = 5
const Q4 = 9
const T4 = 10
const H8 = 12
const L3 = 21
const T6 = 22
const Q8 = 23
const Q9 = 28
const T10 = 24
const H20 = 25

VTKtypemap = Dict{UnionAll,Int}(
    FESetP1 => P1,
    FESetL2 => L2,
    FESetT3 => T3,
    FESetQ4 => Q4,
    FESetT4 => T4,
    FESetH8 => H8,
    FESetQ8 => Q8,
    FESetQ9 => Q9,
    FESetL3 => L3,
    FESetT6 => T6,
    FESetT10 => T10,
    FESetH20 => H20,
)

numnodesmap = Dict{Int,Int}(
    P1 => 1,
    L2 => 2,
    T3 => 3,
    Q4 => 4,
    T4 => 4,
    H8 => 8,
    Q8 => 8,
    Q9 => 9,
    L3 => 3,
    T6 => 6,
    T10 => 10,
    H20 => 20,
)

"""
    vtkexportmesh(theFile::String, fens::FENodeSet, fes::T; opts...) where {T<:AbstractFESet}

Export mesh to a VTK 1.0 file as an unstructured grid.

Arguments:
- `theFile` = file name,
- `fens` = finite element node set,
- `fes` = finite element set,
- `opts` = keyword argument list, where
    + `scalars` = array of tuples, (name, data)
    + `vectors` = array of tuples, (name, data)

For the `scalars`: If `data` is a vector, that data is exported as a single field.
On the other hand, if it is an 2d array, each column is exported  as a separate field.

"""
function vtkexportmesh(
    theFile::String,
    fens::FENodeSet,
    fes::T;
    opts...,
) where {T<:AbstractFESet}
    Cell_type = _match_type(VTKtypemap, typeof(fes))
    Cell_type === nothing && error("Unknown type $(typeof(fes))")
    return vtkexportmesh(theFile, connasarray(fes), fens.xyz, Cell_type; opts...)
end

"""
    vtkexportmesh(theFile::String, Connectivity, Points, Cell_type;
        vectors=nothing, scalars=nothing)

Export mesh to a VTK 1.0 file as an unstructured grid.

Arguments:
- `theFile` = file name,
- `Connectivity` = array of connectivities, one row per element,
- `Points` = array of node coordinates, one row per node,
- `Cell_type` = type of the cell, refer to the predefined constants P1, L2, ..., H20, ...
- `scalars` = array of tuples, (name, data)
- `vectors` = array of tuples, (name, data)

For the `scalars`: If `data` is a vector, that data is exported as a single field.
On the other hand, if it is an 2d array, each column is exported  as a separate field.
"""
function vtkexportmesh(
    theFile::String,
    Connectivity,
    Points,
    Cell_type;
    vectors = nothing,
    scalars = nothing,
)
    X = Points
    if size(Points, 2) < 3
        X = zeros(size(Points, 1), 3)
        for j in axes(Points, 1)
            X[j, 1:size(Points, 2)] = Points[j, :]
        end
    end
    numnodes = get(() -> error("Wrong number of connected nodes!"), numnodesmap, Cell_type)
    if typeof(Connectivity[1]) <: Tuple # Vector of connectivity tuples: convert to an array
        c = fill(0, length(Connectivity), length(Connectivity[1]))
        for i in eachindex(Connectivity)
            c[i, :] = [Connectivity[i]...]
        end
        Connectivity = c
    end
    (numnodes == size(Connectivity, 2)) || error("Wrong number of connected nodes!")

    fid = open(theFile, "w")
    if (fid == -1)
        error(["Could not open " * theFile])
        return nothing
    end
    print(fid, "# vtk DataFile Version 1.0\n")
    print(fid, "FinEtools Export\n")
    print(fid, "ASCII\n")
    print(fid, "\n")
    print(fid, "DATASET UNSTRUCTURED_GRID\n")
    print(fid, "POINTS ", size(X, 1), " double\n")
    for i in axes(X, 1)
        for j = 1:(size(X, 2)-1)
            print(fid, X[i, j], " ")
        end
        print(fid, X[i, end], "\n")
    end
    print(fid, "\n")
    print(fid, "\n")

    print(
        fid,
        "CELLS ",
        size(Connectivity, 1),
        " ",
        (size(Connectivity, 1) * (size(Connectivity, 2) + 1)),
        "\n",
    )
    for i in axes(Connectivity, 1)
        print(fid, size(Connectivity, 2), " ")
        for j = 1:(size(Connectivity, 2)-1)
            print(fid, Connectivity[i, j] - 1, " ")
        end
        print(fid, Connectivity[i, end] - 1, "\n")
    end
    print(fid, "\n")
    print(fid, "\n")
    print(fid, "CELL_TYPES ", size(Connectivity, 1), "\n")
    for i in axes(Connectivity, 1)
        print(fid, Cell_type, "\n")
    end
    print(fid, "\n")
    print(fid, "\n")

    did_point_data = false

    # First try to write point data
    if (scalars !== nothing)
        did_point_data = false
        for sx in eachindex(scalars)
            name, data = scalars[sx]
            if (size(data, 1) == size(Points, 1)) # point data
                if (!did_point_data)
                    print(fid, "POINT_DATA ", size(data, 1), "\n")
                    did_point_data = true
                end
                if size(data, 2) > 1 # there are multiple scalar fields here
                    for i in axes(data, 2)
                        print(fid, "SCALARS ", name * "$(i)", " double\n")
                        print(fid, "LOOKUP_TABLE default\n")
                        for j in axes(data, 1)
                            print(fid, data[j, i], "\n")
                        end
                    end
                else  # there's just one scalar field here
                    print(fid, "SCALARS ", name, " double\n")
                    print(fid, "LOOKUP_TABLE default\n")
                    for j in axes(data, 1)
                        print(fid, data[j], "\n")
                    end
                end
                print(fid, "\n")
            elseif (size(data, 1) == size(Connectivity, 1)) # cell data
                # Handled below, in the section for cell data
            end
        end
    end

    print(fid, "\n")

    if vectors !== nothing
        for vx in eachindex(vectors)
            name, data = vectors[vx]
            if (!did_point_data)
                print(fid, "POINT_DATA ", size(data, 1), "\n")
                did_point_data = true
            end
            print(fid, "VECTORS ", name, " double\n")
            #print(fid,"LOOKUP_TABLE default\n");
            if size(data, 2) < 3
                X = zeros(size(data, 1), 3)
                for j in axes(data, 1)
                    X[j, 1:size(data, 2)] = data[j, :]
                end
            else
                X = data
            end
            for j in axes(X, 1)
                k = 1
                print(fid, X[j, k])
                for k = 2:lastindex(X, 2)
                    print(fid, " ", X[j, k])
                end
                print(fid, "\n")
            end
            print(fid, "\n")
        end
    end
    print(fid, "\n")

    # Are there any cell data?  If so, write out.
    if (scalars !== nothing)
        did_cell_data = false
        for sx in eachindex(scalars)
            name, data = scalars[sx]
            if (size(data, 1) == size(Points, 1)) # point data
            # Handled above for the case of point data
            elseif (size(data, 1) == size(Connectivity, 1)) # cell data
                if (!did_cell_data)
                    print(fid, "CELL_DATA ", size(data, 1), "\n")
                    did_cell_data = true
                end
                if size(data, 2) > 1 # there are multiple scalar fields here
                    for i in axes(data, 2)
                        print(fid, "SCALARS ", name * "$(i)", " double\n")
                        print(fid, "LOOKUP_TABLE default\n")
                        for j in axes(data, 1)
                            print(fid, data[j, i], "\n")
                        end
                    end
                else  # there's just one scalar field here
                    print(fid, "SCALARS ", name, " double\n")
                    print(fid, "LOOKUP_TABLE default\n")
                    for j in axes(data, 1)
                        print(fid, data[j], "\n")
                    end
                end
                print(fid, "\n")
            end
        end
    end

    fid = close(fid)
    return true
end

"""
    vtkexportvectors(theFile::String, Points, vectors)

Export vector data to a VTK 1.0 file.

Arguments:
- `theFile` = file name,
- `Points` = array of collection of coordinates (tuples or vectors), 
- `vectors` = array of tuples, `(name, data)`, where `name` is a string, and
  `data` is array of collection of coordinates (tuples or vectors).

# Example
```
Points = [(1.0, 3.0), (0.0, -1.0)]
vectors = [("v", [(-1.0, -2.0), (1.0, 1.0)])]
vtkexportvectors("theFile.VTK", Points, vectors)
```
will produce file with
```
# vtk DataFile Version 1.0
FinEtools Export
ASCII

DATASET UNSTRUCTURED_GRID
POINTS 2 double
1.0 3.0 0.0
0.0 -1.0 0.0


POINT_DATA 2
VECTORS v double
-1.0 -2.0 0.0
1.0 1.0 0.0
```

!!! note
    The filter "Glyph" must be used within Paraview. Also in the drop-down 
    menu "Glyph mode" select "all points".
"""
function vtkexportvectors(theFile::String, Points, vectors)
    fid = open(theFile, "w")
    if (fid == -1)
        error(["Could not open " * theFile])
        return nothing
    end
    print(fid, "# vtk DataFile Version 1.0\n")
    print(fid, "FinEtools Export\n")
    print(fid, "ASCII\n")
    print(fid, "\n")
    print(fid, "DATASET UNSTRUCTURED_GRID\n")
    print(fid, "POINTS ", length(Points), " double\n")
    pt = fill(0.0, 3)
    for i in eachindex(Points)
        for k in eachindex(Points[i])
            pt[k] = Points[i][k]
        end
        for j = 1:(length(pt)-1)
            print(fid, pt[j], " ")
        end
        print(fid, pt[end], "\n")
    end
    print(fid, "\n")
    print(fid, "\n")

    did_point_data = false

    if vectors !== nothing
        for vx in eachindex(vectors)
            name, data = vectors[vx]
            if (!did_point_data)
                print(fid, "POINT_DATA ", length(data), "\n")
                did_point_data = true
            end
            print(fid, "VECTORS ", name, " double\n")
            #print(fid,"LOOKUP_TABLE default\n");
            for i in eachindex(data)
                for k in eachindex(data[i])
                    pt[k] = data[i][k]
                end
                for j = 1:(length(pt)-1)
                    print(fid, pt[j], " ")
                end
                print(fid, pt[end], "\n")
            end
            print(fid, "\n")
        end
    end
    print(fid, "\n")

    fid = close(fid)
    return true
end

end # VTK

module Abaqus
################################################################################
# Abaqus export
################################################################################
import ...FESetModule:
    AbstractFESet,
    FESetP1,
    FESetL2,
    FESetT3,
    FESetQ4,
    FESetT4,
    FESetH8,
    FESetQ8,
    FESetL3,
    FESetT6,
    FESetT10,
    FESetH20,
    connasarray
import ...FENodeSetModule: FENodeSet
import Base.close
using LinearAlgebra

"""
    AbaqusExporter

Export mesh to Abaqus.
"""
mutable struct AbaqusExporter
    filename::AbstractString
    ios::IO
    element_range::Tuple{Int64,Int64}
    function AbaqusExporter(filename::AbstractString)
        if match(r".*\.inp$", filename) === nothing
            filename = filename * ".inp"
        end
        ios = open(filename, "w+")
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
    println(self.ios, "*HEADING")
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
function NODE(self::AbaqusExporter, xyz::AbstractArray{T,2}) where {T}
    println(self.ios, "*NODE")
    for j in axes(xyz, 1)
        print(self.ios, "$j,$(xyz[j, 1])")
        for ixxxx = 2:1:size(xyz, 2)
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
function ELEMENT(
    self::AbaqusExporter,
    TYPE::AbstractString,
    ELSET::AbstractString,
    start::Integer,
    conn::AbstractArray{T,2},
) where {T<:Integer}
    # Check that start is valid
    (start > 0) || error("The  starting element number must be > 0")
    # Check that the current element number range is disjoint from the
    # range of the elements on input: they must not overlap
    ((self.element_range[2] < start) ||
        (self.element_range[1] > start + size(conn, 1) - 1)) || error("Elements  must be given unique numbers")
    # Update the element range
    self.element_range = (
        min(self.element_range[1], start),
        max(self.element_range[2], start + size(conn, 1) - 1),
    )
    println(self.ios, "*ELEMENT, TYPE =" * TYPE * ", ELSET=" * ELSET)
    for j in axes(conn, 1)
        print(self.ios, "$(j+start-1),")
        for ixxxx = 1:(size(conn, 2)-1)
            if ixxxx > 15
                print(self.ios, "\n")
            end
            print(self.ios, "$(conn[j,ixxxx]),")
        end
        print(self.ios, "$(conn[j,size(conn,2)])\n")
    end
end

function ELEMENT(
    self::AbaqusExporter,
    TYPE::AbstractString,
    ELSET::AbstractString,
    conn::AbstractArray{T,2},
) where {T<:Integer}
    start = self.element_range[2] + 1
    ELEMENT(self, TYPE, ELSET, start, conn)
end

function ELEMENT(
    self::AbaqusExporter,
    TYPE::AbstractString,
    ELSET::AbstractString,
    start::Integer,
    conn::Vector{C},
) where {C<:Tuple}
    c = fill(0, length(conn), length(conn[1]))
    for i in eachindex(conn)
        c[i, :] = [conn[i]...]
    end
    ELEMENT(self, TYPE, ELSET, start, c)
end

function ELEMENT(
    self::AbaqusExporter,
    TYPE::AbstractString,
    ELSET::AbstractString,
    conn::Vector{C},
) where {C<:Tuple}
    start = self.element_range[2] + 1
    ELEMENT(self, TYPE, ELSET, start, conn)
end

"""
    NSET_NSET(self::AbaqusExporter, NSET::AbstractString,
      n::AbstractVector{T}) where {T<:Integer}

Write out the `*NSET` option.

`NSET` = name of the set,
`n` = array of the node numbers
"""
function NSET_NSET(
    self::AbaqusExporter,
    NSET::AbstractString,
    n::AbstractVector{T},
) where {T<:Integer}
    println(self.ios, "*NSET, NSET=" * NSET)
    for j = 1:(length(n)-1)
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
function ELSET_ELSET(
    self::AbaqusExporter,
    ELSET::AbstractString,
    n::AbstractArray{T,1},
) where {T<:Integer}
    println(self.ios, "*ELSET, ELSET=" * ELSET)
    for j = 1:(length(n)-1)
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
function ORIENTATION(
    self::AbaqusExporter,
    ORIENTATION::AbstractString,
    a::AbstractArray{T,1},
    b::AbstractArray{T,1},
) where {T<:Real}
    println(self.ios, "*ORIENTATION,NAME=" * ORIENTATION)
    println(self.ios, "$(a[1]),$(a[2]),$(a[3]),$(b[1]),$(b[2]),$(b[3])")
    println(self.ios, "1,0.0")
end

"""
    MATERIAL(self::AbaqusExporter, MATERIAL::AbstractString)

Write out the `*MATERIAL` option.
"""
function MATERIAL(self::AbaqusExporter, MATERIAL::AbstractString)
    println(self.ios, "*MATERIAL,NAME=" * MATERIAL)
end

"""
    ELASTIC(self::AbaqusExporter, E1::F, E2::F, E3::F, nu12::F, nu13::F, nu23::F,
      G12::F, G13::F, G23::F) where {F}

Write out the `*ELASTIC,TYPE=ENGINEERING CONSTANTS` option.
"""
function ELASTIC(
    self::AbaqusExporter,
    E1::F,
    E2::F,
    E3::F,
    nu12::F,
    nu13::F,
    nu23::F,
    G12::F,
    G13::F,
    G23::F,
) where {F}
    println(self.ios, "*ELASTIC,TYPE=ENGINEERING CONSTANTS")
    println(self.ios, "$E1,$E2,$E3,$nu12,$nu13,$nu23,$G12,$G13,")
    println(self.ios, "$G23")
end

"""
    ELASTIC(self::AbaqusExporter, E::F, nu::F) where {F}

Write out the `*ELASTIC,TYPE=ISOTROPIC` option.
"""
function ELASTIC(self::AbaqusExporter, E::F, nu::F) where {F}
    println(self.ios, "*ELASTIC,TYPE=ISOTROPIC ")
    println(self.ios, "$E,$nu")
end

"""
    EXPANSION(self::AbaqusExporter, CTE::F) where {F}

Write out the `*EXPANSION` option.
"""
function EXPANSION(self::AbaqusExporter, CTE::F) where {F}
    println(self.ios, "*EXPANSION")
    println(self.ios, "$CTE,")
end

"""
    DENSITY(self::AbaqusExporter, rho)

Write out the `*DENSITY` option.
"""
function DENSITY(self::AbaqusExporter, rho)
    println(self.ios, "*DENSITY ")
    println(self.ios, "$rho")
end

"""
    SECTION_CONTROLS(self::AbaqusExporter, NAME::AbstractString,
      OPTIONAL::AbstractString)

Write out the `*SECTION CONTROLS` option.

`OPTIONAL` = string, for instance
                     HOURGLASS=ENHANCED
"""
function SECTION_CONTROLS(
    self::AbaqusExporter,
    NAME::AbstractString,
    OPTIONAL::AbstractString,
)
    println(self.ios, "*SECTION CONTROLS, NAME=" * NAME * "," * OPTIONAL)
end

"""
    SOLID_SECTION(self::AbaqusExporter, MATERIAL::AbstractString,
      ORIENTATION::AbstractString, ELSET::AbstractString,
      CONTROLS::AbstractString)

Write out the `*SOLID SECTION` option.

Level: Part,  Part instance

"""
function SOLID_SECTION(
    self::AbaqusExporter,
    MATERIAL::AbstractString,
    ORIENTATION::AbstractString,
    ELSET::AbstractString,
    CONTROLS::AbstractString,
)
    println(
        self.ios,
        "*SOLID SECTION,MATERIAL=" *
        MATERIAL *
        ",ORIENTATION =" *
        ORIENTATION *
        ",ELSET=" *
        ELSET *
        ", CONTROLS  =" *
        CONTROLS,
    )
end

"""
    SOLID_SECTION(self::AbaqusExporter, MATERIAL::AbstractString,
      ORIENTATION::AbstractString, ELSET::AbstractString)

Write out the `*SOLID SECTION` option.

Level: Part,  Part instance

"""
function SOLID_SECTION(
    self::AbaqusExporter,
    MATERIAL::AbstractString,
    ORIENTATION::AbstractString,
    ELSET::AbstractString,
)
    println(
        self.ios,
        "*SOLID SECTION,MATERIAL=" *
        MATERIAL *
        ",ORIENTATION =" *
        ORIENTATION *
        ",ELSET=" *
        ELSET,
    )
end

"""
    SOLID_SECTION(self::AbaqusExporter, MATERIAL::AbstractString,
      ORIENTATION::AbstractString, ELSET::AbstractString)

Write out the `*SOLID SECTION` option.

Level: Part,  Part instance

"""
function SOLID_SECTION(
    self::AbaqusExporter,
    MATERIAL::AbstractString,
    ORIENTATION::AbstractString,
    ELSET::AbstractString,
    thickness::F,
) where {F}
    SOLID_SECTION(self, MATERIAL, ORIENTATION, ELSET)
    println(self.ios, "$(thickness),")
end

# """
#     HOURGLASS(self::AbaqusExporter, KIND::AbstractString, VALUE::F) where {F}
#
# Write out the `*HOURGLASS` option.
#
# ```
# Example:
# *SOLID SECTION,ELSET=SOLID3,MATERIAL=MAT,CONTROL=Ac
# *HOURGLASS STIFFNESS
# 5.E8
# ```
# """
# function HOURGLASS(self::AbaqusExporter, KIND::AbstractString, VALUE::F) where {F}
#   println(self.ios, "*HOURGLASS " * KIND);
#   println(self.ios, "$( VALUE )");
# end

"""
    SURFACE_SECTION(self::AbaqusExporter, ELSET::AbstractString)

Write out the `*SURFACE SECTION` option.
"""
function SURFACE_SECTION(self::AbaqusExporter, ELSET::AbstractString)
    println(self.ios, "*SURFACE SECTION, ELSET=" * ELSET)
end

"""
    STEP_PERTURBATION_STATIC(self::AbaqusExporter)

Write out the `*STEP,PERTURBATION` option for linear static analysis.
"""
function STEP_PERTURBATION_STATIC(self::AbaqusExporter)
    println(self.ios, "*STEP,PERTURBATION")
    println(self.ios, "*STATIC")
end

"""
    STEP_PERTURBATION_BUCKLE(self::AbaqusExporter, neigv::Integer)

Write out the `*STEP,PERTURBATION` option for linear buckling analysis.
"""
function STEP_PERTURBATION_BUCKLE(self::AbaqusExporter, neigv::Integer)
    println(self.ios, "*STEP, name=Buckling, nlgeom=NO, perturbation")
    println(self.ios, "*BUCKLE")
    println(self.ios, "$(neigv), , , , ")
end

"""
    STEP_FREQUENCY(self::AbaqusExporter, Nmodes::Integer)

Write out the `*STEP,FREQUENCY` option.
"""
function STEP_FREQUENCY(self::AbaqusExporter, Nmodes::Integer)
    println(self.ios, "*STEP")
    println(self.ios, "*FREQUENCY, EIGENSOLVER=LANCZOS")
    println(self.ios, "$(Nmodes), , ,-1.E6 \n")
end

"""
    BOUNDARY(self::AbaqusExporter, mesh, nodes, is_fixed::AbstractArray{B,2}, fixed_value::AbstractArray{F,2}) where {B, F}

Write out the `*BOUNDARY` option. 

The boundary condition is applied to the nodes specified in 
the array `nodes`, in the mesh (or node set) `meshornset`.

`meshornset` = mesh or node set (default is empty)
`nodes` = array of node numbers, the node numbers are attached to the mesh label,
`is_fixed`= array of Boolean flags (true means fixed, or prescribed),  one row per node,
`fixed_value`=array of displacements to which the corresponding displacement components is fixed

# Example

```
BOUNDARY(AE, "ASSEM1.INSTNC1", 1:4, fill(true, 4, 1), reshape([uy(fens.xyz[i, :]...) for i in 1:4], 4, 1))
```
"""
function BOUNDARY(
    self::AbaqusExporter,
    meshornset,
    nodes,
    is_fixed::AbstractArray{B,2},
    fixed_value::AbstractArray{F,2},
) where {B,F}
    println(self.ios, "*BOUNDARY")
    if meshornset == ""
        meshlabel = ""
    else
        meshlabel = meshornset * "."
    end
    for j in nodes
        for k in axes(is_fixed, 2)
            #<node number>, <first dof>, <last dof>, <magnitude of displacement>
            if is_fixed[j, k]
                println(self.ios, "$(meshlabel)$(j),$k,$k,$(fixed_value[j,k])")
            end
        end
    end
end

function BOUNDARY(
    self::AbaqusExporter,
    nodes,
    is_fixed::AbstractVector{B},
    fixed_value::AbstractVector{F},
) where {B,F}
    BOUNDARY(
        self,
        nodes,
        reshape(is_fixed, length(is_fixed), 1),
        reshape(fixed_value, length(fixed_value), 1),
    )
end

"""
    BOUNDARY(self::AbaqusExporter, meshornset, is_fixed::AbstractArray{B,2}, fixed_value::AbstractArray{F,2}) where {B, F}

Write out the `*BOUNDARY` option.

- `meshornset` = name of a mesh or a node set,
- `is_fixed`= array of Boolean flags (true means fixed, or prescribed),  one row per node, as many columns as there are degrees of freedom per node,
- `fixed_value`=array of displacements to which the corresponding displacement components is fixed, as many columns as there are degrees of freedom per node
"""
function BOUNDARY(
    self::AbaqusExporter,
    meshornset,
    is_fixed::AbstractArray{B,2},
    fixed_value::AbstractArray{F,2},
) where {B,F}
    BOUNDARY(self, meshornset, 1:size(is_fixed, 1), is_fixed, fixed_value)
end

"""
    BOUNDARY(self::AbaqusExporter, NSET::AbstractString, dof::Integer,  fixed_value)

Write out the `*BOUNDARY` option.

- `NSET` = name of a node set,
- `is_fixed`= array of Boolean flags (true means fixed, or prescribed),  one row per node,
- `fixed_value`=array of displacements to which the corresponding displacement components is fixed
"""
function BOUNDARY(self::AbaqusExporter, NSET::AbstractString, dof::Integer, fixed_value)
    println(self.ios, "*BOUNDARY")
    println(self.ios, NSET * ",$dof,$dof,$(fixed_value)")
end

"""
    BOUNDARY(self::AbaqusExporter, NSET::AbstractString, dof::Integer)

Write out the `*BOUNDARY` option to fix displacements at zero for a node set.

Invoke at Level: Model,  Step

- `NSET`= node set,
- `dof`=Degree of freedom, 1, 2, 3
"""
function BOUNDARY(self::AbaqusExporter, NSET::AbstractString, dof::Integer)
    BOUNDARY(self, NSET, dof, 0.0)
end

"""
    BOUNDARY(self::AbaqusExporter, NSET::AbstractString, dof::Integer,
      value::F) where {F}

Write out the `*BOUNDARY` option to fix displacements at nonzero value for a
node set.

- `NSET`= node set,
- `dof`=Degree of freedom, 1, 2, 3
- `typ` = DISPLACEMENT
"""
function BOUNDARY(
    self::AbaqusExporter,
    NSET::AbstractString,
    dof::Integer,
    value::F,
    typ::AbstractString,
) where {F}
    println(self.ios, "*BOUNDARY,TYPE=$(typ)")
    println(self.ios, NSET * ",$dof,$dof,$value")
end

"""
    DLOAD(self::AbaqusExporter, ELSET::AbstractString,
      traction::AbstractVector{F}) where {F}

Write out the `*DLOAD` option.
"""
function DLOAD(
    self::AbaqusExporter,
    ELSET::AbstractString,
    traction::AbstractVector{F},
) where {F}
    println(self.ios, "*DLOAD, follower=NO")
    nt = norm(traction)
    print(self.ios, ELSET * ",TRVEC," * "$(nt),")
    println(self.ios, "$(traction[1]/nt),$(traction[2]/nt),$(traction[3]/nt)")
end

"""
    CLOAD(self::AbaqusExporter, NSET::AbstractString, dof::Integer,
      magnitude::F) where {F}

Write out the `*CLOAD` option.

NSET=Node set
dof= 1, 2, 3,
magnitude= signed multiplier
"""
function CLOAD(
    self::AbaqusExporter,
    NSET::AbstractString,
    dof::Integer,
    magnitude::F,
) where {F}
    println(self.ios, "*CLOAD")
    println(self.ios, NSET * ",$dof,$magnitude")
end

"""
    CLOAD(self::AbaqusExporter, nodenumber::Integer, dof::Integer,
      magnitude::F) where {F}

Write out the `*CLOAD` option.

nodenumber=Number of node
dof= 1, 2, 3,
magnitude= signed multiplier
"""
function CLOAD(
    self::AbaqusExporter,
    nodenumber::Integer,
    dof::Integer,
    magnitude::F,
) where {F}
    println(self.ios, "*CLOAD")
    println(self.ios, "$(nodenumber),$dof,$magnitude")
end

"""
    TEMPERATURE(self::AbaqusExporter, nlist::AbstractArray{I, 1},
      tlist::AbstractArray{F, 1}) where {I, F}

Write out the `*TEMPERATURE` option.
"""
function TEMPERATURE(
    self::AbaqusExporter,
    Classifier::AbstractString,
    nlist::AbstractArray{I,1},
    tlist::AbstractArray{F,1},
) where {I,F}
    (length(nlist) == length(tlist)) || error("The number of nodes and temperatures must be the same")
    println(self.ios, "*TEMPERATURE")
    for ixxxx in eachindex(nlist)
        println(self.ios, Classifier * "$(nlist[ixxxx]),$(tlist[ixxxx])")
    end
end

"""
    NODE_PRINT(self::AbaqusExporter, NSET::AbstractString)

Write out the `*NODE PRINT` option.
"""
function NODE_PRINT(self::AbaqusExporter, NSET::AbstractString; KEYS = "U")
    println(self.ios, "*NODE PRINT, NSET=" * NSET)
    println(self.ios, KEYS)
end

"""
    EL_PRINT(self::AbaqusExporter, ELSET::AbstractString, KEYS::AbstractString)

Write out the `*EL PRINT` option.
"""
function EL_PRINT(self::AbaqusExporter, ELSET::AbstractString, KEYS::AbstractString)
    println(
        self.ios,
        "*EL PRINT, ELSET=" * ELSET * ", POSITION=INTEGRATION POINTS, SUMMARY= YES",
    )
    println(self.ios, "")
    println(self.ios, KEYS)
end

"""
    ENERGY_PRINT(self::AbaqusExporter)

Write out the `*ENERGY PRINT` option.
"""
function ENERGY_PRINT(self::AbaqusExporter)
    println(self.ios, "*ENERGY PRINT")
end

"""
    END_STEP(self::AbaqusExporter)

Write out the `*END STEP` option.
"""
function END_STEP(self::AbaqusExporter)
    println(self.ios, "*END STEP")
end

"""
    close(self::AbaqusExporter)

Close  the stream opened by the exporter.
"""
function close(self::AbaqusExporter)
    close(self.ios)
end

end # Abaqus

module CSV
################################################################################
# CSV export
################################################################################
using DataDrop

"""
    savecsv(name::String; kwargs...)

Save arrays as a CSV file.

Example:
```julia
savecsv("ab", a = rand(3), b = rand(3))
```
"""
function savecsv(name::String; kwargs...)
    colnames = Symbol[]
    columns = Any[]
    for apair in pairs(kwargs)
        sy, val = apair
        push!(colnames, sy)
        push!(columns, val)
    end
    ncol = length(colnames)
    nrow = length(columns[1])
    for j = 1:ncol
        (length(columns[j]) == nrow) || error("Columns must have the same number of rows")
    end
    if DataDrop.file_extension(name) == ""
        name = name * ".csv"
    end
    open(name, "w") do fid
        for (j, c) in enumerate(colnames)
            print(fid, "\"", c, "\"")
            j < ncol && print(fid, ",")
        end
        print(fid, "\n")
        for r = 1:nrow
            for j = 1:ncol
                print(fid, columns[j][r])
                j < ncol && print(fid, ",")
            end
            print(fid, "\n")
        end
        print(fid, "\n")
    end
    return true
end

end # CSV

module NASTRAN
################################################################################
# NASTRAN export
################################################################################
import ...FESetModule:
    AbstractFESet,
    FESetP1,
    FESetL2,
    FESetT3,
    FESetQ4,
    FESetT4,
    FESetH8,
    FESetQ8,
    FESetL3,
    FESetT6,
    FESetT10,
    FESetH20,
    connasarray
import ...FENodeSetModule: FENodeSet
import Base.close
using Printf

"""
NASTRANExporter

Exporter of the mesh to NASTRAN.
"""
mutable struct NASTRANExporter
    filename::AbstractString
    ios::IO
    element_range::Tuple{Int64,Int64}
    function NASTRANExporter(filename::AbstractString)
        if match(r".*\.nas$", filename) === nothing
            filename = filename * ".nas"
        end
        ios = open(filename, "w+")
        return new(deepcopy(filename), ios, (typemax(Int64), 0))
    end
end

"""
    CEND(self::NASTRANExporter)

Terminate the Executive Control section.
"""
function CEND(self::NASTRANExporter)
    println(self.ios, "CEND")
end

"""
    BEGIN_BULK(self::NASTRANExporter)

Terminate the Case Control section by starting the bulk section.
"""
function BEGIN_BULK(self::NASTRANExporter)
    println(self.ios, "BEGIN BULK")
end

"""
    ENDDATA(self::NASTRANExporter)

Terminate the bulk section.
"""
function ENDDATA(self::NASTRANExporter)
    println(self.ios, "ENDDATA")
end

"""
    GRID(self::NASTRANExporter, n::Int, xyz)

Write a grid-point statement.
"""
function GRID(self::NASTRANExporter, n::Int, xyz)
    @printf self.ios "GRID,%d,,%g,%g,%g\n" n xyz[1] xyz[2] xyz[3]
end

"""
    PSOLID(self::NASTRANExporter, pid::Int, mid::Int)

Write solid-property statement.
"""
function PSOLID(self::NASTRANExporter, pid::Int, mid::Int)
    @printf self.ios "PSOLID,%d,%d\n" pid mid
end

"""
    MAT1(
        self::NASTRANExporter,
        mid::Int,
        E::T,
        G::T,
        nu::T,
        rho::T,
        A::T,
        TREF::T,
        GE::T,
    ) where {T}

Write a statement for an isotropic elastic material.
"""
function MAT1(
    self::NASTRANExporter,
    mid::Int,
    E::T,
    G::T,
    nu::T,
    rho::T,
    A::T,
    TREF::T,
    GE::T,
) where {T}
    @printf self.ios "MAT1,%d,%g,%g,%g,%g,%g,%g,%g\n" mid E G nu rho A TREF GE
end

"""
    MAT1(
        self::NASTRANExporter,
        mid::Int,
        E::T,
        nu::T,
        rho::T = 0.0,
        A::T = 0.0,
        TREF::T = 0.0,
        GE::T = 0.0,
    ) where {T}

Write a statement for an isotropic elastic material.
"""
function MAT1(
    self::NASTRANExporter,
    mid::Int,
    E::T,
    nu::T,
    rho::T = 0.0,
    A::T = 0.0,
    TREF::T = 0.0,
    GE::T = 0.0,
) where {T}
    @printf self.ios "MAT1,%d,%g,,%g,%g,%g,%g,%g\n" mid E nu rho A TREF GE
end

"""
    CTETRA(self::NASTRANExporter, eid::Int, pid::Int, conn::Vector{Int})

Write a statement for a single tetrahedron element.
"""
function CTETRA(self::NASTRANExporter, eid::Int, pid::Int, conn::Vector{Int})
    nc = length(conn)
    @printf self.ios "CTETRA,%d,%d" eid pid
    for k = 1:nc
        @printf self.ios ",%d" conn[k]
        if k == 6
            @printf self.ios "\n"
        end
    end
    @printf self.ios "\n"
end

"""
    close(self::NASTRANExporter)

Close  the stream opened by the exporter.
"""
function close(self::NASTRANExporter)
    close(self.ios)
end

end # NASTRAN

module STL
################################################################################
# STL export
################################################################################
using Printf
import LinearAlgebra: norm, cross
import Base.close

"""
STLExporter

Exporter of a surface mesh as an STL file.
"""
mutable struct STLExporter
    filename::AbstractString
    ios::IO
    element_range::Tuple{Int64,Int64}
    function STLExporter(filename::AbstractString)
        if match(r".*\.stl$", filename) === nothing
            filename = filename * ".stl"
        end
        ios = open(filename, "w+")
        return new(deepcopy(filename), ios, (typemax(Int64), 0))
    end
end

"""
    solid(self::STLExporter, name::AbstractString = "thesolid")

Write a statement to begin the solid.
"""
function solid(self::STLExporter, name::AbstractString = "thesolid")
    @printf self.ios "solid %s\n" name
end

"""
    facet(self::STLExporter, v1::Vector{T}, v2::Vector{T}, v3::Vector{T}) where {T}

Write a single facet.
"""
function facet(self::STLExporter, v1::Vector{T}, v2::Vector{T}, v3::Vector{T}) where {T}
    V = v2 - v1
    W = v3 - v1
    normal = cross(V, W)
    normal = normal / norm(normal)
    @printf self.ios "facet normal %f %f %f\n" normal[1] normal[2] normal[3]
    @printf self.ios "    outer loop\n"
    @printf self.ios "        vertex %f %f %f\n" v1[1] v1[2] v1[3]
    @printf self.ios "        vertex %f %f %f\n" v2[1] v2[2] v2[3]
    @printf self.ios "        vertex %f %f %f\n" v3[1] v3[2] v3[3]
    @printf self.ios "    endloop\n"
    @printf self.ios "endfacet\n"
end

"""
    endsolid(self::STLExporter, name::AbstractString = "thesolid")

Write statement to end the solid.
"""
function endsolid(self::STLExporter, name::AbstractString = "thesolid")
    @printf self.ios "endsolid %s\n" name
end

"""
    close(self::STLExporter)

Close  the stream opened by the exporter.
"""
function close(self::STLExporter)
    close(self.ios)
end

end # STL

module H2Lib
################################################################################
# H2Lib triangular-surface export
################################################################################
using Printf
import Base.BitSet

mutable struct _HypFace
    n::Int # serial number of the hyperface
    hf::Array{Int} # numbers of vertices defining the hyperface
    a::Int # anchor node number
    o::BitSet # numbers of the other nodes on the hyperface
end

mutable struct _HypFaceContainer
    nv::Int # total number of vertices per hyperface
    d::Dict{Int,Array{_HypFace}} # dictionary
    n::Int # total number of hyperfaces in the container
end

makehypfacedict(nv) = _HypFaceContainer(nv, Dict{Int,Array{_HypFace}}(), 0)

function hypfacen(container, hypf)
    h = sort([i for i in hypf])
    anchor = h[1]
    other = BitSet(h[2:end])
    ca = get(container.d, anchor, _HypFace[])
    for k in eachindex(ca)
        if (ca[k].a == anchor) && (ca[k].o == other)
            return ca[k].n # do have this hyper face: here is its serial number
        end
    end
    return 0 # do not have this hyper face
end

function inserthypface!(container, hypf)
    if hypfacen(container, hypf) == 0
        (length(hypf) == container.nv) || error("Hyperface does not have the right number of vertices")
        h = sort([i for i in hypf])
        anchor = h[1]
        other = BitSet(h[2:end])
        ca = get(container.d, anchor, _HypFace[])
        container.n = container.n + 1
        push!(ca, _HypFace(container.n, hypf, anchor, other))
        container.d[anchor] = ca
    end
    return container
end

function hypfacedicttoarray(container)
    a = fill(0, container.n, container.nv)
    for ca in values(container.d)
        for i in eachindex(ca)
            a[ca[i].n, :] = ca[i].hf
        end
    end
    return a
end

"""
    h2libexporttri(theFile::String, Connectivity, Points)

Write a file in the H2Lib format.
"""
function h2libexporttri(theFile::String, Connectivity, Points)
    (size(Connectivity, 2) == 3) || error("Only triangles accepted")
    adjust = 1 # the library is C-based and expects indexes to be zero-based

    ed = makehypfacedict(2)
    for i in axes(Connectivity, 1)
        inserthypface!(ed, Connectivity[i, [1, 2]])
        inserthypface!(ed, Connectivity[i, [2, 3]])
        inserthypface!(ed, Connectivity[i, [3, 1]])
    end
    ea = hypfacedicttoarray(ed)

    fid = open(theFile, "w")
    if (fid == -1)
        error(["Could not open " * theFile])
        return nothing
    end

    print(fid, size(Points, 1), " ", size(ea, 1), " ", size(Connectivity, 1), "\n")

    for i in axes(Points, 1)
        for j = 1:(size(Points, 2)-1)
            print(fid, Points[i, j], " ")
        end
        print(fid, Points[i, end], "\n")
    end

    for i in axes(ea, 1)
        for j = 1:(size(ea, 2)-1)
            print(fid, ea[i, j] - adjust, " ")
        end
        print(fid, ea[i, end] - adjust, "\n")
    end

    for i in axes(Connectivity, 1)
        s =
            [hypfacen(ed, Connectivity[i, [2, 3]]) hypfacen(ed, Connectivity[i, [3, 1]]) hypfacen(
                ed,
                Connectivity[i, [1, 2]],
            )]
        for j in axes(Connectivity, 2)
            print(fid, Connectivity[i, j] - adjust, " ")
        end
        for j = 1:(length(s)-1)
            print(fid, s[j] - adjust, " ")
        end
        print(fid, s[end] - adjust, "\n")
    end

    fid = close(fid)
    return true
end

end # H2Lib

module MESH
################################################################################
# MESH export
################################################################################
using DelimitedFiles
import ...FESetModule:
    AbstractFESet,
    FESetP1,
    FESetL2,
    FESetT3,
    FESetQ4,
    FESetT4,
    FESetH8,
    FESetQ8,
    FESetQ9,
    FESetL3,
    FESetT6,
    FESetT10,
    FESetH20,
    connasarray
import ...FENodeSetModule: FENodeSet
using ..MeshExportModule: _match_type

MESHtypemap = Dict{UnionAll,String}(
    FESetP1 => "P1",
    FESetL2 => "L2",
    FESetT3 => "T3",
    FESetQ4 => "Q4",
    FESetT4 => "T4",
    FESetH8 => "H8",
    FESetQ8 => "Q8",
    FESetQ9 => "Q9",
    FESetL3 => "L3",
    FESetT6 => "T6",
    FESetT10 => "T10",
    FESetH20 => "H20",
)

"""
    write_MESH(meshfile::String, fens::FENodeSet, fes::T) where {T<:AbstractFESet}

Write the mesh in the MESH format.
"""
function write_MESH(meshfile::String, fens::FENodeSet, fes::T) where {T<:AbstractFESet}
    meshfilebase, ext = splitext(meshfile)
    t = _match_type(MESHtypemap, typeof(fes))
    t === nothing && error("Unknown type $(typeof(fes))")
    dinfo = [meshfilebase * "-xyz.dat", t, meshfilebase * "-conn.dat"]
    # write out a file with the coordinates of the nodes
    open(dinfo[1], "w") do file
        writedlm(file, fens.xyz, ' ')
    end
    # write out a file with the connectivity
    open(dinfo[3], "w") do file
        writedlm(file, connasarray(fes), ' ')
    end
    # if any label is different from the default (0), write out the labels
    if any(i -> i != 0, fes.label)
        push!(dinfo, meshfilebase * "-label.dat")
        # write out a file with the labels
        open(dinfo[4], "w") do file
            writedlm(file, fes.label, ' ')
        end
    end

    # write out a file with the metadata
    open(meshfilebase * ".mesh", "w") do file
        writedlm(file, dinfo)
    end
    return true
end

end

module H5MESH
################################################################################
# HDF5 mesh export
################################################################################
using HDF5
using DataDrop
import ...FESetModule:
    AbstractFESet,
    FESetP1,
    FESetL2,
    FESetT3,
    FESetQ4,
    FESetT4,
    FESetH8,
    FESetQ8,
    FESetQ9,
    FESetL3,
    FESetT6,
    FESetT10,
    FESetH20,
    connasarray
import ...FENodeSetModule: FENodeSet
import ..VTK: VTKtypemap
using ..MeshExportModule: _match_type

"""
    write_H5MESH(meshfile::String, fens::FENodeSet, fes::T) where {T<:AbstractFESet}

Write the mesh in the H5MESH format.

The mesh is stored in a HDF5 file.
"""
function write_H5MESH(meshfile::String, fens::FENodeSet, fes::T) where {T<:AbstractFESet}
    meshfilebase, ext = splitext(meshfile)
    if ext == ""
        ext = ".h5mesh"
    end
    fname = DataDrop.with_extension(meshfile, ext)
    etype = _match_type(VTKtypemap, typeof(fes))
    etype === nothing && error("Unknown VTK type!")
    # If the file exists, delete all contents
    h5open(fname, "w") do fid
    end
    # Store the mesh data
    DataDrop.store_matrix(fname, "xyz", fens.xyz)
    DataDrop.store_matrix(fname, "label", fes.label)
    C = connasarray(fes)
    mn = maximum(C)
    if mn < typemax(Int32)
        C = Int32.(C)
    end
    DataDrop.store_matrix(fname, "conn", C)
    DataDrop.store_value(fname, "etype", etype)
    return true
end

end # H5MESH

module VTKWrite
################################################################################
# VTK export using WriteVTK. Unstructured, binary, file is written.
################################################################################

using WriteVTK
using DataDrop
import ...FESetModule:
    AbstractFESet,
    FESetP1,
    FESetL2,
    FESetT3,
    FESetQ4,
    FESetT4,
    FESetH8,
    FESetQ8,
    FESetQ9,
    FESetL3,
    FESetT6,
    FESetT10,
    FESetH20,
    connasarray
import ...FENodeSetModule: FENodeSet
using ..MeshExportModule: _match_type

using Printf
import LinearAlgebra: norm, cross

const _VTK_TYPE_MAP = Dict{UnionAll,Int}(
    FESetP1 => 1,
    FESetL2 => 3,
    FESetT3 => 5,
    FESetQ4 => 9,
    FESetT4 => 10,
    FESetH8 => 12,
    FESetQ8 => 23,
    FESetQ9 => 28,
    FESetL3 => 21,
    FESetT6 => 22,
    FESetT10 => 24,
    FESetH20 => 25,
)

const P1 = WriteVTK.VTKCellTypes.VTKCellType(_VTK_TYPE_MAP[FESetP1])
const L2 = WriteVTK.VTKCellTypes.VTKCellType(_VTK_TYPE_MAP[FESetL2])
const T3 = WriteVTK.VTKCellTypes.VTKCellType(_VTK_TYPE_MAP[FESetT3])
const Q4 = WriteVTK.VTKCellTypes.VTKCellType(_VTK_TYPE_MAP[FESetQ4])
const T4 = WriteVTK.VTKCellTypes.VTKCellType(_VTK_TYPE_MAP[FESetT4])
const H8 = WriteVTK.VTKCellTypes.VTKCellType(_VTK_TYPE_MAP[FESetH8])
const L3 = WriteVTK.VTKCellTypes.VTKCellType(_VTK_TYPE_MAP[FESetL3])
const T6 = WriteVTK.VTKCellTypes.VTKCellType(_VTK_TYPE_MAP[FESetT6])
const Q8 = WriteVTK.VTKCellTypes.VTKCellType(_VTK_TYPE_MAP[FESetQ8])
const Q9 = WriteVTK.VTKCellTypes.VTKCellType(_VTK_TYPE_MAP[FESetQ9])
const T10 = WriteVTK.VTKCellTypes.VTKCellType(_VTK_TYPE_MAP[FESetT10])
const H20 = WriteVTK.VTKCellTypes.VTKCellType(_VTK_TYPE_MAP[FESetH20])

"""
    vtkwrite(theFile::String, fens::FENodeSet, fes::T; opts...) where {T<:AbstractFESet}

Export mesh to VTK as an unstructured grid (binary file).

Arguments:
- `theFile` = file name,
- `fens` = finite element node set,
- `fes` = finite element set,
- `opts` = keyword argument list, where `scalars` = array of tuples, 
    `(name, data)`, `vectors` = array of tuples, `(name, data)`

For the `scalars`: If `data` is a vector, that data is exported as a single
field. On the other hand, if it is an 2d array, each column is exported  as a
separate field.

"""
function vtkwrite(
    theFile::String,
    fens::FENodeSet,
    fes::T;
    opts...,
) where {T<:AbstractFESet}
    t = _match_type(_VTK_TYPE_MAP, typeof(fes))
    t === nothing && error("Unknown type $(typeof(fes))")
    celltype = WriteVTK.VTKCellTypes.VTKCellType(t)
    return vtkwrite(theFile, connasarray(fes), fens.xyz, celltype; opts...)
end

"""
    vtkwrite(theFile::String, Connectivity, Points, celltype; vectors=nothing, scalars=nothing)

Export mesh to VTK as an unstructured grid (binary format).

Arguments:
- `theFile` = file name,
- `Connectivity` = array of connectivities, one row per element,
- `Points` = array of node coordinates, one row per node,
- `Cell_type` = type of the cell, refer to the predefined 
    constants `WriteVTK.P1`, `WriteVTK.L2`, ..., `WriteVTK.H20``, ...
- `scalars` = array of tuples, (name, data)
- `vectors` = array of tuples, (name, data)

For the `scalars`: If `data` is a vector, that data is exported as a single
field. On the other hand, if it is an 2d array, each column is exported  as a
separate field.

Return the `vtk` file.
"""
function vtkwrite(
    theFile::String,
    Connectivity,
    Points,
    celltype;
    vectors = nothing,
    scalars = nothing,
)
    # The array `Points` has number of rows equal to the number of points,
    # whereas `points` expect this to be the transpose
    points = fill(zero(eltype(Points)), size(Points, 2), size(Points, 1))
    # Copy
    for i in axes(Points, 1)
        points[:, i] .= Points[i, :]
    end
    if typeof(Connectivity[1]) <: Tuple # Vector of connectivity tuples: convert to an array
        c = fill(0, length(Connectivity), length(Connectivity[1]))
        for i in eachindex(Connectivity)
            c[i, :] = [Connectivity[i]...]
        end
        Connectivity = c
    end
    # Prepare an array of the cells
    cells = [MeshCell(celltype, view(Connectivity, i, :)) for i in axes(Connectivity, 1)]
    vtkfile =
        vtk_grid(DataDrop.with_extension(theFile, ".vtu"), points, cells, compress = 3)
    if scalars !== nothing
        for nt in scalars
            an = nt[1] # name of the shape collection
            d = nt[2]
            pdata = fill(0.0, size(d, 2), size(d, 1))
            for j in axes(pdata, 1)
                for i in axes(pdata, 2)
                    pdata[j, i] = d[i, j]
                end
            end
            vtkfile[an] = pdata
        end
    end
    if vectors !== nothing
        for nt in vectors
            an = nt[1] # name of the shape collection
            d = nt[2]
            pdata = fill(0.0, size(d, 2), size(d, 1))
            for j in axes(pdata, 1)
                for i in axes(pdata, 2)
                    pdata[j, i] = d[i, j]
                end
            end
            vtkfile[an] = pdata
        end
    end
    vtk_save(vtkfile)
    return vtkfile
end

"""
    vtkwritecollection(theFile::String, fens::FENodeSet, fes::T, times; opts...) where {T<:AbstractFESet}

Write a collection of VTK files (`.pvd` file).

`times`: array of times

All the other arguments are the same as for `vtkwrite`. If `scalars` or
`vectors` are supplied, they correspond to the times in the `times` array.

See the `vtkwritecollection` methods.
"""
function vtkwritecollection(
    theFile::String,
    fens::FENodeSet,
    fes::T,
    times;
    opts...,
) where {T<:AbstractFESet}
    t = _match_type(_VTK_TYPE_MAP, typeof(fes))
    t === nothing && error("Unknown type $(typeof(fes))")
    celltype = WriteVTK.VTKCellTypes.VTKCellType(t)
    return vtkwritecollection(theFile, connasarray(fes), fens.xyz, celltype, times; opts...)
end

"""
    vtkwritecollection(theFile::String, Connectivity, Points, celltype, times; vectors=nothing, scalars=nothing)

Write a collection of VTK files (`.pvd` file).

`times`: array of times

All the other arguments are the same as for `vtkwrite`. If `scalars` or
`vectors` are supplied, they correspond to the times in the `times` array.

See the `vtkwritecollection` methods.
"""
function vtkwritecollection(
    theFile::String,
    Connectivity,
    Points,
    celltype,
    times;
    vectors = nothing,
    scalars = nothing,
)
    pvd = paraview_collection(DataDrop.with_extension(theFile, "pvd"))
    if scalars !== nothing
        (length(times) == length(scalars)) || error("The number of times and scalars must be the same")
        for (i, t, nt) in zip(1:length(times), times, scalars)
            an = nt[1] # name of the shape collection
            d = nt[2]
            pvd[t] = vtkwrite(
                DataDrop.with_extension(theFile * "-$i", ".vtu"),
                Connectivity,
                Points,
                celltype;
                scalars = [(an, d)],
            )
        end
    end
    if vectors !== nothing
        (length(times) == length(vectors)) || error("The number of times and vectors must be the same")
        for (i, t, nt) in zip(1:length(times), times, vectors)
            an = nt[1] # name of the shape collection
            d = nt[2]
            pvd[t] = vtkwrite(
                DataDrop.with_extension(theFile * "-$i", ".vtu"),
                Connectivity,
                Points,
                celltype;
                vectors = [(an, d)],
            )
        end
    end
    vtk_save(pvd)
    return pvd
end

end # VTKWrite

end
