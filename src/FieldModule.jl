"""
    FieldModule

Module for abstract fields.
"""
module FieldModule

__precompile__(true)

import Base.copyto!

"""
	AbstractField

Abstract field.

Expected  attributes:
  + `values::Array{T,2}`: Array of degree of freedom parameters,  indexed by entity number
  + `dofnums::Array{IT,2}`: Array of degree of freedom numbers, indexed by entity number
  + `kind::Matrix{Int8}`: Array of Boolean flags, indexed by entity number
  + `ranges::Dict(Int8, UnitRange{IT})`: Dictionary of ranges for the degrees of freedom.

See also: [`@add_Field_fields()`](@ref) .
"""
abstract type AbstractField end

"""
    KIND_INT

Constant representing the type of the integer representing the `kind` of a degree of freedom.
"""
const KIND_INT = Int8

"""
Predefined kinds of degrees of freedom.
"""
const DOF_KIND_ALL::KIND_INT = 0
const DOF_KIND_FREE::KIND_INT = 1
const DOF_KIND_DATA::KIND_INT = 2

"""
    add_Field_fields()

Generate the attributes (i. e. fields) of a `Field`. The methods defined for
the abstract type depend on these attributes to be present.
"""
macro add_Field_fields()
    return esc(:(values::Array{T,2};
    dofnums::Array{IT,2};
    kind::Matrix{KIND_INT};
    ranges::Vector{UnitRange{IT}}
    ))
end

"""
    ndofs(self::F)

How many degrees of freedom per entity? 

Ie. number of columns in the `values` array.
"""
ndofs(self::F) where {F<:AbstractField} = size(self.values, 2)

"""
    nents(self::F)

Number of entities associated with the field.
"""
nents(self::F) where {F<:AbstractField} = size(self.values, 1)

"""
    dofrange(self::F, kind) where {F<:AbstractField} 

Return the range of the degrees of freedom of `kind`.
"""
function dofrange(self::F, kind) where {F<:AbstractField} 
    if kind == DOF_KIND_ALL
        return 1:nalldofs(self)
    end
    return self.ranges[kind]
end

"""
    nfreedofs(self::F) where {F<:AbstractField}

Return the number of FREE degrees of freedom (known, data).
"""
nfreedofs(self::F) where {F<:AbstractField} = length(dofrange(self, DOF_KIND_FREE)) 

"""
    nfixeddofs(self::F)

Return the number of FIXED degrees of freedom (known, data).
"""
nfixeddofs(self::F) where {F<:AbstractField} = length(dofrange(self, DOF_KIND_DATA)) 

"""
    nalldofs(self::F) where {F<:AbstractField}

Return the number of ALL degrees of freedom (known, data).
"""
nalldofs(self::F) where {F<:AbstractField} = sum(length(v) for v in values(self.ranges))

"""
    freedofs(self::F) where {F<:AbstractField}

Return range corresponding to the free degrees of freedom.
"""
freedofs(self::F) where {F<:AbstractField} = dofrange(self, DOF_KIND_FREE)

"""
    fixeddofs(self::F) where {F<:AbstractField}

Return range corresponding to the fixed degrees of freedom.
"""
fixeddofs(self::F) where {F<:AbstractField} = dofrange(self, DOF_KIND_DATA)

"""
    copyto!(DEST::F,  SRC::F) where {F<:AbstractField}

Copy data from one field to another.
"""
function copyto!(DEST::F, SRC::F) where {F<:AbstractField}
    copyto!(DEST.values, SRC.values)
    copyto!(DEST.dofnums, SRC.dofnums)
    copyto!(DEST.kind, SRC.kind)
    DEST.ranges = deepcopy(SRC.ranges)
    return DEST
end

"""
    wipe!(self::F) where {F<:AbstractField}

Wipe all the data from the field.

This includes values, prescribed values, degree of freedom numbers, and "is
fixed" flags. The number of free degrees of freedom is set to zero.
"""
function wipe!(self::F) where {F<:AbstractField}
    Zer = zero(eltype(self.values[1]))
    fill!(self.dofnums, 0)
    fill!(self.kind, DOF_KIND_FREE)
    fill!(self.values, Zer)
    self.ranges = UnitRange{eltype(self.dofnums)}[]
    return self
end

"""
    gathersysvec(self::F, kind::KIND_INT = DOF_KIND_FREE) where {F<:AbstractField}

Gather values from the field for the system vector.

# Arguments
- `self`: field;
- `kind`: kind of degrees of freedom to gather; default is `DOF_KIND_FREE`.
"""
function gathersysvec(self::F, kind::KIND_INT = DOF_KIND_FREE) where {F<:AbstractField}
    N = length(dofrange(self, kind))
    vec = zeros(eltype(self.values), N)
    return gathersysvec!(self, vec, kind)
end

"""
    gathersysvec(self::F, kind::Symbol) where {F<:AbstractField}

Gather values from the field for the system vector.

This is a compatibility version, using a symbol.

# Arguments
- `self::F`: The field object.
- `kind::Symbol`: The kind of system vector to gather. 
"""
function gathersysvec(self::F, kind::Symbol) where {F<:AbstractField}
    if kind == :f
        return gathersysvec(self, DOF_KIND_FREE)
    elseif kind == :d
        return gathersysvec(self, DOF_KIND_DATA)
    else
        return gathersysvec(self, DOF_KIND_ALL)
    end
end

"""
    gathersysvec!(self::F, vec::Vector{T}, kind::KIND_INT = DOF_KIND_FREE) where {F<:AbstractField,T}

Gather values from the field for the system vector.

# Arguments
- `self`: field;
- `vec`: destination buffer;
- `kind`: integer, kind of degrees of freedom to gather: default is `DOF_KIND_FREE`.
"""
function gathersysvec!(self::F, vec::Vector{T}, kind::KIND_INT = DOF_KIND_FREE) where {F<:AbstractField,T}
    nents, dim = size(self.values)
    from, upto = first(dofrange(self, kind)), last(dofrange(self, kind))
    length(vec) == length(from:upto) || error("Vector needs to be of length equal to $(length(from:upto))")
    @inbounds for i in 1:nents
        for j in 1:dim
            en = self.dofnums[i, j]
            if from <= en <= upto
                p = en - from + 1
                vec[p] = self.values[i, j]
            end
        end
    end
    return vec
end

"""
    gathervalues_asvec!(
        self::F,
        dest::AbstractArray{T,1},
        conn::CC,
    ) where {F<:AbstractField, T, CC}

Gather values from the field into a vector.

The order is: for each node  in the connectivity, copy into the buffer all the
degrees of freedom,  then the next node and so on.

`dest` = destination buffer: overwritten  inside,  must be preallocated in the
correct size

The order of the loops matters, outer loop goes through the connectivity, inner
loop goes through the degrees of freedom for each entity.
"""
function gathervalues_asvec!(
    self::F,
    dest::AbstractArray{T,1},
    conn::CC,
) where {F<:AbstractField,T,CC}
    length(dest) == length(conn) * ndofs(self) || error("Destination buffer has wrong size")
    en = 1
    @inbounds for i in eachindex(conn)
        for j in axes(self.values, 2)
            dest[en] = self.values[conn[i], j]
            en = en + 1
        end
    end
    return dest
end

"""
    gathervalues_asmat!(
        self::F,
        dest::AbstractArray{T,2},
        conn::CC,
    ) where {F<:AbstractField, T, CC}

Gather values from the field into a two-dimensional array.

The order is: for each node  in the connectivity, copy into the corresponding
row of the buffer all the degrees of freedom,  then the next node into the next
row and so on.

`dest` = destination buffer: overwritten  inside,  must be preallocated
in the correct size

The order of the loops matters, outer loop goes through the connectivity, inner
loop goes through the degrees of freedom for each entity.
"""
function gathervalues_asmat!(
    self::F,
    dest::AbstractArray{T,2},
    conn::CC,
) where {F<:AbstractField,T,CC}
    size(dest) == (length(conn), ndofs(self)) || error("Destination buffer has wrong size")
    @inbounds for i in eachindex(conn)
        for j in axes(self.values, 2)
            dest[i, j] = self.values[conn[i], j]
        end
    end
    return dest
end

"""
    anyfixedvaluenz(self::F, conn::CC) where {F<:AbstractField, CC}

Is any degree of freedom fixed (prescribed) to be non-zero?
"""
function anyfixedvaluenz(self::F, conn::CC) where {F<:AbstractField,CC}
    Zer = zero(eltype(self.values))
    for i in eachindex(conn)
        for j in axes(self.values, 2)
            if self.kind[conn[i], j] == DOF_KIND_DATA # free degree of freedom
                if abs(self.values[conn[i], j]) > 0.0
                    return true
                end
            end
        end
    end
    return false
end

"""
    gatherdofnums!(self::F, dest::A, conn::CC) where {F<:AbstractField, A, CC}

Gather dofnums from the field.

The order is: for each node  in the connectivity, copy into the buffer all the degrees of
freedom for that node,  then the next node  and so on.
"""
function gatherdofnums!(self::F, dest::A, conn::CC) where {F<:AbstractField,A,CC}
    length(dest) == length(conn) * ndofs(self) || error("Destination buffer has wrong size")
    en = 1
    @inbounds for i in eachindex(conn)
        for j in axes(self.dofnums, 2)
            dest[en] = self.dofnums[conn[i], j]
            en = en + 1
        end
    end
    return dest
end

"""
    numberdofs!(self::F) where {F<:AbstractField}

Number the degrees of freedom.

The free components in the field are numbered consecutively, then all the fixed
components are numbered, again consecutively.

No effort is made to optimize the numbering in any way. If you'd like to
optimize the numbering of the degrees of freedom, use a form that sets the
permutation of the degrees of freedom, or the permutation of the nodes.
"""
function numberdofs!(self::F) where {F<:AbstractField}
    return numberdofs!(self, 1:nents(self))
end

"""
    numberdofs!(self::F, entperm) where {F<:AbstractField}

Number the degrees of freedom.

The free components in the field are numbered consecutively, then all the fixed
components are numbered, again consecutively.

The sequence of the entities is given by the `entperm` permutation (array or
range).
"""
function numberdofs!(self::F, entperm) where {F<:AbstractField}
    return numberdofs!(self, entperm, [DOF_KIND_FREE, DOF_KIND_DATA])
end

"""
    numberdofs!(self::F, entperm, kinds) where {F<:AbstractField}

Number the degrees of freedom.

# Arguments
- `self::F`: The field object.
- `entperm`: The permutation of entities.
- `kinds`: The kinds of degrees of freedom. The degrees of freedom are numbered
  in the order in which the kinds are given here.

# Examples
"""
function numberdofs!(self::F, entperm, kinds) where {F<:AbstractField}
    nents, dim = size(self.values)
    self.ranges = fill(UnitRange{eltype(self.dofnums)}(0:0), length(kinds))
    next = 1
    for k in kinds
        from = next
        for i in entperm
            for j  in 1:dim
                if self.kind[i, j] == k
                    self.dofnums[i, j] = next
                    next += 1
                end
            end
        end
        self.ranges[k] = from:(next - 1)
    end
    return self
end

"""
    setebc!(
        self::F,
        fenid::IT1,
        is_fixed::Bool,
        comp::IT2,
        val::T,
    ) where {F<:AbstractField,T<:Number,IT1<:Integer,IT2<:Integer}

Set the EBCs (essential boundary conditions).

`fenids`         - array of N node identifiers
`is_fixed` = scaler Boolean: are the degrees of freedom being fixed (true)
             or released (false),
`comp` = integer, which  degree of freedom (component),
`val` = array of N values of type T

Note:  Any call to `setebc!()` potentially changes the current assignment
which degrees of freedom are free and which are fixed and therefore is
presumed to invalidate the current degree-of-freedom numbering.
"""
function setebc!(
    self::F,
    fenid::IT1,
    is_fixed::Bool,
    comp::IT2,
    val::T,
) where {F<:AbstractField,T<:Number,IT1<:Integer,IT2<:Integer}
    (1 <= comp <= size(self.values, 2)) || error("Requested  nonexistent  degree of freedom")
    (1 <= fenid <= size(self.values, 1)) || error("Requested nonexistent node")
    if is_fixed
        self.kind[fenid, comp] = DOF_KIND_DATA
        self.values[fenid, comp] = val
    else
        self.kind[fenid, comp] = DOF_KIND_FREE
        self.values[fenid, comp] = zero(T)
    end
    self.ranges = eltype(self.ranges)[]
    return self
end

"""
    setebc!(
        self::F,
        fenids::AbstractVector{IT},
        is_fixed::Bool,
        comp::IT,
        val::AbstractVector{T},
    ) where {T<:Number, IT<:Integer}

Set the EBCs (essential boundary conditions).

`fenids`         - array of N node identifiers
`is_fixed` = scaler Boolean: are the degrees of freedom being fixed (true)
             or released (false),
`comp` = integer, which  degree of freedom (component),
`val` = array of N values of type T

Note:  Any call to `setebc!()` potentially changes the current assignment
which degrees of freedom are free and which are fixed and therefore is
presumed to invalidate the current degree-of-freedom numbering.
"""
function setebc!(
    self::F,
    fenids::AbstractVector{IT},
    is_fixed::Bool,
    comp::IT,
    val::AbstractVector{T},
) where {F<:AbstractField,T<:Number,IT<:Integer}
    (comp <= size(self.values, 2)) || error("Requested  nonexistent  degree of freedom")
    (maximum(fenids) <= size(self.values, 1)) || error("Requested nonexistent node")
    (size(fenids) == size(val)) || error("Arrays of mismatched sizes")
    for j in eachindex(fenids)
        setebc!(self, fenids[j], is_fixed, comp, val[j])
    end
    return self
end

"""
    setebc!(
        self::F,
        fenids::AbstractVector{IT},
        is_fixed::Bool,
        comp::IT,
        val::T = 0.0,
    ) where {T<:Number, IT<:Integer}

Set the EBCs (essential boundary conditions).

`fenids`         - array of N node identifiers
`is_fixed` = scaler Boolean: are the degrees of freedom being fixed (true)
             or released (false),
`comp` = integer, which  degree of freedom (component),
`val` = scalar of type T

Note:  Any call to `setebc!()` potentially changes the current assignment
which degrees of freedom are free and which are fixed and therefore is
presumed to invalidate the current degree-of-freedom numbering.
"""
function setebc!(
    self::F,
    fenids::AbstractVector{IT},
    is_fixed::Bool,
    comp::IT,
    val::T = 0.0,
) where {F<:AbstractField,T<:Number,IT<:Integer}
    ((comp >= 1 && comp <= size(self.values, 2))) || error("Requested  nonexistent  degree of freedom")
    (maximum(fenids) <= size(self.values, 1)) || error("Requested nonexistent node")
    (minimum(fenids) >= 1) || error("Requested nonexistent node")
    for j in eachindex(fenids)
        setebc!(self, fenids[j], is_fixed, comp, val)
    end
    return self
end

"""
    setebc!(
        self::F,
        fenids::AbstractVector{IT},
        comp::IT,
        val::AbstractVector{T},
    ) where {T<:Number, IT<:Integer}

Set the EBCs (essential boundary conditions).

`fenids` = array of N node identifiers
`comp` = integer, which  degree of freedom (component),
`val` = array of N values of type `T`

Note:  Any call to `setebc!()` potentially changes the current assignment
which degrees of freedom are free and which are fixed and therefore is
presumed to invalidate the current degree-of-freedom numbering.
"""
function setebc!(
    self::F,
    fenids::AbstractVector{IT},
    comp::IT,
    val::AbstractVector{T},
) where {F<:AbstractField,T<:Number,IT<:Integer}
    return setebc!(self, fenids, true, comp, val)
end

"""
    setebc!(
        self::F,
        fenids::AbstractVector{IT},
        comp::IT,
        val::T = 0.0,
    ) where {T<:Number, IT<:Integer}

Set the EBCs (essential boundary conditions).

`fenids` = array of N node identifiers
`comp` = integer, which  degree of freedom (component),
`val` = scalar of type `T`

Note:  Any call to `setebc!()` potentially changes the current assignment
which degrees of freedom are free and which are fixed and therefore is
presumed to invalidate the current degree-of-freedom numbering.
"""
function setebc!(
    self::F,
    fenids::AbstractVector{IT},
    comp::IT,
    val::T = 0.0,
) where {F<:AbstractField,T<:Number,IT<:Integer}
    return setebc!(self, fenids, true, comp, val)
end

"""
    setebc!(
        self::F,
        fenids::AbstractVector{IT},
        is_fixed::Bool,
        comp::AbstractVector{IT},
        val::T = 0.0,
    ) where {T<:Number, IT<:Integer}

Set the EBCs (essential boundary conditions).

`fenids` = array of N node identifiers
`comp` = integer vector, which degree of freedom (component),
`val` = scalar of type `T`, default is `0.0`

Note:  Any call to `setebc!()` potentially changes the current assignment
which degrees of freedom are free and which are fixed and therefore is
presumed to invalidate the current degree-of-freedom numbering.
"""
function setebc!(
    self::F,
    fenids::AbstractVector{IT},
    is_fixed::Bool,
    comp::AbstractVector{IT},
    val::T = 0.0,
) where {F<:AbstractField,T<:Number,IT<:Integer}
    for j in comp
        setebc!(self, fenids, is_fixed, j, val)
    end
    return self
end

"""
    setebc!(self::F, fenids::AbstractVector{IT})  where {IT<:Integer}

Set the EBCs (essential boundary conditions).

Suppress all degrees of freedom at the given nodes.

`fenids`         - array of N node identifiers

Note:  Any call to `setebc!()` potentially changes the current assignment
which degrees of freedom are free and which are fixed and therefore is
presumed to invalidate the current degree-of-freedom numbering.
"""
function setebc!(self::F, fenids::AbstractVector{IT}) where {F<:AbstractField,IT<:Integer}
    zer = zero(eltype(self.values[1]))
    for comp in axes(self.values, 2)
        setebc!(self, fenids, true, comp, zer)
    end
    return self
end

"""
    setebc!(self::F, fenid::IT) where {IT<:Integer}

Set the EBCs (essential boundary conditions).

Suppress all degrees of freedom at the given node.

`fenid` - One integer as a node identifier

All degrees of freedom at the node are set to zero.

Note:  Any call to `setebc!()` potentially changes the current assignment
which degrees of freedom are free and which are fixed and therefore is
presumed to invalidate the current degree-of-freedom numbering.
"""
function setebc!(self::F, fenid::IT) where {F<:AbstractField,IT<:Integer}
    return setebc!(self, [fenid])
end

"""
    setebc!(self::F)

Set the EBCs (essential boundary conditions).

All essential boundary conditions are CLEARED.

Note:  Any call to `setebc!()` potentially changes the current assignment
which degrees of freedom are free and which are fixed and therefore is
presumed to invalidate the current degree-of-freedom numbering.
"""
function setebc!(self::F) where {F<:AbstractField}
    zer = zero(eltype(self.values[1]))
    for fenid in axes(self.values, 1)
        for comp in axes(self.values, 2)
            setebc!(self, fenid, false, comp, zer)
        end
    end
    return self
end

"""
    applyebc!(self::F) where {F<:AbstractField}

Apply EBCs (essential boundary conditions).
"""
function applyebc!(self::F) where {F<:AbstractField}
    # NO-OP: all has already been done by setebc!
    return self
end

"""
    scattersysvec!(self::F, vec::AbstractVector{T}, kind::KIND_INT = DOF_KIND_FREE) where {F<:AbstractField,T<:Number}

Scatter values to the field from a system vector.

The vector may be for an arbitrary kind of degrees of freedom (default is the free degrees of freedom).
"""
function scattersysvec!(self::F, vec::AbstractVector{T}, kind::KIND_INT = DOF_KIND_FREE) where {F<:AbstractField,T<:Number}
    nents, dim = size(self.values)
    from, upto = first(dofrange(self, kind)), last(dofrange(self, kind))
    length(vec) == length(from:upto) || error("Vector needs to be of length equal to $(length(from:upto))")
    @inbounds for i in 1:nents
        for j in 1:dim
            en = self.dofnums[i, j]
            if from <= en <= upto
                p = en - from + 1
                self.values[i, j] = vec[p] 
            end
        end
    end
    return self
end

"""
    incrscattersysvec!(self::F, vec::AbstractVector{T}) where {F<:AbstractField, T<:Number}

Increment values of the field by scattering a system vector.

The vector may be either for just the free degrees of freedom, or for all the
degrees of freedom.
"""
function incrscattersysvec!(
    self::F,
    vec::AbstractVector{T},
) where {F<:AbstractField,T<:Number}
    nents, dim = size(self.values)
    nve = length(vec)
    for i = 1:nents
        for j = 1:dim
            dn = self.dofnums[i, j]
            if (dn > 0) && (dn <= nve)
                self.values[i, j] += vec[dn]
            end
        end
    end
    return self
end

"""
    prescribeddofs(uebc::F1, u::F2) where {F1<:AbstractField,  F2<:AbstractField}

Find which degrees of freedom are prescribed.
`uebc` = field which defines the constraints (is the dof fixed and to which value),
`u` = field which does not have the constraints applied, and serves as the source of equation numbers;
`uebc` and `u` may be one and the same field.

"""
function prescribeddofs(uebc::F1, u::F2) where {F1<:AbstractField,F2<:AbstractField}
    dofnums = eltype(u.dofnums)[]
    prescribedvalues = eltype(uebc.values)[]
    nents, dim = size(uebc.values)
    (size(uebc.values) == size(u.values)) || error("Fields of different sizes")
    for i = 1:nents
        for j = 1:dim
            if uebc.kind[i, j] == DOF_KIND_DATA
                dn = u.dofnums[i, j]
                push!(prescribedvalues, uebc.values[i, j])
                push!(dofnums, dn)
            end
        end
    end
    return dofnums, prescribedvalues
end

end
