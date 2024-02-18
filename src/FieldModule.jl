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
  + `is_fixed::Matrix{Bool}`: Array of Boolean flags, indexed by entity number
  + `_nfreedofs::IT`: Total number of free degrees of freedom

See also: [`@add_Field_fields()`](@ref) .
"""
abstract type AbstractField end

"""
    add_Field_fields()

Generate the attributes (i. e. fields) of a `Field`. The methods defined for
the abstract type depend on these attributes to be present.
"""
macro add_Field_fields()
    return esc(:(values::Array{T,2};
    dofnums::Array{IT,2};
    is_fixed::Matrix{Bool};
    _nfreedofs::IT))
end

"""
    ndofs(self::F)

Dimension of the degree of freedom parameters (i. e. how many degrees of
freedom per entity).
"""
ndofs(self::F) where {F<:AbstractField} = size(self.values, 2)

"""
    nents(self::F)

Number of nodes associated with the field.
"""
nents(self::F) where {F<:AbstractField} = size(self.values, 1)

"""
    nfreedofs(self::F)

Return to number of FREE degrees of freedom.
"""
nfreedofs(self::F) where {F<:AbstractField} = self._nfreedofs

"""
    nalldofs(self::F)

Return to number of ALL degrees of freedom (total number of degrees of freedom,
which is equal to the number of degrees of freedom per entity times the number
of entities).
"""
nalldofs(self::F) where {F<:AbstractField} = prod(size(self.values))

"""
    nfixeddofs(self::F)

Return to number of FIXED degrees of freedom (known, data).
"""
nfixeddofs(self::F) where {F<:AbstractField} = nalldofs(self) - nfreedofs(self)

"""
    freedofs(self::F) where {F<:AbstractField}

Return range corresponding to the free degrees of freedom.
"""
freedofs(self::F) where {F<:AbstractField} = 1:nfreedofs(self)

"""
    fixeddofs(self::F) where {F<:AbstractField}

Return range corresponding to the fixed degrees of freedom.
"""
fixeddofs(self::F) where {F<:AbstractField} = (nfreedofs(self)+1):nalldofs(self)

"""
    copyto!(DEST::F,  SRC::F) where {F<:AbstractField}

Copy data from one field to another.
"""
function copyto!(DEST::F, SRC::F) where {F<:AbstractField}
    copyto!(DEST.values, SRC.values)
    copyto!(DEST.dofnums, SRC.dofnums)
    copyto!(DEST.is_fixed, SRC.is_fixed)
    DEST._nfreedofs = SRC._nfreedofs
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
    self._nfreedofs = 0
    fill!(self.dofnums, 0)
    fill!(self.is_fixed, false)
    fill!(self.values, Zer)
    return self
end

"""
    gathersysvec(self::F, which = :f) where {F<:AbstractField}

Gather values from the field for the system vector.

# Arguments
- `self`: field;
- `which`:
    + `:f` - Collect a vector that includes the free degrees of freedom (default);
    + `:d` - Collect a vector that includes the fixed (data) degrees of freedom;
    + `:a` - Collect a vector that includes all the degrees of freedom, free and fixed.

The system vector consists of two parts: the first part, from `1` to `nfreedofs
(self)` are the free degrees of freedom, the second part from `nfreedofs
(self)+1` to `nalldofs(self)` are the fixed degrees of freedom.

This function returns either the entire vector, or one of the parts.
"""
function gathersysvec(self::F, which = :f) where {F<:AbstractField}
    nents, dim = size(self.values)
    N = 0
    if which == :f
        N = nfreedofs(self)
    else
        if which == :a
            N = nalldofs(self)
        else
            N = nfixeddofs(self)
        end
    end
    vec = zeros(eltype(self.values), N)
    return gathersysvec!(self, vec, which)
end

"""
    gathersysvec!(self::F,
        vec::Vector{T}, which = :f) where {F<:AbstractField, T}

Gather values from the field for the system vector.

# Arguments
- `self`: field;
- `which`:
    + `:f` - Collect a vector that includes the free degrees of freedom (default);
    + `:d` - Collect a vector that includes the fixed (data) degrees of freedom;
    + `:a` - Collect a vector that includes all the degrees of freedom, free and fixed.

The system vector consists of two parts: the first part, from `1` to `nfreedofs
(self)` are the free degrees of freedom, the second part from `nfreedofs
(self)+1` to `nalldofs(self)` are the fixed degrees of freedom.

This function gathers either the entire vector, or one of the parts. The length
of the supplied buffer `vec` must be correct, either `nfreedofs
(self)`, `nalldofs(self)`, or `nfixeddofs(self)`.
"""
function gathersysvec!(self::F, vec::Vector{T}, which = :f) where {F<:AbstractField,T}
    nents, dim = size(self.values)
    upto = length(vec)
    if which == :f
        N = nfreedofs(self)
        (upto == N) || error("Vector needs to be of length equal to $N")
        for j = 1:dim
            for i = 1:nents
                en = self.dofnums[i, j]
                if 1 <= en <= nfreedofs(self)
                    vec[en] = self.values[i, j]
                end
            end
        end
    else
        if which == :a
            N = nalldofs(self)
            (upto == N) || error("Vector needs to be of length equal to $N")
            for j = 1:dim
                for i = 1:nents
                    en = self.dofnums[i, j]
                    if 1 <= en <= nalldofs(self)
                        vec[en] = self.values[i, j]
                    end
                end
            end
        else
            N = nfixeddofs(self)
            (upto == N) || error("Vector needs to be of length equal to $N")
            for j = 1:dim
                for i = 1:nents
                    en = self.dofnums[i, j]
                    if en > nfreedofs(self)
                        vec[en-nfreedofs(self)] = self.values[i, j]
                    end
                end
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

`dest` = destination buffer: overwritten  inside,  must be preallocated
in the correct size
"""
function gathervalues_asvec!(
    self::F,
    dest::AbstractArray{T,1},
    conn::CC,
) where {F<:AbstractField,T,CC}
    # The order of the loops matters, first i, then j
    en = 1
    for i in eachindex(conn)
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
"""
function gathervalues_asmat!(
    self::F,
    dest::AbstractArray{T,2},
    conn::CC,
) where {F<:AbstractField,T,CC}
    @inbounds for j in axes(self.values, 2)
        for i in eachindex(conn)
            dest[i, j] = self.values[conn[i], j]
        end
    end
    return dest
end

"""
    gatherfixedvalues_asvec!(
        self::F,
        dest::AbstractArray{T,1},
        conn::CC,
    ) where {F<:AbstractField, T, CC}

Gather FIXED values from the field into a vector.

The order is: for each node  in the connectivity, copy into the buffer all the
fixed degrees of freedom,  then the next node and so on. If a degree of freedom
is NOT fixed, the corresponding entry is  set to zero.

`dest` = destination buffer: overwritten  inside,  must be preallocated
in the correct size
"""
function gatherfixedvalues_asvec!(
    self::F,
    dest::AbstractArray{T,1},
    conn::CC,
) where {F<:AbstractField,T,CC}
    Zer = zero(eltype(self.values))
    # The order of the loops matters here! It must be i, j
    en = 1
    for i in eachindex(conn)
        for j in axes(self.values, 2)
            if self.is_fixed[conn[i], j] # free degree of freedom
                dest[en] = self.values[conn[i], j]
            else
                dest[en] = Zer
            end
            en = en + 1
        end
    end
    return dest
end

"""
    gatherfixedvalues_asmat!(
        self::F,
        dest::AbstractArray{T,2},
        conn::CC,
    ) where {F<:AbstractField, T, CC}

Gather FIXED values from the field into a two-dimensional array.

The order is: for each node  in the connectivity, copy into the corresponding
row of the buffer all the degrees of freedom,  then the next node into the next
row and so on.  If a degree of freedom
is NOT fixed, the corresponding entry is set to zero.

`dest` = destination buffer: overwritten  inside,  must be preallocated
in the correct size
"""
function gatherfixedvalues_asmat!(
    self::F,
    dest::AbstractArray{T,2},
    conn::CC,
) where {F<:AbstractField,T,CC}
    Zer = zero(eltype(self.values))
    for j in axes(self.values, 2)
        for i in eachindex(conn)
            if self.is_fixed[conn[i], j] # fixed degree of freedom
                dest[i, j] = self.values[conn[i], j]
            else
                dest[i, j] = Zer
            end
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
            if self.is_fixed[conn[i], j] # free degree of freedom
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
    en = 1
    for i in eachindex(conn)
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
    nents, dim = size(self.values)
    self._nfreedofs = 0
    # First free
    for i in entperm
        for j = 1:dim
            if !self.is_fixed[i, j] # free degree of freedom
                self._nfreedofs = self._nfreedofs + 1
                self.dofnums[i, j] = self._nfreedofs
            end
        end
    end
    # Then fixed
    nfixeddofs = self._nfreedofs
    for i in entperm
        for j = 1:dim
            if self.is_fixed[i, j] # free degree of freedom
                nfixeddofs = nfixeddofs + 1
                self.dofnums[i, j] = nfixeddofs
            end
        end
    end
    return self
end

function _setebc!(
    self::F,
    fenid::IT,
    is_fixed::Bool,
    comp::IT,
    val::T,
) where {F<:AbstractField,T<:Number,IT<:Integer}
    self.is_fixed[fenid, comp] = is_fixed
    if self.is_fixed[fenid, comp]
        self.values[fenid, comp] = val
    else
        self.values[fenid, comp] = zero(T)
    end
    return self
end

"""
    setebc!(
        self::F,
        fenid::IT,
        is_fixed::Bool,
        comp::IT,
        val::T,
    ) where {T<:Number, IT<:Integer}

Set the EBCs (essential boundary conditions).

`fenids`         - array of N node identifiers
`is_fixed` = scaler Boolean: are the degrees of freedom being fixed (true)
             or released (false),
`comp` = integer, which  degree of freedom (component),
`val` = array of N values of type T

Note:  Any call to `setebc!()` potentially changes the current assignment
which degrees of freedom are free and which are fixed and therefore is
presumed to invalidate the current degree-of-freedom numbering. In such a case
this method sets `_nfreedofs = 0`; and  `dofnums=0`.
"""
function setebc!(
    self::F,
    fenid::IT,
    is_fixed::Bool,
    comp::IT,
    val::T,
) where {F<:AbstractField,T<:Number,IT<:Integer}
    @assert 1 <= comp <= size(self.values, 2) "Requested  nonexistent  degree of freedom"
    @assert 1 <= fenid <= size(self.values, 1) "Requested nonexistent node"
    _setebc!(self, fenid, is_fixed, comp, val)
    self._nfreedofs = 0
    fill!(self.dofnums, 0)
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
presumed to invalidate the current degree-of-freedom numbering. In such a case
this method sets `_nfreedofs = 0`; and  `dofnums=0`.
"""
function setebc!(
    self::F,
    fenids::AbstractVector{IT},
    is_fixed::Bool,
    comp::IT,
    val::AbstractVector{T},
) where {F<:AbstractField,T<:Number,IT<:Integer}
    @assert comp <= size(self.values, 2) "Requested  nonexistent  degree of freedom"
    @assert maximum(fenids) <= size(self.values, 1) "Requested nonexistent node"
    @assert size(fenids) == size(val) "Arrays of mismatched sizes"
    for j in eachindex(fenids)
        _setebc!(self, fenids[j], is_fixed, comp, val[j])
    end
    self._nfreedofs = 0
    fill!(self.dofnums, 0)
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
presumed to invalidate the current degree-of-freedom numbering. In such a case
this method sets `_nfreedofs = 0`; and  `dofnums=0`.
"""
function setebc!(
    self::F,
    fenids::AbstractVector{IT},
    is_fixed::Bool,
    comp::IT,
    val::T = 0.0,
) where {F<:AbstractField,T<:Number,IT<:Integer}
    @assert (comp >= 1 && comp <= size(self.values, 2)) "Requested  nonexistent  degree of freedom"
    @assert maximum(fenids) <= size(self.values, 1) "Requested nonexistent node"
    @assert minimum(fenids) >= 1 "Requested nonexistent node"
    for j in eachindex(fenids)
        _setebc!(self, fenids[j], is_fixed, comp, val)
    end
    self._nfreedofs = 0
    fill!(self.dofnums, 0)
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
presumed to invalidate the current degree-of-freedom numbering. In such a case
this method sets `_nfreedofs = 0`; and  `dofnums=0`.
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
presumed to invalidate the current degree-of-freedom numbering. In such a case
this method sets `_nfreedofs = 0`; and  `dofnums=0`.
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

Note:  Any call to `setebc!()` potentially changes the current assignment which
degrees of freedom are free and which are fixed and therefore is presumed to
invalidate the current degree-of-freedom numbering. In such a case this method
sets `_nfreedofs = 0`; and  `dofnums=0`.
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
presumed to invalidate the current degree-of-freedom numbering. In such a case
this method sets `_nfreedofs = 0`; and  `dofnums=0`.
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

`fenid`         - One integer as a node identifier

Note:  Any call to setebc!() potentially changes the current assignment
which degrees of freedom are free and which are fixed
and therefore is presumed to invalidate the
current degree-of-freedom numbering. In such a case this method sets
`_nfreedofs = 0`; and  `dofnums=0`.
"""
function setebc!(self::F, fenid::IT) where {F<:AbstractField,IT<:Integer}
    return setebc!(self, [fenid])
end

"""
    setebc!(self::F)

Set the EBCs (essential boundary conditions).

All essential boundary conditions are CLEARED.

Note:  Any call to setebc!() potentially changes the current assignment
which degrees of freedom are free and which are fixed
and therefore is presumed to invalidate the
current degree-of-freedom numbering. In such a case this method sets
`_nfreedofs = 0`; and  `dofnums=0`.
"""
function setebc!(self::F) where {F<:AbstractField}
    self._nfreedofs = 0
    fill!(self.dofnums, 0)
    fill!(self.is_fixed, false)
    fill!(self.values, zero(eltype(self.values[1])))
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
    scattersysvec!(self::F, vec::AbstractVector{T}) where {F<:AbstractField, T<:Number}

Scatter values to the field from a system vector.

The vector may be either for just the free degrees of freedom, or for all the
degrees of freedom.
"""
function scattersysvec!(self::F, vec::AbstractVector{T}) where {F<:AbstractField,T<:Number}
    nents, dim = size(self.values)
    nve = length(vec)
    for i = 1:nents
        for j = 1:dim
            dn = self.dofnums[i, j]
            if (dn > 0) && (dn <= nve)
                self.values[i, j] = vec[dn]
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
    @assert size(uebc.values) == size(u.values)
    for i = 1:nents
        for j = 1:dim
            if uebc.is_fixed[i, j]
                dn = u.dofnums[i, j]
                push!(prescribedvalues, uebc.values[i, j])
                push!(dofnums, dn)
            end
        end
    end
    return dofnums, prescribedvalues
end

end
