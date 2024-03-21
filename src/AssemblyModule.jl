"""
    AssemblyModule

Module for assemblers  of system matrices and vectors.
"""
module AssemblyModule

__precompile__(true)

using SparseArrays: sparse, spzeros, SparseMatrixCSC
using LinearAlgebra: diag
import Base: eltype
using ..MatrixUtilityModule: matrix_blocked_ff, vector_blocked_f, _zeros_via_calloc

"""
    AbstractSysmatAssembler

Abstract type of system-matrix assembler.
"""
abstract type AbstractSysmatAssembler end

"""
    eltype(a::A) where {A <: AbstractSysmatAssembler}

What is the type of the matrix buffer entries?
"""
eltype(a::A) where {A<:AbstractSysmatAssembler} = eltype(a._matbuffer)

setnomatrixresult(a::A, nomatrixresult) where {A<:AbstractSysmatAssembler} =
    let
        a._nomatrixresult = nomatrixresult
        a
    end

setforceinit(a::A, force_init) where {A<:AbstractSysmatAssembler} =
    let
        a._force_init = force_init
        a
    end

"""
    expectedntriples(
        a::A,
        elem_mat_nrows::IT,
        elem_mat_ncols::IT,
        n_elem_mats::IT,
    ) where {A<:AbstractSysmatAssembler, IT}

How many triples (I, J, V) does the assembler expect to store?

Default is: the product of the size of the element matrices times the
number of matrices. 
"""
expectedntriples(
    a::A,
    elem_mat_nrows::IT,
    elem_mat_ncols::IT,
    n_elem_mats::IT,
) where {A<:AbstractSysmatAssembler,IT} = elem_mat_nrows * elem_mat_ncols * n_elem_mats

function _resize_buffers(self, chunk)
    new_buffer_length = self._buffer_length + chunk
    resize!(self._rowbuffer, new_buffer_length)
    resize!(self._colbuffer, new_buffer_length)
    resize!(self._matbuffer, new_buffer_length)
    if self._force_init
        setvectorentries!(@view(self._rowbuffer[(self._buffer_length+1):end]), 1)
        setvectorentries!(@view(self._colbuffer[(self._buffer_length+1):end]), 1)
        setvectorentries!(
            @view(self._matbuffer[(self._buffer_length+1):end]),
            zero(eltype(self._matbuffer)),
        )
    end
    self._buffer_length = new_buffer_length
    return self
end

"""
    SysmatAssemblerSparse{IT, MBT, IBT} <: AbstractSysmatAssembler

Type for assembling a sparse global matrix from elementwise matrices.

!!! note

    All fields of the datatype are private. The type is manipulated by the
    functions `startassembly!`, `assemble!`, and `makematrix!`.
"""
mutable struct SysmatAssemblerSparse{IT,T,MBT<:AbstractVector{T},IJT<:AbstractVector{IT}} <:
               AbstractSysmatAssembler
    _buffer_length::IT
    _matbuffer::MBT
    _rowbuffer::IJT
    _colbuffer::IJT
    _buffer_pointer::IT
    _row_nalldofs::IT
    _col_nalldofs::IT
    _nomatrixresult::Bool
    _force_init::Bool
end

"""
    SysmatAssemblerSparse(z = zero(T), nomatrixresult = false) where {T}

Construct a sparse system matrix assembler.

The matrix entries are of type `T`. The assembler either produces a sparse
matrix (when `nomatrixresult = true`), or does not (when `nomatrixresult =
false`). When the assembler does not produce the sparse matrix when
`makematrix!` is called, it still can be constructed from the buffers stored in
the assembler, until they are cleared when the assembler is destroyed.

# Example

This is how a sparse matrix is assembled from two rectangular dense matrices.
```
    a = SysmatAssemblerSparse(0.0)
    startassembly!(a, 5, 5, 3, 7, 7)
    m = [0.24406   0.599773    0.833404  0.0420141
        0.786024  0.00206713  0.995379  0.780298
        0.845816  0.198459    0.355149  0.224996]
    assemble!(a, m, [1 7 5], [5 2 1 4])
    m = [0.146618  0.53471   0.614342    0.737833
         0.479719  0.41354   0.00760941  0.836455
         0.254868  0.476189  0.460794    0.00919633
         0.159064  0.261821  0.317078    0.77646
         0.643538  0.429817  0.59788     0.958909]
    assemble!(a, m, [2 3 1 7 5], [6 7 3 4])
    A = makematrix!(a)
```
Here `A` is a sparse matrix of the size 7x7.

When the `nomatrixresult` is set as true, no matrix is produced.
```
    a = SysmatAssemblerSparse(0.0, true)
    startassembly!(a, 5, 5, 3, 7, 7)
    m = [0.24406   0.599773    0.833404  0.0420141
        0.786024  0.00206713  0.995379  0.780298
        0.845816  0.198459    0.355149  0.224996]
    assemble!(a, m, [1 7 5], [5 2 1 4])
    m = [0.146618  0.53471   0.614342    0.737833
         0.479719  0.41354   0.00760941  0.836455
         0.254868  0.476189  0.460794    0.00919633
         0.159064  0.261821  0.317078    0.77646
         0.643538  0.429817  0.59788     0.958909]
    assemble!(a, m, [2 3 1 7 5], [6 7 3 4])
    A = makematrix!(a)
```

Here `A` is a named tuple of four sparse zero matrices. To construct the correct matrix is still
possible, for instance like this:
```
    a.nomatrixresult = false
    A = makematrix!(a)
```
At this point all the buffers of the assembler have potentially been cleared,
and `makematrix!(a) ` is no longer possible.

"""
function SysmatAssemblerSparse(z::T, nomatrixresult = false) where {T}
    return SysmatAssemblerSparse(0, T[z], Int[0], Int[0], 0, 0, 0, nomatrixresult, false)
end

function SysmatAssemblerSparse()
    return SysmatAssemblerSparse(zero(Float64))
end

"""
    startassembly!(self::SysmatAssemblerSparse{T},
        elem_mat_nrows::IT,
        elem_mat_ncols::IT,
        n_elem_mats::IT,
        row_nalldofs::IT,
        col_nalldofs::IT;
        force_init = false
        ) where {T, IT<:Integer}

Start the assembly of a global matrix.

The method makes buffers for matrix assembly. It must be called before
the first call to the method `assemble!`.

# Arguments
- `elem_mat_nrows` = row dimension of the element matrix;
- `elem_mat_ncols` = column dimension of the element matrix;
- `n_elem_mats` = number of element matrices;
- `row_nalldofs`= The total number of rows as a tuple;
- `col_nalldofs`= The total number of columns as a tuple.

If the `buffer_pointer` field of the assembler is 0, which is the case after
that assembler was created, the buffers are resized appropriately given the
dimensions on input. Otherwise, the buffers are left completely untouched. The
`buffer_pointer` is set to the beginning of the buffer, namely 1.

# Returns
- `self`: the modified assembler.

!!! note

    The buffers are initially not filled with anything meaningful.
    After the assembly, only the `(self._buffer_pointer - 1)` entries
    are meaningful numbers. Beware!

!!! note

    The buffers may be automatically resized if more numbers are assembled
    then initially intended. However, this operation will not necessarily
    be efficient and fast.
"""
function startassembly!(
    self::SysmatAssemblerSparse{T},
    elem_mat_nrows::IT,
    elem_mat_ncols::IT,
    n_elem_mats::IT,
    row_nalldofs::IT,
    col_nalldofs::IT;
    force_init = false,
) where {T,IT<:Integer}
    expected_ntriples = expectedntriples(self, elem_mat_nrows, elem_mat_ncols, n_elem_mats)
    # Only resize the buffers if the pointer is less than 1; otherwise the
    # buffers are already initialized and in use.
    if self._buffer_pointer < 1
        self._buffer_length = expected_ntriples
        resize!(self._rowbuffer, self._buffer_length)
        resize!(self._colbuffer, self._buffer_length)
        resize!(self._matbuffer, self._buffer_length)
        self._buffer_pointer = 1
        self._row_nalldofs = row_nalldofs
        self._col_nalldofs = col_nalldofs
    end
    # Leave the buffers uninitialized, unless the user requests otherwise
    self._force_init = force_init
    if self._force_init
        setvectorentries!(self._rowbuffer, 1)
        setvectorentries!(self._colbuffer, 1)
        setvectorentries!(self._matbuffer, zero(eltype(self._matbuffer)))
    end
    return self
end

"""
    assemble!(
        self::SysmatAssemblerSparse,
        mat::MT,
        dofnums_row::IT,
        dofnums_col::IT,
    ) where {MT, IT}

Assemble a rectangular matrix.
"""
function assemble!(
    self::SysmatAssemblerSparse,
    mat::MBT,
    dofnums_row::CIT,
    dofnums_col::CIT,
) where {MBT,CIT}
    # Assembly of a rectangular matrix.
    # The method assembles a rectangular matrix using the two vectors of
    # equation numbers for the rows and columns.
    nrows = length(dofnums_row)
    ncolumns = length(dofnums_col)
    p = self._buffer_pointer
    if self._buffer_length < p + ncolumns * nrows - 1
        self = _resize_buffers(self, ncolumns * nrows * 1000)
    end
    size(mat) == (nrows, ncolumns) || error("Wrong size of matrix")
    @inbounds for j = 1:ncolumns
        dj = dofnums_col[j]
        dj < 1 && error("Column degree of freedom < 1")
        dj > self._col_nalldofs && error("Column degree of freedom > size")
        for i = 1:nrows
            di = dofnums_row[i]
            di < 1 && error("Row degree of freedom < 1")
            di > self._row_nalldofs && error("Row degree of freedom > size")
            self._matbuffer[p] = mat[i, j] # serialized matrix
            self._rowbuffer[p] = di
            self._colbuffer[p] = dj
            p = p + 1
        end
    end
    self._buffer_pointer = p
    return self
end

"""
    makematrix!(self::SysmatAssemblerSparse)

Make a sparse matrix.

The sparse matrix is returned.

!!! note

    If `nomatrixresult` is set to true, dummy (zero) sparse matrices are
    returned. The entire result of the assembly is preserved in the assembler
    buffers. The ends of the buffers are filled with illegal
    (ignorable) values.

!!! note

    When the matrix is constructed (`nomatrixresult` is false), the buffers
    are not deallocated, and the `buffer_pointer` is set to 1. It is then
    possible to immediately start assembling another matrix.
"""
function makematrix!(self::SysmatAssemblerSparse)
    @assert length(self._rowbuffer) >= self._buffer_pointer - 1
    @assert length(self._colbuffer) >= self._buffer_pointer - 1
    # We have the option of retaining the assembled results, but not
    # constructing the sparse matrix.
    if self._nomatrixresult
        # Dummy (zero) sparse matrix is returned. The entire result of the
        # assembly is preserved in the assembler buffers. The ends of the
        # buffers are filled with legal (ignorable) values.
        self._rowbuffer[(self._buffer_pointer):end] .= 1
        self._colbuffer[(self._buffer_pointer):end] .= 1
        self._matbuffer[(self._buffer_pointer):end] .= 0.0
        return spzeros(self._row_nalldofs, self._col_nalldofs)
    end
    # The sparse matrix is constructed and returned.
    S = sparse(
        view(self._rowbuffer, 1:(self._buffer_pointer-1)),
        view(self._colbuffer, 1:(self._buffer_pointer-1)),
        view(self._matbuffer, 1:(self._buffer_pointer-1)),
        self._row_nalldofs,
        self._col_nalldofs,
    )
    # Get ready for more assembling
    self._buffer_pointer = 1
    return S
end

"""
    SysmatAssemblerSparseSymm{IT, MBT, IBT} <: AbstractSysmatAssembler

Assembler for a **symmetric square** matrix  assembled from symmetric square
matrices.

!!! note

    All fields of the datatype are private. The type is manipulated by the
    functions `startassembly!`, `assemble!`, and `makematrix!`.
"""
mutable struct SysmatAssemblerSparseSymm{
    IT,
    T,
    MBT<:AbstractVector{T},
    IJT<:AbstractVector{IT},
} <: AbstractSysmatAssembler
    _buffer_length::IT
    _matbuffer::MBT
    _rowbuffer::IJT
    _colbuffer::IJT
    _buffer_pointer::IT
    _row_nalldofs::IT
    _col_nalldofs::IT
    _nomatrixresult::Bool
    _force_init::Bool
end

expectedntriples(
    a::SysmatAssemblerSparseSymm,
    elem_mat_nrows::IT,
    elem_mat_ncols::IT,
    n_elem_mats::IT,
) where {IT} = Int((elem_mat_nrows * elem_mat_ncols + elem_mat_nrows) / 2 * n_elem_mats)

"""
    SysmatAssemblerSparseSymm(z::T, nomatrixresult = false) where {T}

Construct blank system matrix assembler for symmetric matrices. The matrix
entries are of type `T`.

# Example

This is how a symmetric sparse matrix is assembled from two square dense matrices.
```
    a = SysmatAssemblerSparseSymm(0.0)
    startassembly!(a, 5, 5, 3, 7, 7)
    m = [0.24406   0.599773    0.833404  0.0420141
        0.786024  0.00206713  0.995379  0.780298
        0.845816  0.198459    0.355149  0.224996]
    assemble!(a, m'*m, [5 2 1 4], [5 2 1 4])
    m = [0.146618  0.53471   0.614342    0.737833
         0.479719  0.41354   0.00760941  0.836455
         0.254868  0.476189  0.460794    0.00919633
         0.159064  0.261821  0.317078    0.77646
         0.643538  0.429817  0.59788     0.958909]
    assemble!(a, m'*m, [2 3 1 5], [2 3 1 5])
    A = makematrix!(a)
```
# See also
SysmatAssemblerSparse
"""
function SysmatAssemblerSparseSymm(z::T, nomatrixresult = false) where {T}
    return SysmatAssemblerSparseSymm(
        0,
        T[z],
        Int[0],
        Int[0],
        0,
        0,
        0,
        nomatrixresult,
        false,
    )
end

function SysmatAssemblerSparseSymm()
    return SysmatAssemblerSparseSymm(zero(Float64))
end

"""
    startassembly!(self::SysmatAssemblerSparseSymm{T},
        elem_mat_nrows::IT,
        elem_mat_ncols::IT,
        n_elem_mats::IT,
        row_nalldofs::IT,
        col_nalldofs::IT;
        force_init = false
        ) where {T, IT<:Integer}

Start the assembly of a global matrix.

The method makes buffers for matrix assembly. It must be called before
the first call to the method `assemble!`.

# Arguments
- `elem_mat_nrows` = row dimension of the element matrix;
- `elem_mat_ncols` = column dimension of the element matrix;
- `n_elem_mats` = number of element matrices;
- `row_nalldofs`= The total number of rows as a tuple;
- `col_nalldofs`= The total number of columns as a tuple.

If the `buffer_pointer` field of the assembler is 0, which is the case after
that assembler was created, the buffers are resized appropriately given the
dimensions on input. Otherwise, the buffers are left completely untouched. The
`buffer_pointer` is set to the beginning of the buffer, namely 1.

# Returns
- `self`: the modified assembler.

!!! note

    The buffers may be automatically resized if more numbers are assembled
    then initially intended. However, this operation will not necessarily
    be efficient and fast.

!!! note

    The buffers are initially not filled with anything meaningful.
    After the assembly, only the `(self._buffer_pointer - 1)` entries
    are meaningful numbers. Beware!
"""
function startassembly!(
    self::SysmatAssemblerSparseSymm{T},
    elem_mat_nrows::IT,
    elem_mat_ncols::IT,
    n_elem_mats::IT,
    row_nalldofs::IT,
    col_nalldofs::IT;
    force_init = false,
) where {T,IT<:Integer}
    @assert elem_mat_nrows == elem_mat_ncols "Symmetric sparse matrix is assumed to be assembled from square matrices"
    expected_ntriples = expectedntriples(self, elem_mat_nrows, elem_mat_ncols, n_elem_mats)
    # Only resize the buffers if the pointer is less than 1; otherwise the
    # buffers are already initialized and in use.
    if self._buffer_pointer < 1
        self._buffer_length = expected_ntriples
        resize!(self._rowbuffer, self._buffer_length)
        resize!(self._colbuffer, self._buffer_length)
        resize!(self._matbuffer, self._buffer_length)
        self._buffer_pointer = 1
        row_nalldofs == col_nalldofs || error("Row and column info do not agree")
        self._row_nalldofs = row_nalldofs
        self._col_nalldofs = col_nalldofs
    end
    # Leave the buffers uninitialized, unless the user requests otherwise
    self._force_init = force_init
    if force_init
        setvectorentries!(self._rowbuffer, 1)
        setvectorentries!(self._colbuffer, 1)
        setvectorentries!(self._matbuffer, zero(eltype(self._matbuffer)))
    end
    return self
end

"""
    assemble!(
        self::SysmatAssemblerSparseSymm,
        mat::MT,
        dofnums::IT,
        ignore
    ) where {MT, IT}

Assemble a square symmetric matrix.

`dofnums` are the row degree of freedom numbers, the column degree of freedom
number input is ignored (the row and column numbers are assumed to be the same).
"""
function assemble!(
    self::SysmatAssemblerSparseSymm,
    mat::MBT,
    dofnums_row::CIT,
    dofnums_col::CIT,
) where {MBT,CIT}
    # Assembly of a square symmetric matrix.
    # The method assembles the lower triangle of the square symmetric matrix using the two vectors of
    # equation numbers for the rows and columns.
    nrows = length(dofnums_row)
    ncolumns = length(dofnums_col)
    @assert nrows == ncolumns
    p = self._buffer_pointer
    if self._buffer_length < p + (ncolumns * nrows + nrows) / 2 - 1
        self = _resize_buffers(self, ncolumns * nrows * 1000)
    end
    @assert size(mat) == (nrows, ncolumns)
    @inbounds for j = 1:ncolumns
        dj = dofnums_col[j]
        dj < 1 && error("Column degree of freedom < 1")
        dj > self._col_nalldofs && error("Column degree of freedom > size")
        for i = j:nrows
            di = dofnums_row[i]
            di < 1 && error("Row degree of freedom < 1")
            di > self._row_nalldofs && error("Row degree of freedom > size")
            self._matbuffer[p] = mat[i, j] # serialized matrix
            self._rowbuffer[p] = dofnums_row[i]
            self._colbuffer[p] = dofnums_col[j]
            p = p + 1
        end
    end
    self._buffer_pointer = p
    return self
end

"""
    makematrix!(self::SysmatAssemblerSparseSymm)

Make a sparse symmetric square matrix.

!!! note

    If `nomatrixresult` is set to true, dummy (zero) sparse matrix is returned. The entire result of the assembly is preserved in the assembler buffers.
    The ends of the buffers are filled with illegal (ignorable) values.

!!! note

    When the matrix is constructed (`nomatrixresult` is false), the buffers
    are not deallocated, and the `buffer_pointer` is set to 1. It is then
    possible to immediately start assembling another matrix.
"""
function makematrix!(self::SysmatAssemblerSparseSymm)
    @assert length(self._rowbuffer) >= self._buffer_pointer - 1
    @assert length(self._colbuffer) >= self._buffer_pointer - 1
    # We have the option of retaining the assembled results, but not
    # constructing the sparse matrix.
    if self._nomatrixresult
        # Dummy (zero) sparse matrix is returned. The entire result of the
        # assembly is preserved in the assembler buffers. The ends of the
        # buffers are filled with legal (ignorable) values.
        self._rowbuffer[(self._buffer_pointer):end] .= 1
        self._colbuffer[(self._buffer_pointer):end] .= 1
        self._matbuffer[(self._buffer_pointer):end] .= 0.0
        return spzeros(self._row_nalldofs, self._col_nalldofs)
    end
    # The sparse matrix is constructed and returned. The  buffers used for
    # the assembly are cleared.
    S = sparse(
        view(self._rowbuffer, 1:(self._buffer_pointer-1)),
        view(self._colbuffer, 1:(self._buffer_pointer-1)),
        view(self._matbuffer, 1:(self._buffer_pointer-1)),
        self._row_nalldofs,
        self._col_nalldofs,
    )
    #  Now we need to construct the other triangle of the matrix. The diagonal
    #  will be duplicated.
    S = S + transpose(S)
    @inbounds for j in axes(S, 1)
        S[j, j] *= 0.5      # the diagonal is there twice; fix it;
    end
    # Get ready for more assembling
    self._buffer_pointer = 1
    return S
end

"""
    SysmatAssemblerSparseDiag{T<:Number} <: AbstractSysmatAssembler

Assembler for a **symmetric square diagonal** matrix  assembled from symmetric
square diagonal matrices.

Warning: off-diagonal elements of the elementwise matrices will be ignored
during assembly!

!!! note

    All fields of the datatype are private. The type is manipulated by the
    functions `startassembly!`, `assemble!`, and `makematrix!`.
"""
mutable struct SysmatAssemblerSparseDiag{
    IT,
    T,
    MBT<:AbstractVector{T},
    IJT<:AbstractVector{IT},
} <: AbstractSysmatAssembler
    # Type for assembling of a sparse global matrix from elementwise matrices.
    _buffer_length::IT
    _matbuffer::MBT
    _rowbuffer::IJT
    _colbuffer::IJT
    _buffer_pointer::IT
    _row_nalldofs::IT
    _col_nalldofs::IT
    _nomatrixresult::Bool
    _force_init::Bool
end

expectedntriples(
    a::SysmatAssemblerSparseDiag,
    elem_mat_nrows::IT,
    elem_mat_ncols::IT,
    n_elem_mats::IT,
) where {IT} = (max(elem_mat_nrows, elem_mat_ncols) * n_elem_mats)

"""
    SysmatAssemblerSparseDiag(z::T, nomatrixresult = false) where {T}

Construct blank system matrix assembler for square diagonal matrices. The matrix
entries are of type `T`.

"""
function SysmatAssemblerSparseDiag(z::T, nomatrixresult = false) where {T}
    return SysmatAssemblerSparseDiag(
        0,
        T[z],
        Int[0],
        Int[0],
        0,
        0,
        0,
        nomatrixresult,
        false,
    )
end

function SysmatAssemblerSparseDiag()
    return SysmatAssemblerSparseDiag(zero(Float64))
end

"""
    startassembly!(self::SysmatAssemblerSparseDiag{T},
        elem_mat_nrows::IT,
        elem_mat_ncols::IT,
        n_elem_mats::IT,
        row_nalldofs::IT,
        col_nalldofs::IT;
        force_init = false
        ) where {T, IT<:Integer}

Start the assembly of a symmetric square diagonal matrix.

The method makes buffers for matrix assembly. It must be called before
the first call to the method `assemble!`.

# Arguments
- `elem_mat_nrows` = row dimension of the element matrix;
- `elem_mat_ncols` = column dimension of the element matrix;
- `n_elem_mats` = number of element matrices;
- `row_nalldofs`= The total number of rows as a tuple;
- `col_nalldofs`= The total number of columns as a tuple.

The values stored in the buffers are initially undefined!

# Returns
- `self`: the modified assembler.
"""
function startassembly!(
    self::SysmatAssemblerSparseDiag{T},
    elem_mat_nrows::IT,
    elem_mat_ncols::IT,
    n_elem_mats::IT,
    row_nalldofs::IT,
    col_nalldofs::IT;
    force_init = false,
) where {T,IT<:Integer}
    @assert elem_mat_nrows == elem_mat_ncols "Diagonal sparse matrix is assumed to be assembled from square matrices"
    expected_ntriples = expectedntriples(self, elem_mat_nrows, elem_mat_ncols, n_elem_mats)
    # Only resize the buffers if the pointer is less than 1; otherwise the
    # buffers are already initialized and in use.
    if self._buffer_pointer < 1
        self._buffer_length = expected_ntriples
        resize!(self._rowbuffer, self._buffer_length)
        resize!(self._colbuffer, self._buffer_length)
        resize!(self._matbuffer, self._buffer_length)
        self._buffer_pointer = 1
        row_nalldofs == col_nalldofs || error("Row and column info do not agree")
        self._row_nalldofs = row_nalldofs
        self._col_nalldofs = col_nalldofs
    end
    # Leave the buffers uninitialized, unless the user requests otherwise
    self._force_init = force_init
    if force_init
        setvectorentries!(self._rowbuffer, 1)
        setvectorentries!(self._colbuffer, 1)
        setvectorentries!(self._matbuffer, zero(eltype(self._matbuffer)))
    end
    return self
end

"""
    assemble!(
        self::SysmatAssemblerSparseDiag,
        mat::MT,
        dofnums_row::IV,
        dofnums_col::IV,
    ) where {MT, IV}

Assemble a square symmetric diagonal matrix.

- `dofnums` = the row degree of freedom numbers, the column degree of freedom
number input is ignored (the row and column numbers are assumed to be the same).
- `mat` = **diagonal** square matrix
"""
function assemble!(
    self::SysmatAssemblerSparseDiag,
    mat::MBT,
    dofnums_row::CIT,
    dofnums_col::CIT,
) where {MBT,CIT}
    # Assembly of a square symmetric matrix.
    # The method assembles the lower triangle of the square symmetric matrix using the two vectors of
    # equation numbers for the rows and columns.
    nrows = length(dofnums_row)
    ncolumns = length(dofnums_col)
    @assert nrows == ncolumns
    p = self._buffer_pointer
    if self._buffer_length < p + ncolumns * nrows - 1
        self = _resize_buffers(self, ncolumns * nrows * 1000)
    end
    @assert size(mat) == (nrows, ncolumns)
    @inbounds for j = 1:ncolumns
        dj = dofnums_col[j]
        dj < 1 && error("Column degree of freedom < 1")
        dj > self._col_nalldofs && error("Column degree of freedom > size")
        self._matbuffer[p] = mat[j, j] # serialized matrix
        self._rowbuffer[p] = dj
        self._colbuffer[p] = dj
        p = p + 1
    end
    self._buffer_pointer = p
    return self
end

"""
    makematrix!(self::SysmatAssemblerSparseDiag)

Make a sparse symmetric square diagonal matrix.

!!! note

    If `nomatrixresult` is set to true, dummy (zero) sparse matrix is returned. The entire result of the assembly is preserved in the assembler buffers.
    The ends of the buffers are filled with illegal (ignorable) values.

!!! note

    When the matrix is constructed (`nomatrixresult` is false), the buffers
    are not deallocated, and the `buffer_pointer` is set to 1. It is then
    possible to immediately start assembling another matrix.
"""
function makematrix!(self::SysmatAssemblerSparseDiag)
    @assert length(self._rowbuffer) >= self._buffer_pointer - 1
    @assert length(self._colbuffer) >= self._buffer_pointer - 1
    # We have the option of retaining the assembled results, but not
    # constructing the sparse matrix.
    if self._nomatrixresult
        # Dummy (zero) sparse matrix is returned. The entire result of the
        # assembly is preserved in the assembler buffers. The ends of the
        # buffers are filled with illegal (ignorable) values.
        self._rowbuffer[(self._buffer_pointer):end] .= 0
        self._colbuffer[(self._buffer_pointer):end] .= 0
        self._matbuffer[(self._buffer_pointer):end] .= 0.0
        return spzeros(self._row_nalldofs, self._col_nalldofs)
    end
    # The sparse matrix is constructed and returned. The  buffers used for
    # the assembly are cleared.
    S = sparse(
        view(self._rowbuffer, 1:(self._buffer_pointer-1)),
        view(self._colbuffer, 1:(self._buffer_pointer-1)),
        view(self._matbuffer, 1:(self._buffer_pointer-1)),
        self._row_nalldofs,
        self._col_nalldofs,
    )
    # Get ready for more assembling
    self._buffer_pointer = 1
    return S
end

# =============================================================================
# =============================================================================
# =============================================================================
# =============================================================================

"""
    AbstractSysvecAssembler

Abstract type of system vector assembler.
"""
abstract type AbstractSysvecAssembler end

eltype(a::A) where {A<:AbstractSysvecAssembler} = eltype(a._buffer)

"""
    startassembly!(self::SysvecAssembler{T}, ndofs_row::FInt) where {T<:Number}

Start assembly.

The method makes the buffer for the vector assembly. It must be called before
the first call to the method assemble.

- `elem_mat_nmatrices` = number of element matrices expected to be processed
  during the assembly.
- `ndofs_row`= Total number of degrees of freedom.

# Returns
- `self`: the modified assembler.
"""
function startassembly!(
    self::SV,
    row_nalldofs::Tuple{IT,IT},
) where {SV<:AbstractSysvecAssembler,IT} end

"""
    assemble!(self::SysvecAssembler{T}, vec::MV,
      dofnums::D) where {T<:Number, MV<:AbstractArray{T}, D<:AbstractArray{FInt}}

Assemble an elementwise vector.

The method assembles a column element vector using the vector of degree of
freedom numbers for the rows.
"""
function assemble!(self::SV, vec::MV, dofnums::IV) where {SV<:AbstractSysvecAssembler,MV,IV} end

"""
    makevector!(self::SysvecAssembler)

Make the global vector.
"""
function makevector!(self::SV) where {SV<:AbstractSysvecAssembler} end

"""
    SysvecAssembler

Assembler for the system vector.
"""
mutable struct SysvecAssembler{VBT,IT} <: AbstractSysvecAssembler
    _buffer::VBT
    _row_nalldofs::IT
end

"""
    SysvecAssembler(z::T) where {T}

Construct blank system vector assembler. The vector entries are of type `T`
determined by the zero value.
"""
function SysvecAssembler(z::T) where {T}
    return SysvecAssembler(T[z], 1)
end

function SysvecAssembler()
    return SysvecAssembler(zero(Float64))
end

"""
    startassembly!(self::SysvecAssembler, ndofs_row)

Start assembly.

The method makes the buffer for the vector assembly. It must be called before
the first call to the method assemble.

`ndofs_row`= Total number of degrees of freedom.

# Returns
- `self`: the modified assembler.
"""
function startassembly!(self::SysvecAssembler, row_nalldofs::IT) where {IT<:Integer}
    self._row_nalldofs = row_nalldofs
    resize!(self._buffer, self._row_nalldofs)
    self._buffer .= 0
    return self
end

"""
    assemble!(self::SysvecAssembler{T}, vec::MV,
      dofnums::D) where {T<:Number, MV<:AbstractArray{T}, D<:AbstractArray{FInt}}

Assemble an elementwise vector.

The method assembles a column element vector using the vector of degree of
freedom numbers for the rows.
"""
function assemble!(self::SysvecAssembler, vec::MV, dofnums::IV) where {MV,IV}
    @inbounds for i in eachindex(dofnums)
        gi = dofnums[i]
        gi < 1 && error("Row degree of freedom < 1")
        gi > self._row_nalldofs && error("Row degree of freedom > size")
        self._buffer[gi] += vec[i]
    end
end

"""
    makevector!(self::SysvecAssembler)

Make the global vector.
"""
function makevector!(self::SysvecAssembler)
    return deepcopy(self._buffer)
end

"""
    SysmatAssemblerSparseHRZLumpingSymm{IT, MBT, IBT} <: AbstractSysmatAssembler

Assembler for a **symmetric lumped square** matrix  assembled from  **symmetric square**
matrices.

Reference: A note on mass lumping and related processes in the finite element method,
E. Hinton, T. Rock, O. C. Zienkiewicz, Earthquake Engineering & Structural Dynamics,
volume 4, number 3, 245--249, 1976.


!!! note

    All fields of the datatype are private. The type is manipulated by the
    functions `startassembly!`, `assemble!`, and `makematrix!`.

!!! note

    This assembler can compute and assemble diagonalized mass matrices.
    However, if the meaning of the entries of the mass matrix  differs
    (translation versus rotation), the mass matrices will not be computed
    correctly. Put bluntly: it can only be used for homogeneous mass matrices,
    all translation degrees of freedom, for instance.
"""
mutable struct SysmatAssemblerSparseHRZLumpingSymm{
    IT,
    T,
    MBT<:AbstractVector{T},
    IJT<:AbstractVector{IT},
} <: AbstractSysmatAssembler
    _buffer_length::IT
    _matbuffer::MBT
    _rowbuffer::IJT
    _colbuffer::IJT
    _buffer_pointer::IT
    _row_nalldofs::IT
    _col_nalldofs::IT
    _nomatrixresult::Bool
    _force_init::Bool
end

expectedntriples(
    a::SysmatAssemblerSparseHRZLumpingSymm,
    elem_mat_nrows::IT,
    elem_mat_ncols::IT,
    n_elem_mats::IT,
) where {IT} = (max(elem_mat_nrows, elem_mat_ncols) * n_elem_mats)

"""
    SysmatAssemblerSparseHRZLumpingSymm(z::T, nomatrixresult = false) where {T}

Construct blank system matrix assembler. The matrix entries are of type `T`.

"""
function SysmatAssemblerSparseHRZLumpingSymm(z::T, nomatrixresult = false) where {T}
    return SysmatAssemblerSparseHRZLumpingSymm(
        0,
        T[z],
        Int[0],
        Int[0],
        0,
        0,
        0,
        nomatrixresult,
        false,
    )
end

function SysmatAssemblerSparseHRZLumpingSymm()
    return SysmatAssemblerSparseHRZLumpingSymm(zero(Float64))
end

"""
    startassembly!(self::SysmatAssemblerSparseHRZLumpingSymm{T},
            elem_mat_nrows::IT,
            elem_mat_ncols::IT,
            n_elem_mats::IT,
            row_nalldofs::IT,
            col_nalldofs::IT;
            force_init = false
            ) where {T, IT<:Integer}

Start the assembly of a symmetric lumped diagonal square global matrix.

The method makes buffers for matrix assembly. It must be called before
the first call to the method `assemble!`.

# Arguments
- `elem_mat_nrows` = row dimension of the element matrix;
- `elem_mat_ncols` = column dimension of the element matrix;
- `n_elem_mats` = number of element matrices;
- `row_nalldofs`= The total number of rows as a tuple;
- `col_nalldofs`= The total number of columns as a tuple.

The values stored in the buffers are initially undefined!

# Returns
- `self`: the modified assembler.
"""
function startassembly!(
    self::SysmatAssemblerSparseHRZLumpingSymm{T},
    elem_mat_nrows::IT,
    elem_mat_ncols::IT,
    n_elem_mats::IT,
    row_nalldofs::IT,
    col_nalldofs::IT;
    force_init = false,
) where {T,IT<:Integer}
    @assert elem_mat_nrows == elem_mat_ncols "Diagonal sparse matrix is assumed to be assembled from square matrices"
    expected_ntriples = expectedntriples(self, elem_mat_nrows, elem_mat_ncols, n_elem_mats)
    # Only resize the buffers if the pointer is less than 1; otherwise the
    # buffers are already initialized and in use.
    if self._buffer_pointer < 1
        self._buffer_length = expected_ntriples
        resize!(self._rowbuffer, self._buffer_length)
        resize!(self._colbuffer, self._buffer_length)
        resize!(self._matbuffer, self._buffer_length)
        self._buffer_pointer = 1
        row_nalldofs == col_nalldofs || error("Row and column info do not agree")
        self._row_nalldofs = row_nalldofs
        self._col_nalldofs = col_nalldofs
    end
    # Leave the buffers uninitialized, unless the user requests otherwise
    self._force_init = force_init
    if force_init
        setvectorentries!(self._rowbuffer, 1)
        setvectorentries!(self._colbuffer, 1)
        setvectorentries!(self._matbuffer, zero(eltype(self._matbuffer)))
    end
    return self
end

"""
    assemble!(
        self::SysmatAssemblerSparseHRZLumpingSymm,
        mat::MT,
        dofnums::IV,
        ignore::IV,
    ) where {MT, IV}

Assemble a HRZ-lumped square symmetric matrix.

Assembly of a HRZ-lumped square symmetric matrix. The method assembles the
scaled diagonal of the square symmetric matrix using the two vectors of
equation numbers for the rows and columns.
"""
function assemble!(
    self::SysmatAssemblerSparseHRZLumpingSymm,
    mat::MBT,
    dofnums_row::CIT,
    dofnums_col::CIT,
) where {MBT,CIT}
    nrows = length(dofnums_row)
    ncolumns = length(dofnums_col)
    @assert nrows == ncolumns
    p = self._buffer_pointer
    if self._buffer_length < p + ncolumns * nrows - 1
        self = _resize_buffers(self, ncolumns * nrows * 1000)
    end
    @assert size(mat) == (nrows, ncolumns)
    # Now comes the lumping procedure
    em2 = sum(sum(mat, dims = 1)) # total mass times the number of space dimensions
    dem2 = zero(eltype(mat)) # total mass on the diagonal times the number of space dimensions
    @inbounds for i = 1:nrows
        dem2 += mat[i, i]
    end
    ffactor = em2 / dem2 # total-element-mass compensation factor
    @inbounds for j = 1:ncolumns
        dj = dofnums_col[j]
        dj < 1 && error("Column degree of freedom < 1")
        dj > self._col_nalldofs && error("Column degree of freedom > size")
        self._matbuffer[p] = mat[j, j] * ffactor # serialized matrix
        self._rowbuffer[p] = dj
        self._colbuffer[p] = dj
        p = p + 1
    end
    self._buffer_pointer = p
    return self
end

"""
    makematrix!(self::SysmatAssemblerSparseHRZLumpingSymm)

Make a sparse HRZ-lumped **symmetric square**  matrix.

!!! note

    If `nomatrixresult` is set to true, dummy (zero) sparse matrix is returned. The entire result of the assembly is preserved in the assembler buffers.
    The ends of the buffers are filled with illegal (ignorable) values.

!!! note

    When the matrix is constructed (`nomatrixresult` is false), the buffers
    are not deallocated, and the `buffer_pointer` is set to 1. It is then
    possible to immediately start assembling another matrix.
"""
function makematrix!(self::SysmatAssemblerSparseHRZLumpingSymm)
    @assert length(self._rowbuffer) >= self._buffer_pointer - 1
    @assert length(self._colbuffer) >= self._buffer_pointer - 1
    # We have the option of retaining the assembled results, but not
    # constructing the sparse matrix.
    if self._nomatrixresult
        # Dummy (zero) sparse matrix is returned. The entire result of the
        # assembly is preserved in the assembler buffers. The ends of the
        # buffers are filled with illegal (ignorable) values.
        self._rowbuffer[(self._buffer_pointer):end] .= 0
        self._colbuffer[(self._buffer_pointer):end] .= 0
        self._matbuffer[(self._buffer_pointer):end] .= 0.0
        return spzeros(self._row_nalldofs, self._col_nalldofs)
    end
    # The sparse matrix is constructed and returned. The  buffers used for
    # the assembly are cleared.
    S = sparse(
        view(self._rowbuffer, 1:(self._buffer_pointer-1)),
        view(self._colbuffer, 1:(self._buffer_pointer-1)),
        view(self._matbuffer, 1:(self._buffer_pointer-1)),
        self._row_nalldofs,
        self._col_nalldofs,
    )
    # Get ready for more assembling
    self._buffer_pointer = 1
    return S
end

"""
    SysmatAssemblerFFBlock{A<:AbstractSysmatAssembler, IT} <: AbstractSysmatAssembler

Type for extracting a free-free matrix, delegating the actual assembly to a
different assembler.
"""
mutable struct SysmatAssemblerFFBlock{A<:AbstractSysmatAssembler,IT} <:
               AbstractSysmatAssembler
    _a::A
    _row_nfreedofs::IT
    _col_nfreedofs::IT
end

eltype(self::A) where {A<:SysmatAssemblerFFBlock} = eltype(self._a)

"""
    SysmatAssemblerFFBlock(row_nfreedofs::IT, col_nfreedofs = row_nfreedofs) where {IT<:Integer}

Constructor, where the wrapped assembler is for general sparse matrices.

Supply the number of free degrees of freedom.
"""
function SysmatAssemblerFFBlock(
    row_nfreedofs::IT,
    col_nfreedofs = row_nfreedofs,
) where {IT<:Integer}
    return SysmatAssemblerFFBlock(SysmatAssemblerSparse(), row_nfreedofs, col_nfreedofs)
end

"""
    startassembly!(self::SysmatAssemblerFFBlock,
        elem_mat_nrows::IT,
        elem_mat_ncols::IT,
        n_elem_mats::IT,
        row_nalldofs::IT,
        col_nalldofs::IT;
        force_init = false) where {IT <: Integer}

Start assembly, delegate to the wrapped assembler.
"""
function startassembly!(
    self::SysmatAssemblerFFBlock,
    elem_mat_nrows::IT,
    elem_mat_ncols::IT,
    n_elem_mats::IT,
    row_nalldofs::IT,
    col_nalldofs::IT;
    force_init = false,
) where {IT<:Integer}
    self._a = startassembly!(
        self._a,
        elem_mat_nrows,
        elem_mat_ncols,
        n_elem_mats,
        row_nalldofs,
        col_nalldofs;
        force_init,
    )
    return self
end

"""
    assemble!(self::SysmatAssemblerFFBlock,
        mat::MBT,
        dofnums_row::CIT,
        dofnums_col::CIT) where {MBT, CIT}

Assemble a matrix, v
"""
function assemble!(
    self::SysmatAssemblerFFBlock,
    mat::MBT,
    dofnums_row::CIT,
    dofnums_col::CIT,
) where {MBT,CIT}
    assemble!(self._a, mat, dofnums_row, dofnums_col)
    return self
end

"""
    makematrix!(self::SysmatAssemblerFFBlock)

Make an assembled matrix. Delegate the construction of the matrix to the wrapped
assembler. Then extract the left upper corner block of the matrix(the free-free
matrix).
"""
function makematrix!(self::SysmatAssemblerFFBlock)
    return matrix_blocked_ff(makematrix!(self._a), self._row_nfreedofs, self._col_nfreedofs)
end

"""
    SysvecAssemblerFBlock

Assembler for the system vector, which extracts the free vector.
"""
mutable struct SysvecAssemblerFBlock{A<:AbstractSysvecAssembler,IT} <:
               AbstractSysvecAssembler
    _a::A
    _row_nfreedofs::IT
end

eltype(self::A) where {A<:SysvecAssemblerFBlock} = eltype(self._a)

"""
    SysvecAssemblerFBlock(row_nfreedofs::IT) where {IT}

Constructor of the free block assembler.
"""
function SysvecAssemblerFBlock(row_nfreedofs::IT) where {IT}
    return SysvecAssemblerFBlock(SysvecAssembler(zero(Float64)), row_nfreedofs)
end

"""
    startassembly!(self::SysvecAssemblerFBlock, row_nalldofs::IT) where {IT <: Integer}

Start assembly.
"""
function startassembly!(self::SysvecAssemblerFBlock, row_nalldofs::IT) where {IT<:Integer}
    startassembly!(self._a, row_nalldofs)
    return self
end

"""
    assemble!(self::SysvecAssemblerFBlock,
        vec::MV,
        dofnums::IV) where {MV, IV}

Assemble an elementwise vector.
"""
function assemble!(self::SysvecAssemblerFBlock, vec::MV, dofnums::IV) where {MV,IV}
    assemble!(self._a, vec, dofnums)
    return self
end

"""
    makevector!(self::SysvecAssemblerFBlock)

Make the global "free" vector.
"""
function makevector!(self::SysvecAssemblerFBlock)
    return vector_blocked_f(makevector!(self._a), self._row_nfreedofs)
end


end # end of module
