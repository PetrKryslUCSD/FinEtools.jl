"""
    AssemblyModule

Module for assemblers  of system matrices and vectors.
"""
module AssemblyModule

__precompile__(true)

using SparseArrays: sparse, spzeros, SparseMatrixCSC
using LinearAlgebra: diag

"""
    AbstractSysmatAssembler

Abstract type of system-matrix assembler.
"""
abstract type AbstractSysmatAssembler end

"""
    matrix_blocked(S, row_nfreedofs, col_nfreedofs)

Partition matrix into blocks.

The function returns the sparse matrix as a named tuple of its constituent
blocks. The matrix is composed out of four blocks
```
A = [A_ff A_fd
     A_df A_dd]
```
which are returned as a named tuple `(ff = A_ff, fd = A_fd, df = A_df, dd = A_dd)`.
Here `f` stands for free, and `d` stands for data (i.e. fixed, prescribed, ...).
The size of the `ff` block is `row_nfreedofs, col_nfreedofs`.
"""
function matrix_blocked(S, row_nfreedofs, col_nfreedofs)
    row_nalldofs, col_nalldofs = size(S)
    row_nfreedofs < row_nalldofs || error("The ff block has too many rows")
    col_nfreedofs < col_nalldofs || error("The ff block has too many columns")
    row_f_dim = row_nfreedofs
    row_d_dim = (row_nfreedofs < row_nalldofs ? row_nalldofs - row_nfreedofs - 1 : 0)
    col_f_dim = col_nfreedofs
    col_d_dim = (col_nfreedofs < col_nalldofs ? col_nalldofs - col_nfreedofs - 1 : 0)

    if row_f_dim > 0 &&  col_f_dim > 0
        S_ff = S[1:row_nfreedofs, 1:col_nfreedofs]
    else
        S_ff = spzeros(row_f_dim, col_f_dim)
    end
    if row_f_dim > 0 && col_d_dim > 0
        S_fd = S[1:row_nfreedofs, col_nfreedofs+1:end]
    else
        S_fd = spzeros(row_f_dim, col_d_dim)
    end
    if row_d_dim > 0 && col_f_dim > 0
        S_df = S[row_nfreedofs+1:end, 1:col_nfreedofs]
    else
        S_df = spzeros(row_d_dim, col_f_dim)
    end
    if row_d_dim > 0 && col_d_dim > 0
        S_dd = S[row_nfreedofs+1:end, col_nfreedofs+1:end]
    else
        S_dd = spzeros(row_d_dim, col_d_dim)
    end
    return (ff = S_ff, fd = S_fd, df = S_df, dd = S_dd)
end

"""
    vector_blocked(V, row_nfreedofs)

Partition vector into two pieces.

The vector is composed out of two blocks
```
V = [V_f
     V_d]
```
which are returned as a named tuple `(f = V_f, d = V_d)`.
"""
function vector_blocked(V, row_nfreedofs)
    row_nalldofs = length(V)
    row_nfreedofs < row_nalldofs || error("The f block has too many rows")
    row_f_dim = row_nfreedofs
    row_d_dim = (row_nfreedofs < row_nalldofs ? row_nalldofs - row_nfreedofs - 1 : 0)
    if row_f_dim > 0
        V_f = V[1:row_nfreedofs]
    else
        V_f = eltype(V)[]
    end
    if row_d_dim > 0 && col_f_dim > 0
        V_d = V[row_nfreedofs+1:end]
    else
        V_d = eltype(V)[]
    end
    return (f = V_f, d = V_d)
end

"""
    SysmatAssemblerSparse{IT, MBT, IBT} <: AbstractSysmatAssembler

Type for assembling a sparse global matrix from elementwise matrices.

!!! note

    All fields of the datatype are private. The type is manipulated by the
    functions `startassembly!`, `assemble!`, and `makematrix!`.
"""
mutable struct SysmatAssemblerSparse{IT, MBT, IBT} <: AbstractSysmatAssembler
    buffer_length::IT
    matbuffer::MBT
    rowbuffer::IBT
    colbuffer::IBT
    buffer_pointer::IT
    row_nalldofs::IT
    col_nalldofs::IT
    nomatrixresult::Bool
    force_init::Bool
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
    startassembly!(a, 5*5*3, 7, 7)
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
    startassembly!(a, 5*5*3, 7, 7)
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
```:1
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
        expected_ntriples::IT,
        row_nalldofs::IT,
        col_nalldofs::IT;
        force_init = false
        ) where {T, IT<:Integer}

Start the assembly of a global matrix.

The method makes buffers for matrix assembly. It must be called before
the first call to the method `assemble!`.

# Arguments
- `expected_ntriples`= expected number of coordinate triples, i.e. the number of
  terms in the `I`, `J`, `V` vectors;
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
    After the assembly, only the `(self.buffer_pointer - 1)` entries
    are meaningful numbers. Beware!

!!! note

    The buffers may be automatically resized if more numbers are assembled
    then initially intended. However, this operation will not necessarily
    be efficient and fast.
"""
function startassembly!(self::SysmatAssemblerSparse{T},
    expected_ntriples::IT,
    row_nalldofs::IT,
    col_nalldofs::IT;
    force_init = false
    ) where {T, IT<:Integer}
    # Only resize the buffers if the pointer is less than 1; otherwise the
    # buffers are already initialized and in use.
    if self.buffer_pointer < 1
        self.buffer_length = expected_ntriples
        resize!(self.rowbuffer, self.buffer_length)
        resize!(self.colbuffer, self.buffer_length)
        resize!(self.matbuffer, self.buffer_length)
        self.buffer_pointer = 1
        self.row_nalldofs = row_nalldofs
        self.col_nalldofs = col_nalldofs
    end
    # Leave the buffers uninitialized, unless the user requests otherwise
    self.force_init = force_init
    if self.force_init
        self.rowbuffer .= 1
        self.colbuffer .= 1
        self.matbuffer .= zero(T)
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
    mat::MT,
    dofnums_row::IT,
    dofnums_col::IT,
) where {MT, IT}
    # Assembly of a rectangular matrix.
    # The method assembles a rectangular matrix using the two vectors of
    # equation numbers for the rows and columns.
    nrows = length(dofnums_row)
    ncolumns = length(dofnums_col)
    p = self.buffer_pointer
    if p + ncolumns * nrows >= self.buffer_length
        new_buffer_length = self.buffer_length + ncolumns * nrows * 1000
        resize!(self.rowbuffer, new_buffer_length)
        resize!(self.colbuffer, new_buffer_length)
        resize!(self.matbuffer, new_buffer_length)
        if self.force_init
            self.rowbuffer[self.buffer_length+1:end] .= 1
            self.colbuffer[self.buffer_length+1:end] .= 1
            self.matbuffer[self.buffer_length+1:end] .= 0
        end
        self.buffer_length = new_buffer_length
    end
    @assert size(mat) == (nrows, ncolumns)
    @inbounds for j in 1:ncolumns
        dj = dofnums_col[j]
        dj < 1 && error("Column degree of freedom < 1")
        dj > self.col_nalldofs && error("Column degree of freedom > size")
        for i in 1:nrows
            di = dofnums_row[i]
            di < 1 && error("Row degree of freedom < 1")
            di > self.row_nalldofs && error("Row degree of freedom > size")
            self.matbuffer[p] = mat[i, j] # serialized matrix
            self.rowbuffer[p] = dofnums_row[i]
            self.colbuffer[p] = dofnums_col[j]
            p = p + 1
        end
    end
    self.buffer_pointer = p
    return self
end

"""
    makematrix!(self::SysmatAssemblerSparse)

Make a sparse matrix.

A named tuple of the four matrix blocks is returned, as described in the
documentation for the constructor.

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
    @assert length(self.rowbuffer) >= self.buffer_pointer - 1
    @assert length(self.colbuffer) >= self.buffer_pointer - 1
    # We have the option of retaining the assembled results, but not
    # constructing the sparse matrix.
    if self.nomatrixresult
        # Dummy (zero) sparse matrix is returned. The entire result of the
        # assembly is preserved in the assembler buffers. The ends of the
        # buffers are filled with legal (ignorable) values.
        self.rowbuffer[self.buffer_pointer:end] .= 1
        self.colbuffer[self.buffer_pointer:end] .= 1
        self.matbuffer[self.buffer_pointer:end] .= 0.0
        return spzeros(self.row_nalldofs, self.col_nalldofs)
    end
    # The sparse matrix is constructed and returned. The  buffers used for
    # the assembly are cleared.
    S = sparse(
        view(self.rowbuffer, 1:self.buffer_pointer-1),
        view(self.colbuffer, 1:self.buffer_pointer-1),
        view(self.matbuffer, 1:self.buffer_pointer-1),
        self.row_nalldofs,
        self.col_nalldofs,
        )
    # Get ready for more assembling
    self.buffer_pointer = 1
    # Construct the blocks of the matrix
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
mutable struct SysmatAssemblerSparseSymm{IT, MBT, IBT} <: AbstractSysmatAssembler
    buffer_length::IT
    matbuffer::MBT
    rowbuffer::IBT
    colbuffer::IBT
    buffer_pointer::IT
    row_nalldofs::IT
    col_nalldofs::IT
    nomatrixresult::Bool
    force_init::Bool
end

"""
    SysmatAssemblerSparseSymm(z::T, nomatrixresult = false) where {T}

Construct blank system matrix assembler for symmetric matrices. The matrix
entries are of type `T`.

# Example

This is how a symmetric sparse matrix is assembled from two square dense matrices.
```
    a = SysmatAssemblerSparseSymm(0.0)
    startassembly!(a, 5*5*3, 7, 7)
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
    return SysmatAssemblerSparseSymm(0, T[z], Int[0], Int[0], 0, 0, 0, nomatrixresult, false)
end

function SysmatAssemblerSparseSymm()
    return SysmatAssemblerSparseSymm(zero(Float64))
end

"""
    startassembly!(self::SysmatAssemblerSparseSymm{T},
        expected_ntriples::IT,
        row_nalldofs::IT,
        col_nalldofs::IT;
        force_init = false
        ) where {T, IT<:Integer}

Start the assembly of a global matrix.

The method makes buffers for matrix assembly. It must be called before
the first call to the method `assemble!`.

# Arguments
- `expected_ntriples`= expected number of coordinate triples, i.e. the number of
  terms in the `I`, `J`, `V` vectors;
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
    After the assembly, only the `(self.buffer_pointer - 1)` entries
    are meaningful numbers. Beware!
"""
function startassembly!(self::SysmatAssemblerSparseSymm{T},
    expected_ntriples::IT,
    row_nalldofs::IT,
    col_nalldofs::IT;
    force_init = false
    ) where {T, IT<:Integer}
    # Only resize the buffers if the pointer is less than 1; otherwise the
    # buffers are already initialized and in use.
    if self.buffer_pointer < 1
        self.buffer_length = expected_ntriples
        resize!(self.rowbuffer, self.buffer_length)
        resize!(self.colbuffer, self.buffer_length)
        resize!(self.matbuffer, self.buffer_length)
        self.buffer_pointer = 1
        row_nalldofs == col_nalldofs || error("Row and column info do not agree")
        self.row_nalldofs = row_nalldofs
        self.col_nalldofs = col_nalldofs
    end
    # Leave the buffers uninitialized, unless the user requests otherwise
    if force_init
        self.rowbuffer .= 1
        self.colbuffer .= 1
        self.matbuffer .= 0
        self.force_init = force_init
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
    mat::MT,
    dofnums_row::IT,
    dofnums_col::IT,
) where {MT, IT}
    # Assembly of a square symmetric matrix.
    # The method assembles the lower triangle of the square symmetric matrix using the two vectors of
    # equation numbers for the rows and columns.
    nrows = length(dofnums_row)
    ncolumns = length(dofnums_col)
    @assert nrows == ncolumns
    p = self.buffer_pointer
    if p + ncolumns * nrows >= self.buffer_length
        new_buffer_length = self.buffer_length + ncolumns * nrows * 1000
        resize!(self.rowbuffer, new_buffer_length)
        resize!(self.colbuffer, new_buffer_length)
        resize!(self.matbuffer, new_buffer_length)
        if self.force_init
            self.rowbuffer[self.buffer_length+1:end] .= 1
            self.colbuffer[self.buffer_length+1:end] .= 1
            self.matbuffer[self.buffer_length+1:end] .= 0
        end
        self.buffer_length = new_buffer_length
    end
    @assert size(mat) == (nrows, ncolumns)
    @inbounds for j in 1:ncolumns
        dj = dofnums_col[j]
        dj < 1 && error("Column degree of freedom < 1")
        dj > self.col_nalldofs && error("Column degree of freedom > size")
        for i in j:nrows
            di = dofnums_row[i]
            di < 1 && error("Row degree of freedom < 1")
            di > self.row_nalldofs && error("Row degree of freedom > size")
            self.matbuffer[p] = mat[i, j] # serialized matrix
            self.rowbuffer[p] = dofnums_row[i]
            self.colbuffer[p] = dofnums_col[j]
            p = p + 1
        end
    end
    self.buffer_pointer = p
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
    @assert length(self.rowbuffer) >= self.buffer_pointer - 1
    @assert length(self.colbuffer) >= self.buffer_pointer - 1
    # We have the option of retaining the assembled results, but not
    # constructing the sparse matrix.
    if self.nomatrixresult
        # Dummy (zero) sparse matrix is returned. The entire result of the
        # assembly is preserved in the assembler buffers. The ends of the
        # buffers are filled with legal (ignorable) values.
        self.rowbuffer[self.buffer_pointer:end] .= 1
        self.colbuffer[self.buffer_pointer:end] .= 1
        self.matbuffer[self.buffer_pointer:end] .= 0.0
        return spzeros(self.row_nalldofs, self.col_nalldofs)
    end
    # The sparse matrix is constructed and returned. The  buffers used for
    # the assembly are cleared.
    S = sparse(
        view(self.rowbuffer, 1:self.buffer_pointer-1),
        view(self.colbuffer, 1:self.buffer_pointer-1),
        view(self.matbuffer, 1:self.buffer_pointer-1),
        self.row_nalldofs,
        self.col_nalldofs,
    )
    #  Now we need to construct the other triangle of the matrix. The diagonal
    #  will be duplicated.
    S = S + transpose(S)
    @inbounds for j in axes(S, 1)
        S[j, j] *= 0.5      # the diagonal is there twice; fix it;
    end
    # Get ready for more assembling
    self.buffer_pointer = 1
    # Construct the blocks of the matrix
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
mutable struct SysmatAssemblerSparseDiag{IT, MBT, IBT} <: AbstractSysmatAssembler
    # Type for assembling of a sparse global matrix from elementwise matrices.
    buffer_length::IT
    matbuffer::MBT
    rowbuffer::IBT
    colbuffer::IBT
    buffer_pointer::IT
    row_nalldofs::IT
    col_nalldofs::IT
    nomatrixresult::Bool
    force_init::Bool
end

"""
    SysmatAssemblerSparseDiag(z::T, nomatrixresult = false) where {T}

Construct blank system matrix assembler for square diagonal matrices. The matrix
entries are of type `T`.

"""
function SysmatAssemblerSparseDiag(z::T, nomatrixresult = false) where {T}
    return SysmatAssemblerSparseDiag(0, T[z], Int[0], Int[0], 0, 0, 0, nomatrixresult, false)
end

function SysmatAssemblerSparseDiag()
    return SysmatAssemblerSparseDiag(zero(Float64))
end

"""
    startassembly!(self::SysmatAssemblerSparseDiag{T},
        expected_ntriples::IT,
        row_nalldofs::IT,
        col_nalldofs::IT;
        force_init = false
        ) where {T, IT<:Integer}

Start the assembly of a symmetric square diagonal matrix.

The method makes buffers for matrix assembly. It must be called before
the first call to the method `assemble!`.

# Arguments
- `expected_ntriples`= expected number of coordinate triples, i.e. the number of
  terms in the `I`, `J`, `V` vectors;
- `row_nalldofs`= The total number of rows as a tuple;
- `col_nalldofs`= The total number of columns as a tuple.

The values stored in the buffers are initially undefined!

# Returns
- `self`: the modified assembler.
"""
function startassembly!(self::SysmatAssemblerSparseDiag{T},
    expected_ntriples::IT,
    row_nalldofs::IT,
    col_nalldofs::IT;
    force_init = false
    ) where {T, IT<:Integer}
    # Only resize the buffers if the pointer is less than 1; otherwise the
    # buffers are already initialized and in use.
    if self.buffer_pointer < 1
        self.buffer_length = expected_ntriples
        resize!(self.rowbuffer, self.buffer_length)
        resize!(self.colbuffer, self.buffer_length)
        resize!(self.matbuffer, self.buffer_length)
        self.buffer_pointer = 1
        row_nalldofs == col_nalldofs || error("Row and column info do not agree")
        self.row_nalldofs = row_nalldofs
        self.col_nalldofs = col_nalldofs
    end
    # Leave the buffers uninitialized, unless the user requests otherwise
    if force_init
        self.rowbuffer .= 1
        self.colbuffer .= 1
        self.matbuffer .= 0
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
    mat::MT,
    dofnums_row::IV,
    dofnums_col::IV,
) where {MT, IV}
    # Assembly of a square symmetric matrix.
    # The method assembles the lower triangle of the square symmetric matrix using the two vectors of
    # equation numbers for the rows and columns.
    nrows = length(dofnums_row)
    ncolumns = length(dofnums_col)
    @assert nrows == ncolumns
    p = self.buffer_pointer
    @assert p + ncolumns <= self.buffer_length + 1
    @assert size(mat) == (nrows, ncolumns)
    @inbounds for j in 1:ncolumns
        dj = dofnums_col[j]
        dj < 1 && error("Column degree of freedom < 1")
        dj > self.col_nalldofs && error("Column degree of freedom > size")
        self.matbuffer[p] = mat[j, j] # serialized matrix
        self.rowbuffer[p] = dj
        self.colbuffer[p] = dj
        p = p + 1
    end
    self.buffer_pointer = p
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
    @assert length(self.rowbuffer) >= self.buffer_pointer - 1
    @assert length(self.colbuffer) >= self.buffer_pointer - 1
    # We have the option of retaining the assembled results, but not
    # constructing the sparse matrix.
    if self.nomatrixresult
        # Dummy (zero) sparse matrix is returned. The entire result of the
        # assembly is preserved in the assembler buffers. The ends of the
        # buffers are filled with illegal (ignorable) values.
        self.rowbuffer[self.buffer_pointer:end] .= 0
        self.colbuffer[self.buffer_pointer:end] .= 0
        self.matbuffer[self.buffer_pointer:end] .= 0.0
        return spzeros(self.row_nalldofs, self.col_nalldofs)
    end
    # The sparse matrix is constructed and returned. The  buffers used for
    # the assembly are cleared.
    S = sparse(
        view(self.rowbuffer, 1:self.buffer_pointer-1),
        view(self.colbuffer, 1:self.buffer_pointer-1),
        view(self.matbuffer, 1:self.buffer_pointer-1),
        self.row_nalldofs,
        self.col_nalldofs,
    )
    # Get ready for more assembling
    self.buffer_pointer = 1
    # Construct the blocks of the matrix
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
function startassembly!(self::SV, row_nalldofs::Tuple{IT, IT}) where {SV<:AbstractSysvecAssembler, IT} end

"""
    assemble!(self::SysvecAssembler{T}, vec::MV,
      dofnums::D) where {T<:Number, MV<:AbstractArray{T}, D<:AbstractArray{FInt}}

Assemble an elementwise vector.

The method assembles a column element vector using the vector of degree of
freedom numbers for the rows.
"""
function assemble!(
    self::SV,
    vec::MV,
    dofnums::IV,
) where {SV<:AbstractSysvecAssembler,MV,IV} end

"""
    makevector!(self::SysvecAssembler)

Make the global vector.
"""
function makevector!(self::SV) where {SV<:AbstractSysvecAssembler} end

"""
    SysvecAssembler

Assembler for the system vector.
"""
mutable struct SysvecAssembler{VBT, IT} <: AbstractSysvecAssembler
    F_buffer::VBT
    row_nalldofs::IT
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
    self.row_nalldofs = row_nalldofs
    resize!(self.F_buffer, self.row_nalldofs)
    self.F_buffer .= 0
    return self
end

"""
    assemble!(self::SysvecAssembler{T}, vec::MV,
      dofnums::D) where {T<:Number, MV<:AbstractArray{T}, D<:AbstractArray{FInt}}

Assemble an elementwise vector.

The method assembles a column element vector using the vector of degree of
freedom numbers for the rows.
"""
function assemble!(
    self::SysvecAssembler,
    vec::MV,
    dofnums::IV,
) where {MV, IV}
    @inbounds for i in eachindex(dofnums)
        gi = dofnums[i]
        gi < 1 && error("Row degree of freedom < 1")
        gi > self.row_nalldofs && error("Row degree of freedom > size")
        self.F_buffer[gi] += vec[i]
    end
end

"""
    makevector!(self::SysvecAssembler)

Make the global vector.
"""
function makevector!(self::SysvecAssembler)
    return deepcopy(self.F_buffer)
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
mutable struct SysmatAssemblerSparseHRZLumpingSymm{IT, MBT, IBT} <: AbstractSysmatAssembler
    # Type for assembling of a sparse global matrix from elementwise matrices.
    buffer_length::IT
    matbuffer::MBT
    rowbuffer::IBT
    colbuffer::IBT
    buffer_pointer::IT
    row_nalldofs::IT
    col_nalldofs::IT
    nomatrixresult::Bool
    force_init::Bool
end

"""
    SysmatAssemblerSparseHRZLumpingSymm(z::T, nomatrixresult = false) where {T}

Construct blank system matrix assembler. The matrix entries are of type `T`.

"""
function SysmatAssemblerSparseHRZLumpingSymm(z::T, nomatrixresult = false) where {T}
    return SysmatAssemblerSparseHRZLumpingSymm(0, T[z], Int[0], Int[0], 0, 0, 0, nomatrixresult, false)
end

function SysmatAssemblerSparseHRZLumpingSymm()
    return SysmatAssemblerSparseHRZLumpingSymm(zero(Float64))
end

"""
    startassembly!(self::SysmatAssemblerSparseHRZLumpingSymm{T},
            expected_ntriples::IT,
            row_nalldofs::IT,
            col_nalldofs::IT;
            force_init = false
            ) where {T, IT<:Integer}

Start the assembly of a symmetric lumped diagonal square global matrix.

The method makes buffers for matrix assembly. It must be called before
the first call to the method `assemble!`.

# Arguments
- `expected_ntriples`= expected number of coordinate triples, i.e. the number of
  terms in the `I`, `J`, `V` vectors;
- `row_nalldofs`= The total number of rows as a tuple;
- `col_nalldofs`= The total number of columns as a tuple.

The values stored in the buffers are initially undefined!

# Returns
- `self`: the modified assembler.
"""
function startassembly!(self::SysmatAssemblerSparseHRZLumpingSymm{T},
        expected_ntriples::IT,
        row_nalldofs::IT,
        col_nalldofs::IT;
        force_init = false
        ) where {T, IT<:Integer}
    # Only resize the buffers if the pointer is less than 1; otherwise the
    # buffers are already initialized and in use.
    if self.buffer_pointer < 1
        self.buffer_length = expected_ntriples
        resize!(self.rowbuffer, self.buffer_length)
        resize!(self.colbuffer, self.buffer_length)
        resize!(self.matbuffer, self.buffer_length)
        self.buffer_pointer = 1
        row_nalldofs == col_nalldofs || error("Row and column info do not agree")
        self.row_nalldofs = row_nalldofs
        self.col_nalldofs = col_nalldofs
    end
    # Leave the buffers uninitialized, unless the user requests otherwise
    if force_init
        self.rowbuffer .= 1
        self.colbuffer .= 1
        self.matbuffer .= 0
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
    mat::MT,
    dofnums_row::IV,
    dofnums_col::IV,
) where {MT, IV}
    nrows = length(dofnums_row)
    ncolumns = length(dofnums_col)
    @assert nrows == ncolumns
    p = self.buffer_pointer
    if p + ncolumns * nrows >= self.buffer_length
        new_buffer_length = self.buffer_length + ncolumns * nrows * 1000
        resize!(self.rowbuffer, new_buffer_length)
        resize!(self.colbuffer, new_buffer_length)
        resize!(self.matbuffer, new_buffer_length)
        if self.force_init
            self.rowbuffer[self.buffer_length+1:end] .= 1
            self.colbuffer[self.buffer_length+1:end] .= 1
            self.matbuffer[self.buffer_length+1:end] .= 0
        end
        self.buffer_length = new_buffer_length
    end
    @assert size(mat) == (nrows, ncolumns)
    # Now comes the lumping procedure
    em2 = sum(sum(mat, dims = 1)) # total mass times the number of space dimensions
    dem2 = zero(eltype(mat)) # total mass on the diagonal times the number of space dimensions
    @inbounds for i in 1:nrows
        dem2 += mat[i, i]
    end
    ffactor = em2 / dem2 # total-element-mass compensation factor
    @inbounds for j in 1:ncolumns
        dj = dofnums_col[j]
        dj < 1 && error("Column degree of freedom < 1")
        dj > self.col_nalldofs && error("Column degree of freedom > size")
        self.matbuffer[p] = mat[j, j] * ffactor # serialized matrix
        self.rowbuffer[p] = dj
        self.colbuffer[p] = dj
        p = p + 1
    end
    self.buffer_pointer = p
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
    @assert length(self.rowbuffer) >= self.buffer_pointer - 1
    @assert length(self.colbuffer) >= self.buffer_pointer - 1
    # We have the option of retaining the assembled results, but not
    # constructing the sparse matrix.
    if self.nomatrixresult
        # Dummy (zero) sparse matrix is returned. The entire result of the
        # assembly is preserved in the assembler buffers. The ends of the
        # buffers are filled with illegal (ignorable) values.
        self.rowbuffer[self.buffer_pointer:end] .= 0
        self.colbuffer[self.buffer_pointer:end] .= 0
        self.matbuffer[self.buffer_pointer:end] .= 0.0
        return spzeros(self.row_nalldofs, self.col_nalldofs)
    end
    # The sparse matrix is constructed and returned. The  buffers used for
    # the assembly are cleared.
    S = sparse(
        view(self.rowbuffer, 1:self.buffer_pointer-1),
        view(self.colbuffer, 1:self.buffer_pointer-1),
        view(self.matbuffer, 1:self.buffer_pointer-1),
        self.row_nalldofs,
        self.col_nalldofs,
    )
    # Get ready for more assembling
    self.buffer_pointer = 1
    # Construct the blocks of the matrix
    return S
end


"""
    SysmatAssemblerReduced{T<:Number} <: AbstractSysmatAssembler

Type for assembling a sparse global matrix from elementwise matrices.

!!! note

    All fields of the datatype are private. The type is manipulated by the
    functions `startassembly!`, `assemble!`, and `makematrix!`.
"""
mutable struct SysmatAssemblerReduced{MT, TMT, IT} <: AbstractSysmatAssembler
    m::MT # reduced system matrix
    ndofs_row::IT
    ndofs_col::IT
    red_ndofs_row::IT
    red_ndofs_col::IT
    t::Matrix{TMT} # transformation matrix
    nomatrixresult::Bool
end

function SysmatAssemblerReduced(t::TMT, z = zero(Float64), nomatrixresult = false) where {TMT}
    ndofs_row = ndofs_col = size(t, 1)
    red_ndofs_row = red_ndofs_col = size(t, 2)
    m = fill(z, red_ndofs_row, red_ndofs_col)
    return SysmatAssemblerReduced(
        m,
        ndofs_row,
        ndofs_col,
        red_ndofs_row,
        red_ndofs_col,
        t,
        nomatrixresult
    )
end


function startassembly!(
    self::SysmatAssemblerReduced,
    elem_mat_nrows,
    elem_mat_ncols,
    elem_mat_nmatrices,
    ndofs_row,
    ndofs_col,
)
    @assert self.ndofs_row == ndofs_row
    @assert self.ndofs_col == ndofs_col
    self.m .= zero(eltype(self.m))
    return self
end

function assemble!(
    self::SysmatAssemblerReduced,
    mat::MT,
    dofnums_row::IV,
    dofnums_col::IV,
)  where {MT, IV}
    R, C = size(mat)
    T = eltype(mat)
    @assert R == C "The matrix must be square"
    @assert dofnums_row == dofnums_col "The degree of freedom numbers must be the same for rows and columns"
    for i in 1:R
        gi = dofnums_row[i]
        if gi < 1
            mat[i, :] .= zero(T)
            dofnums_row[i] = 1
        end
    end
    for i in 1:C
        gi = dofnums_col[i]
        if gi < 1
            mat[:, i] .= zero(T)
            dofnums_col[i] = 1
        end
    end
    lt = self.t[dofnums_row, :]
    self.m .+= lt' * mat * lt
    return self
end

"""
    makematrix!(self::SysmatAssemblerReduced)

Make a sparse matrix.
"""
function makematrix!(self::SysmatAssemblerReduced)
    return self.m
end


end # end of module
