module ParallelAssemblyModule

using ..FEMMBaseModule: AbstractFEMM, finite_elements, iselementbased
using ..FieldModule: ndofs, nalldofs
using ..NodalFieldModule: NodalField
using ..AssemblyModule:
    AbstractSysmatAssembler,
    startassembly!
using ..FESetModule: nodesperelem

function _update_buffer_range(buffer_length, iend)
    istart = iend + 1
    iend = iend + buffer_length
    buffer_range = istart:iend
    return buffer_range, iend
end

function _task_local_assembler(a::AT, buffer_range) where {AT}
    buffer_length = maximum(buffer_range) - minimum(buffer_range) + 1
    matbuffer = view(a.matbuffer, buffer_range)
    rowbuffer = view(a.rowbuffer, buffer_range)
    colbuffer = view(a.colbuffer, buffer_range)
    @show typeof(matbuffer)
    buffer_pointer = 1
    nomatrixresult = true
    force_init = false
    _a =  AT(
        buffer_length,
        matbuffer,
        rowbuffer,
        colbuffer,
        buffer_pointer,
        a.row_nalldofs,
        a.col_nalldofs,
        nomatrixresult,
        force_init,
    )
    @show typeof(_a)
    return _a
end

"""
    parallel_tasks_matrix_assembly(
        femms::AbstractArray{FEMM, 1},
        assembler::AT,
        uf::NodalField{UFT},
        matrixcomputation!::F
    ) where {FEMM<:AbstractFEMM, AT<:AbstractSysmatAssembler, F<:Function}


"""
function parallel_matrix_assembly(
    femms::AbstractArray{FEMM, 1},
    assembler::AT,
    uf::NodalField{UFT},
    matrixcomputation!::F
) where {FEMM<:AbstractFEMM, AT<:AbstractSysmatAssembler, UFT<:Number, F<:Function}
    # The number of tasks is equal to the number of the FEMMs
    ntasks = length(femms)
    ntasks >= 1 || error("Number of tasks must be >= 1")
    elem_mat_nmatrices = count.(finite_elements.(femms))
    # Check compatibility of the FEMMs: they need to be the same "type" -- all
    # conductivity, for instance, and all using the same type of finite element 
    for j in eachindex(femms)
        # Check that the femms are element-based
        # TO DO
        iselementbased(femms[j]) || error("FEMM is not element-based")
        (supertype(typeof(femms[j])) === supertype(typeof(femms[1]))) || 
            error("All FEMMs must be of the same type")
        (typeof(finite_elements(femms[j])) === typeof(finite_elements(femms[1]))) || 
            error("All finite elements must be of the same type")
    end
    # Count the total number of COO triples
    nne = nodesperelem(finite_elements(femms[1]))
    ndn = ndofs(uf)
    elem_mat_nrows = nne * ndn
    elem_mat_ncols = nne * ndn
    triples = elem_mat_nrows * elem_mat_ncols .* elem_mat_nmatrices
    totaltriples = sum(triples)
    ndofs_row = nalldofs(uf)
    ndofs_col = nalldofs(uf)
    # Started the assembly
    startassembly!(assembler, totaltriples, ndofs_row, ndofs_col,)
    # Parallel assembly loop
    iend = 0
    Threads.@sync begin
        for j in eachindex(femms)
            Threads.@spawn let 
                buffer_range, iend = _update_buffer_range(triples[j], iend)
                _a = _task_local_assembler(assembler, buffer_range)
                @show typeof(_a)
                matrixcomputation!(femms[j], _a, uf)
            end
        end
    end
    # The buffer should now be completely full
    iend == assembler.buffer_length || error("Buffer pointer mismatch") 
    assembler.buffer_pointer = iend
    # Return the updated assembler; makematrix! can be called on that to retrieve the matrix
    assembler.nomatrixresult = !true
    return assembler
end

end # module