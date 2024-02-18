module ParallelAssemblyModule

using ..FEMMBaseModule: AbstractFEMM, finite_elements, iselementbased
using ..FieldModule: ndofs, nalldofs
using ..NodalFieldModule: NodalField
using ..AssemblyModule:
    AbstractSysmatAssembler,
    startassembly!, makematrix!
using ..FESetModule: nodesperelem

function _update_buffer_range(buffer_length, iend)
    istart = iend + 1
    iend = iend + buffer_length
    buffer_range = istart:iend
    return buffer_range, iend
end

function _task_local_assembler(AssemblerType::Type{AT}, buffer_range, matbuffer, rowbuffer, colbuffer, row_nalldofs, col_nalldofs) where {AT}
    buffer_length = maximum(buffer_range) - minimum(buffer_range) + 1
    buffer_pointer = 1
    nomatrixresult = true
    force_init = false
    return AssemblerType(
        buffer_length,
        view(matbuffer, buffer_range),
        view(rowbuffer, buffer_range),
        view(colbuffer, buffer_range),
        buffer_pointer,
        row_nalldofs,
        col_nalldofs,
        nomatrixresult,
        force_init,
    )
end

function _check_femm_compatibility(femms::AbstractArray{FEMM, 1}) where {FEMM<:AbstractFEMM}
    for j in eachindex(femms)
        iselementbased(femms[j]) || error("FEMM is not element-based")
        (nameof(typeof(femms[j])) ===
         nameof(typeof(femms[1]))) ||
            error("All FEMMs must be of the same type")
        (nameof(typeof(finite_elements(femms[j]))) ===
         nameof(typeof(finite_elements(femms[1])))) ||
            error("All finite elements must be of the same type")
    end
    return true
end

"""
    make_assembler(
        femms::AbstractArray{FEMM, 1},
        AssemblerType::Type{AT},
        uf::NodalField{UFT},
    ) where {FEMM<:AbstractFEMM, AT<:AbstractSysmatAssembler, UFT<:Number}
    
Make the global assembler.    

- `femms`: array of finite element machines, one for each subdomain
- `AssemblerType`: type of the global assembler, unparametrized
- `uf`: nodal field of the unknowns
"""
function make_assembler(
    femms::AbstractArray{FEMM, 1},
    AssemblerType::Type{AT},
    uf::NodalField{UFT},
) where {FEMM<:AbstractFEMM, AT<:AbstractSysmatAssembler, UFT<:Number}
    _check_femm_compatibility(femms)
    elem_mat_nmatrices = count.(finite_elements.(femms))
    # Count the total number of COO triples
    nne = nodesperelem(finite_elements(femms[1]))
    ndn = ndofs(uf)
    elem_mat_nrows = nne * ndn
    elem_mat_ncols = nne * ndn
    triples = elem_mat_nrows * elem_mat_ncols .* elem_mat_nmatrices
    totaltriples = sum(triples)
    ndofs_row = nalldofs(uf)
    ndofs_col = nalldofs(uf)
    # Start the assembly
    assembler = startassembly!(AssemblerType(zero(UFT)), totaltriples, ndofs_row, ndofs_col,)
    assembler.nomatrixresult = true
    return assembler
end

"""
    start_assembler!(
        assembler::AT
    ) where {AT<:AbstractSysmatAssembler}

Start assembly on the global assembler.

The assembler is informed that no matrix result should be computed. All the data
should be preserved in the assembly buffers.
"""
function start_assembler!(
    assembler::AT
) where {AT<:AbstractSysmatAssembler}
    assembler = startassembly!(assembler,
        assembler.buffer_length,
        assembler.row_nalldofs,
        assembler.col_nalldofs,)
    assembler.nomatrixresult = true
    return assembler
end


"""
    make_matrix!(
        assembler::AT
    ) where {AT<:AbstractSysmatAssembler}

Make the global matrix, using the data collected in the assembler.
"""
function make_matrix!(
    assembler::AT
) where {AT<:AbstractSysmatAssembler}
    assembler.buffer_pointer = assembler.buffer_length + 1
    assembler.nomatrixresult = !true
    return makematrix!(assembler)
end

"""
    make_task_assemblers(
        femms::AbstractArray{FEMM, 1},
        assembler::AT,
        AssemblerType::Type{AT},
        uf::NodalField{UFT}
    ) where {FEMM<:AbstractFEMM, AT<:AbstractSysmatAssembler, UFT<:Number}

Make task assemblers for the subdomains.

- `femms`: array of finite element machines, one for each subdomain
- `assembler`: the global assembler 
- `AssemblerType`: type of the global assembler, unparametrized
- `uf`: nodal field of the unknowns

Return an array of the assemblers. Each of the assemblers only refers to a view
of the buffers of the global assembler, hence they are reasonably cheap to
construct.
"""
function make_task_assemblers(
    femms::AbstractArray{FEMM, 1},
    assembler::AT,
    AssemblerType::Type{AT},
    uf::NodalField{UFT}
) where {FEMM<:AbstractFEMM, AT<:AbstractSysmatAssembler, UFT<:Number}
    _check_femm_compatibility(femms)
    # The number of tasks is equal to the number of the FEMMs
    ntasks = length(femms)
    ntasks >= 1 || error("Number of tasks must be >= 1")
    elem_mat_nmatrices = count.(finite_elements.(femms))
    # Check compatibility of the FEMMs: they need to be the same "type" -- all
    # conductivity, for instance, and all using the same type of finite element 
    elem_mat_nmatrices = count.(finite_elements.(femms))
    # Count the total number of COO triples
    nne = nodesperelem(finite_elements(femms[1]))
    ndn = ndofs(uf)
    elem_mat_nrows = nne * ndn
    elem_mat_ncols = nne * ndn
    triples = elem_mat_nrows * elem_mat_ncols .* elem_mat_nmatrices
    # Make there task assemblers
    assemblers = AssemblerType[]
    iend = 0
    for j in eachindex(femms)
        buffer_range, iend = _update_buffer_range(triples[j], iend)
        push!(assemblers, _task_local_assembler(
            AssemblerType, buffer_range, 
            assembler.matbuffer, assembler.rowbuffer, assembler.colbuffer, 
            assembler.row_nalldofs, assembler.col_nalldofs)
        )
    end
    return assemblers
end

"""
    parallel_matrix_assembly(
        femms::AbstractArray{FEMM, 1},
        assemblers::AbstractVector{AT},
        matrixcomputation!::F
    ) where {FEMM<:AbstractFEMM, AT<:AbstractSysmatAssembler, F<:Function}

Execute the assembly in parallel.
"""
function parallel_matrix_assembly(
    femms::AbstractArray{FEMM, 1},
    assemblers::AbstractVector{AT},
    matrixcomputation!::F
) where {FEMM<:AbstractFEMM, AT<:AbstractSysmatAssembler, F<:Function}
    for j in eachindex(assemblers)
        assemblers[j].buffer_pointer = 1
        assemblers[j].nomatrixresult = true
    end
    Threads.@threads for j in eachindex(femms)
        matrixcomputation!(femms[j], assemblers[j])
    end
    # Threads.@sync begin
    #     for j in eachindex(femms)
    #         assemblers[j].buffer_pointer = 1
    #         assemblers[j].nomatrixresult = true
    #         Threads.@spawn let 
    #             matrixcomputation!(femms[j], assemblers[j])
    #         end
    #     end
    # end
    return true
end

end # module