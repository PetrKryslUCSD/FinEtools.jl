using Glob

"""
    catexamples(output, pattern)

Generate a module with examples  on a common theme.  
"""
function catexamples(output, pattern)
    n,ext = splitext(output)
    if ext == ".jl"
        outputname = n
    else 
        outputname = output
        output = output * ".jl"
    end 
    List =  Glob.glob(pattern * ".jl")
    Functions = []
    open(output, "w") do o
        @printf o "module %s\n" outputname
        for ex1 in List
            n,ext = splitext(ex1)
            if ext == ".jl" && n != outputname 
                println("Example $ex1 in $(pwd())\n")
                open(ex1, "r") do f
                    text = read(f, String)
                    @printf o "\n" 
                    @printf o "function %s()\n" n
                    @printf o "%s\n" text
                    @printf o "end # %s\n\n " n
                end
                push!(Functions, n)
            end
        end
        @printf o "function allrun()\n" 
        for foo in Functions 
            @printf o """println("#####################################################") \n"""
            @printf o """println("# %s ")\n""" foo
            @printf o "    %s()\n" foo
        end
        @printf o "    return true\n" 
        @printf o "end # function allrun\n\n" 
        @printf o "end # module %s\n" outputname
    end
    println("Wrote $output in $(pwd())\n")
end
