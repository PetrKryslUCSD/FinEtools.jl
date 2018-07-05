"""
    runexamples()

Run all *_examples.jl files in the working folder.  
"""
function runexamples()
    origcwd = pwd()
    
    pkgdir(pkg::String) = abspath(joinpath(dirname(Base.find_package(pkg)), ".."))
   
    println("\nExamples in folder $(pwd())\n")

    for ex1 in readdir(".")
        if occursin(r".*_examples.jl", ex1)
            println("\nExample $ex1 in $(pwd())\n")
            include(ex1);
            n, ext = splitext(ex1)
            Meta.eval(Meta.parse("Main." * n * ".allrun()"))
        end
    end

    cd(origcwd)
end
