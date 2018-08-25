using Literate

for tut1 in readdir(".")
    if occursin(r".*-tutorial.jl", tut1)
        # println("\nExample $ex1 in $(pwd())\n")
        Literate.markdown(tut1, "."; documenter=false);
    end
end
