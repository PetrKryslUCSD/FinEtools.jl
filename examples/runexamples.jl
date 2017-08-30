function runexamples()
    origcwd = pwd()
    packagedir = Pkg.dir("FinEtools")
    cd(packagedir)
    exds = readdir("./examples")

    function doex(ex1)
        include(ex1)
    end

    for exd1 in exds
        cd("./examples/" * exd1);
        println("\nExamples in folder $(pwd())\n")
        exss = readdir(".")
        for ex1 in exss
            n,ext = splitext(ex1)
            if ext==".jl"
                println("\nExample $ex1 in $(pwd()) ($n $ext)\n")
                doex(ex1)
            end
        end
        cd(packagedir)
    end

    cd(origcwd)
end
