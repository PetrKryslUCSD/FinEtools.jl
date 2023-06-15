# using Pkg;
# Pkg.activate()
# Pkg.instantiate()

using Documenter, FinEtools

makedocs(
    modules = [FinEtools],
    doctest = false,
    clean = true,
    format = Documenter.HTML(prettyurls = false),
    authors = "Petr Krysl",
    sitename = "FinEtools.jl",
    pages = Any[
        "Home"=>"index.md",
        "How to"=>"howto/howto.md",
        "Tutorials"=>"tutorials/tutorials.md",
        "Concepts"=>"concepts/concepts.md",
        "Reference"=>"man/man.md",
    ],
)

deploydocs(repo = "github.com/PetrKryslUCSD/FinEtools.jl.git", devbranch = "experimental")
