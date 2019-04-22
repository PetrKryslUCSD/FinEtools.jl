using Documenter, FinEtools

makedocs(
    format = :html,
    sitename = "FinEtools.jl",
    modules = [FinEtools],
    pages = [
        "index.md",
    ]
)

deploydocs(
    repo = "github.com/JuliaPDE/FinEtools.jl.git",
    target = "build",
    julia  = "nightly",
    deps = nothing,
    make = nothing
)
