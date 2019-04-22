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
    repo = "https://github.com/PetrKryslUCSD/FinEtools.jl.git",
    target = "build",
    julia  = "nightly",
    deps = nothing,
    make = nothing
)
