
Issues:

-- Define the labels in the finite elements that to be a vector no matter what.
```
labels = fes.label; # individual element labels
if length(fes.label) == 0
    labels = zeros(FInt, nfes);
end
```
Implemented  05/07/2017

-- Get rid of the property modules: replace with definition in the materials.
Needs careful deliberation.
Initiated 05/22/2017. Implementation started 05/24/2017.

-- In algorithms, assume that  the finite element model machine  is submitted as
input instead of the material parameters and the finite elements. Similarly to
SPyFE.
Initiated 05/02/2017. It is working in  heat diffusion 05//23/2017.

-- There are some loops  in the finite element model machines which
are not instrumented with @inbounds end @fastmath.

-- Implement coordinate system module. Replace MaterialOrientationModule
with it.
Implemented 05/05/2017.

-- How to work with anisotropic materials?  Strains in local coordinate system,
ddisplacements in the global coordinate system?  Or both in the global
coordinate system, and a transformation for the material property matrix?

-- The  base FEMM should create a coordinate system  that
allows for  pre-computation: identity matrix.

-- Special FEMM  heat transfer  on the surface.
And special FEMM for  the boundary condition on the surface?
We already have one for  the Winkler  boundary condition?
Implemented 05/22/2017..

-- The Poisson  heat conduction problem:
Conductivity
  6.043066 seconds (2.00 M allocations: 1.296 GB, 13.11% gc time)
Why are there 2 million allocations  (this is the number of triangles)?
Profiling? Could be array slicing?
Solution: It was the type instability due to the keyword arguments  used  for
the finite element sets. Keyword arguments  should be removed everywhere.

-- What is __precompile__() for?

-- What is @propagate_inbounds for?

-- 
```
julia> A = rand(3,4);

julia> B = A .- mean(A,1)
3×4 Array{Float64,2}:
  0.343976   0.427378  -0.503356  -0.00448691
 -0.210096  -0.531489   0.168928  -0.128212
 -0.13388    0.104111   0.334428   0.132699

julia> mean(B,1)
1×4 Array{Float64,2}:
 0.0  0.0  0.0  0.0
```

 -- Use return type annotations. Status: Still not sure if I want to do this in 2023.
 It might help with understanding the code from the signature?

 -- Replace explicit loops with in-place matrix mult functions? A_mult_B! and so on...
 Notes: Tested for  the heat conduction model. Around 20% slower: more
 memory allocations.
 ```
 A_mul_Bt!(kappa_bargradNT, (Jac*w[j])*kappa_bar, gradN); # intermediate result
 A_mul_B!(Ke1, gradN, kappa_bargradNT);
 Ke[:,:] += Ke1[:,:]
 ```
 Conductivity
  8.454714 seconds (10.00 M allocations: 2.398 GB, 10.42% gc time)
versus
```
A_mul_Bt!(kappa_bargradNT, (Jac*w[j])*kappa_bar, gradN); # intermediate result
A_mul_B!(Ke1, gradN, kappa_bargradNT);
@inbounds for nx = 1:Kedim # complete the lower triangle
  @inbounds for mx = 1:Kedim
    Ke[mx, nx] = Ke[mx, nx] + Ke1[nx, mx]
  end
end
```
Conductivity
  9.182616 seconds (4.00 M allocations: 1.504 GB, 10.83% gc time)
versus
All-loops
Conductivity
    6.209980 seconds (2.05 M allocations: 1.298 GB, 12.54% gc time)

-- add_nhnt!(): test and include  in the code
Implemented.

-- WARNING: both MeshTriangleModule and MeshQuadrilateralModule export "Q4toT3";
Conflict resolved.

-- For symmetric elementwise matrices compute only the upper triangle in the
integration point loop; complete the matrix only below the loop.
Alternatively don't bother completing the lower triangle if the assembler
works only on the upper triangle.
Implemented 05/1/17: compute only upper triangle in the integration group,
complete the full matrix underneath the loop. Not completing the triangle
would not work if the assembler expected the full matrix.

-- Test the speed  with the material orientation  matrix different from identity.

--
Other workarounds include "types as tags" (search; NumericFuns.jl is one example
among many), FastAnonymous.jl, and in julia 0.4 the new call() overloading.

-- Move the code for axial symmetry and other dimension
into the FEMM modules: type  parameter?

-- SymPy.jl or SymEngine.jl.

-- I think Julia is about giving you the tools to build these things.
If you want "the tool for auto-optimizing scientific computing code",
that would ParallelAccelerator.jl

-- Jacobian functions for zero-dimensional manifolds.
Implemented 05/20/2017

-- export setRm!
export getRm!  Not used?
Implemented 05/20/2017

--
@code_warntype Jacobiansurface(FEMMBase, FFltMat, FFltMat, FIntMat, FFltMat)

-- replace error() with assert in FEMMBaseModule
Implemented 05/01/2017

-- return fens,nfes            # I think I should not be overwriting the input!
What is that about?

-- FEMMBaseModule:  redesign  integration methods  by parameterising with
the function type.
Implementation completed.

-- VTK export: export  cell data.
Implemented  July 2017.

-- Implement  element  fields. (Analogy to  nodal fields.)
Implemented partially: so far only elementwise constant field has been
implemented (or rather used  in code)..

-- To get type stability  (to eliminate that pesky allocation during
calls in a loop), one has to  define:
    femmbase::FEMMBase{S, F} # base finite element modeling machine
Implemented

-- Type function cross(theta, v) and others in the rotation module.
It would seem that  cross() is no longer needed in the rotation module.
Implemented

-- NodalField constructors!  No keywords.
Implemented 05/17/2017

-- ; context...
Replace with  positional arguments.

--  thermalstress{MR<:DeformationModelReduction3D}(::Type{MR},
                        self::MaterialDeformationLinear, dT= 0.0)
    and also update!() should probably use buffer passed in and filled.
    In order to avoid allocations...
Implemented  05/22/2017: removed thermal stress. Stress and strain vectors are
now passed in as buffers.

-- stochastic lifestyle website

-- `update!()` for the material could have  another method  with deformation
gradient going in instead of strain.

--- Right now the tangent moduli are computed  just once  for the
small-deformation elastic model.
`self.material.tangentmoduli!(self.material, D, 0.0, 0.0, loc, 0)`
That is fine when the material is homogeneous within  the  domain,
but it needs to be made more general for heterogeneous structures.

-- Conversion between  matrices and vectors:  must decide which is which
and stick to it.
```
const a = reshape(collect(1:1000),1000,1);
@time vec(a);  # This is MUCH faster
@time dropdims(a, dims = 2);
```

-- Acoustic ABC: surface FEMM?
Implemented 05/30/2017

-- otherdimensionunity(loc::FFltMat, conn::CC, N::FFltMat)::FFlt where {CC}
and elsewhere: CC a union of  vector and column or row matrix?

-- gatherdofnums!(self::Field, dest::A, conn::FIntVec) where {A}
or
gatherdofnums!(self::Field, dest::A, conn::CC) where {A, CC}
Perhaps CC<:AbstractArray
Implemented August 2017

-- Packages:  symbol-gen, compare-files

--- 
```
function inspectintegpoints(self::FEMMDeforLinear,
              geom::NodalField{FFlt},  u::NodalField{T},
              dT::NodalField{FFlt},
              felist::FIntVec,
              inspector::Function,  idat; context...) where {T<:Number}
```
change to
```
function inspectintegpoints(self::FEMMDeforLinear,
              geom::NodalField{FFlt},  u::NodalField{T},
              dT::NodalField{FFlt},
              felist::FIntVec,
              inspector::F,  idat::D; context...) where {T<:Number, F<:Function, D}
```

-- Can I do this for all finite elements in one shot?
```
for i = 1:count(fb.fes)
  getconn!(fb.fes, conn, i);
  gathervalues_asmat!(geom, x, conn);# retrieve element coordinates
  idat = inspectintegpoints(self,  geom,  u,  dT,
            [i],  idi_inspector,  idat; output = output,  context...);
end
```
Implemented 06/03/2017: the inspector function now does not need to refer to
data defined within a loop  over the finite elements, and the inspector can be
called just once for all the finite elements.

-- `inspectintegpoints`:  the inspector function has been redefined.  The
arguments are different.  The description needs to be updated and all references
need to be updated. Implemented: July 2017.

--- Methods that  fill a buffer should also  return the buffer as a convenience.
Apparently, that is the Julia convention.
Implemented  06 17 2017

-- Instead of
```
function test_vectorized(a, a_old)
    @. a[:] = a_old[:] * 2
end
```
use
```
@views function test_vectorized3(a, a_old)
    @. a[:] = a_old[:] * 2
end
```

-- Instead of "restart" use  https://github.com/timholy/Revise.jl

-- Should algorithms print out anything (for instance warnings)?
Option: return an error status as part of the model data.

-- Use static arrays (stack-allocated) in the material update methods.
Status: under deliberation. 10/08/2018

-- rotatestiffness!   Should probably be redesigned for efficiency
On hold.  07/22/2017

-- Atom packages
Julia-client
language-Julia
minimap
auto-indent
split-diff
simple-drag-drop-text

-- Atom themes:
"behave" is an excellent dark theme:  very easy on the eyes

-- Notebooks
```
Pkg.clone("https://github.com/JuliaDiffEq/DiffEqTutorials.jl")
using IJulia
notebook(dir=Pkg.dir("DiffEqTutorials"))
using IJulia;
ENV["IJULIA_DEBUG"] = true
notebook(detached=true)
```


-- When  Atom would not start: run this
  atom --clear-window-state
  in a command window

-- Start notebook:
```
function examples()
  path = joinpath(@__DIR__,"..","examples")
  @eval using IJulia
  @eval notebook(dir=$path)
end
```

-- Mean-strain extrapolation technology: need to consider the thermal properties
of the stabilization material. It needs to produce  thermal strains  just as
the real material. Status: OPEN.

-- Search for JFinEALE: some examples have not been converted  yet. Done.

-- For the QR factorization use qrfact!(A, pivot=Val(false))

-- Return additional data from Abaqus  reader  as a dictionary.
Implemented: 09/09/2017

-- This is a bit of a puzzle: How to transform (rotate) stresses for
  lower-dimensional  finite elements embedded in higher dimensional spaces?
  The transformation matrix is rectangular, not square.

-- in  keymap.cson
'atom-text-editor':
  'ctrl-space': 'autocomplete-plus:activate'

'.platform-win32 .find-and-replace, .platform-linux .find-and-replace':
  'alt-shift-ctrl-r': 'find-and-replace:replace-next'

-- Mesh import: return dictionary of output data
Implemented: 09/09/2017

-- Running Julia in  the git bash terminal (mintty) is an excellent alternative:
the shell mode works.

-- VS code:
Extensions:
Bookmarks, Julia Language Support 0.8, Partial Diff

Setup of the executable:
```
{
    "julia.executablePath": "C:\\Users\\PetrKrysl\\AppData\\Local\\Julia-0.7.0-DEV\\bin\\julia.exe",
}
```

Key bindings:
```
// Place your key bindings in this file to overwrite the defaults
[
{ "key": "ctrl+1", "command": "extension.partialDiff.markSection1",
        "when": "editorTextFocus" },
{ "key": "ctrl+2", "command": "extension.partialDiff.markSection2AndTakeDiff",
        "when": "editorTextFocus" },
{ "key": "ctrl+3", "command": "extension.partialDiff.diffSelectionWithClipboard",
        "when": "editorTextFocus" },
{ "key": "ctrl+shift+i", "command": "editor.action.reindentlines",
        "when": "editorTextFocus" },
]
```
User settings:
```
{
    "julia.executablePath": "C:\\Users\\PetrKrysl\\AppData\\Local\\Julia-0.6.0\\bin\\julia.exe",
    "terminal.integrated.shell.windows": "C:\\Program Files\\Git\\bin\\bash.exe",
    "workbench.colorTheme": "Tomorrow Night Blue",
    "julia.runlinter": false,
    "window.zoomLevel": 0,
    // When opening a file, `editor.tabSize` and `editor.insertSpaces` will be detected based on the file contents.
    "editor.detectIndentation": true,
    // Controls if the editor should automatically adjust the indentation when users type, paste or move lines. Indentation rules of the language must be available.
    "editor.autoIndent": true,
    "editor.wordWrap": "on",
    "editor.wrappingIndent": "indent",
}
```

-- Speed up the elementwise-field transfer the same way it's been done for the nodal fields.
(Using a search structure based on partitioning and sub- meshes.)

-- Atom platformio terminal:  set Auto run command to
/c/Users/PetrKrysl/AppData/Local/Julia-0.7.0-DEV/bin/julia.exe
set Show override to
C:\Program Files\Git\bin\bash.exe

-- Key bindings for Atom

'atom-text-editor':
  'ctrl-space': 'autocomplete-plus:activate'

'.platform-win32 .find-and-replace, .platform-linux .find-and-replace':
  'alt-shift-ctrl-r': 'find-and-replace:replace-next'

'.platform-linux atom-workspace, .platform-win32 atom-workspace':
  'alt-ctrl-shift-enter': 'platformio-ide-terminal:insert-selected-text'

-- Failures to indent code properly in Atom and VS code for lines with
continuation are not specific to  Julia.
Both fail for instance for  C-language code.
Solution: don't bother with splitting lines, rely on auto-wrap.

-- Configuring PyPlot to use Anaconda
```
ENV["PYTHON"]="C:\\ProgramData\\Anaconda3\\python.exe"
ENV["PYTHON"]="C:\\Users\\PetrKrysl\\Documents\\WinPython-64bit-3.6.3.0Qt5\\python-3.6.3.amd64\\python.exe"
Pkg.build("PyCall")
```

-- VS code: To insert a symbol, simply execute the Unicode: Insert Math Symbol command and type in the name of your desired symbol or select it from the drop-down list. (Unicode Latex extension..)

-- SIMD optimizations are not turned on by default, because they only speed up very particular kinds of code, and turning them on everywhere would slow down the compiler too much.
```
function mysum3(A)
    s = zero(eltype(A))
    @simd for a in A
        s += a
    end
    return s
end
```
-- Implement
fen2fe  = FENodeToFEMap(fes, nnodes(geom))
or in other words use the tuple-based connectivity. Status: implemented in November 2017.

-- Do not export ANYTHING? Then  one knows where methods come from because
they need to be explicitly brought in. Implemented: no source files export
functions, only the top-level module exports the public interface.

-- Instead of a geometry field, pass the finite element nodes? Then the node
data could be stored as immutable, with possible performance gains.

-- Have  multiple concrete types for fields? Some mutable, some immutable?

-- When fusing nodes, how about finding the intersection of the bounding
boxes, then finding which nodes fall into those boxes, and then only
checking for proximity for those nodes? Implemented: search speeded up by
    a factor of around 3.

-- Why are there so many allocations in fusenodes()?

-- Number all degrees of freedom (instead of just the free ones). This will
then allow the assembler to assemble the entire matrix, not just
the free-free block. How to supply the number of degrees of freedom to the
assembler when it is getting constructed? Actually, all methods allow
for the assembler to be constructed outside of the method. In that case,
we can construct that assembler not just from the number of free degrees
of freedom, but from the total number of  degrees of freedom.

-- Naming of formal parameters that are modified inside a function:

```
module mmmmmmmmmmmm
function f2!(a!)
  a![1] = 2;
end
end
B=[5]
using .mmmmmmmmmmmm: f2!
f2!(B)
println("B = $(B)")
```b
I like that a! is telling me inside the function that I’m changing something passed in as argument.

-- https://github.com/tonsky/FiraCode/wiki/VS-Code-Instructions
// Controls the font size in pixels.
    "editor.fontSize": 14,
    // Controls the font family.
    "editor.fontFamily": "Fira Code",
    "editor.fontLigatures": true

– Consider MMatrix as a buffer (https://discourse.julialang.org/t/inbounds-is-the-compiler-now-so-smart-that-this-is-no-longer-necessary/12357/17?u=petrkryslucsd)

- uview fast? https://discourse.julialang.org/t/views-performance-in-julia-1-0-dev/12871/2

- Would it make sense to use named tuples instead of data dictionary?

- VS Code JSON settings (user)
```
{
    "terminal.integrated.shell.windows": "C:\\Program Files\\Git\\bin\\bash.exe",
    "terminal.integrated.fontSize": 11,
    "terminal.integrated.fontFamily": "Fira Code",
    "julia.runlinter": true,
    "window.zoomLevel": 1,
    // When opening a file, `editor.tabSize` and `editor.insertSpaces` will be detected based on the file contents.
    "editor.detectIndentation": true,
    // Controls if the editor should automatically adjust the indentation when users type, paste or move lines. Indentation rules of the language must be available.
    "editor.autoIndent": true,
    "editor.wordWrap": "on",
    "editor.wrappingIndent": "indent",
    "workbench.colorTheme": "United Ubuntu",
    "workbench.iconTheme": "vs-minimal",
    "workbench.startupEditor": "newUntitledFile",
    // Controls the font size in pixels.
    "editor.fontSize": 11,
    // Controls the font family.
    // "editor.fontFamily": "Fira Code",
    "editor.fontFamily": "Fira Code",
    "editor.fontLigatures": true,
    "julia.enableTelemetry": true,
    "editor.renderWhitespace": "boundary",
    "window.menuBarVisibility": "visible",
}
```

- Change the name IntegData to IntegDomain. A much better name for the concept of an object which enables integration over some geometrical domain. Status: Implemented 11/08/2018.

-- Rename abstract types for consistency: Change
```
abstract type AbstractSysmatAssembler end;
abstract type AbstractSysvecAssembler end;
abstract type IntegRule end
abstract type DeforModelRed end;
abstract type FEMMAbstractBase; end
abstract type AbstractFEMMDeforLinear <: FEMMAbstractBase end
abstract type AbstractFEMMDeforLinearMS <: AbstractFEMMDeforLinear end
abstract type AbstractFEMMDeforLinearNICE <: AbstractFEMMDeforLinear end
abstract type Field end
abstract type FESet3Manifold{NODESPERELEM} <: FESet{NODESPERELEM} end
abstract type MatDefor end
```
into
```
abstract type AbstractSysmatAssembler end;
abstract type AbstractSysvecAssembler end;
abstract type AbstractIntegRule end
abstract type AbstractDeforModelRed end;

abstract type AbstractFEMMBase; end
abstract type AbstractFEMMDeforLinear <: AbstractFEMMBase end
abstract type AbstractFEMMDeforLinearMS <: AbstractFEMMDeforLinear end
abstract type AbstractFEMMDeforLinearNICE <: AbstractFEMMDeforLinear end
abstract type AbstractFESet3Manifold{NODESPERELEM} <: AbstractFESet{NODESPERELEM} end

abstract type AbstractField end

abstract type AbstractMatDefor end
```
Status: Implemented 04/27/2019.

-- Do we still need functions stored in the structs for materials?
Status: 04/27/2019. It appears we do: the issue is how to make the module functions call the correct functions defined for concrete types without having to define the methods for the types inside the module.

-- Implement a true growing vector for allocating entries in the Assembly module.
This will eliminate the need to know in advance the size, while avoiding the overhead penalty with push!().

-- Color schemes for ST 3: Make STproject.sublime-project project file with
```
{
	"folders":
	[
		{
			"path": "."
		}
	],
	"settings": {
        // "color_scheme": "Monokai.sublime-color-scheme",
                 // "color_scheme": "Packages/Color Scheme - Dusk/Dusk.tmTheme",
                 "color_scheme": "Packages/Color Scheme - JustBeforeDawn/JustBeforeDawn.tmTheme",
                 "color_scheme": "Packages/Coffee Color Scheme/Coffee.tmTheme",
                 "color_scheme": "Packages/Theme - Asphalt/Asphalt.tmTheme",
    },
}
```

-- `updatecsmat!()`, `update!()`: last argument should be the number of the entity, not just the number of an element.
In order to accommodate nodally integrated finite elements.

-- `inspectintegpoints()`: For nodally integrated elements this could be optimized to sweep through the nodes just once, when possible.

-- Double averaging for NICE?

-- divmat!() (and the gradient matrix) do not need to take a buffer. Efficiency is not a concern for these functions. Status: 05-17-2019 done.

-- Package registration:
https://github.com/JuliaRegistries/Registrator.jl
1. Set the Project.toml version field in your repository to your new desired version.
2. comment on branch
`@JuliaRegistrator register`

-- innerproduct() is not being tested. Status:

-- How to implement the "gather operation" efficiently?

```
module mgather1
using Test
using Random
function test()
  N = 1_000_000 # Number of nodes in the mesh
  nen = 20 # Number of nodes per element
  nloops = 2*N
  indexes = randperm(N)[1:nen]
  buffnen3 = rand(nen, 3)
  buff3nen = rand(3, nen)
  dataN3 = rand(N, 3)
  data3N = Matrix(transpose(dataN3))

  t1 = @elapsed for l in 1:nloops
    for i in 1:nen
      ii = indexes[i]
      for j in 1:3
        buffnen3[i, j] = dataN3[ii, j]
      end
    end
  end

  t2 = @elapsed for l in 1:nloops
    for j in 1:3
      for i in 1:nen
        ii = indexes[i]
        buffnen3[i, j] = dataN3[ii, j]
      end
    end
  end

  t3 = @elapsed for l in 1:nloops
    for i in 1:nen
      ii = indexes[i]
      for j in 1:3
        buff3nen[j, i] = dataN3[ii, j]
      end
    end
  end

  t4 = @elapsed for l in 1:nloops
    for j in 1:3
      for i in 1:nen
        ii = indexes[i]
        buff3nen[j, i] = dataN3[ii, j]
      end
    end
  end


  t5 = @elapsed for l in 1:nloops
    for i in 1:nen
      ii = indexes[i]
      for j in 1:3
        buffnen3[i, j] = data3N[j, ii]
      end
    end
  end

  t6 = @elapsed for l in 1:nloops
    for j in 1:3
      for i in 1:nen
        ii = indexes[i]
        buffnen3[i, j] = data3N[j, ii]
      end
    end
  end

  t7 = @elapsed for l in 1:nloops
    for i in 1:nen
      ii = indexes[i]
      for j in 1:3
        buff3nen[j, i] = data3N[j, ii]
      end
    end
  end

  t8 = @elapsed for l in 1:nloops
    for j in 1:3
      for i in 1:nen
        ii = indexes[i]
        buff3nen[j, i] = data3N[j, ii]
      end
    end
  end

  [t1, t2, t3, t4, t5, t6, t7, t8] ./ nloops .* 1e6 # In microseconds
end
end
using Main.mgather1
ts = [0.0 for i in 1:8]
ntries = 100
for i in 1:ntries
  ts .+= mgather1.test()
end
ts ./= ntries
ts = Float32.(ts)

println("Mesh data N x 3, Element buffer nen x 3, Loop i, j: Time $(ts[1]) [mus]")
println("Mesh data N x 3, Element buffer nen x 3, Loop j, i: Time $(ts[2]) [mus]")
println("Mesh data N x 3, Element buffer 3 x nen, Loop i, j: Time $(ts[3]) [mus]")
println("Mesh data N x 3, Element buffer 3 x nen, Loop j, i: Time $(ts[4]) [mus]")
println("Mesh data 3 x N, Element buffer nen x 3, Loop i, j: Time $(ts[5]) [mus]")
println("Mesh data 3 x N, Element buffer nen x 3, Loop j, i: Time $(ts[6]) [mus]")
println("Mesh data 3 x N, Element buffer 3 x nen, Loop i, j: Time $(ts[7]) [mus]")
println("Mesh data 3 x N, Element buffer 3 x nen, Loop j, i: Time $(ts[8]) [mus]")

Results:

Mesh data N x 3, Element buffer nen x 3, Loop i, j: Time 0.09680472 [mus]                                                   
Mesh data N x 3, Element buffer nen x 3, Loop j, i: Time 0.04763569 [mus]          <====== FASTEST                                         
Mesh data N x 3, Element buffer 3 x nen, Loop i, j: Time 0.09928007 [mus]                                                   
Mesh data N x 3, Element buffer 3 x nen, Loop j, i: Time 0.052601 [mus]                                                     
Mesh data 3 x N, Element buffer nen x 3, Loop i, j: Time 0.06879589 [mus]                                                   
Mesh data 3 x N, Element buffer nen x 3, Loop j, i: Time 0.052395895 [mus]                                                  
Mesh data 3 x N, Element buffer 3 x nen, Loop i, j: Time 0.07161312 [mus]                                                   
Mesh data 3 x N, Element buffer 3 x nen, Loop j, i: Time 0.056804158 [mus]         
```

– What about accessing arrays? Here is an example with static arrays and regular arrays:

```
using StaticArrays

struct Point3D{T} <: FieldVector{3, T}
    x::T
    y::T
    z::T
end

mutable struct MPoint3D{T} <: FieldVector{3, T}
    x::T
    y::T
    z::T
end

function f1!(v)
    for i in eachindex(v)
        x, y, z = v[i]
        v[i] = Point3D(2x, x*y, x*z)
    end
end

function f2!(v)
    for i in eachindex(v)
        x, y, z = v[i]
        v[i] .= (2x, x*y, x*z)
    end
end

function f5!(v)
    for i in axes((v, 1)
        x, y, z = v[i, 1], v[i, 2], v[i, 3]
        v[i, 1], v[i, 2], v[i, 3] = 2x, x*y, x*z
    end
end

function f6!(v)
    for i in axes((v, 2)
        x, y, z = v[1, i], v[2, i], v[3, i]
        v[1, i], v[2, i], v[3, i] = 2x, x*y, x*z
    end
end

using BenchmarkTools
v1 = [@SVector(rand(3)) for i = 1:1000000];
v2 = [@MVector(rand(3)) for i = 1:1000000];
v3 = [Point3D(rand(3))  for i = 1:1000000];
v4 = [MPoint3D(rand(3)) for i = 1:1000000];
v5 = rand(1000000, 3);
v6 = rand(3, 1000000);

@btime f1!($v1)   # 1.693 ms (0 allocations: 0 bytes)  
@btime f2!($v2)  # 2.821 ms (0 allocations: 0 bytes)  

@btime f1!($v3)  # 1.728 ms (0 allocations: 0 bytes)     
@btime f2!($v4)  # 57.962 ms (3000000 allocations: 45.78 MiB)   

@btime f5!($v5)  # 1.518 ms (0 allocations: 0 bytes)  <========== FASTEST
@btime f6!($v6)    # 1.739 ms (0 allocations: 0 bytes)
```

– For Loops 2.0: Index Notation And The Future Of Tensor Compilers | Peter Ahrens

-- Why is R .= slower?    
```                             
using BenchmarkTools                                          
R = rand(3,3)    
a = rand(3)                                                         
@btime $R .= $a*$a'                                                                                                                         
  101.931 ns (2 allocations: 176 bytes)                                                                                                            
3×3 Array{Float64,2}:                                                                                                                              
 0.0638104  0.0338781  0.198319                                                                                                                    
 0.0338781  0.0179865  0.105291                                                                                                                    
 0.198319   0.105291   0.616362                                                                                                                    
                                                                                                                                                   
@btime $R = $a*$a'                                                                                                                          
  85.166 ns (2 allocations: 176 bytes)                                                                                                             
3×3 Array{Float64,2}:                                                                                                                              
 0.0638104  0.0338781  0.198319                                                                                                                    
 0.0338781  0.0179865  0.105291                                                                                                                    
 0.198319   0.105291   0.616362     
```

- Registry update. It can be speeded up with
```
using Pkg
pkg"registry rm General"
ENV["JULIA_PKG_SERVER"] = ""
pkg"registry add General"
```
 
-- DOCUMENTER_KEY
```                                    
using Pkg; Pkg.add("DocumenterTools");                                 
using DocumenterTools                                                  
DocumenterTools.genkeys(user="PetrKryslUCSD", repo="git@github.com:PetrKryslUCSD/FinEtools.jl.git")                                                  
using Pkg; Pkg.rm("DocumenterTools");  
```

-- Go-to-error regex

Examples of err messages:

     @ Main.sphere_mode_examples C:\Users\pkonl\Documents\00WIP\FinEtoolsAcoustics.jl\examples\modal\sphere_mode_examples.jl:72 
   [2] allrun()        
   @ Main.sphere_mode_examples C:\Users\pkonl\sphere_mode_examples.jl:72                                                                                                          
     @ Main.sphere_mode_examples C:\Users\pkonl\Documents\00WIP\FinEtoolsAcoustics.jl\examples\modal\sphere_mode_examples.jl:156                                                                                                                                


     julia> include("C:\\Users\\pkonl\\Documents\\00WIP\\FinEtools.jl\\play.jl")
     ERROR: LoadError: UndefVarError: drdrdd not defined
     Stacktrace:                                                                                                                                                                        
        [1] top-level scope                                                                                                                                                             
          @ C:\Users\pkonl\Documents\00WIP\FinEtools.jl\play.jl:1                                                                                                                       
        [2] include(fname::String)                                                                                                                                                      
          @ Base.MainInclude .\client.jl:476                                                                                                                                            
     in expression starting at C:\Users\pkonl\Documents\00WIP\FinEtools.jl\play.jl:1                                                                                                    
                                                                                                                                                                                        
     julia> 
     
{ // This build system simply opens a new interactive Julia REPL 
    // 2023 Petr Krysl
    // Depends on AbbreviatedStackTraces
    "title": "Julia REPL",

    "target": "terminus_open",

    "auto_close": false,
    "shell_cmd": "julia -t4", // Assume that PATH is set before invoking the editor
    "cwd": "${file_path:${folder}}",
    // "selector": "source.julia",
    "focus": false,
   
    "file_regex": "(?:[@](?:\\s\\S[\\w.]*[^@])??|in expression starting at)\\s(\\S+[.]jl):([0-9]+)\\s.*"
}

"file_regex": "(?:[@][\s][\s\S.*]*|in expression starting at)\s(\S.+[.jl]):([0-9]+)\s*"  




"""
    SysmatAssemblerSparseThr{T<:Number} <: AbstractSysmatAssembler

Type for multithreaded assembling of a sparse global matrix from elementwise matrices.

!!! note

    All fields of the datatype are private. No need to access them directly.

# Example

```

nth = Base.Threads.nthreads()

function ass(a, assembly_line, r)
    for i in r
        AM.assemble!(a, assembly_line[i]...)
    end
    return a
end

AM.startassembly!(a, 5, 5, length(assembly_line), N, N)
Base.Threads.@threads :static for th in 1:Base.Threads.nthreads()
    a, r = AM.startassembly!(a, 5, 5, length(assembly_line), N, N)
    a = ass(a, assembly_line, r)
    A = makematrix!(a)
end
A = makematrix!(a)
```
The assembler `a` contains a thread-private storage. Notice the `:static`
argument of the macro `@threads`: this will guarantee that each call of
`AM.assemble!(a, assembly_line[i]...)` will be executed by a fixed thread
(i.e. calls to `threadid()` will yield the same value within the `assemble!`
function).
"""
mutable struct SysmatAssemblerSparseThr{T<:Number} <: AbstractSysmatAssembler
    matbuffer::Vector{T}
    rowbuffer::Vector{FInt}
    colbuffer::Vector{FInt}
    buffer_length::FInt
    thread_done::Vector{Bool}
    thread_begin::Vector{FInt}
    thread_end::Vector{FInt}
    thread_pointer::Vector{FInt}
    ndofs_row::FInt
    ndofs_col::FInt
    nomatrixresult::Bool
end

"""
    SysmatAssemblerSparseThr(zero::T=0.0, nomatrixresult = false) where {T<:Number}

Construct blank system matrix assembler.

This assembler is like `SysmatAssemblerSparse`, except it can be used with
multiple threads.
"""
function SysmatAssemblerSparseThr(zero::T = 0.0, nomatrixresult = false) where {T<:Number}
    nth = Base.Threads.nthreads()
    return SysmatAssemblerSparseThr{T}(
        Vector{T}[],
        Vector{FInt}[],
        Vector{FInt}[],
        0,
        [false for _  in 1:nth],
        [0 for _  in 1:nth],
        [0 for _  in 1:nth],
        [0 for _  in 1:nth],
        0,
        0,
        nomatrixresult,
        )
end

"""
    startassembly!(self::SysmatAssemblerSparseThr{T},
      elem_mat_nrows::FInt, elem_mat_ncols::FInt, elem_mat_nmatrices::FInt,
      ndofs_row::FInt, ndofs_col::FInt) where {T<:Number}

Start the assembly of a global matrix.

This assembler is like `SysmatAssemblerSparse`, except it can be used with
multiple threads.

The first call is expected to occur *before* the threaded loop.
Therefore, the data on input indicates the total load.

# Returns
- `self`: the modified assembler.
"""
function startassembly!(
    self::SysmatAssemblerSparseThr{T},
    elem_mat_nrows::FInt,
    elem_mat_ncols::FInt,
    elem_mat_nmatrices::FInt,
    ndofs_row::FInt,
    ndofs_col::FInt,
) where {T<:Number}
    nth = Base.Threads.nthreads()
    # The first call is carried out by thread 1, and the buffer is empty at that
    # point.
    th = Base.Threads.threadid()
    elem_per_thread = Int(ceil(elem_mat_nmatrices / nth))
    first_call = (th == 1) && (self.buffer_length == 0)
    if first_call # the assembly buffer has not been initialized yet
        # Increase the allocated length a little bit to allow for irregularities
        self.buffer_length = Int(round(1.1 * elem_mat_nmatrices * elem_mat_nrows * elem_mat_ncols))
        # Initialize the rows and columns to 1, and the coefficients to 0.0.
        # This will make it possible to assemble even entries not set by
        # the threads.
        self.rowbuffer = fill(1, self.buffer_length)
        self.colbuffer = fill(1, self.buffer_length)
        self.matbuffer = fill(zero(T), self.buffer_length)
        self.ndofs_row = ndofs_row
        self.ndofs_col = ndofs_col
    else
        # Now in this case the function is called from each individual thread.
        chunk = Int64(round(self.buffer_length / nth)) + 1 # Number of elements per thread
        self.thread_begin[th] = chunk*(th-1)+1
        self.thread_end[th] = chunk*(th)+1-1
        @assert chunk >= elem_per_thread * elem_mat_nrows * elem_mat_ncols
        if self.thread_end[th] > self.buffer_length
            self.thread_end[th] = self.buffer_length
        end
        self.thread_pointer[th] = self.thread_begin[th]
    end

    return self
end

"""
    assemble!(self::SysmatAssemblerSparseThr{T}, mat::FMat{T}, dofnums_row::Union{FIntVec, FIntMat}, dofnums_col::Union{FIntVec, FIntMat}) where {T<:Number}

Assemble a rectangular matrix.

This assembler is like `SysmatAssemblerSparse`, except it can be used with
multiple threads.
"""
function assemble!(
    self::SysmatAssemblerSparseThr{T},
    mat::FMat{T},
    dofnums_row::Union{FIntVec,FIntMat},
    dofnums_col::Union{FIntVec,FIntMat},
) where {T<:Number}
    th = Base.Threads.threadid()
    nrows = length(dofnums_row)
    ncolumns = length(dofnums_col)
    p = self.thread_pointer[th]
    @assert p + ncolumns * nrows <= self.thread_end[th] + 1
    @assert size(mat) == (nrows, ncolumns)
    @inbounds for j in 1:ncolumns,  i in 1:nrows
        self.matbuffer[p] = mat[i, j] # serialized matrix
        self.rowbuffer[p] = dofnums_row[i]
        self.colbuffer[p] = dofnums_col[j]
        p = p + 1
    end
    self.thread_pointer[th] = p
    return self
end

"""
    makematrix!(self::SysmatAssemblerSparseThr)

Make a sparse matrix.

This assembler is like `SysmatAssemblerSparse`, except it can be used with
multiple threads.
"""
# function makematrix!(self::SysmatAssemblerSparseThr{T}) where {T<:Number}
#     # Here we will go through the rows and columns, and whenever the row or
#     # the column refer to indexes outside of the limits of the matrix, the
#     # corresponding value will be set to 0 and assembled to row and column 1.

#     th = Base.Threads.threadid()
#     all_done = true
#     for k in eachindex(self.thread_done)
#         all_done = all_done && (
#             (self.thread_done[k] || (
#                 isempty(self.rowbuffer[k]) &&
#                 isempty(self.colbuffer[k]) &&
#                 isempty(self.matbuffer[k]))
#             )
#             )
#     end
#     if self.nomatrixresult || (!all_done)
#         # No actual sparse matrix is returned. The entire result of the assembly
#         # is preserved in the assembler buffers.
#         # @inbounds for j in 1:length(self.rowbuffer[th])
#         #     if (self.rowbuffer[th][j] > self.ndofs_row) || (self.rowbuffer[th][j] <= 0)
#         #         self.rowbuffer[th][j] = 1
#         #         self.matbuffer[th][j] = 0.0
#         #     end
#         #     if (self.colbuffer[th][j] > self.ndofs_col) || (self.colbuffer[th][j] <= 0)
#         #         self.colbuffer[th][j] = 1
#         #         self.matbuffer[th][j] = 0.0
#         #     end
#         # end
#         self.thread_done[th] = true
#         return spzeros(self.ndofs_row, self.ndofs_col)
#     elseif all_done && th == 1
#         # The sparse matrix is constructed and returned. The  buffers used for
#         # the assembly are cleared.
#         S = sparse(
#             cat(self.rowbuffer..., dims=1),
#             cat(self.colbuffer..., dims=1),
#             cat(self.matbuffer..., dims=1),
#             self.ndofs_row,
#             self.ndofs_col,
#         )
#         self = SysmatAssemblerSparseThr(zero(T)) # get rid of the buffers
#         return S
#     end
# end

function makematrix!(self::SysmatAssemblerSparseThr{T}) where {T<:Number}
    th = Base.Threads.threadid()
    all_done = !any(x -> !x, self.thread_done)
    self.thread_done[th] = true

    # Here we will go through the rows and columns, and whenever the row or
    # the column refer to indexes outside of the limits of the matrix, the
    # corresponding value will be set to 0 and assembled to row and column 1.
    if all_done && th == 1
        #@inbounds
        for j in eachindex(self.rowbuffer)
            if (self.rowbuffer[j] > self.ndofs_row) || (self.rowbuffer[j] <= 0)
                self.rowbuffer[j] = 1
                self.matbuffer[j] = 0.0
            end
            if (self.colbuffer[j] > self.ndofs_col) || (self.colbuffer[j] <= 0)
                self.colbuffer[j] = 1
                self.matbuffer[j] = 0.0
            end
        end
        @assert th == 1 # This is supposed to happen only in serial mode
        if self.nomatrixresult
            # No actual sparse matrix is returned. The entire result of the assembly
            # is preserved in the assembler buffers.
            return spzeros(self.ndofs_row, self.ndofs_col)
        end
        S = sparse(
            self.rowbuffer,
            self.colbuffer,
            self.matbuffer,
            self.ndofs_row,
            self.ndofs_col,
            )
        # The buffers are cleared
        self.rowbuffer = []
        self.colbuffer = []
        self.matbuffer = []
        self.buffer_length = 0

        return S
    else
        return spzeros(self.ndofs_row, self.ndofs_col)
    end
end


module mutil_th001

using Test
using LinearAlgebra
using FinEtools
import FinEtools.AssemblyModule as SM

function _test()
    # @show Base.Threads.nthreads()
    m = [
    0.24406   0.599773    0.833404  0.0420141
    0.786024  0.00206713  0.995379  0.780298
    0.845816  0.198459    0.355149  0.224996
    ]
    m1 = m'*m; r1 = [5 2 1 4]; c1 = r1
    m = [
    0.146618  0.53471   0.614342    0.737833
    0.479719  0.41354   0.00760941  0.836455
    0.254868  0.476189  0.460794    0.00919633
    0.159064  0.261821  0.317078    0.77646
    0.643538  0.429817  0.59788     0.958909
    ]
    m2 = m'*m; r2 = [2 3 1 5]; c2 = r2
    m = [
    -0.146618  0.53471   0.614342    0.737833
    0.479719  0.41354   0.00760941  -0.836455
    0.254868  -0.476189  0.460794    0.00919633
    0.159064  0.261821  0.317078    0.77646
    ]
    m3 = m; r3 = [2 3 1 5]; c3 = [7 6 2 4]
    m4 = Matrix(m'); r4 = [7 6 2 4]; c4 = [2 3 1 5];
    assembly_line = [
    (m1, r1, c1),
    (m1, r1, c1),
    (m2, r2, c2),
    (m2, r2, c2),
    (m3, r3, c3),
    (m4, r4, c4),
    (m1, r1, c1),
    (m1, r1, c1),
    (m2, r2, c2),
    (m2, r2, c2),
    (m4, r4, c4),
    (m3, r3, c3),
    (m3, r3, c3),
    (m1, r1, c1),
    (m1, r1, c1),
    (m2, r2, c2),
    (m2, r2, c2),
    (m3, r3, c3),
    (m4, r4, c4),
    (m1, r1, c1),
    (m1, r1, c1),
    (m2, r2, c2),
    (m2, r2, c2),
    (m4, r4, c4),
    (m3, r3, c3),
    (m3, r3, c3),
    (m1, r1, c1),
    (m1, r1, c1),
    (m2, r2, c2),
    (m2, r2, c2),
    (m3, r3, c3),
    (m4, r4, c4),
    (m1, r1, c1),
    (m1, r1, c1),
    (m2, r2, c2),
    (m2, r2, c2),
    (m4, r4, c4),
    (m3, r3, c3),
    (m3, r3, c3),
    (m1, r1, c1),
    (m1, r1, c1),
    (m2, r2, c2),
    (m2, r2, c2),
    (m3, r3, c3),
    (m4, r4, c4),
    (m1, r1, c1),
    (m1, r1, c1),
    (m2, r2, c2),
    (m2, r2, c2),
    (m4, r4, c4),
    (m3, r3, c3),
    (m3, r3, c3),
    ]
    a = SM.SysmatAssemblerSparse(0.0)
    SM.startassembly!(a, 5, 5, length(assembly_line), 7, 7)
    for i in eachindex(assembly_line)
        SM.assemble!(a, assembly_line[i]...)
    end
    refA = SM.makematrix!(a)

    # @show a

    a = SM.SysmatAssemblerSparseThr(0.0)

    nth = Base.Threads.nthreads()

    SM.startassembly!(a, 5, 5, length(assembly_line), 7, 7)
    Base.Threads.@threads :static for th in 1:Base.Threads.nthreads()
        a, r = SM.startassembly!(a, 5, 5, length(assembly_line), 7, 7)
        for i in r
            SM.assemble!(a, assembly_line[i]...)
        end
        A = makematrix!(a)
    end
    # @info "Final"
    A = makematrix!(a)
    @test norm(A - refA) / norm(refA) < 1.0e-9
    return true
end

_test()

end # module


module mutil_th002

using Test
using LinearAlgebra
using FinEtools
import FinEtools.AssemblyModule as AM
using Random


function _test()
    Random.seed!(1234);
    # @show Base.Threads.nthreads()
    m = [
    0.24406   0.599773    0.833404  0.0420141
    0.786024  0.00206713  0.995379  0.780298
    0.845816  0.198459    0.355149  0.224996
    ]
    m1 = m'*m;
    m = [
    0.146618  -0.53471   0.614342    0.737833
    0.479719  0.41354   -0.00760941  0.836455
    0.254868  0.476189  0.460794    0.00919633
    0.159064  0.261821  0.317078    0.77646
    -0.643538  0.429817  0.59788     0.958909
    ]
    m2 = m'*m;
    m = [
    -0.146618  0.53471   0.614342    0.737833
    0.479719  0.41354   0.00760941  -0.836455
    0.254868  -0.476189  0.460794    0.00919633
    0.159064  0.261821  0.317078    0.77646
    ]
    m3 = m;
    m4 = Matrix(m');
    ms = [m1, m2, m3, m4, ]
    N = 16000
    p = randperm(N)
    assembly_line = []
    for i in 1:N
        for k in eachindex(ms)
            if length(p) < sum(size(m))
                p = randperm(N)
            end
            m = ms[k]
            r = [popfirst!(p) for _ in 1:size(m, 1)]
            c = [popfirst!(p) for _ in 1:size(m, 2)]
            push!(assembly_line, (m, r, c))
        end
    end

    start = time()
    a = AM.SysmatAssemblerSparse(0.0)
    AM.startassembly!(a, 5, 5, length(assembly_line), N, N)
    for i in eachindex(assembly_line)
        AM.assemble!(a, assembly_line[i]...)
    end
    refA = AM.makematrix!(a)
    # @show time() - start

    # display(spy(refA))
    # @show a

    a = AM.SysmatAssemblerSparseThr(0.0)

    nth = Base.Threads.nthreads()

    function ass(a, assembly_line, r)
        for i in r
            AM.assemble!(a, assembly_line[i]...)
        end
        return a
    end

    start = time()
    AM.startassembly!(a, 5, 5, length(assembly_line), N, N)
    Base.Threads.@threads :static for th in 1:Base.Threads.nthreads()
        a, r = AM.startassembly!(a, 5, 5, length(assembly_line), N, N)
        a = ass(a, assembly_line, r)
        A = makematrix!(a)
    end
    # @info "Final"
    A = makematrix!(a)
    # @show time() - start

    @test norm(A - refA) / norm(refA) < 1.0e-9
    return true
end

_test()

end # module

-- Updating of coordinate systems. updatecsmat!() perhaps should receive the id
   of the finite element instead of its label?
   Done. Also the quadrature point id.
   
-- The delegation mechanism needs to be described. 09/29/23

-- The finite element model machine documentation needs to be refreshed. 09/29/23

-- Describe bilinear form and linear form functions.

-- JuliaFormatter: annotate_untyped_fields_with_any=true


-- Formatting: JuliaFormatter
using JuliaFormatter
format("./src", SciMLStyle(), annotate_untyped_fields_with_any=true)   

-- integratefieldfunction: should there be a variant which takes something like
   the gradients of the basis functions and similar quantities?
   

-- Assembler for Reduced Models. Not used at the moment.

```

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

function SysmatAssemblerReduced(t::TMT,
    z = zero(Float64),
    nomatrixresult = false) where {TMT}
    ndofs_row = ndofs_col = size(t, 1)
    red_ndofs_row = red_ndofs_col = size(t, 2)
    m = fill(z, red_ndofs_row, red_ndofs_col)
    return SysmatAssemblerReduced(m,
        ndofs_row,
        ndofs_col,
        red_ndofs_row,
        red_ndofs_col,
        t,
        nomatrixresult)
end

function startassembly!(self::SysmatAssemblerReduced,
    elem_mat_nrows,
    elem_mat_ncols,
    elem_mat_nmatrices,
    ndofs_row,
    ndofs_col)
    @assert self.ndofs_row == ndofs_row
    @assert self.ndofs_col == ndofs_col
    self.m .= zero(eltype(self.m))
    return self
end

function assemble!(self::SysmatAssemblerReduced,
    mat::MT,
    dofnums_row::IV,
    dofnums_col::IV) where {MT, IV}
    R, C = size(mat)
    T = eltype(mat)
    @assert R==C "The matrix must be square"
    @assert dofnums_row==dofnums_col "The degree of freedom numbers must be the same for rows and columns"
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
```

-- Why is DeforModelRedModule in FinEtools? Should it (could it) be moved to the linear deformation package? 
