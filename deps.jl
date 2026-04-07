using Pkg

Pkg.add(Pkg.PackageSpec(url="https://github.com/Hydro-Umass/Sad.jl"))
Pkg.add("ArgParse")
Pkg.add("PackageCompiler")
Pkg.add("Distributions")
Pkg.add("NCDatasets")
Pkg.add("JSON")
Pkg.add("PyCall")
Pkg.precompile()
