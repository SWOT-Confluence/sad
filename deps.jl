using Pkg

# Pkg.add(Pkg.PackageSpec(url="https://gitlab.com/kandread/Sad.jl.git"))
Pkg.add(Pkg.PackageSpec(url="https://gitlab.com/nikki-t/Sad.jl.git"))
Pkg.add("PackageCompiler")
Pkg.add("Distributions")
Pkg.add("NCDatasets")
Pkg.add("JSON")
Pkg.precompile()