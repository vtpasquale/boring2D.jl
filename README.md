# boring2D

[![Build Status](https://github.com/vtpasquale/boring2D.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/vtpasquale/boring2D.jl/actions/workflows/CI.yml?query=branch%3Amain)



## Setup Julia packages and Python
```
# In julia terminal

# Add PyCall to base with package manager
julia>] 
(@v1.10) pkg> add PyCall

# Force Julia to use its own Python distribution
julia> ENV["PYTHON"]="" 

# Rebuild PyCall
julia> using Pkg
julia> using PyCall
julia> Pkg.build("PyCall")

# Add Conda to base with package manager
julia>] 
(@v1.10) pkg> add Conda

# Add meshio Python libray
julia> using Conda
julia> Conda.add("meshio")

] # Dev boring2D from package manager
julia>] 
(@v1.10) pkg> dev https://github.com/vtpasquale/boring2D.jl.git


# activate base package when needed
(boring2D) pkg> activate 
(@v1.10) pkg>

# activate boring2D package when needed
(@v1.10) pkg> activate boring2D
(boring2D) pkg> 

```
