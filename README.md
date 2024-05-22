# boring2D

[![Build Status](https://github.com/vtpasquale/boring2D.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/vtpasquale/boring2D.jl/actions/workflows/CI.yml?query=branch%3Amain)



## Set Python Distribution
To force Julia to use its own Python distribution, via Conda, simply set `ENV["PYTHON"]` to the empty string `""` and re-run `Pkg.build("PyCall")`.
```
using Pkg
using PyCall
ENV["PYTHON"]=""
Pkg.build("PyCall")
```
