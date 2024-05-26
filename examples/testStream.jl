
# using Revise
using boring2D

cd("examples")

# Cl = solveStream("n0012.toml")
# println("Cl = $(Cl)")

Cl,dCl = solveStreamLiftAdjoint("n0012.toml")

println("Cl = $(Cl)")
println("dCl = $(dCl)")






