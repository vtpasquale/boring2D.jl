
# using Revise
using boring2D

cd("examples")

Cl, ∂Cl∂α = solveStream("n0012.toml")
println("Cl = $(Cl)")
println("Cl*180/pi = $(Cl*180/pi)")
println("∂Cl∂α     = $(∂Cl∂α)")