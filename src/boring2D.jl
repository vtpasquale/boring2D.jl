module boring2D

export solveConvection, solveStream

include("mesh2D.jl")
include("triangleElements.jl")
include("closedBoundary2D.jl")
include("sparseTriplet.jl")
include("solveConvection.jl")
include("solveStream.jl")

end