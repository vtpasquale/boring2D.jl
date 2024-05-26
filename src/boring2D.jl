module boring2D

export solveConvection, solveStream, solveStreamLiftAdjoint

struct Mesh2D
    "[nNodes,2 Float64] Node locations (x,y)."
    nodes::Matrix{Float64}

    "[nNodes,3 Int32] Edge nodes numbers and boundary ID number."
    edges::Matrix{Int32}

    "[nTriangles,3 Int32] Triangle element node numbers."
    triangles::Matrix{Int32}
    
end

include("meshIO.jl")
include("triangleElements.jl")
include("sparseTriplet.jl")
include("solveConvection.jl")
include("solveStream.jl")
include("solveStreamLiftAdjoint.jl")

end