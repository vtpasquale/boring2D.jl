module boring2D

export Mesh2D, readMesh, writeMesh, writeSolution, solveConvection, TriangleElements, assembleConvectionStiffness, solvePotentialFlow, solveStream

struct Mesh2D
    "[nNodes,2 Float64] Node locations (x,y)."
    nodes::Matrix{Float64}

    "[nNodes,3 Int32] Edge nodes numbers and boundary ID number."
    edges::Matrix{Int32}

    "[nTriangles,3 Int32] Triangle element node numbers."
    triangles::Matrix{Int32}

    # "Internal constructor"
    # Gmf() = new()
end

include("meshIO.jl")
include("triangleElements.jl")
include("sparseTriplet.jl")
include("solveConvection.jl")
include("solvePotentialFlow.jl")
include("solveStream.jl")

end