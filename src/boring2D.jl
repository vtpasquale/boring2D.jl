module boring2D

export Mesh2D, readMesh, writeMesh

struct Mesh2D
    " [nNodes,2 Float64] Node locations (x,y)."
    nodes::Matrix{Float64}

    "[nNodes,3 Int32] Edge nodes numbers and boundary ID number."
    edges::Matrix{Int32}

    "[nTriangles,3 Int32] Triangle element node numbers."
    triangles::Matrix{Int32}

    # "Internal constructor"
    # Gmf() = new()
end

include("meshIO.jl")

end