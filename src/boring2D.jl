module boring2D

export Mesh2D, readMesh, writeMesh

struct Mesh2D
    " [nNodes,2 Float64] Node locations (X,Y)."
    nodes::Matrix{Float64}

    "[nNodes,3 Int32] Edge nodes numbers and boundary number."
    edges::Matrix{Int32}

    "[nSurfTrias,3 Int32] Node indices for triangular boundary surface faces."
    triangles::Matrix{Int32}

    # "Internal constructor"
    # Gmf() = new()
end

include("meshIO.jl")

end