using PyCall

function readSU2(fileName::AbstractString)
    py"""
    import meshio
    """

    m = py"meshio.read"(fileName)
    nodes      = m.points
    triangles  = m.get_cells_type("triangle") .+ 1
    edgeNodes  = m.get_cells_type("line") .+ 1
    edgeBCs    = m.get_cell_data("su2:tag","line")
    edges = hcat(edgeNodes,edgeBCs)

    return Mesh2D(nodes,edges,triangles)
end

