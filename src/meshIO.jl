using PyCall

# Read mesh with format infered from file extension (most of the time)
function readMesh(meshFileName::AbstractString)
    extension = getFileExtension(meshFileName)
    return readMeshFormat(meshFileName::AbstractString,extension::AbstractString)
end

# Read mesh with specified format (less often)
function readMesh(meshFileName::AbstractString,extension::AbstractString)
    return readMeshFormat(meshFileName::AbstractString,extension::AbstractString)
end

function readMeshFormat(meshFileName::AbstractString,format::AbstractString)
    # Check that mesh file exists
    if !isfile(meshFileName)
        println("""Mesh file "$meshFileName" not found.""")
        return 1
    end

    # Read format
    if format=="su2"
        return readSU2(meshFileName::AbstractString)
    else
        println("""Mesh format "$format" not supported for reading.""")
    end
    
end

function getFileExtension(fileName::AbstractString)
    splitFileName = split(fileName,".")
    extension = splitFileName[end]
    return lowercase(extension)
end

function readSU2(meshFileName::AbstractString)
    py"""
    import meshio
    """

    m = py"meshio.read"(meshFileName)
    nodes      = m.points
    triangles  = m.get_cells_type("triangle") .+ 1
    edgeNodes  = m.get_cells_type("line") .+ 1
    edgeBCs    = m.get_cell_data("su2:tag","line")
    edges = hcat(edgeNodes,edgeBCs)

    return Mesh2D(nodes,edges,triangles)
end

