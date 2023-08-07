using PyCall

# --------------- READ --------------------------


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
        return 1
    end
    
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

# --------------- WRITE --------------------------

# write mesh with format infered from file extension (most of the time)
function writeMesh(meshFileName::AbstractString,mesh::Mesh2D)
    extension = getFileExtension(meshFileName)
    return writeMeshFormat(meshFileName::AbstractString,extension::AbstractString,mesh::Mesh2D)
end

# write mesh with specified format (less often)
function writeMesh(meshFileName::AbstractString,extension::AbstractString,mesh::Mesh2D)
    return writeMeshFormat(meshFileName::AbstractString,extension::AbstractString,mesh::Mesh2D)
end

function writeMeshFormat(meshFileName::AbstractString,format::AbstractString,mesh::Mesh2D)

    # write format
    if format=="su2"
        writeMeshio(meshFileName::AbstractString,mesh::Mesh2D)
        return 0
    elseif format=="vtu"
        writeMeshio(meshFileName::AbstractString,mesh::Mesh2D)
        return 0
    elseif format=="vtk"
        writeMeshio(meshFileName::AbstractString,mesh::Mesh2D)
        return 0
    elseif format=="ugrid"
        writeMeshio(meshFileName::AbstractString,mesh::Mesh2D)
        return 0
    elseif format=="nas"
        writeMeshio(meshFileName::AbstractString,mesh::Mesh2D)
        return 0
    else
        println("""Mesh format "$format" not supported for writing.""")
        return 0
    end
    
end

function writeMeshio(meshFileName::AbstractString,mesh::Mesh2D)
    points = mesh.nodes
    cells = [("triangle",mesh.triangles.-1)]

    py"""
    import meshio
    """
    m = py"meshio.Mesh"(points,cells)
    m.write(meshFileName)

    return 0
end

# --------------- HELPER --------------------------


function getFileExtension(fileName::AbstractString)
    splitFileName = split(fileName,".")
    extension = splitFileName[end]
    return lowercase(extension)
end