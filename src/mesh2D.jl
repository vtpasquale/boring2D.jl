using PyCall


struct Mesh2D
    "[nNodes,2 Float64] Node locations (x,y)."
    nodes::Matrix{Float64}

    "[nNodes,3 Int32] Edge nodes numbers and boundary ID number."
    edges::Matrix{Int32}

    "[nTriangles,3 Int32] Triangle element node numbers."
    triangles::Matrix{Int32}
    
end

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
    elseif format=="plt"
        return readPLT(meshFileName::AbstractString)
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


function readPLT(filename::String)
    # Construct from CBSFlow ascii file
    fid = open(filename,"r")
                
    # Process header
    getLine = readline(fid)
    splitLine = split(strip(getLine))  
    nTriangles = parse(Int32,splitLine[1])
    nNodes = parse(Int32,splitLine[2])
    nEdges = parse(Int32,splitLine[3])

    # Read triangles
    tri = Array{Int32}(undef,nTriangles,3)
    for i = 1:nTriangles
        splitLine = split(strip(readline(fid)))
        tri[i,1] = parse(Int32,splitLine[2])
        tri[i,2] = parse(Int32,splitLine[3])
        tri[i,3] = parse(Int32,splitLine[4])
    end
    
    # Read nodes
    nodes = Array{Float64}(undef,nNodes,2)
    for i = 1:nNodes
        splitLine = split(strip(readline(fid)))
        nodes[i,1] = parse(Float64,splitLine[2])
        nodes[i,2] = parse(Float64,splitLine[3])
    end

    # Read edges
    edges = Array{Int32}(undef,nEdges,3)
    for i = 1:nEdges
        splitLine = split(strip(readline(fid)))
        edges[i,1] = parse(Int32,splitLine[1])
        edges[i,2] = parse(Int32,splitLine[2])
        edges[i,3] = parse(Int32,splitLine[4])
    end
    
    close(fid)
    return Mesh2D(nodes,edges,tri)
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
    elseif format=="nas"
        writeMeshio(meshFileName::AbstractString,mesh::Mesh2D)
        return 0
    else
        println("""Mesh format "$format" not supported for writing.""")
        return 0
    end
    
end

function writeMeshio(meshFileName::AbstractString,mesh::Mesh2D)
    py"""
    import numpy as np
    import meshio
    def pyWriteMeshIo(meshFileName,mesh):
        nTriangles = np.size(mesh.triangles,0)
        points = mesh.nodes
        cells = [("triangle",mesh.triangles-1),("line",mesh.edges[:,0:2]-1)]
        cell_data_out = {"su2:tag": [np.zeros(nTriangles,dtype=int),mesh.edges[:,2]]}
        m = meshio.Mesh(points,cells,cell_data=cell_data_out)
        m.write(meshFileName)
        return 0
    """
    py"pyWriteMeshIo"(meshFileName,mesh)
    return 0
end

function writeSolution(solutionFileName::AbstractString,mesh::Mesh2D,pointOutput::Dict,cellOutput::Dict)
    # Write solution to output file with meshio. VTU is the suggested output format. pointOutput and cellOutput can be empty dictionaries.
    py"""
    import numpy as np
    import meshio
    def pyWriteSolution(meshFileName,mesh,pointOutputDict,cellOutputDict):
        nTriangles = np.size(mesh.triangles,0)
        points = mesh.nodes
        cells = [("triangle",mesh.triangles-1)]
        m = meshio.Mesh(points,cells,point_data=pointOutputDict,cell_data=cellOutputDict)
        m.write(meshFileName)
        return 0
    """
    py"pyWriteSolution"(solutionFileName,mesh,pointOutput,cellOutput)
    return 0
end

# --------------- HELPER --------------------------


function getFileExtension(fileName::AbstractString)
    splitFileName = split(fileName,".")
    extension = splitFileName[end]
    return lowercase(extension)
end