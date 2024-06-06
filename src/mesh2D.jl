using PyCall


struct Mesh2D
    "[nNodes,2 Float64] Node locations (x,y)."
    nodes::AbstractArray # Matrix{Float64}

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

function getNextSu2Line(fid)
    getLine = readline(fid)
    splitLine = split(strip(getLine))
    if startswith(splitLine[1],"%")
        nextLine = getNextSu2Line(fid)
    else
        return splitLine
    end
end

function readSU2(meshFileName::AbstractString)

    # meshio is not reading meshes from pointwise correctly
    # py"""
    # import meshio
    # """
    # m = py"meshio.read"(meshFileName)
    # nodes      = m.points
    # triangles  = m.get_cells_type("triangle") .+ 1
    # edgeNodes  = m.get_cells_type("line") .+ 1
    # edgeBCs    = m.get_cell_data("su2:tag","line")
    # edges = hcat(edgeNodes,edgeBCs)

    fid = open(meshFileName,"r")

    # Process header
    nextLine = getNextSu2Line(fid)
    if cmp(uppercase(nextLine[1]),"NDIME=") != 0
        error("Error processing NDIME in .su2 file.")
    end
    ndime = parse(Int32,nextLine[2])
    if ndime != 2
        error("NDIME should be 2 in .su2 file.")
    end
    
    # Read triangles
    nextLine = getNextSu2Line(fid)
    if cmp(uppercase(nextLine[1]),"NELEM=") != 0
        error("Error processing NELEM in .su2 file.")
    end
    nTriangles = parse(Int32,nextLine[2])
    tri = Array{Int32}(undef,nTriangles,3)
    for i = 1:nTriangles
        nextLine = getNextSu2Line(fid)
        if parse(Int32,nextLine[1]) != 5
            error("Only triangle elements supported.")
        end
        tri[i,1] = parse(Int32,nextLine[2])
        tri[i,2] = parse(Int32,nextLine[3])
        tri[i,3] = parse(Int32,nextLine[4])
    end
    tri .+= 1 # convert base 0 to base 1
    
    # Read nodes
    nextLine = getNextSu2Line(fid)
    if cmp(uppercase(nextLine[1]),"NPOIN=") != 0
        error("Error processing NPOIN in .su2 file.")
    end
    nNodes = parse(Int32,nextLine[2])
    nodes = Array{Float64}(undef,nNodes,2)
    for i = 1:nNodes
        nextLine = getNextSu2Line(fid)
        nodes[i,1] = parse(Float64,nextLine[1])
        nodes[i,2] = parse(Float64,nextLine[2])
    end
    
    # Read edges
    Int321 = Int32(1)
    edges = []
    nextLine = getNextSu2Line(fid)
    if cmp(uppercase(nextLine[1]),"NMARK=") != 0
        error("Error processing NMARK in .su2 file.")
    end
    nBoundary = parse(Int32,nextLine[2])
    for i = 1:nBoundary
        nextLine = getNextSu2Line(fid)
        if cmp(uppercase(nextLine[1]),"MARKER_TAG=") != 0
            error("Error processing MARKER_TAG in .su2 file.")
        end
        nextLine = getNextSu2Line(fid)
        if cmp(uppercase(nextLine[1]),"MARKER_ELEMS=") != 0
            error("Error processing MARKER_ELEMS in .su2 file.")
        end
        nEdge = parse(Int32,nextLine[2])
        for j = 1:nEdge
            nextLine = getNextSu2Line(fid)
            push!(edges,[parse(Int32,nextLine[2])+Int321 parse(Int32,nextLine[3])+Int321 i])
        end
    end
    nEdges_ = size(edges,1)
    edges_ = Array{Int32}(undef,nEdges_,3)
    for i = 1:nEdges_
        edges_[i,:] = edges[i]
    end
    

    return Mesh2D(nodes,edges_,tri)
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