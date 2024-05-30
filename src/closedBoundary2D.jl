
struct ClosedBoundary2D

    "[nEdges] Edge lengths"
    edgeLength::Vector{Float64}

    "[2,nEdges] Surface tangent vector at each edge."
    edgeTangent::Matrix{Float64}

    "[2,nEdges] Surface normal vector at each edge."
    edgeNormal::Matrix{Float64}

    "[nNodes] node ID numbers along boundary"
    nodeIDs::Vector{Int32}

    "[nNodes] Edge length attributed to each node"
    nodeLength::Vector{Float64}

    "[2,nNodes] Surface tangent vector at each node."
    nodeTangent::Matrix{Float64}

    "[2,nNodes] Surface normal vector at each node."
    nodeNormal::Matrix{Float64}

    "[nNodes,2] Position vector at each node."
    nodeLocation::Matrix{Float64}

    "Boundary ID number in mesh fil."
    ID::Int64
end

function processMeshBoundaries(mesh::Mesh2D,inputData::Dict)

    # Boundary ID logic
    meshBoundaryIDs = sort(unique(mesh.edges[:,3]))
    nBoundaries = size(meshBoundaryIDs,1)
    if nBoundaries != 2
        error("This process is set up for one airfoil and one freestream boundary.")
    end

    # All boundary conditions are essential. Otherwise neeed additional logic | boundaryEdges = mesh.edges[mesh.edges[:,3].==boundaryIDs[i],1:2]
    boundaryEdgeNodes = sort(unique(mesh.edges[:,1:2]))
    nDof = size(mesh.nodes,1)
    f = trues(nDof) # free nodes
    f[boundaryEdgeNodes] .= false # free dof
    s = .!f  # fixed dof

    # Airfoil boundary
    airfoilBoundaryID = inputData["mesh"]["airfoilBoundaryID"]
    locb = indexin(airfoilBoundaryID,meshBoundaryIDs)
    if isnothing(locb)
        error("The airfoil boundary ID specified in the input file was not found in the mesh.")
    end
    airfoilBoundary = processClosedBoundary2D(mesh,airfoilBoundaryID)

    # Farfield boundary
    farfieldBoundaryID = inputData["mesh"]["farfieldBoundaryID"]
    locb = indexin(farfieldBoundaryID,meshBoundaryIDs)
    if isnothing(locb)
        error("The farfield boundary ID specified in the input file was not found in the mesh.")
    end
    farfieldBoundary = processClosedBoundary2D(mesh,farfieldBoundaryID)

    return airfoilBoundary, farfieldBoundary, f, s
end

function processClosedBoundary2D(mesh::Mesh2D,boundaryID::Int64)

    # using boring2D
    # cd("examples")
    # inputFileName = "n0012.toml"
    # using TOML
    # inputData = TOML.parsefile(inputFileName)
    # mesh = boring2D.readMesh(inputData["mesh"]["fileName"])
    # boundaryID = inputData["mesh"]["airfoilBoundaryID"]

    # get boundary edges
    edges = mesh.edges[mesh.edges[:,3].==boundaryID,1:2]
    nEdges = size(edges,1)

    # Check closed loop
    if edges[1,1] != edges[end,2]
        error("The boundary should be a closed loop.")
    end
    if edges[2:end,1] != edges[1:end-1,2]
        error("The boundary should be a closed loop.")
    end

    # compute edge lengths
    xe = mesh.nodes[edges,1]
    ye = mesh.nodes[edges,2]
    dx = xe[:,2].-xe[:,1]
    dy = ye[:,2].-ye[:,1]
    edgeLength = sqrt.(dx.^2+dy.^2)

    # Compute angles
    theta = atan.(dy,dx)
    Rz90 =  [0.0 -1.0
            1.0 0.0]

    # edge tangent and normal vectors
    edgeTangent = zeros(2,nEdges)
    edgeTangent[1,:] = cos.(theta)
    edgeTangent[2,:] = sin.(theta)
    edgeNormal = Rz90*edgeTangent

    # node quantities
    nodeIDs = edges[:,1]
    nodeLength = 0.5*[edgeLength[1]+edgeLength[end];edgeLength[1:end-1]+edgeLength[2:end]]
    nodeTangentSum = [edgeTangent[:,1]+edgeTangent[:,end] edgeTangent[:,1:end-1].+edgeTangent[:,2:end]]

    nodeTangentSumNorm = sqrt.(nodeTangentSum[1,:].^2 + nodeTangentSum[2,:].^2)
    nodeTangent = nodeTangentSum./ transpose(nodeTangentSumNorm)
    nodeNormal = Rz90*nodeTangent

    edgeTangent = zeros(2,nEdges)
    edgeTangent[1,:] = cos.(theta)
    edgeTangent[2,:] = sin.(theta)

    nodeLocation = mesh.nodes[nodeIDs,:]

    # Checks
    if abs( sum(edgeLength) - sum(nodeLength)) > 1e-12
        error("Summed edge length and node lengths should be equal across the boundary.")
    end

    # output data in structure
    return ClosedBoundary2D(edgeLength,edgeTangent,edgeNormal,nodeIDs,nodeLength,nodeTangent,nodeNormal,nodeLocation,boundaryID)
end

function getAdjacentElementIDs(mesh::Mesh2D,boundary::ClosedBoundary2D)

    # find element adjacent to edge by finding matching nodes
    boundaryNodes = mesh.edges[mesh.edges[:,3].==boundary.ID,1:2]

    nEdges = size(boundary.edgeLength,1)
    elementID = zeros(Int32,nEdges,1)
    for i = 1:nEdges
        edge1 = mesh.triangles .== boundaryNodes[i,1]
        edge2 = mesh.triangles .== boundaryNodes[i,2]
        bothEdge = edge1 .+ edge2
        sumBothEdge = sum(bothEdge,dims=2)
        elementIndex = findfirst(sumBothEdge.==2)
        elementID[i] = elementIndex[1]
    end

    return elementID
end

function computeGxγ(mesh::Mesh2D,triangleElements::Vector{TriangleElements},boundary::ClosedBoundary2D,f::BitVector)
    nDof = size(mesh.nodes,1)
    nEdges = size(boundary.edgeLength,1)
    
    adjacentElement = boring2D.getAdjacentElementIDs(mesh,boundary)
    
    gψγ = zeros(nDof)
    for i = 1:nEdges
        nodeDof = mesh.triangles[adjacentElement[i],:]
        
        velYCoeff =     transpose(triangleElements[adjacentElement[i]].dNdX[1,:]) # * psi[dof]
        velXCoeff = -1* transpose(triangleElements[adjacentElement[i]].dNdX[2,:]) # * psi[dof]
    
        vTanCoeff = velXCoeff.*boundary.edgeTangent[1,i] .+ velYCoeff.*boundary.edgeTangent[2,i]
        circulationCoeff = boundary.edgeLength[i].*vTanCoeff
    
        gψγ[nodeDof] = gψγ[nodeDof] .+ transpose(circulationCoeff)
    end

    # x only includes free dof | this approch will have issues with natural boundary conditions
    gxγ = transpose( [gψγ[f]; 0.0] )

    return gxγ
end

