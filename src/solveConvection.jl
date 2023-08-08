using SparseArrays

function solveConvection(meshFileName::AbstractString,k::Float64,boundaryTemps::Vector{Float64})
    # k is conductivity

    mesh = readMesh(meshFileName)
    triangleElements = TriangleElements(mesh)
    K = assembleConvectionStiffness(mesh,triangleElements,k)

    # Determine fixed and free dof
    nDof = size(mesh.nodes,1)
    nEdges = size(mesh.edges,1)
    edgeNodes = sort(unique(mesh.edges[:,1:2])); # Fixed dof
    f = trues(nDof);
    f[edgeNodes] .= false;
    s = .!f

    # Assemble boundary temperatures
    u = zeros(nDof,1);
    for i = 1:nEdges
        # this will overwrite temperatures if nodes are on more than one edge
        u[mesh.edges[i,1:2]] .= boundaryTemps[mesh.edges[i,3]]; 
    end

    # Partition and solve
    Kff = K[f,f]
    Kfs = K[f,s]
    us = u[s]
    rhs = -Kfs*us
    # uf = Kff\rhs
    # u[f] = uf
    # # u[f] = K[f,f]\(-K[f,s]*u[s]);

    # # Create solution dictionaries
    # pointOutput = Dict("Temp"=>u)
    # cellOutput = Dict()

    # # Write output data to file
    # boring2D.writeSolution("out.vtu",mesh,pointOutput,cellOutput)

    return (Kff,rhs)
end

function assembleConvectionStiffness(mesh::Mesh2D,triangleElements::Vector{TriangleElements},k::Float64)

    # Number of elements and nodes
    nElements = size(mesh.triangles,1);
    nNodes = size(mesh.nodes,1);

    # One Temperature DOF per node
    nDof = nNodes;

    # # Global matrices (slow method)
    # K = spzeros(nDof,nDof);

    # Global matrices with triplets
    Ks = SparseTriplet(nNodes);

    # Assemble stiffness matrix
    for i = 1:nElements
        dof = mesh.triangles[i,:]
        dNdX = triangleElements[i].dNdX
        area = triangleElements[i].area
        k_e = transpose(dNdX)*k*dNdX*area
        addMatrix!(Ks,k_e,dof)
    end
    K = convertToSparseMatrix(Ks,nDof,nDof);
    return K
end
