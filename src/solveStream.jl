using SparseArrays

function solveStream(meshFileName::AbstractString,boundaryConditionType::Vector{Int64},freestreamV::Float64,freestreamAlpha::Float64,essentialPhi::Vector{Float64},zeroVyEle::Int64)

    mesh = readMesh(meshFileName)

    u0 = freestreamV*cos(freestreamAlpha)
    v0 = freestreamV*sin(freestreamAlpha)
    
    triangleElements = TriangleElements(mesh)
    
    k = 1.0 # for potential flow equations
    K = assembleConvectionStiffness(mesh,triangleElements,k)
    
    # Determine fixed and free dof
    nDof = size(mesh.nodes,1)
    nEle = size(mesh.triangles,1)   

    boundaryIDs = sort(unique(mesh.edges[:,3]))
    nBoundaries = size(boundaryIDs,1)    

    f = trues(nDof) # free nodes
    phi = zeros(nDof,1)
    
    for i = 1:nBoundaries
        if boundaryConditionType[i] == 0
            # natural BC
        else 
            boundaryEdges = mesh.edges[mesh.edges[:,3].==boundaryIDs[i],1:2]
            boundaryEdgeNodes = sort(unique(boundaryEdges[:,1:2]))
            f[boundaryEdgeNodes] .= false
    
            if boundaryConditionType[i] == 1
                phi[boundaryEdgeNodes] .= essentialPhi[i]
    
            elseif boundaryConditionType[i] == 2
                xs = mesh.nodes[boundaryEdgeNodes,1]
                ys = mesh.nodes[boundaryEdgeNodes,2]
                phi[boundaryEdgeNodes] .= v0.*xs - u0.*ys
    
            else
                error("Unknown boundary condition")
    
            end
        end
    end
    
    # fixed dof
    s = .!f
    
    # # ---------------------------
    # # Farfield vortex (results in zero phi)
    # ds = ( xs.^2 + ys.^2 ).^(.5)
    # thetas = atan.(ys,xs)
    # u1 =  (1.0 ./ds).*sin.(thetas)
    # u2 = -(1.0 ./ds).*cos.(thetas)
    
    # using Plots
    # plot(thetas,phi[s] )
    # plot(xs,ys)
    
    # # ---------------------------
    # # create source load vector
    # w3 = 1.0/3.0
    # p = zeros(nDof,1)
    # for i = 1:nEle
    #     weightedArea = w3*triangleElements[i].area
    #     p[mesh.triangles[i,1]] += weightedArea
    #     p[mesh.triangles[i,2]] += weightedArea
    #     p[mesh.triangles[i,3]] += weightedArea
    # end

    # # create constraint equation 
    # constEqu = spzeros(nDof)
    # dof = mesh.triangles[zeroVyEle,:]
    # coefficents = triangleElements[zeroVyEle].dNdX[2,:]
    # constEqu[dof] = coefficents

    # # Solve equations with constraint and source strength
    # A = [K[f,f] -p[f]; transpose(constEqu[f]) 0.0]
    # b =  [- K[f,s]*phi[s]; 0.0] 
    # x = A\b
    # phi[f] = x[1:end-1]

    # ---------------------------
    # Partition and solve
    # Kff = K[f,f]
    # Kfs = K[f,s]
    # us = u[s]
    # rhs = -Kfs*us
    # uf = Kff\rhs
    # u[f] = uf
    phi[f] = K[f,f]\(-K[f,s]*phi[s])

    # ---------------------------
    # Recover velocities
    velX = zeros(nEle,1)
    velY = zeros(nEle,1)
    for i = 1:nEle
        dof = mesh.triangles[i,:]
        velY[i] =     transpose(triangleElements[i].dNdX[1,:]) * phi[dof]
        velX[i] = -1* transpose(triangleElements[i].dNdX[2,:]) * phi[dof]
    end
    # Create solution dictionaries
    pointOutput = Dict("phi"=>phi)
    cellOutput = Dict("vX"=>[velX],"vY"=>[velY])
    
    # Write output data to file
    # boring2D.writeSolution("potentialFlow.vtu",mesh,pointOutput,cellOutput)   

    return pointOutput,cellOutput
end