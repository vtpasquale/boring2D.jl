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
    farfieldOnes = zeros(nDof,1) # Vector of ones at fairfield dof

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
                farfieldOnes[boundaryEdgeNodes] .= 1.0
            else
                error("Unknown boundary condition")
    
            end
        end
    end
    
    # fixed dof
    s = .!f
    
    # Farfield constant term
    kfshat = K[f,s]*farfieldOnes[s]

    # create constraint equation 
    constEqu = spzeros(nDof)
    dof = mesh.triangles[zeroVyEle,:]
    coefficents = triangleElements[zeroVyEle].dNdX[1,:]
    constEqu[dof] = coefficents

    # Solve equations unknown boundary constant + Kutta-condition constraint
    A = [K[f,f] kfshat; transpose(constEqu[f]) 0.0]
    b =  [- K[f,s]*phi[s]; 0.0] 
    x = A\b
    phi[f] = x[1:end-1]

    # ---------------------------
    # Partition and solve
    # Kff = K[f,f]
    # Kfs = K[f,s]
    # us = u[s]
    # rhs = -Kfs*us
    # uf = Kff\rhs
    # u[f] = uf
    # phi[f] = K[f,f]\(-K[f,s]*phi[s])

     # ---------------------------
     # Approximate nodal velocities
     nAdjacentElements = zeros(Int8,nDof,1)
     sumAdjacentVx = zeros(nDof,1)
     sumAdjacentVy = zeros(nDof,1)

    # ---------------------------
    # Recover element velocities
    velX = zeros(nEle,1)
    velY = zeros(nEle,1)
    for i = 1:nEle
        dof = mesh.triangles[i,:]
        velY[i] =     transpose(triangleElements[i].dNdX[1,:]) * phi[dof]
        velX[i] = -1* transpose(triangleElements[i].dNdX[2,:]) * phi[dof]

        nAdjacentElements[dof] .+= 1
        sumAdjacentVy[dof] .+= velY[i]
        sumAdjacentVx[dof] .+= velX[i]
    end
    nodeVy = sumAdjacentVy ./ nAdjacentElements
    nodeVx = sumAdjacentVx ./ nAdjacentElements
    nodeVmag = sqrt.(nodeVx.^2 + nodeVy.^2)
    Cp = 1.0 .- (nodeVmag ./ freestreamV).^2

    # Create solution dictionaries
    pointOutput = Dict("psi"=>phi,"vX"=>nodeVx,"vY"=>nodeVy,"vMag"=>nodeVmag,"Cp"=>Cp)
    cellOutput = Dict("vX"=>[velX],"vY"=>[velY])
    
    # Write output data to file
    boring2D.writeSolution("streamSolution.vtu",mesh,pointOutput,cellOutput)   


    for i = 1:nBoundaries
        if boundaryConditionType[i] == 1
            # integrate airfoil circulations - update if logic to generalize

            boundaryEdges = mesh.edges[mesh.edges[:,3].==boundaryIDs[i],1:2]
            xe = mesh.nodes[boundaryEdges,1]
            ye = mesh.nodes[boundaryEdges,2]
            dx = xe[:,2].-xe[:,1]
            dy = ye[:,2].-ye[:,1]
            dl= sqrt.(dx.^2+dy.^2)
            
            Cpe = 0.5*(Cp[boundaryEdges[:,1]].+Cp[boundaryEdges[:,2]])
            # Cfx=sum(Cpe.*dy)
            # Cfy=sum(Cpe.*dx)

            theta = atan.(dy,dx)
            # OutNormalAngle = theta .+ pi/2
            inNormalAngle = theta .- pi/2
            Cfx = sum( Cpe.*dl.*cos.(inNormalAngle) )
            Cfy = sum( Cpe.*dl.*sin.(inNormalAngle) )

            CL = Cfx*sin(freestreamAlpha) + Cfy*cos(freestreamAlpha) 
            CD = Cfx*cos(freestreamAlpha) + Cfy*sin(freestreamAlpha) 

            println("CL=$(CL)")
            println("CD=$(CD)")

            # vXe = 0.5*(nodeVx[boundaryEdges,1].+nodeVx[boundaryEdges,2])
            # vYe = 0.5*(nodeVy[boundaryEdges,1].+nodeVy[boundaryEdges,2])
            # Vt = vXe.*cos.(theta) .+ vYe.*sin.(theta)
            # Vn = vXe.*sin.(theta) .+ vYe.*cos.(theta)
            # dl= sqrt.(dx.^2+dy.^2)
            # ut = u.*cos(theta)+v.*sin(theta)
            # Circ = -sum(Vt.*dl);
        end
    end


    return 0
end