using SparseArrays
using TOML

struct SteamFunctionLinearSystem
    "[nFreePsi+1,nFreePsi+1] Coefficent matrix."
    A::SparseMatrixCSC{Float64}

    "[nFreePsi+1,1] Right-hand side."
    b::Vector{Float64}

    "[nPsi] Free psi dof."
    f::BitVector
    
    Kfs::SparseMatrixCSC{Float64}

end

function getFxfromF(sfls::SteamFunctionLinearSystem)
    nDof = size(sfls.f,1)
    fx = falses(nDof+1)
    fx[1:end-1] .= sfls.f
    fx[end] = true
    return fx
end

struct FarfieldEdgeData
    "[nFarfieldEdgeNodes] Node ID numbers."
    farfieldEdgeNodes::Vector{Int32}

    "[nFarfieldEdgeNodes] x position of farfield edge nodes."
    xs::Vector{Float64}

    "[nFarfieldEdgeNodes] y position of farfield edge nodes."
    ys::Vector{Float64}

    "Freestream velocity"
    Vinf::Float64
end

function solveStream(inputFileName::AbstractString)

    # inputFileName = "n0012.toml"
    # using SparseArrays
    # using TOML
    
    # Read input file
    inputData = TOML.parsefile(inputFileName)

    # Echo input file
    TOML.print(inputData)

    # Process mesh
    mesh = boring2D.readMesh(inputData["mesh"]["fileName"])
    nDof = size(mesh.nodes,1)
    # nEle = size(mesh.triangles,1)
    triangleElements = boring2D.TriangleElements(mesh)

    # Assemble and solve linear system
    sfls, psi, farfieldEdgeData = boring2D.assembleStreamFunctionLinearSystem(mesh,inputData,triangleElements)
    A = sfls.A
    b = sfls.b
    f = sfls.f
    x = A\b
    psi[f] = x[1:end-1]

    # Recover lift coefficent
    fx = boring2D.getFxfromF(sfls)
    xToCl = boring2D.compute_xToCl(mesh,inputData,triangleElements,fx)
    Cl = xToCl*x

    # Lift coefficent adjoint solution
    λ = transpose(A)\transpose(xToCl)
    lamdaOut = zeros(nDof+1,1)
    lamdaOut[fx] = λ

    # Lift derivative with respect to α
    ∂b∂α = boring2D.compute_∂b∂α(inputData,sfls,farfieldEdgeData)
    ∂Cl∂α = transpose(λ)*∂b∂α

    # Create solution dictionaries
    pointOutput = Dict("psi"=>psi,"lambda"=>lamdaOut[1:end-1])
    cellOutput = Dict()

    # Write output data to file
    boring2D.writeSolution("streamSolution.vtu",mesh,pointOutput,cellOutput)

    return Cl, ∂Cl∂α
end

function assembleStreamFunctionLinearSystem(mesh::Mesh2D,inputData::Dict,triangleElements::Vector{TriangleElements})
    # Assemble domain and boundary matricies
    K = assembleConvectionStiffness(mesh,triangleElements,1.0)
    f, s, psi, esa, farfieldEdgeData = streamBoundaryConditions(mesh,inputData)
    
    # free set in x vector
    nDof = size(mesh.nodes,1)
    fx = falses(nDof+1)
    fx[1:end-1] .= f
    fx[end] = true

    # Farfield constant term
    keas = K[f,s]*esa[s]

    # create constraint equation
    zeroVyEle = inputData["mesh"]["kuttaConditionElementID"]
    g = spzeros(nDof)
    dof = mesh.triangles[zeroVyEle,:]
    coefficents = triangleElements[zeroVyEle].dNdX[1,:]
    g[dof] = coefficents

    # Solve equations with unknown boundary constant + Kutta-condition constraint
    A = [K[f,f] keas; transpose(g[f]) 0.0]
    b =  [-K[f,s]*psi[s]; 0.0] 

    Kfs = K[f,s]
    sfls = SteamFunctionLinearSystem(A,b,f,Kfs)    

    return sfls, psi, farfieldEdgeData
end

function compute_∂b∂α(inputData::Dict,sfls::SteamFunctionLinearSystem,farfieldEdgeData::FarfieldEdgeData)
    
    Kfs = sfls.Kfs
    f = sfls.f
    s = .!f
    nDof = size(f,1)
    dpsi = zeros(nDof,1)
    
    Vinf = inputData["freestream"]["velocity"] # Moved to input argument to enable sensitivities
    alphaDeg = inputData["freestream"]["alphaDeg"]
    du0 = -Vinf*sind(alphaDeg) # This is the derivative wrt alpha in radians, not degres
    dv0 =  Vinf*cosd(alphaDeg)
    dpsi[farfieldEdgeData.farfieldEdgeNodes] .= dv0.*farfieldEdgeData.xs - du0.*farfieldEdgeData.ys

    ∂b∂α = [-Kfs*dpsi[s]; 0.0] 

    return ∂b∂α
end

function streamBoundaryConditions(mesh::Mesh2D,inputData::Dict)       
    # Velocities
    Vinf = inputData["freestream"]["velocity"] # Moved to input argument to enable sensitivities
    alphaDeg = inputData["freestream"]["alphaDeg"]
    u0 = Vinf*cosd(alphaDeg)
    v0 = Vinf*sind(alphaDeg)

    # Preallocate arrays
    nDof = size(mesh.nodes,1)
    f = trues(nDof) # free nodes
    psi = zeros(nDof,1)
    esa = zeros(nDof,1) # Vector of ones at fairfield dof

    # Boundary ID logic
    meshBoundaryIDs = sort(unique(mesh.edges[:,3]))
    nBoundaries = size(meshBoundaryIDs,1)
    if nBoundaries != 2
        error("This process is set up for one airfoil and one freestream boundary.")
    end

    # All boundary conditions are essential. Otherwise neeed additional logic | boundaryEdges = mesh.edges[mesh.edges[:,3].==boundaryIDs[i],1:2]
    boundaryEdgeNodes = sort(unique(mesh.edges[:,1:2]))
    f[boundaryEdgeNodes] .= false # free dof
    s = .!f  # fixed dof

    # Airfoil boundary
    airfoilBoundaryID = inputData["mesh"]["airfoilBoundaryID"]
    locb = indexin(airfoilBoundaryID,meshBoundaryIDs)
    if isnothing(locb)
        error("The airfoil boundary ID specified in the input file was not found in the mesh.")
    end
    # psi is preallocated to zero, which is the correct boundary value at the airfoil (c2=0)

    # Farfield boundary
    farfieldBoundaryID = inputData["mesh"]["farfieldBoundaryID"]
    locb = indexin(farfieldBoundaryID,meshBoundaryIDs)
    if isnothing(locb)
        error("The farfield boundary ID specified in the input file was not found in the mesh.")
    end
    farfieldEdges = mesh.edges[mesh.edges[:,3].==farfieldBoundaryID,1:2]
    farfieldEdgeNodes = sort(unique(farfieldEdges[:,1:2]))
    xs = mesh.nodes[farfieldEdgeNodes,1]
    ys = mesh.nodes[farfieldEdgeNodes,2]
    psi[farfieldEdgeNodes] .= v0.*xs - u0.*ys
    esa[farfieldEdgeNodes] .= 1.0

    farfieldEdgeData = FarfieldEdgeData(farfieldEdgeNodes,xs,ys,Vinf)
    return f, s, psi, esa, farfieldEdgeData
end

function compute_xToCl(mesh::Mesh2D,inputData::Dict,triangleElements::Vector{TriangleElements},fx::BitVector)
    nDof = size(mesh.nodes,1)

    # get airfoil edges
    airfoilBoundaryID = inputData["mesh"]["airfoilBoundaryID"]
    airfoilEdges = mesh.edges[mesh.edges[:,3].==airfoilBoundaryID,1:2]

    # compute edge lengths
    xe = mesh.nodes[airfoilEdges,1]
    ye = mesh.nodes[airfoilEdges,2]
    dx = xe[:,2].-xe[:,1]
    dy = ye[:,2].-ye[:,1]
    dl = sqrt.(dx.^2+dy.^2)

    # Compute angles
    theta = atan.(dy,dx)

    # find element adjacent to edge by finding matching nodes
    nEdges = size(airfoilEdges,1)
    elementID = zeros(Int32,nEdges,1)
    for i = 1:nEdges
        edge1 = mesh.triangles .== airfoilEdges[i,1]
        edge2 = mesh.triangles .== airfoilEdges[i,2]
        bothEdge = edge1 .+ edge2
        sumBothEdge = sum(bothEdge,dims=2)
        elementIndex = findfirst(sumBothEdge.==2)
        elementID[i] = elementIndex[1]
    end

    # Matrix that maps x=[psi;0] to circulation
    xToCirculationMap = zeros(1,nDof+1)
    for i = 1:nEdges
        nodeDof = mesh.triangles[elementID[i],:]
        
        velYCoeff =     transpose(triangleElements[elementID[i]].dNdX[1,:]) # * psi[dof]
        velXCoeff = -1* transpose(triangleElements[elementID[i]].dNdX[2,:]) # * psi[dof]

        vTanCoeff = velXCoeff.*cos.(theta[i]) .+ velYCoeff.*sin.(theta[i])
        circulationCoeff = dl[i].*vTanCoeff

        xToCirculationMap[1,nodeDof] = xToCirculationMap[1,nodeDof] .+ transpose(circulationCoeff)
    end

    # Kutta-Joukowski with chord=1
    xToCl = (2/inputData["freestream"]["velocity"])* 
            transpose(xToCirculationMap[fx])  # indexing converts array to ID so transpose required

    return xToCl
end

# function postprocess()
#      # ---------------------------
#      # Approximate nodal velocities
#      nAdjacentElements = zeros(Int8,nDof,1)
#      sumAdjacentVx = zeros(nDof,1)
#      sumAdjacentVy = zeros(nDof,1)

#     # ---------------------------
#     # Recover element velocities
#     velX = zeros(nEle,1)
#     velY = zeros(nEle,1)
#     for i = 1:nEle
#         dof = mesh.triangles[i,:]
#         velY[i] =     transpose(triangleElements[i].dNdX[1,:]) * psi[dof]
#         velX[i] = -1* transpose(triangleElements[i].dNdX[2,:]) * psi[dof]

#         nAdjacentElements[dof] .+= 1
#         sumAdjacentVy[dof] .+= velY[i]
#         sumAdjacentVx[dof] .+= velX[i]
#     end
#     nodeVy = sumAdjacentVy ./ nAdjacentElements
#     nodeVx = sumAdjacentVx ./ nAdjacentElements
#     nodeVmag = sqrt.(nodeVx.^2 + nodeVy.^2)
#     Cp = 1.0 .- (nodeVmag ./ freestreamV).^2

#     # Create solution dictionaries
#     pointOutput = Dict("psi"=>psi,"vX"=>nodeVx,"vY"=>nodeVy,"vMag"=>nodeVmag,"Cp"=>Cp)
#     cellOutput = Dict("vX"=>[velX],"vY"=>[velY])
    
#     # Write output data to file
#     boring2D.writeSolution("streamSolution.vtu",mesh,pointOutput,cellOutput)   


#     for i = 1:nBoundaries
#         if boundaryConditionType[i] == 1
#             # integrate airfoil circulations - update if logic to generalize

#             boundaryEdges = mesh.edges[mesh.edges[:,3].==boundaryIDs[i],1:2]
#             xe = mesh.nodes[boundaryEdges,1]
#             ye = mesh.nodes[boundaryEdges,2]
#             dx = xe[:,2].-xe[:,1]
#             dy = ye[:,2].-ye[:,1]
#             dl= sqrt.(dx.^2+dy.^2)
            
#             Cpe = 0.5*(Cp[boundaryEdges[:,1]].+Cp[boundaryEdges[:,2]])
#             # Cfx=sum(Cpe.*dy)
#             # Cfy=sum(Cpe.*dx)

#             theta = atan.(dy,dx)
#             # OutNormalAngle = theta .+ pi/2
#             inNormalAngle = theta .- pi/2
#             Cfx = sum( Cpe.*dl.*cos.(inNormalAngle) )
#             Cfy = sum( Cpe.*dl.*sin.(inNormalAngle) )

#             CL = -Cfx*sin(freestreamAlpha) + Cfy*cos(freestreamAlpha) 
#             CD =  Cfx*cos(freestreamAlpha) + Cfy*sin(freestreamAlpha) 

#             println("CL=$(CL)")
#             println("CD=$(CD)")

#             # println("Cfx*sin(freestreamAlpha)=$(Cfx*sin(freestreamAlpha))")
#             # println("Cfy*cos(freestreamAlpha)=$(Cfy*cos(freestreamAlpha))")
#             # println("Cfx*cos(freestreamAlpha)=$(Cfx*cos(freestreamAlpha))")
#             # println("Cfy*sin(freestreamAlpha)=$(Cfy*sin(freestreamAlpha))")

#             vXe = 0.5*(nodeVx[boundaryEdges[:,1]].+nodeVx[boundaryEdges[:,2]])
#             vYe = 0.5*(nodeVy[boundaryEdges[:,1]].+nodeVy[boundaryEdges[:,2]])
#             vTan = vXe.*cos.(theta) .+ vYe.*sin.(theta)
#             circulation = sum(vTan .*dl)
#             CLc = 2*circulation/freestreamV

#             println("CLc=$(CLc)")

#             # Vn = vXe.*sin.(theta) .+ vYe.*cos.(theta)
#             # dl= sqrt.(dx.^2+dy.^2)
#             # ut = u.*cos(theta)+v.*sin(theta)
#             # Circ = -sum(Vt.*dl);
#         end
#     end
#     return 0
# end