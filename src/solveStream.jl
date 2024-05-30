using SparseArrays
using TOML

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
    nEle = size(mesh.triangles,1)
    triangleElements = boring2D.TriangleElements(mesh)
    airfoilBoundary, farfieldBoundary, f, s = boring2D.processMeshBoundaries(mesh,inputData)
    gxγ = boring2D.computeGxγ(mesh,triangleElements,airfoilBoundary,f)

    # Assemble coefficent matrix
    K = boring2D.assembleConvectionStiffness(mesh,triangleElements,1.0)
    g = boring2D.assembleKuttaConditionConstraint(mesh,inputData,triangleElements)
    esa = boring2D.assembleEsa(nDof,farfieldBoundary)
    A = boring2D.assembleA(K,g,esa,f,s)

    # Assemble right-hand side
    ψsk = boring2D.assembleψsk(nDof,farfieldBoundary,inputData)
    b = boring2D.assembleB(K,ψsk,f,s)

    # Stream function solution
    x = A\Array(b)

    # Stream contour output
    ψ = boring2D.recoverψ(x,ψsk,esa,f,s)

    # Lift coefficent
    Cl = (2/inputData["freestream"]["velocity"])* gxγ*x

    # Adjoint solution
    ∂Cl∂x = (2/inputData["freestream"]["velocity"])* gxγ
    λ = transpose(A)\transpose(∂Cl∂x)

    # Adjoint contour output
    fx = boring2D.getFxfromF(f)
    lamdaOut = zeros(nDof+1,1)
    lamdaOut[fx] = λ

    # Sensitivity wrt α
    ∂ψsk_∂α = boring2D.assemble_∂ψsk_∂α(nDof,farfieldBoundary,inputData)
    ∂b_∂α = boring2D.assemble_∂b_∂α(K,∂ψsk_∂α,f,s)
    ∂Cl_∂α = transpose(λ)*∂b_∂α # per radian | ∂Cl_∂α_deg = ∂Cl_∂α*pi/180

    pointOutput = Dict("psi"=>ψ,"lamda"=>lamdaOut[1:end-1])
    cellOutput = Dict()

    # Write output data to file
    boring2D.writeSolution("streamSolution.vtu",mesh,pointOutput,cellOutput)

    #  # ---------------------------
    #  # Approximate nodal velocities
    #  nAdjacentElements = zeros(Int8,nDof,1)
    #  sumAdjacentVx = zeros(nDof,1)
    #  sumAdjacentVy = zeros(nDof,1)
    #  sumAdjacentdVx = zeros(nDof,1)
    #  sumAdjacentdVy = zeros(nDof,1)
    # # ---------------------------
    # # Recover element velocities and adjoint velocities
    # velX = zeros(nEle,1)
    # velY = zeros(nEle,1)
    # dVelX = zeros(nEle,1)
    # dVelY = zeros(nEle,1)
    # for i = 1:nEle
    #     dof = mesh.triangles[i,:]
    #     velY[i] =     transpose(triangleElements[i].dNdX[1,:]) * psi[dof]
    #     velX[i] = -1* transpose(triangleElements[i].dNdX[2,:]) * psi[dof]
    #     dVelX[i] =    transpose(triangleElements[i].dNdX[1,:]) * lamdaOut[dof]
    #     dVelY[i] = -1*transpose(triangleElements[i].dNdX[2,:]) * lamdaOut[dof]

    #     nAdjacentElements[dof] .+= 1
    #     sumAdjacentVy[dof] .+= velY[i]
    #     sumAdjacentVx[dof] .+= velX[i]
    #     sumAdjacentdVy[dof] .+= dVelX[i]
    #     sumAdjacentdVx[dof] .+= dVelY[i]
    # end
    # nodeVy = sumAdjacentVy ./ nAdjacentElements
    # nodeVx = sumAdjacentVx ./ nAdjacentElements
    # nodeDVy = sumAdjacentdVy ./ nAdjacentElements
    # nodeDVx = sumAdjacentdVx ./ nAdjacentElements
    # nodeVmag = sqrt.(nodeVx.^2 + nodeVy.^2)
    # Cp = 1.0 .- (nodeVmag ./ inputData["freestream"]["velocity"]).^2

    # # Create solution dictionaries
    # pointOutput = Dict("psi"=>psi,"lambda"=>lamdaOut[1:end-1],"vX"=>nodeVx,"vY"=>nodeVy,"vMag"=>nodeVmag,"Cp"=>Cp,"nodeDVy"=>nodeDVy,"nodeDVx"=>nodeDVx)
    # cellOutput = Dict("vX"=>[velX],"vY"=>[velY],"dVelX"=>[dVelX],"dVelY"=>[dVelY])

    # # Write output data to file
    # boring2D.writeSolution("streamSolution.vtu",mesh,pointOutput,cellOutput)

    return Cl, ∂Cl_∂α
end

function assembleKuttaConditionConstraint(mesh::Mesh2D,inputData::Dict,triangleElements::Vector{TriangleElements})
    nDof = size(mesh.nodes,1)
    zeroVyEle = inputData["mesh"]["kuttaConditionElementID"]
    g = spzeros(nDof)
    dof = mesh.triangles[zeroVyEle,:]
    coefficents = triangleElements[zeroVyEle].dNdX[1,:]
    g[dof] = coefficents
    return g
end

function assembleEsa(nDof::Int64,farfieldBoundary::ClosedBoundary2D)
    esa = spzeros(nDof) # Vector of ones at fairfield dof
    esa[farfieldBoundary.nodeIDs] .= 1.0
    return esa
end

function assembleA(K::SparseMatrixCSC{Float64},g::SparseVector{Float64},esa::SparseVector{Float64},f::BitVector,s::BitVector)
    keas = K[f,s]*esa[s]
    A = [K[f,f] keas; transpose(g[f]) 0.0]
    return A
end

function assembleψsk(nDof::Int64,farfieldBoundary::ClosedBoundary2D,inputData::Dict)     
    Vinf = inputData["freestream"]["velocity"] # Moved to input argument to enable sensitivities
    alphaDeg = inputData["freestream"]["alphaDeg"]
    u0 = Vinf*cosd(alphaDeg)
    v0 = Vinf*sind(alphaDeg)
    ψsk = spzeros(nDof) 
    xs = farfieldBoundary.nodeLocation[:,1]
    ys = farfieldBoundary.nodeLocation[:,2]
    ψsk[farfieldBoundary.nodeIDs] .= v0.*xs - u0.*ys
    return ψsk
end

function assemble_∂ψsk_∂α(nDof::Int64,farfieldBoundary::ClosedBoundary2D,inputData::Dict)     
    Vinf = inputData["freestream"]["velocity"] # Moved to input argument to enable sensitivities
    alphaDeg = inputData["freestream"]["alphaDeg"]
    du0 = -Vinf*sind(alphaDeg) # This is the derivative wrt alpha in radians, not degres
    dv0 =  Vinf*cosd(alphaDeg)
    ∂ψsk_∂α = spzeros(nDof) 
    xs = farfieldBoundary.nodeLocation[:,1]
    ys = farfieldBoundary.nodeLocation[:,2]
    ∂ψsk_∂α[farfieldBoundary.nodeIDs] .= dv0.*xs - du0.*ys
    return ∂ψsk_∂α
end

function assembleB(K::SparseMatrixCSC{Float64},ψsk::SparseVector{Float64},f::BitVector,s::BitVector)
    b =  [-K[f,s]*ψsk[s]; 0.0] 
    return b
end

function assemble_∂b_∂α(K::SparseMatrixCSC{Float64},∂ψsk_∂α::SparseVector{Float64},f::BitVector,s::BitVector)
    ∂b_∂α =  [-K[f,s]*∂ψsk_∂α[s]; 0.0] 
    return ∂b_∂α
end

function recoverψ(x::Vector{Float64},ψsk::SparseVector{Float64},esa::SparseVector{Float64},f::BitVector,s::BitVector)
    nDof = size(f,1)
    ψ = Array{Float64}(undef, nDof)
    ψ[f] = x[1:end-1]
    c3 = x[end]
    ψ[s] = ψsk[s] .+ c3*esa[s]
    return ψ
end

function getFxfromF(f::BitVector)
    nDof = size(f,1)
    fx = falses(nDof+1)
    fx[1:end-1] .= f
    fx[end] = true
    return fx
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