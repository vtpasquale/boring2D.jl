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
    lamdaOut = zeros(nDof+1)
    lamdaOut[fx] = λ

    # Sensitivity wrt α
    ∂ψsk_∂α = boring2D.assemble_∂ψsk_∂α(nDof,farfieldBoundary,inputData)
    ∂r_∂α = boring2D.assemble_∂r_∂α(K,∂ψsk_∂α,f,s)
    ∂Cl_∂α = transpose(λ)*∂r_∂α # per radian | ∂Cl_∂α_deg = ∂Cl_∂α*pi/180
    
    # Shape sensitivities
    ∂g_∂b, ∂r_∂b = boring2D.computeShapePartials(mesh,x,inputData,size(airfoilBoundary.nodeIDs,1))

    # This sign    v     matches finite difference result. The math needs to be resolved
    ∂Cl_∂b = ∂g_∂b - transpose(λ)*∂r_∂b

    boring2D.surfaceOutput2Vtk("adjointSurfaceSensitivity.vtk",airfoilBoundary,∂Cl_∂b[:]./airfoilBoundary.nodeLength)

    # Element Output
    velX, velY = recoverVelocities(ψ,mesh,triangleElements)
    dVelX, dVelY = recoverVelocities(lamdaOut,mesh,triangleElements)
    velMag = sqrt.(velX.^2 + velY.^2)
    Cp = 1.0 .- (velMag ./ inputData["freestream"]["velocity"]).^2

    # Print Kutta–Joukowski Cl results
    println("Lift coefficient from Kutta-Joukowski")
    println("Cl = $Cl")

    # Print force coefficients from integrated pressure coefficient
    boring2D.computeCfFromCp(Cp,mesh,airfoilBoundary,inputData)

    # Output dictionaries for output file
    pointOutput = Dict("psi"=>ψ,"lamda"=>lamdaOut[1:end-1])
    cellOutput = Dict("vX"=>[velX],"vY"=>[velY],"dVelX"=>[dVelX],"dVelY"=>[dVelY],"velMag"=>[velMag],"Cp"=>[Cp])

    # Write output data to file
    boring2D.writeSolution("streamSolution.vtu",mesh,pointOutput,cellOutput)

    return Cl, ∂Cl_∂α
end

function computeCirculationMatrixAndResidual(shapeVar::Vector,mesh0::Mesh2D,x::Vector{Float64},inputData::Dict)
    # This function is intended for shape sensitivity analysis only. 

    # mesh = boring2D.readMesh(inputData["mesh"]["fileName"])
    mesh = deepcopy(mesh0)
    nDof = size(mesh.nodes,1)

    # Placeholder boundary data for design update
    airfoilBoundary1, farfieldBoundary1, f1, s1 = boring2D.processMeshBoundaries(mesh,inputData)

    # Deform mesh normals based on design variables
    if size(shapeVar,1) != size(airfoilBoundary1.nodeIDs,1)
        error("Check design variable definition. size(shapeVar,1) != size(airfoilBoundary.nodeIDs,1)")
    end
    mesh.nodes[airfoilBoundary1.nodeIDs,1] += shapeVar.*airfoilBoundary1.nodeNormal[1,:]
    mesh.nodes[airfoilBoundary1.nodeIDs,2] += shapeVar.*airfoilBoundary1.nodeNormal[2,:]

    # Process mesh
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

    # Residual
    r = A*x-b

    return gxγ, r
end

function computeShapePartials(mesh0::Mesh2D,x::Vector{Float64},inputData::Dict,nVars::Int64)
    # Assemble partial derivatives w.r.t. shape design varaiables using finite difference

    # Check residual for zero-value design variables
    # nVars = size(airfoilBoundary.nodeIDs,1)
    shapeVar = zeros(nVars)
    gxγ0, r0 = computeCirculationMatrixAndResidual(shapeVar,mesh0,x,inputData)
    if maximum(abs.(r0)) > 1e-10
        error("The residual is higher than expected. There could be a design varaible or processing issue.")
    end

    # Compute partials using finite difference - automatic differntiation has unresolved issues
    nR = size(r0,1)
    ∂r_∂b = zeros(nR,nVars)
    ∂g_∂b = zeros(1, nVars)

    Δ = 1e-6
    for i = 1:nVars
        shapeVar = zeros(nVars)
        shapeVar[i] = Δ
        gxγ, r = boring2D.computeCirculationMatrixAndResidual(shapeVar,mesh0,x,inputData)
        
        ∂gxγ_∂b = (gxγ.-gxγ0) ./ Δ
        ∂g_∂bi = (2/inputData["freestream"]["velocity"]) * ( ∂gxγ_∂b*x )
        ∂g_∂b[i] = ∂g_∂bi[1] 

        ∂r_∂b[:,i] = (r.-r0) ./ Δ
    end

    return ∂g_∂b, ∂r_∂b
end

function computeRotatedResidual(δα::Real,mesh::Mesh2D,x::Vector{Float64},inputData0::Dict)
    # This function is for checking the analytic α sensitivy only

    inputData = deepcopy(inputData0)
    inputData["freestream"]["alphaDeg"] += 180.0/pi*δα

    # Process mesh
    # mesh = boring2D.readMesh(inputData["mesh"]["fileName"])
    nDof = size(mesh.nodes,1)
    airfoilBoundary, farfieldBoundary, f, s = boring2D.processMeshBoundaries(mesh,inputData)

    # Process elements
    triangleElements = boring2D.TriangleElements(mesh)

    # Assemble coefficent matrix
    K = boring2D.assembleConvectionStiffness(mesh,triangleElements,1.0)
    g = boring2D.assembleKuttaConditionConstraint(mesh,inputData,triangleElements)
    esa = boring2D.assembleEsa(nDof,farfieldBoundary)
    A = boring2D.assembleA(K,g,esa,f,s)

    # Assemble right-hand side
    ψsk = boring2D.assembleψsk(nDof,farfieldBoundary,inputData) 
    b = boring2D.assembleB(K,ψsk,f,s)

    r = Array( A*x-b )

    return r
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

function assemble_∂r_∂α(K::SparseMatrixCSC{Float64},∂ψsk_∂α::SparseVector{Float64},f::BitVector,s::BitVector)
    ∂r_∂α =  [-K[f,s]*∂ψsk_∂α[s]; 0.0] 
    return ∂r_∂α
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

function recoverVelocities(ψ::Vector{Float64},mesh::Mesh2D,triangleElements::Vector{TriangleElements})
    nEle = size(mesh.triangles,1)
    velX = zeros(nEle)
    velY = zeros(nEle)
    for i = 1:nEle
        dof = mesh.triangles[i,:]
        velY[i] =     transpose(triangleElements[i].dNdX[1,:]) * ψ[dof]
        velX[i] = -1* transpose(triangleElements[i].dNdX[2,:]) * ψ[dof]
    end
    return velX, velY
end

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
