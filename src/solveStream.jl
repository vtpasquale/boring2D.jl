using SparseArrays
using TOML
using ForwardDiff

function solveStream(inputFileName::AbstractString)
    # inputFileName = "n0012_input.toml"
    # using SparseArrays
    # using TOML
    
    # Read input file
    inputData = TOML.parsefile(inputFileName)
    
    # Echo input file
    TOML.print(inputData)
    
    # Output dictionaries
    textOutput =  Dict() # Dict("KuttaJoukowski" => KuttaJoukowskiTextOut,"IntegratedPressure" => IntegratedPressureTextOut, "Adjoint" => AdjointTextOut)
    # pointOutput = Dict() # Dict("psi"=>ψ,"lamda"=>lamdaOut[1:end-1])
    # cellOutput = Dict() # Dict("vX"=>[velX],"vY"=>[velY],"dVelX"=>[dVelX],"dVelY"=>[dVelY],"velMag"=>[velMag],"Cp"=>[Cp])

    # Process mesh
    mesh = boring2D.readMesh(inputData["mesh"]["fileName"])
    nDof = size(mesh.nodes,1)
    nEle = size(mesh.triangles,1)
    triangleElements = boring2D.TriangleElements(mesh)
    airfoilBoundary, farfieldBoundary, f, s = boring2D.processMeshBoundaries(mesh,inputData)
    gxγ = boring2D.computeGxγ(mesh,triangleElements,airfoilBoundary,f)
    
    # Assemble coefficient matrix
    K = boring2D.assembleConvectionStiffness(mesh,triangleElements,1.0)
    g = boring2D.assembleKuttaConditionConstraint(mesh,inputData,triangleElements)
    esa = boring2D.assembleEsa(nDof,farfieldBoundary)
    A = boring2D.assembleA(K,g,esa,f,s)
    
    # Assemble right-hand side
    ψsk = boring2D.assembleψsk(nDof,farfieldBoundary,inputData)
    b = boring2D.assembleB(K,ψsk,f,s)
    
    # Stream function solution
    x = convert(Matrix{Float64},A)\convert(Vector{Float64},b)
    
    # Stream function values
    ψ = boring2D.recoverψ(x,ψsk,esa,f,s)
    
    # Lift coefficient from Kutta-Joukowski
    Cl = (2/inputData["freestream"]["velocity"])* gxγ*x
    textOutput["KuttaJoukowski"] = Dict("Cl" => Cl)

    # Velocity and Cp
    velX, velY = boring2D.recoverVelocities(ψ,mesh,triangleElements)
    velMag = sqrt.(velX.^2 + velY.^2)
    Cp = 1.0 .- (velMag ./ inputData["freestream"]["velocity"]).^2

    # Integrate pressure coefficients
    textOutput["IntegratedPressure"] = boring2D.computeCfFromCp(Cp,mesh,airfoilBoundary,inputData)

    # Volume output dictionaries
    pointOutput = Dict("psi"=>ψ)
    cellOutput = Dict("vX"=>[velX],"vY"=>[velY],"velMag"=>[velMag],"Cp"=>[Cp])

    # Adjoint solution
    ∂Cl∂x = (2/inputData["freestream"]["velocity"])* gxγ
    λ = transpose(A)\transpose(∂Cl∂x)
    
    # Adjoint contour output
    fx = boring2D.getFxfromF(f)
    lamdaOut = zeros(nDof+1)
    lamdaOut[fx] = λ
    cellOutput["λ"] = [lamdaOut]
    
    # Sensitivity wrt α
    ∂ψsk_∂α = boring2D.assemble_∂ψsk_∂α(nDof,farfieldBoundary,inputData)
    ∂r_∂α = boring2D.assemble_∂r_∂α(K,∂ψsk_∂α,f,s)
    ∂Cl_∂α = transpose(λ)*∂r_∂α # per radian | ∂Cl_∂α_deg = ∂Cl_∂α*pi/180
    textOutput["Adjoint"] = Dict("∂Cl_∂α_rad" => ∂Cl_∂α,"∂Cl_∂α_deg" => ∂Cl_∂α*pi/180)
    
    # Shape sensitivities
    ∂g_∂b, ∂r_∂b = boring2D.computeShapePartials(mesh,x,inputData,size(airfoilBoundary.nodeIDs,1))
    ∂Cl_∂b = ∂g_∂b - transpose(λ)*∂r_∂b
    
    # Length-normalized shape sensitivities to surface output file
    boring2D.surfaceOutput2Vtk(inputData["output"]["surfaceFile"],airfoilBoundary,∂Cl_∂b[:]./airfoilBoundary.nodeLength)
    
    # Recover element output
    dVelX, dVelY = boring2D.recoverVelocities(lamdaOut,mesh,triangleElements)
    cellOutput["dVelX"] = [dVelX]
    cellOutput["dVelY"] = [dVelY]

    # Text output dictionary and output file
    TOML.print(textOutput)
    open(inputData["output"]["tomlFile"], "w") do io
        TOML.print(io, textOutput)
    end
    
    # Volume output dictionaries and output file
    boring2D.writeSolution(inputData["output"]["volumeFile"],mesh,pointOutput,cellOutput)

    return 0
end

function computeCirculationMatrixAndResidual(shapeVar::Vector,mesh0::Mesh2D,x::Vector{Float64},inputData::Dict)
    # This function is intended for shape sensitivity analysis only. 

    # Covert mesh nodes array to high-level type supported by dual
    nDof = size(mesh0.nodes,1)
    meshNodes = Array{Number}(undef, (nDof,2))
    meshNodes[:,1] = mesh0.nodes[:,1] 
    meshNodes[:,2] = mesh0.nodes[:,2]

    airfoilBoundary1, farfieldBoundary1, f1, s1 = boring2D.processMeshBoundaries(mesh0,inputData)

    # Deform mesh normals based on design variables
    if size(shapeVar,1) != size(airfoilBoundary1.nodeIDs,1)
        error("Check design variable definition. size(shapeVar,1) != size(airfoilBoundary.nodeIDs,1)")
    end
    meshNodes[airfoilBoundary1.nodeIDs,1] += shapeVar.*airfoilBoundary1.nodeNormal[1,:]
    meshNodes[airfoilBoundary1.nodeIDs,2] += shapeVar.*airfoilBoundary1.nodeNormal[2,:]

    # New dual-number-compatible mesh with design perturbations
    mesh = boring2D.Mesh2D(meshNodes,mesh0.edges,mesh0.triangles)

    triangleElements = boring2D.TriangleElements(mesh)
    airfoilBoundary, farfieldBoundary, f, s = boring2D.processMeshBoundaries(mesh,inputData);
    
    # cirulation matrix
    gxγ = boring2D.computeGxγ(mesh,triangleElements,airfoilBoundary,f);

    # Assemble coefficent matrix
    K = boring2D.assembleConvectionStiffness(mesh,triangleElements,1.0)
    g = boring2D.assembleKuttaConditionConstraint(mesh,inputData,triangleElements)
    esa = boring2D.assembleEsa(nDof,farfieldBoundary)
    A = boring2D.assembleA(K,g,esa,f,s)

    # Assemble right-hand side
    ψsk = boring2D.assembleψsk(nDof,farfieldBoundary,inputData)
    b = boring2D.assembleB(K,ψsk,f,s)

    # Residual
    r = Array(A*x-b)

    return gxγ, r
end

function computeCirculationMatrix(shapeVar::Vector,mesh0::Mesh2D,inputData::Dict)
    # This function is intended for shape sensitivity analysis only. 

    # Covert mesh nodes array to high-level type supported by dual
    nDof = size(mesh0.nodes,1)
    meshNodes = Array{Number}(undef, (nDof,2))
    meshNodes[:,1] = mesh0.nodes[:,1] 
    meshNodes[:,2] = mesh0.nodes[:,2]

    airfoilBoundary1, farfieldBoundary1, f1, s1 = boring2D.processMeshBoundaries(mesh0,inputData)

    # Deform mesh normals based on design variables
    if size(shapeVar,1) != size(airfoilBoundary1.nodeIDs,1)
        error("Check design variable definition. size(shapeVar,1) != size(airfoilBoundary.nodeIDs,1)")
    end
    meshNodes[airfoilBoundary1.nodeIDs,1] += shapeVar.*airfoilBoundary1.nodeNormal[1,:]
    meshNodes[airfoilBoundary1.nodeIDs,2] += shapeVar.*airfoilBoundary1.nodeNormal[2,:]

    # New dual-number-compatible mesh with design perturbations
    mesh = boring2D.Mesh2D(meshNodes,mesh0.edges,mesh0.triangles)

    triangleElements = boring2D.TriangleElements(mesh)
    airfoilBoundary, farfieldBoundary, f, s = boring2D.processMeshBoundaries(mesh,inputData);
    
    # cirulation matrix
    gxγ = boring2D.computeGxγ(mesh,triangleElements,airfoilBoundary,f)

    return gxγ
end

function computeResidual(shapeVar::Vector,mesh0::Mesh2D,x::Vector{Float64},inputData::Dict)
    # This function is intended for shape sensitivity analysis only. 

    # Covert mesh nodes array to high-level type supported by dual
    nDof = size(mesh0.nodes,1)
    meshNodes = Array{Number}(undef, (nDof,2))
    meshNodes[:,1] = mesh0.nodes[:,1] 
    meshNodes[:,2] = mesh0.nodes[:,2]

    airfoilBoundary1, farfieldBoundary1, f1, s1 = boring2D.processMeshBoundaries(mesh0,inputData)

    # Deform mesh normals based on design variables
    if size(shapeVar,1) != size(airfoilBoundary1.nodeIDs,1)
        error("Check design variable definition. size(shapeVar,1) != size(airfoilBoundary.nodeIDs,1)")
    end
    meshNodes[airfoilBoundary1.nodeIDs,1] += shapeVar.*airfoilBoundary1.nodeNormal[1,:]
    meshNodes[airfoilBoundary1.nodeIDs,2] += shapeVar.*airfoilBoundary1.nodeNormal[2,:]

    # New dual-number-compatible mesh with design perturbations
    mesh = boring2D.Mesh2D(meshNodes,mesh0.edges,mesh0.triangles)

    triangleElements = boring2D.TriangleElements(mesh)
    airfoilBoundary, farfieldBoundary, f, s = boring2D.processMeshBoundaries(mesh,inputData);

    # Assemble coefficent matrix
    # Converting sparse to dense for dual-number compatiblity - not a great solution
    K = Array( boring2D.assembleConvectionStiffness(mesh,triangleElements,1.0) )
    g = Array( boring2D.assembleKuttaConditionConstraint(mesh,inputData,triangleElements) )
    esa = Array( boring2D.assembleEsa(nDof,farfieldBoundary) )
    A = boring2D.assembleA(K,g,esa,f,s)

    # Assemble right-hand side
    ψsk = boring2D.assembleψsk(nDof,farfieldBoundary,inputData)
    b = boring2D.assembleB(K,ψsk,f,s)

    # Residual
    r = Array(A*x-b)

    return r
end

function computeShapePartialsFD(mesh0::Mesh2D,x::Vector{Float64},inputData::Dict,nVars::Int64)
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


function computeShapePartials(mesh0::Mesh2D,x::Vector{Float64},inputData::Dict,nVars::Int64)
    # Assemble partial derivatives w.r.t. shape design varaiables using automatic differentiation

    # Create unary functions for derivative calculation
    gxγ(shapeVar::Vector) = computeCirculationMatrix(shapeVar,mesh0,inputData)
    r(shapeVar::Vector)   = computeResidual(shapeVar,mesh0,x,inputData)

    # Automatic differentiation
    shapeVars = zeros(nVars)
    ∂gxγ_∂b = ForwardDiff.jacobian(gxγ,shapeVars)
    ∂r_∂b =   ForwardDiff.jacobian(r,  shapeVars)

    ∂g_∂b = (2/inputData["freestream"]["velocity"]) * ( transpose(∂gxγ_∂b) *x )

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
    g = Vector{Number}(undef,nDof)
    g[:] .= 0.0
    dof = mesh.triangles[zeroVyEle,:]
    coefficents = triangleElements[zeroVyEle].dNdX[1,:]
    g[dof] = coefficents
    return g
end

function assembleEsa(nDof::Int64,farfieldBoundary::ClosedBoundary2D)
    esa = Vector{Number}(undef,nDof)
    esa[:] .= 0.0
    esa[farfieldBoundary.nodeIDs] .= 1.0
    return esa
end

function assembleA(K::AbstractMatrix,g::AbstractVector,esa::AbstractVector,f::BitVector,s::BitVector)
    keas = K[f,s]*esa[s]
    A = [K[f,f] keas; transpose(g[f]) 0.0]
    return A
end

function assembleψsk(nDof::Int64,farfieldBoundary::ClosedBoundary2D,inputData::Dict)     
    Vinf = inputData["freestream"]["velocity"] # Moved to input argument to enable sensitivities
    alphaDeg = inputData["freestream"]["alphaDeg"]
    u0 = Vinf*cosd(alphaDeg)
    v0 = Vinf*sind(alphaDeg)
    ψsk = Vector{Number}(undef,nDof)
    ψsk[:] .= 0.0
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
    ∂ψsk_∂α = Vector{Number}(undef,nDof)
    ∂ψsk_∂α[:] .= 0.0
    xs = farfieldBoundary.nodeLocation[:,1]
    ys = farfieldBoundary.nodeLocation[:,2]
    ∂ψsk_∂α[farfieldBoundary.nodeIDs] .= dv0.*xs - du0.*ys
    return ∂ψsk_∂α
end

function assembleB(K::AbstractArray,ψsk::AbstractVector,f::BitVector,s::BitVector)
    b =  [-K[f,s]*ψsk[s]; 0.0] 
    return b
end

function assemble_∂r_∂α(K::AbstractMatrix,∂ψsk_∂α::AbstractVector,f::BitVector,s::BitVector)
    ∂r_∂α =  [-K[f,s]*∂ψsk_∂α[s]; 0.0] 
    return ∂r_∂α
end

function recoverψ(x::Vector{Float64},ψsk::AbstractArray,esa::AbstractArray,f::BitVector,s::BitVector)
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
