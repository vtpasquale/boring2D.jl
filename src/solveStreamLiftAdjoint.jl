using SparseArrays
using TOML

function solveStreamLiftAdjoint(inputFileName::AbstractString)
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
    
     # Assemble domain and boundary matricies
     K = boring2D.assembleConvectionStiffness(mesh,triangleElements,1.0)
     f, s, psi, esa, dpsi = boring2D.streamBoundaryConditionsDiffAlpha(mesh,inputData)
     
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
     db = [-K[f,s]*dpsi[s]; 0.0] 
    
     x = A\b
     psi[f] = x[1:end-1]
    
     # compute lift coefficent from psi
     psiToCirculationMap = boring2D.computePsiToCirculationMap(inputData["mesh"]["airfoilBoundaryID"],mesh,triangleElements)
     circulation = psiToCirculationMap*psi
     Cl = 2*circulation[1]/inputData["freestream"]["velocity"] # Kutta-Joukowski with chord=1
    
     # Adjoint Solution
     xToCirculationMap = zeros(nDof+1,1);
     xToCirculationMap[1:end-1]=psiToCirculationMap
     xToLift = (2/inputData["freestream"]["velocity"])* xToCirculationMap  # Kutta-Joukowski with chord=1
     fa = falses(nDof+1)
     fa[1:end-1] .= f
     fa[end] = true
     lam = transpose(A)\xToLift[fa]
     dCl = transpose(lam)*db
     dCl_deg = dCl *pi/180
    
     println(Cl)
     println(dCl_deg)
    
    lamdaOut = zeros(nDof+1,1)
    lamdaOut[fa] = lam
    
     # Create solution dictionaries
     pointOutput = Dict("psi"=>psi,"lambda"=>lamdaOut[1:end-1])
     cellOutput = Dict()
     
     # Write output data to file
     boring2D.writeSolution("streamSolution.vtu",mesh,pointOutput,cellOutput)

    return Cl, dCl_deg
end

function streamBoundaryConditionsDiffAlpha(mesh::Mesh2D,inputData::Dict)    
    
    # Velocities
    Vinf = inputData["freestream"]["velocity"] # Moved to input argument to enable sensitivities
    alphaDeg = inputData["freestream"]["alphaDeg"]
    u0 = Vinf*cosd(alphaDeg)
    v0 = Vinf*sind(alphaDeg)

    du0 = -Vinf*sind(alphaDeg) # I think this derivative wrt alpha in radians, not degres
    dv0 =  Vinf*cosd(alphaDeg)

    # Preallocate arrays
    nDof = size(mesh.nodes,1)
    f = trues(nDof) # free nodes
    psi = zeros(nDof,1)
    esa = zeros(nDof,1) # Vector of ones at fairfield dof

    dpsi = zeros(nDof,1)

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
    dpsi[farfieldEdgeNodes] .= dv0.*xs - du0.*ys
    esa[farfieldEdgeNodes] .= 1.0

    return f, s, psi, esa, dpsi
end