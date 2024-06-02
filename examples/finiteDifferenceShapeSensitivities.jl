using Revise
using boring2D
using SparseArrays
using TOML

cd("examples")

function computeClForFiniteDifference(shapeVar::Vector,mesh0::boring2D.Mesh2D,inputData::Dict)

    # Process mesh
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

    return Cl, x, A
end

#%% --------------------------------------

inputFileName = "n0012.toml"

# Read input file
inputData = TOML.parsefile(inputFileName)

# Echo input file
TOML.print(inputData)

# Process mesh
mesh = boring2D.readMesh(inputData["mesh"]["fileName"])
airfoilBoundary, farfieldBoundary, f, s = boring2D.processMeshBoundaries(mesh,inputData)

# Check residual for zero-value design variables
nVars = size(airfoilBoundary.nodeIDs,1)
shapeVar = zeros(nVars)
Cl0, x0, A0 = computeClForFiniteDifference(shapeVar,mesh,inputData)

# Compute ∂Cl_∂b using finite difference - automatic differntiation has unresolved issues
∂Cl_∂b = zeros(nVars)
∂r_∂b = zeros(size(x0,1),nVars)

Δ = (1.0/20.0)*minimum( airfoilBoundary.nodeLength )
for i = 1:nVars
    shapeVar = zeros(nVars)
    shapeVar[i] = Δ
    ClFD,xFD = computeClForFiniteDifference(shapeVar,mesh,inputData)
    ∂Cl_∂b[i] = (ClFD-Cl0) ./ Δ
    ∂x_∂b = (xFD-x0) ./ Δ
    ∂r_∂b[:,i] = A0*∂x_∂b
end

boring2D.surfaceOutput2Vtk("FiniteDiffSurfaceSensitivity.vtk",airfoilBoundary,∂Cl_∂b[:]./airfoilBoundary.nodeLength)